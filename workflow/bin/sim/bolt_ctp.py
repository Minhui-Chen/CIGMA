import shutil
import sys, os, re, tempfile, subprocess
import numpy as np, pandas as pd
from cigma import util


def main():
    #
    data = np.load(sys.argv[1], allow_pickle=True).item()
    outs = []

    tmpdir = tempfile.mkdtemp()
    tmp = os.path.join(tmpdir, 'tmp')
    # tmp = os.path.join('bolt_test', 'tmp')

    # identifiers
    chromosome = 1
    ref_allele = "A"
    alt_allele = "G"

    for key in data.keys():
        C = data[key]['P'].shape[1]
        cts = ['CT' + str(c+1) for c in range(C)]
        n_ind = data[key]['G'].shape[0]
        n_snp = data[key]['G'].shape[1]

        # identifiers
        individual_ids = ["ind" + str(i+1) for i in range(n_ind)]
        snp_ids = ["rs" + str(i+1) for i in range(n_snp)]

        # file name
        prefix = tmp + f'.{key}'
        pheno_f = prefix + '.pheno'
        ped_f = prefix + '.ped'
        map_f = prefix + '.map'
        # covar_f = prefix + '.covar'

        # prepare pheno files
        y = data[key]['Y']
        with open(pheno_f, "w") as phenef:
            head = '\t'.join(['FID', 'IID'] + cts) + '\n'
            phenef.write(head)

            for i in range(n_ind):
                phenef.write(f"{individual_ids[i]}\t{individual_ids[i]}")
                for c in range(C):
                    phenef.write(f"\t{y[i,c]}")
                phenef.write("\n")


        # prepare geno files
        rawG = data[key]['rawG']
        with open(ped_f, "w") as pedf:
            for i, ind in enumerate(individual_ids):
                pedf.write(f"{ind}\t{ind}\t0\t0\t0\t-9\t")
                alleles = []
                for g in rawG[i]:
                    if g == 0:
                        alleles += [ref_allele, ref_allele]
                    elif g == 1:
                        alleles += [ref_allele, alt_allele]
                    elif g == 2:
                        alleles += [alt_allele, alt_allele]
                    else:
                        alleles += ["0", "0"]
                pedf.write(" ".join(alleles) + "\n")


        with open(map_f, "w") as mapf:
            for i, snp in enumerate(snp_ids):
                bp_pos = i + 1  # if no real positions
                mapf.write(f"{chromosome}\t{snp}\t0\t{bp_pos}\n")

        util.subprocess_popen(['plink', '--file', prefix, '--make-bed', '--out', prefix])

        # # prepare qcovar file
        # P = data[key]['P']
        # with open(covar_f, "w") as covarf:
        #     head = '\t'.join(['FID', 'IID'] + cts[:-1]) + '\n'
        #     covarf.write(head)

        #     for i in range(n_ind):
        #         covarf.write(f"{individual_ids[i]}\t{individual_ids[i]}")
        #         for c in range(C-1):
        #             covarf.write(f"\t{P[i,c]}")
        #         covarf.write("\n")


        # bolt command
        cmd = ['bolt', '--bfile=' + prefix, '--phenoFile=' + pheno_f]
        cmd += [f'--phenoCol=CT{c}' for c in range(1, C+1)]
        # cmd += ['--covarFile=' + covar_f, '--qCovarCol=CT{1:' + f'{C-1}' + '}']
        cmd += ['--reml']

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            print("=== STDOUT ===")
            print(result.stdout)
            print("=== STDERR ===")
            print(result.stderr)
        except subprocess.CalledProcessError as e:
            stderr = e.stderr or ""
            print("=== STDOUT (before error) ===")
            print(e.stdout or "")
            print("=== STDERR ===")
            print(stderr)
            if "ERROR: Heritability estimate is close to 0;" in stderr or "std::runtime_error" in stderr or 'boost::numeric::ublas::internal_logic' in stderr:
                print(f"Skipping gene {key}.")
                outs.append({'gene': key, 'shared_h2': np.nan, 'specific_h2': np.nan,
                             'sigma_g2': np.full((C, C), fill_value=np.nan),
                             'sigma_e2': np.full((C, C), fill_value=np.nan),
                             'hom_g2': np.nan, 'vbar': np.nan,
                             'hom_e2': np.nan, 'wbar': np.nan,
                             })
                continue
            else:
                print(f"Error running BOLT-LMM for gene {key}: ")
                print(e)
                raise

        # extract results
        sigma_e2 = np.full((C, C), fill_value=np.nan)
        sigma_g2 = np.full((C, C), fill_value=np.nan)
        sigma = np.full(C, fill_value=np.nan)
         
        for line in result.stdout.splitlines():
            if re.search(rf"^Variance component 0: \[{C},{C}\]", line.strip()):
                content = re.search(r'\(\((.*)\)\)', line).group(1)
                rows = [list(map(float, row.split(','))) for row in content.split('),(')]
                sigma_e2 = np.array(rows)
            elif re.search(rf"^Variance component 1: \[{C},{C}\]", line.strip()):
                content = re.search(r'\(\((.*)\)\)', line).group(1)
                rows = [list(map(float, row.split(','))) for row in content.split('),(')]
                sigma_g2 = np.array(rows)
            for i in range(C):
                if re.search(f"^Phenotype {i+1} variance sigma2:", line.strip()):
                    sigma[i] = float(line.strip().split()[-2])
        assert not np.isnan(sigma_e2).any()
        assert not np.isnan(sigma_g2).any()
        assert not np.isnan(sigma).any()

        # scale variance to raw scale
        sigma_e2 = sigma_e2 * np.sqrt(np.outer(sigma, sigma))
        sigma_g2 = sigma_g2 * np.sqrt(np.outer(sigma, sigma))

        # compute shared vs specific variance
        hom_g2 = sigma_g2[np.triu_indices(C, k=1)].mean()
        vbar = np.diag(sigma_g2).mean() - hom_g2
        hom_e2 = sigma_e2[np.triu_indices(C, k=1)].mean()
        wbar = np.diag(sigma_e2).mean() - hom_e2

        # calculate h2
        t = hom_g2 + vbar + hom_e2 + wbar
        shared_h2 = hom_g2 / t
        specific_h2 = vbar / t


        # save
        outs.append({'gene': key, 'shared_h2': shared_h2, 'specific_h2': specific_h2,
                     'sigma_g2': sigma_g2, 'sigma_e2': sigma_e2,
                     'hom_g2': hom_g2, 'vbar': vbar,
                     'hom_e2': hom_e2, 'wbar': wbar,
                     })

    shutil.rmtree(tmpdir)
    
    np.save(sys.argv[2], outs)


if __name__ == '__main__':
    main()
