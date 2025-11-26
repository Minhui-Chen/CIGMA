import gzip, re, os
import subprocess
import pandas as pd


def main():
    for bim, annot_all_f, all_f in zip(snakemake.input.bim, snakemake.output.annot_all, snakemake.output.all):
        tmp_all = re.sub(r'\.gz$', '.tmp.gz', annot_all_f)
        tmp_bim = re.sub(r'\.gz$', '.tmp.bim', annot_all_f)
        with open(tmp_bim, 'w') as f:
            f.write('CHR\tBP\tSNP\tCM\n')
            for line in open(bim):
                f.write('\t'.join([line.split()[0], line.split()[3], line.split()[1], line.split()[2]]) + '\n')

        proc = subprocess.Popen(['python', 'ldsc/make_annot.py', '--gene-set-file', snakemake.input.all,
                '--gene-coord-file', snakemake.input.gene_coord, '--windowsize', snakemake.wildcards.window,
                '--bimfile', bim, '--annot-file', tmp_all], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True)
        stdout, stderr = proc.communicate()
        print("STDOUT:\n", stdout)
        print("STDERR:\n", stderr)
        print("Return code:", proc.returncode)

        # add snp info
        with gzip.open(annot_all_f, 'wt') as f_out:
            f_ins = gzip.open(tmp_all, 'rt').readlines()
            f_bims = open(tmp_bim, 'r').readlines()
            for line1, line2 in zip(f_bims, f_ins):
                f_out.write(line1.strip() + '\t' + line2.strip() + '\n')
        
        os.remove(tmp_all)

        proc = subprocess.Popen(['python', 'ldsc/ldsc.py', '--l2', '--bfile', os.path.splitext(bim)[0],
                '--ld-wind-cm', '1', '--annot', annot_all_f, '--out', re.sub(r'\.l2\.ldscore\.gz$', '', all_f),
                '--print-snps', snakemake.input.hapmap3], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True)
        stdout, stderr = proc.communicate()
        print("STDOUT:\n", stdout)
        print("STDERR:\n", stderr)
        print("Return code:", proc.returncode)

    genes = pd.read_table(snakemake.input.genes)
    cts = genes['ct'].unique()
    with open(snakemake.output.ldcts, 'w') as f:
        for ct in cts:
            ct_genes = genes.loc[genes['ct'] == ct, 'gene']
            ct = ct.replace(" ", "")
            ct_genes.to_csv(snakemake.output.ct_genes, sep='\t', index=False, header=False)

            for chr, bim, annot_all_f, annot_genes_f, genes_f in zip(snakemake.params.chrs, snakemake.input.bim, snakemake.output.annot_all,
                                                        snakemake.params.annot_genes, snakemake.params.genes):
                tmp_genes = re.sub(r'\.gz$', '.tmp.gz', annot_genes_f)
                tmp_bim = re.sub(r'\.gz$', '.tmp.bim', annot_all_f)

                proc = subprocess.Popen(['python', 'ldsc/make_annot.py', '--gene-set-file', snakemake.output.ct_genes,
                        '--gene-coord-file', snakemake.input.gene_coord, '--windowsize', snakemake.wildcards.window,
                        '--bimfile', bim, '--annot-file', tmp_genes], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True)
                stdout, stderr = proc.communicate()
                print("STDOUT:\n", stdout)
                print("STDERR:\n", stderr)
                print("Return code:", proc.returncode)

                # add snp info
                with gzip.open(annot_genes_f, 'wt') as f_out:
                    f_ins = gzip.open(tmp_genes, 'rt').readlines()
                    f_bims = open(tmp_bim, 'r').readlines()
                    for line1, line2 in zip(f_bims, f_ins):
                        f_out.write(line1.strip() + '\t' + line2.strip() + '\n')

                os.remove(tmp_genes)

                proc = subprocess.Popen(['python', 'ldsc/ldsc.py', '--l2', '--bfile', os.path.splitext(bim)[0],
                        '--ld-wind-cm', '1', '--annot', annot_genes_f, '--out', re.sub(str(chr) + r'\.l2\.ldscore\.gz$', f'{ct}.{chr}', genes_f),
                        '--print-snps', snakemake.input.hapmap3], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True)
                stdout, stderr = proc.communicate()
                print("STDOUT:\n", stdout)
                print("STDERR:\n", stderr)
                print("Return code:", proc.returncode)

            f.write('{}\t{},{}\n'.format(ct, re.sub('1.l2.ldscore.gz$', f'{ct}.', snakemake.params.genes[0]),
                                        re.sub('1.l2.ldscore.gz$', '', snakemake.output.all[0])))


if __name__ == '__main__':
    main()