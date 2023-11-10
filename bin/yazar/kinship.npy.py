import sys, os, re, shutil, gzip
import numpy as np, pandas as pd
from gxctmm import util

def main():
    # 
    genes_f, P_f = sys.argv[1], sys.argv[2]
    r = int(sys.argv[3])
    bed_f = sys.argv[4]
    chr = int(sys.argv[5])
    kinship_f = sys.argv[6]

    # read
    genes = pd.read_table(genes_f)
    P = pd.read_table(P_f, index_col=0)
    bfile = os.path.splitext(bed_f)[0]

    # filter variants and inds
    tmpfn = util.generate_tmpfn()
    chr_prefix = f'{tmpfn}.chr{chr}'
    util.extract_vcf(bfile, samples=P.index.tolist(), maf='0.05', output_bfile=chr_prefix, update_bim=False)

    kinship = {'gene': [], 'K': [], 'nsnp': []}
    # make kinship for each gene
    for _, row in genes.loc[genes['chr']==chr].iterrows():
        gene, start, end = row['feature'], row['start'], row['end']
        kinship_prefix = f'{chr_prefix}.{gene}'
        
        # grm
        nsnp = util.grm(chr_prefix, chr, start, end, r, kinship_prefix, 'plink', format='gz')

        # collect grm
        if nsnp > 0:
            K = []
            for line in gzip.open(f'{kinship_prefix}.rel.gz', 'rt'):
                K += line.strip().split()

            kinship['gene'].append(gene)
            kinship['K'].append(np.array(K).astype('float'))
            kinship['nsnp'].append(nsnp)

            # read id
            ids = np.array([line.strip().split()[0] for line in open(f'{kinship_prefix}.rel.id')])
            if 'ids' not in kinship.keys():
                kinship['ids'] = ids
            else:
                if np.any(kinship['ids'] != ids):
                    sys.exit('IDs not matching!')

    np.save(kinship_f, kinship)


if __name__ == '__main__':
    main()
