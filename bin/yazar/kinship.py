import sys, os
import numpy as np, pandas as pd
from gxctmm import util

def main():
    # 
    genes_f, P_f = sys.argv[1], sys.argv[2]
    r = int(sys.argv[3])
    bed_f = sys.argv[4]
    chr = int(sys.argv[5])
    kinship_tmp = sys.argv[6]
    kinship_f = sys.argv[7]

    # read
    genes = pd.read_table(genes_f)
    P = pd.read_table(P_f, index_col=0)
    bfile = os.path.splitext(bed_f)[0]

    # filter variants and inds
    chr_prefix = os.path.join(os.path.dirname(kinship_tmp), f'chr{chr}')
    util.extract_vcf(bfile, samples=P.index.tolist(), maf='0.05', output_bfile=chr_prefix, update_bim=False)

    with open(kinship_f, 'w') as f:
        f.write('gene\tK\tsnps\n')
    
        # make kinship for each gene
        bim = pd.read_csv(chr_prefix+'.bim', sep='\s+', names=['chr', 'snp', 'cm', 'bp','a1','a2'])
        for _, row in genes.loc[genes['chr']==chr].iterrows():
            gene, start, end = row['feature'], row['start'], row['end']
            kinship_prefix = os.path.join(os.path.dirname(kinship_tmp),gene)
            nsnp = util.grm(chr_prefix, chr, start, end, r, kinship_prefix)
            f.write(f'{gene}\t{kinship_prefix}.rel.bin\t{nsnp}\n')

    os.remove(chr_prefix+'.bed')
    os.remove(chr_prefix+'.bim')
    os.remove(chr_prefix+'.fam')

if __name__ == '__main__':
    main()
