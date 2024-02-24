import sys, os, re, shutil, gzip
import numpy as np, pandas as pd
from gxctmm import util


def cre_region(data, gene, chr):
        # promoter grm
        data = data.loc[(data['feature'] == gene) & (data['chr'] == chr)]
        if data.shape[0] == 1:
            cre_start, cre_end = data['cre_start'].tolist()[0], data['cre_end'].tolist()[0]
        else:
            cre_start, cre_end = data['cre_start'].tolist(), data['cre_end'].tolist()
        return cre_start, cre_end



def main():
    # 
    genes_f, ctp_f = sys.argv[1], sys.argv[2]
    r = int(sys.argv[3])
    promoter_f, enhancer_f = sys.argv[4], sys.argv[5]
    bed_f = sys.argv[6]
    chr = int(sys.argv[7])
    promoter_kinship_f, enhancer_kinship_f = sys.argv[8], sys.argv[9]

    # read
    genes = pd.read_table(genes_f, usecols=['feature', 'GeneSymbol', 'chr', 'start', 'end'])
    genes = genes.drop_duplicates()
    ctp = pd.read_table(ctp_f, index_col=(0, 1))
    bfile = os.path.splitext(bed_f)[0]
    promoters = pd.read_table(promoter_f)
    enhancers = pd.read_table(enhancer_f, names=['chr', 'cre_start', 'cre_end', 'feature', 'a', 'b'])

    # TODO: promoter 1kb to tss
    # TODO: is there large promoters that mid point is 1kb away from tss
    dis_cut = 1000
    promoters['cre_mid'] = (promoters['cre_start'] + promoters['cre_end']) / 2
    promoters = promoters.loc[(promoters['cre_mid'] - promoters['tss']).abs() < dis_cut]

    # enhancers: autosome
    enhancers = enhancers.loc[enhancers['chr'] != 'chrX']
    enhancers['chr'] = enhancers['chr'].str[3:].astype(int)
    enhancers['cre_start'] = enhancers['cre_start'] + 1
    # TODO: enhancers within r
    enhancers = enhancers.merge(genes, on=['chr', 'feature'])
    enhancers = enhancers.loc[(enhancers['cre_start'] > (enhancers['start'] - r)) &
                              (enhancers['cre_end'] < (enhancers['end'] + r))]

    # exclude genes not in ctp
    genes = genes.loc[genes['feature'].isin(ctp.columns)]
    print(genes.shape)
    # genes with promoter and enhancers
    genes = genes.loc[genes['feature'].isin(enhancers['feature'])]
    print(genes.shape)
    genes = genes.loc[genes['feature'].isin(promoters['feature'])]
    print(genes.shape)

    # filter variants and inds
    tmpfn = util.generate_tmpfn()
    chr_prefix = f'{tmpfn}.chr{chr}'
    # maf 1%
    maf = '0.01'
    util.extract_vcf(bfile, samples=ctp.index.get_level_values(0).tolist(), maf=maf, 
                     output_bfile=chr_prefix, update_bim=False)

    promoter_kinship = {'gene': [], 'K': [], 'nsnp': []}
    enhancer_kinship = {'gene': [], 'K': [], 'nsnp': []}
    # make kinship for each gene
    for _, row in genes.loc[genes['chr']==chr].iterrows():
        gene, start, end = row['feature'], row['start'], row['end']
        promoter_kinship_prefix = f'{chr_prefix}.{gene}.promoter'
        enhancer_kinship_prefix = f'{chr_prefix}.{gene}.enhancer'
        
        # promoter grm
        cre_start, cre_end = cre_region(promoters, gene, chr)
        promoter_nsnp = util.grm(chr_prefix, promoter_kinship_prefix, chr=chr, start=cre_start, end=cre_end, 
                        r=0, tool='plink', format='gz')

        # enhancer grm
        cre_start, cre_end = cre_region(enhancers, gene, chr)
        enhancer_nsnp = util.grm(chr_prefix, enhancer_kinship_prefix, chr=chr, start=cre_start, end=cre_end, 
                        r=0, tool='plink', format='gz')

        # collect grm
        if promoter_nsnp > 0 and enhancer_nsnp > 0:
            promoter_kinship['gene'].append(gene)
            enhancer_kinship['gene'].append(gene)

            # promoter
            K = []
            for line in gzip.open(f'{promoter_kinship_prefix}.rel.gz', 'rt'):
                K += line.strip().split()

            promoter_kinship['K'].append(np.array(K).astype('float'))
            promoter_kinship['nsnp'].append(promoter_nsnp)

            # enhancer
            K = []
            for line in gzip.open(f'{enhancer_kinship_prefix}.rel.gz', 'rt'):
                K += line.strip().split()

            enhancer_kinship['K'].append(np.array(K).astype('float'))
            enhancer_kinship['nsnp'].append(enhancer_nsnp)

            # read id
            ids = np.array([line.strip().split()[0] for line in open(f'{promoter_kinship_prefix}.rel.id')])
            enhancer_ids = np.array([line.strip().split()[0] for line in open(f'{enhancer_kinship_prefix}.rel.id')])

            if np.any(ids != enhancer_ids):
                sys.exit('IDs not matching')

            if 'ids' not in promoter_kinship.keys():
                promoter_kinship['ids'] = ids
            else:
                if np.any(promoter_kinship['ids'] != ids):
                    sys.exit('IDs not matching!')

            if 'ids' not in enhancer_kinship.keys():
                enhancer_kinship['ids'] = ids

    np.save(promoter_kinship_f, promoter_kinship)
    np.save(enhancer_kinship_f, enhancer_kinship)


if __name__ == '__main__':
    main()
