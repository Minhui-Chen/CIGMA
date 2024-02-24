import sys
import numpy as np
import pandas as pd


def find_closest_gene(mid, genes):
    # TODO: what about promoters within genes
    # TODO: what about overlapping genes that more than one genes are close to the promoter
    data = genes.copy()
    data['dis'] = (data['tss'] - mid).abs()

    # one promoter can have multiple closest genes
    return data.loc[data['dis'] == data['dis'].min(), ['feature', 'start', 'end', 'tss']]


def main():
    bed = pd.read_table(snakemake.input.bed, names=['chr', 'cre_start', 'cre_end', 'cre_id', 'id2', 'type'])
    # bed = bed.loc[bed['type'].str.contains('PLS') | bed['type'].str.contains('DNase-H3K4me3')]  # TODO: including DNase-H3K4me3, which is considered promoter in Moore. from Mostafavi Table S8, seems only PLS is used
    bed = bed.loc[bed['type'].str.contains('PLS')]  # TODO: including DNase-H3K4me3, which is considered promoter in Moore. from Mostafavi Table S8, seems only PLS is used
    bed = bed.drop(columns=['id2', 'type'])

    # autosome
    bed['chr'] = bed['chr'].str[3:]
    bed = bed.loc[bed['chr'].isin([str(x) for x in range(1, 23)])]
    bed['chr'] = bed['chr'].astype(int)
    bed['cre_start'] = bed['cre_start'] + 1

    # mid point of promoter
    bed['cre_mid'] = (bed['cre_start'] + bed['cre_end']) / 2

    # 
    genes = pd.read_table(snakemake.input.genes, usecols=['feature', 'chr', 'start', 'end', 'tss'])

    #
    data = []
    for chr in bed['chr'].unique():
        chr_bed = bed.loc[bed['chr'] == chr]
        chr_genes = genes.loc[genes['chr'] == chr]

        for index, row in chr_bed.iterrows():
            closest_gene = find_closest_gene(row['cre_mid'], chr_genes)
            closest_gene[row.index] = row.to_numpy()
            data.append(closest_gene)

    # 
    data = pd.concat(data, axis=0, ignore_index=True)
    data = data.drop(columns='cre_mid')
    data.to_csv(snakemake.output.promoter, sep='\t', index=False)


if __name__ == '__main__':
    main()