import numpy as np, pandas as pd


def main():
    #
    genes = pd.read_table(snakemake.input.genes)
    gene_info = pd.read_table(snakemake.input.gene_info)
    gene_info = gene_info.rename(columns={'chromosome_name':'chr', 'start_position':'start', 'end_position':'end'})
    gene_info['tss'] = gene_info['transcript_start']
    gene_info.loc[gene_info['strand'] == -1, 'tss'] = gene_info.loc[gene_info['strand'] == -1, 'transcript_end']

    # merge
    genes = genes.merge(gene_info, left_on=['feature', 'GeneSymbol'], right_on=['ensembl_gene_id', 'hgnc_symbol'])
    genes = genes.drop(columns=['features', 'ensembl_gene_id', 'hgnc_symbol'])

    # drop 6 genes with identical location
    grouped = genes.groupby(['chr', 'start', 'end'])
    dup_genes = grouped.filter(lambda x: x['feature'].nunique() > 1)
    print(dup_genes)
    genes = genes.loc[~genes['feature'].isin(dup_genes['feature'])]

    #
    genes = genes.sort_values(by=['chr', 'start', 'end'])
    genes.to_csv(snakemake.output.meta, sep='\t', index=False)
    genes = genes.drop(columns=['entrezgene', 'transcript_start', 'transcript_end', 'strand', 'tss'])  # ENSG00000225880 has two entrezgene ids
    genes = genes.drop_duplicates()
    genes.to_csv(snakemake.output.location, sep='\t', index=False)


if __name__ == '__main__':
    main()

