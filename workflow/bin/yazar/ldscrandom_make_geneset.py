import sys
import numpy as np
import pandas as pd

from cigma import util


def main():
    out = np.load(snakemake.input.out, allow_pickle=True).item()
    df = util.read_out(out, ['gene', 'p:free:V', 'p:free:hom_g2'])
    df['var_beta'] = np.var(out['free']['ct_beta'], axis=1)

    # gcta
    gcta = pd.read_table(snakemake.input.gcta)
    gcta = gcta.rename(columns={'p': 'gcta_p'})
    df = df.merge(gcta[['gene', 'gcta_p']], on='gene')

    # remove hla genes: https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37
    gene_df = pd.read_table(snakemake.input.location)
    gene_df = gene_df.loc[~((gene_df['chr'] == snakemake.params.chr) & (gene_df['end'] > snakemake.params.start) 
                            & (gene_df['start'] < snakemake.params.end))]
    df = df.loc[df['gene'].isin(gene_df['feature'])]
    gene_df = gene_df.loc[gene_df['feature'].isin(df['gene'])]

    # control
    np.savetxt(snakemake.output.all, df['gene'], fmt='%s')

    # number of genes to select
    ngene = int(snakemake.wildcards.ngene)
    if ngene == 0:
        ngene = (df['p:free:V'] < (0.05 / len(out['gene']))).sum()

    # find top genes
    if snakemake.wildcards.gset == 'var':
        target_genes = df.sort_values('p:free:V')['gene'].values[:ngene]
    elif snakemake.wildcards.gset == 'shared':
        target_genes = df.sort_values('p:free:hom_g2')['gene'].values[:ngene]
    elif snakemake.wildcards.gset == 'mean':
        target_genes = df.sort_values('var_beta', ascending=False)['gene'].values[:ngene]
    elif snakemake.wildcards.gset == 'gcta':
        target_genes = df.sort_values('gcta_p')['gene'].values[:ngene]
    else:
        sys.exit(f'Unknown gset: {snakemake.wildcards.gset}')
    np.savetxt(snakemake.output.target, np.array(target_genes), fmt='%s')

    # op
    ctp = pd.read_table(snakemake.input.ctp, index_col=(0,1))
    P = pd.read_table(snakemake.input.P, index_col=0)
    ops = []
    for gene in df['gene'].values:
        gene_ctp = ctp[gene].unstack()
        assert all(gene_ctp.index == P.index)
        assert all(gene_ctp.columns == P.columns)
        gene_op = (gene_ctp * P).sum(axis=1)
        gene_op.name = gene
        ops.append(gene_op)
    ops = pd.concat(ops, axis=1)
    mean_expr = ops.mean(axis=0)
    mean_expr.index.name = 'feature'
    mean_expr.name = 'mean_expr'
    mean_expr = mean_expr.reset_index()

    # gene length 
    gene_df = gene_df.merge(mean_expr, on='feature')
    gene_df['mean_expr_rank'] = gene_df['mean_expr'].rank(method='dense')
    gene_df['gene_length'] = gene_df['end'] - gene_df['start']
    assert all(gene_df['gene_length'] > 0)
    gene_df['gene_length_rank'] = gene_df['gene_length'].rank(method='dense')
    gene_df.to_csv(snakemake.output.genes, sep='\t', index=False)

    # nearby genes
    nearby_genes = []
    window = int(float(snakemake.wildcards.nearbywindow))
    for gene in target_genes:
        chr = gene_df.loc[gene_df['feature'] == gene, 'chr'].values[0]
        start = gene_df.loc[gene_df['feature'] == gene, 'start'].values[0]
        end = gene_df.loc[gene_df['feature'] == gene, 'end'].values[0]
        window_start = start - window
        window_end = end + window
        # overlap
        filter = (gene_df['chr'] == chr) & \
                    (((gene_df['start'] >= window_start) & (gene_df['start'] <= window_end)) | \
                    ((gene_df['end'] >= window_start) & (gene_df['end'] <= window_end)) | \
                    ((gene_df['start'] <= window_start) & (gene_df['end'] >= window_end)))
        nearby = gene_df.loc[filter, 'feature'].values
        nearby_genes.append(nearby)
    nearby_genes = np.unique(np.concatenate(nearby_genes))

    # random
    rng = np.random.default_rng(seed=snakemake.params.seed)
    k = len(snakemake.output.random)
    bin_size = int(snakemake.wildcards.binsize)
    randoms = []
    for gene in target_genes:
        rank = gene_df.loc[gene_df['feature'] == gene, 'mean_expr_rank'].values[0]
        length = gene_df.loc[gene_df['feature'] == gene, 'gene_length_rank'].values[0]
        filter = (gene_df['mean_expr_rank'] >= rank - bin_size) & (gene_df['mean_expr_rank'] <= rank + bin_size) & \
                    (gene_df['gene_length_rank'] >= length - bin_size) & (gene_df['gene_length_rank'] <= length + bin_size) & \
                    (~gene_df['feature'].isin(nearby_genes))
        candidates = gene_df.loc[filter, 'feature'].values
        if len(candidates) == 0:
            sys.exit(f'No candidates found for gene {gene} (rank: {rank}, length: {length})')
        selected = rng.choice(candidates, size=300 * k, replace=True)
        randoms.append(selected)
    randoms = np.array(randoms).T  # k x ngene

    # remove gene sets with duplicates
    filters = [len(set(randoms[i])) == len(randoms[i]) for i in range(randoms.shape[0])]
    randoms = randoms[filters]
    if randoms.shape[0] < len(snakemake.output.random):
        print(randoms.shape)
        sys.exit(f'Warning: only {randoms.shape[0]} random sets generated, need {len(snakemake.output.random)}')

    for i, f in enumerate(snakemake.output.random):
        np.savetxt(f, randoms[i], fmt='%s')


if __name__ == '__main__':
    main()