import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from cigma import preprocess


def main():
    # par
    max_cells = snakemake.wildcards.get('max_cells', None)
    prop_reads = snakemake.wildcards.get('prop_reads', None)
    transform = 'logp_cp10k'
    if snakemake.params.get('transform'):
        transform = snakemake.params.transform
    if snakemake.wildcards.get('transform'):
        transform = snakemake.wildcards.transform
    if snakemake.wildcards.get('L'):
        transform = snakemake.wildcards.L

    # filter genes
    var = pd.read_table(snakemake.input.var, index_col=0)
    if 'feature_is_filtered' in var.columns:
        genes = var.loc[~var['feature_is_filtered']].index.to_numpy()
    else:
        genes = var.index.to_numpy()

    if 'subset_gene' in snakemake.params.keys():
        # NOTE: just used for test since it would impact transformation
        # random select genes
        rng = np.random.default_rng(seed=int(snakemake.params.seed))
        genes = rng.choice(genes, snakemake.params.subset_gene, replace=False)

    # read
    ann = sc.read_h5ad(snakemake.input.h5ad)

    # exclude replicates
    obs = pd.read_table(snakemake.input.obs, index_col=0)
    ind_pool = np.unique(obs[snakemake.params.ind_col].astype('str')+'+'+obs[snakemake.params.pool_col].astype('str'))
    cells = ((~ann.obs[snakemake.params.ind_col].isna())
            & (~ann.obs[snakemake.params.ct_col].isna())
            & (ann.obs[snakemake.params.ind_col].astype('str')+'+'+ann.obs[snakemake.params.pool_col].astype('str')).isin(ind_pool))

    if max_cells is not None:
        # Select up to max_cells for ind-ct pairs with more than max_cells
        max_cells = int(max_cells)

        grouped = ann.obs[cells].groupby([snakemake.params.ind_col, snakemake.params.ct_col])
        filtered = grouped.filter(lambda x: len(x) > max_cells)
        sampled = filtered.groupby([snakemake.params.ind_col, snakemake.params.ct_col]).sample(n=max_cells, random_state=0) # NOTE: random state

        # For groups with 100 or fewer rows, keep all rows
        small_groups = grouped.filter(lambda x: len(x) <= max_cells)
        kept = pd.concat([sampled, small_groups])
        
        cells = cells & (ann.obs.index.isin(kept.index))

    data = ann[cells, genes]

    # transform (sparse matrix X only)
    X = data.X
    if prop_reads is not None:
        # proportion of reads to sample
        prop_reads = float(prop_reads)
        X = preprocess.sample_reads(X=X, p_reads=prop_reads)
    X = preprocess.transform(X=X, transform=transform)

    # save transformed data
    sparse.save_npz(snakemake.output.X, X)
    data.obs.rename_axis('cell').to_csv(snakemake.output.obs, sep='\t')
    data.var.rename_axis('feature').to_csv(snakemake.output.var, sep='\t')


if __name__ == '__main__':
    main()