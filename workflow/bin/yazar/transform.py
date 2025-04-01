import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from cigma import preprocess


def main():
    # par
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
    data = ann[cells, genes]

    # transform
    if sparse.issparse(ann.X):
        X = data.X
    else:
        # when X is dense, data.X give error: Only one indexing vector or array is currently allowed for fancy indexing
        X = ann[:, genes].X[cells]  # NOTE: transform fun may not work for dense matrix
    if transform == 'pearson':
        X = preprocess.transform(ann=data, transform=transform)
    else:
        X = preprocess.transform(X=X, transform=transform)

    # save transformed data
    sparse.save_npz(snakemake.output.X, X)
    data.obs.rename_axis('cell').to_csv(snakemake.output.obs, sep='\t')
    data.var.rename_axis('feature').to_csv(snakemake.output.var, sep='\t')


if __name__ == '__main__':
    main()