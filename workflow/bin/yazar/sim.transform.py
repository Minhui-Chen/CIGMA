import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def main():
    # ind, cell types, genes
    ctp = pd.read_table(snakemake.input.ctp, index_col=(0, 1))
    cts = np.unique(ctp.index.get_level_values(1))
    inds = np.unique(ctp.index.get_level_values(0))
    genes = ctp.columns.to_numpy()

    obs = pd.read_table(snakemake.input.obs, index_col=0)
    ind_pool = np.unique(obs[snakemake.params.ind_col].astype('str')+'+'+obs[snakemake.params.pool_col].astype('str'))

    ann = sc.read_h5ad(snakemake.input.h5ad, backed='r')
    # ann = sc.read_h5ad(snakemake.input.h5ad)
    print(ann.shape)

    # filter cells
    cells = ((ann.obs[snakemake.params.ind_col].isin(inds))
            & (ann.obs[snakemake.params.ct_col].isin(cts))
            & (ann.obs[snakemake.params.ind_col].astype('str')+'+'+ann.obs[snakemake.params.pool_col].astype('str')).isin(ind_pool))
    data = ann[cells, genes]
    
    # extract matrix
    if sparse.issparse(ann.X):
        X = data.X
    else:
        # when X is dense, data.X give error: Only one indexing vector or array is currently allowed for fancy indexing
        X = ann[:, genes].X[cells]

    # save matrix
    print(X.shape)
    sparse.save_npz(snakemake.output.X, X)
    data.var.rename_axis('feature').to_csv(snakemake.output.var, sep='\t')

    # permute cells across cell types
    obs = data.obs.copy()
    if snakemake.params.seed is not None:
        rng = np.random.default_rng(seed=int(snakemake.params.seed))
        obs_grouped = obs.groupby(snakemake.params.ind_col, observed=True)
        obs[snakemake.params.ct_col] = obs_grouped[snakemake.params.ct_col].transform(
            lambda x: rng.permutation(x))
    else:
        pass
    obs.rename_axis('cell').to_csv(snakemake.output.obs, sep='\t')


if __name__ == '__main__':
    main()