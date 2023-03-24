import numpy as np, pandas as pd
import scanpy as sc
from scipy import stats

def cal_ctnu(data: scipy.sparse._csr.csr_matrix, axis: int=0) -> np.array:
    '''
    compute ctnu 

    Parameters:
        data:   gene expression data for cells (rows) and genes (columns)
        axis:   the axis across cells
    Returns:
        ctnu
    '''
    ctp = np.array(data.mean(axis=axis))[0]
    ctp2 = np.array(data.power(2).mean(axis=axis))[0]
    ctnu = (ctp2 - ctp**2) / data.shape[0]
    return( ctnu )

def main():
    # 
    ann = sc.read_h5ad(snakemake.input.h5ad, backed='r')
    obs = ann.obs
    obs = obs.rename(columns={snakemake.params.ind_col:'ind', snakemake.params.ct_col:'ct'})
    genes = pd.read_table(snakemake.input.genes, sep='\t')
    genes = genes.loc[~genes['feature_is_filtered'], 'feature'].to_numpy()

    # random select genes
    rng = np.random.default_rng(seed=snakemake.params.seed)
    genes = rng.choice(genes, snakemake.params.gene_no, replace=False)

    # pairs of ind and ct
    ind_ct = obs.loc[(~obs['ind'].isna()) & (~obs['ct'].isna()), ['ind', 'ct']].drop_duplicates()

    # bootstrap
    var_ctnu = {'ind':[], 'ct':[]}
    for index, row in ind_ct.iterrows():
        ind, ct = row['ind'], row['ct']
        var_ctnu['ind'].append( ind )
        var_ctnu['ct'].append( ct )
        data = ann[(obs['ind']==ind) & (obs['ct']==ct), genes].X 

        var_ctnu = stats.bootstrap((data,), cal_ctnu).standard_error**2

