from typing import Tuple, Optional, Union

import re, sys
import numpy as np, pandas as pd
from scipy import stats, sparse
from anndata import AnnData
import scanpy as sc

from . import log


def normalize(X: Union[sparse.csr_matrix, np.ndarray], l=1e6) -> Union[sparse.csr_matrix, np.ndarray]:
    """
    Normalize read counts to CPM: x/s, where s = total reads / 10^6

    Parameters:
        X:  ann.X (raw read counts), sparse or dense
        l:  normalized total counts

    Returns:
        Normalized X (sparse input is returned as CSR, dense input as ndarray)
    """
    if sparse.issparse(X):
        s = X.sum(axis=1).A1
        X = X.multiply(l / s[:, np.newaxis]).tocsr()
    else:
        s = X.sum(axis=1)
        X = X * (l / s[:, np.newaxis])

    return X


def transform(X: Optional[Union[sparse.csr_matrix, np.ndarray]] = None,
              transform: str='raw_counts',
              ann: Optional[AnnData] = None, ) -> Union[sparse.csr_matrix, np.ndarray]:
    """
    Transform read counts
    Check more transformation methods in sc.pp and sc.experimental.pp,
    and https://www.sc-best-practices.org/preprocessing_visualization/normalization.html

    Parameters:
        X:  ann.X (raw read counts), in the shape of cell x gene, sparse or dense
        transform:  transformation methods
                    'raw_counts': no transformation
                    'logp_cpm': ln(CPM + 1)
                    'logp_cp10k': ln(cp10k + 1)

    Returns:
        Transformed X
    """
    if X is not None:
        if transform == 'raw_counts':
            X = X.copy()
        elif transform == 'logp_cpm':
            l = 1e6
            X = normalize(X, l)
            X = np.log1p(X)
        elif transform == 'logp_cp10k':
            l = 1e4
            X = normalize(X, l)
            X = np.log1p(X)
        else:
            raise ValueError(f"Unknown transform '{transform}'! Supported: 'raw_counts', 'logp_cpm', 'logp_cp10k'.")

    elif ann is not None:
        if transform == 'raw_counts':
            X = ann.X.copy()
        elif transform == 'logp_cpm' or transform == 'logp_cp10k':
            l = 1e6 if transform == 'logp_cpm' else 1e4
            scales_counts = sc.pp.normalize_total(ann, target_sum=l, inplace=False)
            X = sparse.csr_matrix(sc.pp.log1p(scales_counts["X"]))
        else:
            raise ValueError(f"Unknown transform '{transform}'! Supported: 'raw_counts', 'logp_cpm', 'logp_cp10k'.")

    else:
        raise ValueError('Either X or ann must be provided!')

    return X


def calculate_expected_cell_to_cell_variance(X: sparse.csr_matrix, 
                            indicator: pd.DataFrame) -> np.ndarray:
    """
    Calculate expected cell to cell variance for each pair of individual 
    and cell type

    Parameters:
        X:  ann.X
        indicators:  a dataframe with columns of 'ind+ct', 
                    indicating the assignment of cells to cell types

    Returns:
        A dataframe with columns of 'ind+ct',
        of cell-to-cell variances for each pair of individual and cell type
    """
    # calculate cell numbers
    cell_num = indicator.to_numpy().sum(axis=0)[:, np.newaxis]

    # convert to sparse matrix
    indicator = sparse.csr_matrix(indicator.to_numpy())

    # calculate ctp
    indicator_inv = indicator.multiply(1 / cell_num.T)
    ctp = indicator_inv.T @ X

    # calculate cell-to-cell variance
    ctp2 = indicator_inv.T @ X.power(2)
    cell_var = ctp2 - ctp.power(2)

    return cell_var.toarray()


def sample_cells(x: pd.Series, rng: np.random.Generator, n: Optional[int]=None, 
                 prop: Optional[float]=None) -> pd.Series:
    """
    Sample n cells or proportion of cells for the pair of 
    individual and cell type

    Parameters:
        x:  a column of cell identity
        n:  number of cells to sample
        prop: proportion of cells to sample

    Returns:
        A series of sampled cells
    """
    x = x.copy()
    n_real_cells = int(x.sum())
    if n is None and prop is not None:
        n = int(prop * n_real_cells)
    sampled_cells = rng.choice(n_real_cells, n, replace=True)
    cells, counts = np.unique(sampled_cells, return_counts=True)
    all_counts = np.zeros(n_real_cells)
    all_counts[cells] = counts
    # sampled_cells = rng.integers(0, n + 1, n_real_cells - 1)
    # sampled_cells = np.append(sampled_cells, [0, n])
    # sampled_cells.sort()
    # sampled_cells = np.diff(sampled_cells)
    # rng.shuffle(sampled_cells)
    x.loc[x == 1] = all_counts

    return x


def sample_reads(X: sparse.csr_matrix, n_reads: Optional[int]=None, p_reads: Optional[float]=None,
                 rng: np.random.Generator=np.random.default_rng()) -> sparse.csr_matrix:
    """
    Sampling reads per cell: Row-wise multinomial sampling from a CSR sparse matrix.

    Parameters
    ----------
    X : csr_matrix
        Sparse matrix, each row is a cell containing read counts across genes.
    n_reads : int or None
        If int, fixed number of reads per cell.
    p_reads : float or None
        If float, proportion of reads to sample per cell.
    rng : np.random.Generator

    Returns
    -------
    Y : csr_matrix
        Multinomial resampled counts in CSR format.
    """

    data_out, indices_out, indptr_out = [], [], [0]

    for i in range(X.shape[0]):
        start, end = X.indptr[i], X.indptr[i+1]
        counts = X.data[start:end]
        cols  = X.indices[start:end]

        if counts.sum() == 0:
            # empty row
            indptr_out.append(len(data_out))
            continue

        # normalize to probabilities
        probs = counts / counts.sum()

        # total reads to sample
        n = 0
        if n_reads is not None:
            n = n_reads
        elif p_reads is not None:
            n = int(p_reads * counts.sum())

        draw = rng.multinomial(n, probs)

        # keep only nonzeros
        mask = draw > 0
        data_out.extend(draw[mask])
        indices_out.extend(cols[mask])
        indptr_out.append(len(data_out))

    Y = sparse.csr_matrix(
        (np.array(data_out, dtype=int),
        np.array(indices_out, dtype=int),
        np.array(indptr_out, dtype=int)),
        shape=X.shape
    )

    return Y


def pseudobulk(counts: Optional[pd.DataFrame] = None, 
               meta: Optional[pd.DataFrame] = None, 
               ann: Optional[AnnData] = None, 
               X: Optional[sparse.csr_matrix] = None, 
               obs: Optional[pd.DataFrame] = None, 
               var: Optional[pd.DataFrame] = None,
               ind_col: str = 'ind', ct_col: str = 'ct', cell_col: str = 'cell', 
               ind_cut: int = 0, ct_cut: int = 1, 
               sim_cell_count: Optional[int] = None, 
               sim_cell_prop: Optional[float] = None, 
               output_expected_ctnu: Optional[bool] = False,
               seed: Optional[int] = None,
               ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Compute Cell Type-Specific Pseudobulk and Cell Type-specific noise variance
    and remove individuals with # of cells <= ind_cut
    and for paris of individual and cell type with # of cell <= ct_cut, set CTP and CT noise variance to missing

    Parameters:
        counts: contains processed gene expression level across all cells (from all cell types and individuals).
                Each row corresponds to a single genes and each column corresponds to a single cell.
                Use gene names as dataframe row INDEX and cell IDs as dataframe COLUMNS label.
        meta:   contains cell meta informations.
                It contains at least three columns: cell, ind, ct.
                cell contains cell IDs, corresponding to column labels in counts.
                ind contains individual IDs.
                ct contains cell type, indicating the assignment of cells to cell types. 
        ann:    AnnData object (e.g. data from scanpy)
                each row corresponds to a cell with cell id, 
                and each column corresponds to a gene with a gene id.
                each cell have additional metadata, including cell type (column name: ct) 
                and donor information (column name: ind) for each cell
        X:  ann.X (normalized)
        obs:    ann.obs, cell info, with cell id as index
        var:    ann.var, gene info, with gene id as index
        ind_col:    column name for individual in meta or ann; it's renamed as ind in the output 
        ct_col: column name for cell type in meta or ann; it's renamed as ct in the output
        cell_col:   column name for cell in meta 
        ind_cut:    exclude individuals with # of cells <= ind_cut
        ct_cut: set pairs of individual and cell type with # of cells <= ct_cut to missing
        sim_cell_count: number of cells to simulate for each (ind, ct) pair, 
                        if None, no simulation is performed; only used in ann
        sim_cell_prop: proportion of cells to simulate for each (ind, ct) pair, 
                        if None, no simulation is performed; only used in ann
        output_expected_ctnu: if True, compute expected CT noise variance instead of observed CT noise variance
        seed: random seed for simulation

    Returns:
        A tuple of 
            #. Cell Type-specific Pseudobulk of index: (ind, ct) * columns: genes
            #. Cell Type-specific noise variance of idnex: (ind, ct) * columns: genes
            #. cell type proportion matrix
            #. number of cells in each (ind, ct) before filtering
    """

    if ann is None and X is None:
        # sanity check
        ## missing values
        if np.any(pd.isna(counts)):
            sys.exit('Missing values in gene expression!\n')

        # change column name
        meta = meta.rename(columns={ind_col: 'ind', ct_col: 'ct', cell_col: 'cell'})

        ## identical unique cell ids
        if (len(np.unique(meta['cell'])) != meta.shape[0]):
            sys.exit('Duplicate cells!\n')
        if (meta.shape[0] != counts.shape[1]) or (len(np.setdiff1d(meta['cell'], counts.columns)) != 0):
            sys.exit('Cells not matching in counts and meta!\n')

        # collect genes and cell types
        genes = counts.index.tolist()
        # cts = np.unique(meta['ct'].to_numpy())

        # transfrom gene * cell to cell * gene in counts
        counts = counts.transpose()

        # merge with meta
        data = counts.merge(meta, left_index=True, right_on='cell')  # cell * (gene, ind, ct)

        # group by ind and ct
        data_grouped = data.groupby(['ind', 'ct'])
        # exclude groups with only one cell
        b = len(data_grouped)
        data_grouped = data_grouped.filter(lambda x: len(x) > 1).groupby(['ind', 'ct'])
        a = len(data_grouped)
        if a != b:
            print(f'Exclude {b - a} individual-cell type pairs with only one cell')

        # compute ctp
        ctp = data_grouped[genes].aggregate(np.mean)

        # compute ctnu
        ctnu = data_grouped[genes].aggregate(stats.sem)  # need to rm ind-ct pair with only 1 cell
        ctnu = ctnu ** 2

    else:
        if ann is not None:
            obs = ann.obs
            X = ann.X if sparse.issparse(ann.X) else ann.X[:, :]  # convert sparse dataset to sparse matrix
            var = ann.var
        obs = obs.rename(columns={ind_col: 'ind', ct_col: 'ct'})

        # paires of ind and ct
        print('Finding ind-ct pairs')
        sep = '_+_'
        while np.any(obs['ind'].str.contains(sep)) or np.any(obs['ct'].str.contains(sep)):
            sep = '_' + sep + '_'
        ind_ct = obs['ind'].astype(str) + sep + obs['ct'].astype(str)

        if sim_cell_count is not None and sim_cell_prop is not None:
            raise ValueError('Both sim_cell_count and sim_cell_prop are provided, only one of them is allowed')
        simulating = sim_cell_count is not None or sim_cell_prop is not None

        if simulating or output_expected_ctnu:
            # dense one-hot needed here: sample_cells() operates per (ind,ct) column as a
            # pandas Series, and calculate_expected_cell_to_cell_variance() expects a
            # DataFrame -- for the common case below (neither simulating nor computing
            # expected ctnu), the indicator is built directly as sparse instead, since a
            # dense (n_cells x n_ind_ct_pairs) one-hot matrix is >99% zeros for typical
            # cell/individual/cell-type counts.
            indicator = pd.get_dummies(ind_ct, dtype=float)  # don't use dtype int8

            # calculate cell-to-cell variance
            if output_expected_ctnu:
                cell_var = calculate_expected_cell_to_cell_variance(X, indicator)

            if simulating:
                log.logger.info('Simulating cells')
                rng = np.random.default_rng(seed)
                if sim_cell_count is not None:
                    # simulate cells for each (ind, ct) pair
                    indicator = indicator.apply(lambda x: sample_cells(x, rng,
                                                                n=sim_cell_count),
                                                axis=0)
                    assert np.all(indicator.sum(axis=0) == sim_cell_count)
                else:
                    # simulate cells for each (ind, ct) pair
                    indicator = indicator.apply(lambda x: sample_cells(x, rng,
                                                                prop=sim_cell_prop),
                                                axis=0)

            cell_num = indicator.to_numpy().sum(axis=0)[:, np.newaxis]
            ind_ct = indicator.columns
            indicator = sparse.csr_matrix(indicator.to_numpy())
        else:
            # build the indicator directly as sparse -- avoids ever materializing an
            # (n_cells x n_ind_ct_pairs) dense one-hot matrix
            codes, categories = pd.factorize(ind_ct, sort=True)
            n_cells, n_pairs = len(codes), len(categories)
            indicator = sparse.csr_matrix(
                (np.ones(n_cells), (np.arange(n_cells), codes)), shape=(n_cells, n_pairs))
            cell_num = np.asarray(indicator.sum(axis=0)).reshape(-1, 1)
            ind_ct = categories

        # inverse indicator matrix
        indicator_inv = indicator.multiply(1 / cell_num.T)
        # print(indicator_inv.getrow(0)[0, 0])

        # compute ctp
        log.logger.info('Computing CTP')
        ctp = indicator_inv.T @ X

        # compute ctnu
        log.logger.info('Computing CTNU')
        adj_cell_num = cell_num - 1
        adj_cell_num[adj_cell_num == 0] = 1e-12  # avoid divide by zero  # array from np matrix doesn't work (bug?)
        # print(adj_cell_num)
        # print(np.any(adj_cell_num == 0))
        # np.savetxt('adj_cell_num.txt', adj_cell_num)

        if sparse.issparse(X): 
            ctp2 = indicator_inv.T @ X.power(2)
            ctnu = (ctp2 - ctp.power(2)).multiply(1 / adj_cell_num)

            ctp = ctp.toarray()
            ctnu = ctnu.toarray()
            # print(ctp[:5, :5], ctp2.toarray()[:5, :5], ctnu[:5, :5], adj_cell_num[:5, 0])

        else:
            ctp2 = indicator_inv.T @ (X ** 2)
            ctnu = (ctp2 - ctp ** 2) * (1 / adj_cell_num)
            # print(ctp[:5, :5], ctp2[:5, :5], ctnu[:5, :5], adj_cell_num[0, :5])

        if output_expected_ctnu and sim_cell_count is not None:
            log.logger.info('Outputting expected CTNU')
            ctnu = cell_var / sim_cell_count

        # convert to df
        ctp_index = pd.MultiIndex.from_frame(ind_ct.to_series().str.split(sep, expand=True, regex=False),
                                            names=['ind', 'ct'])
        ctp = pd.DataFrame(ctp, index=ctp_index, columns=var.index)
        ctnu = pd.DataFrame(ctnu, index=ctp_index, columns=var.index)

        # group by ind-ct to compute cell numbers
        data_grouped = obs.reset_index(drop=False, names='cell').groupby(['ind', 'ct'])


    # compute cell numbers
    cell_counts = data_grouped['cell'].count().unstack(fill_value=0).astype('int')
    # sort  
    cell_counts.index = cell_counts.index.astype(str)
    cell_counts = cell_counts.sort_index(axis=0)
    if sim_cell_count is not None:
        # update cell numbers
        cell_counts = cell_counts / cell_counts
        cell_counts = cell_counts * sim_cell_count
        cell_counts = cell_counts.fillna(0)
        cell_counts = cell_counts.astype('int')

    # filter individuals
    cell_counts = cell_counts[cell_counts.sum(axis=1) >= ind_cut]

    # convert cell counts to proportions
    P = cell_counts.divide(cell_counts.sum(axis=1), axis=0)

    # filter cts
    cell_counts2 = cell_counts.stack()
    cell_counts2 = cell_counts2.loc[cell_counts2 >= ct_cut]
    ctp = ctp.loc[ctp.index.isin(cell_counts2.index)]
    ctnu = ctnu.loc[ctnu.index.isin(cell_counts2.index)]

    # print(ctp.shape, ctnu.shape, P.shape)

    return (ctp, ctnu, P, cell_counts)


def std(ctp: pd.DataFrame, ctnu: pd.DataFrame, P: pd.DataFrame, return_all: bool = False
        ) -> Union[Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame], 
                   Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]]:
    """
    For each Gene, stardardize Overall Pseudobulk (OP) to mean 0 and std 1, and scale ctp,
    overall noise variance (nu), ct-specific noise variance (ctnu) correspondingly.

    Parameters:
        ctp:    imputed ctp with shape index: (ind, ct) * columns: genes
        ctnu:   imputed ctnu with shape index: (ind, ct) * columns: genes
        P:  cell type proportions matrix of shape ind * ct
        return_all: return ctnu for OP and nu for CTP
    Returns:
        A tuple of
            #. op
            #. nu
            #. ctnu_op (optional)
            #. ctp
            #. nu_ctp (optional)
            #. ctnu
    """
    genes = ctp.columns.tolist()

    # sanity check that all genes are expressed in all cell types (vectorized across all genes)
    all_zero_by_ct = (ctp == 0).groupby(level='ct').transform('all')
    if all_zero_by_ct.any().any():
        bad_genes = all_zero_by_ct.columns[all_zero_by_ct.any(axis=0)]
        for gene in bad_genes:
            bad_cts = all_zero_by_ct.index.get_level_values('ct')[all_zero_by_ct[gene]].unique().tolist()
            raise ValueError(f'Gene {gene} is not expressed in cell types {bad_cts}!\nCIGMA assumes all genes are expressed in all cell types, please check your data and consider removing genes that are not expressed in some cell types.')

    # sanity check order of inds and cts in ctp, ctnu, and P (checked once, not per gene)
    inds = P.index
    cts = P.columns
    gene0_ctp = ctp[genes[0]].unstack()
    gene0_ctnu = ctnu[genes[0]].unstack()
    if not (gene0_ctnu.index.equals(inds) and gene0_ctp.index.equals(inds)):
        sys.exit('Individuals not matching!')
    if not (gene0_ctnu.columns.equals(cts) and gene0_ctp.columns.equals(cts)):
        sys.exit('Cell types not matching!')

    # broadcast P (ind * ct) against ctp/ctnu (index: (ind, ct) * columns: genes)
    P_flat = P.stack()
    P_flat.index = P_flat.index.set_names(ctp.index.names)

    # compute op and nu (all genes at once)
    op = ctp.mul(P_flat, axis=0).groupby(level='ind').sum()
    neg = ctnu < 0
    ctnu_op = ctnu.where(~neg, 0)  # set negative ctnu to 0 for OP data
    nu_op = ctnu_op.mul(P_flat ** 2, axis=0).groupby(level='ind').sum()

    # set negative ctnu to max (per gene, per cell type, over individuals) for CTP data
    max_by_ct = ctnu.groupby(level='ct').transform('max')
    ctnu_ctp = ctnu.where(~neg, max_by_ct)
    nu_ctp = ctnu_ctp.mul(P_flat ** 2, axis=0).groupby(level='ind').sum()

    # standardize op (per gene mean/std/var, computed for all genes at once)
    mean, sd, var = op.mean(axis=0), op.std(axis=0), op.var(axis=0)

    std_op = op.sub(mean, axis=1).div(sd, axis=1)
    std_ctp = ctp.sub(mean, axis=1).div(sd, axis=1)
    nu_op = nu_op.div(var, axis=1)
    nu_ctp = nu_ctp.div(var, axis=1)
    ctnu_op = ctnu_op.div(var, axis=1)
    ctnu_ctp = ctnu_ctp.div(var, axis=1)

    # santity check order of inds and cts doesn't change after standardization
    # assert np.all(std_op.index.equals(inds) and std_op.columns.equals(cts)
    if not (std_ctp.index.equals(ctp.index) and std_ctp.columns.equals(ctp.columns)):
        if not std_ctp.index.equals(ctp.index):
            index1 = std_ctp.index
            index2 = ctp.index
            print(std_ctp.index)
            print(ctp.index)
            print(index1.difference(index2))  
            print(index2.difference(index1))
            print(index1.dtype == index2.dtype)
            print(index1.name == index2.name)
            sys.exit('Orders of individuals changed after standardization!')
        elif not std_ctp.columns.equals(ctp.columns):
            # print(std_ctp.columns)
            # print(ctp.columns)
            # print(std_ctp.columns[std_ctp.columns != ctp.columns])
            # print(ctp.columns[std_ctp.columns != ctp.columns])
            sys.exit('Orders of cell types changed after standardization!')

    if return_all:
        return std_op, nu_op, ctnu_op, std_ctp, nu_ctp, ctnu_ctp
    else:
        return std_op, nu_op, std_ctp, ctnu_ctp
