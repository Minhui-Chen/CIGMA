import os, math, re, sys
import numpy as np
import pandas as pd

from scipy import stats
from scipy import linalg



def add_fixed(levels, ss, rng):
    ''' Add a test fixed effect'''
    a = rng.integers(levels, size=ss)
    while len(np.unique(a)) < levels:
        a = rng.integers(levels, size=ss)
    X = pd.get_dummies(a, dtype='int8')

    b = X.columns.to_numpy() / np.std(a) # to shrink
    # centralize b
    b = b - np.mean(X @ b)

    return X.to_numpy(), b


def add_random(levels, ss, rng):
    ''' Add a test random effect'''
    b = rng.normal(0, 1, levels)
    a = rng.choice(b, ss)
    while len(np.unique(a)) < levels:
        a = rng.choice(b, ss)

    X = pd.get_dummies(a, dtype='int8')
    b = X.columns.to_numpy()

    return X.to_numpy(), b


def adjust_min_value(arr, cut=0.01):
    for i in range(arr.shape[0]):
        row = arr[i]
        deficit = 0
        
        # Identify elements less than 0.01
        mask = row < cut
        deficit = np.sum(cut - row[mask])
        
        # Set the minimum value to 0.01 for elements that are less than 0.01
        row[mask] = cut
        
        # Redistribute the deficit proportionally to elements greater than 0.01
        row[~mask] -= deficit * (row[~mask] / np.sum(row[~mask]))
    
    return arr


def prune_correlated_snps(G_std, n_keep, r2_threshold=0.8, rng=None):
    """
    Select SNPs that are relatively independent (r² < threshold).
    
    Parameters:
    -----------
    G_std : ndarray
        Standardized genotype matrix (samples x SNPs)
    n_keep : int
        Number of SNPs to keep
    r2_threshold : float
        Maximum r² allowed between selected SNPs
    rng : Generator
        Random number generator
        
    Returns:
    --------
    selected_indices : ndarray
        Indices of selected SNPs
    """
    n_snps = G_std.shape[1]
    
    # Shuffle SNP order for random selection
    if rng is not None:
        snp_order = rng.permutation(n_snps)
    else:
        snp_order = np.arange(n_snps)
    
    selected = []
    
    for idx in snp_order:
        if len(selected) >= n_keep:
            break
            
        # Check correlation with already selected SNPs
        if len(selected) == 0:
            selected.append(idx)
        else:
            # Compute r² with all selected SNPs
            correlations = np.abs(np.corrcoef(G_std[:, idx], G_std[:, selected].T)[0, 1:])
            r2 = correlations ** 2
            
            # Add if not highly correlated with any selected SNP
            if np.all(r2 < r2_threshold):
                selected.append(idx)
    
    # If we don't have enough SNPs, lower threshold iteratively
    # if len(selected) < n_keep:
    #     print(f"Warning: Could not find {n_keep} SNPs with r2<{r2_threshold}. Found {len(selected)}. Relaxing threshold...")
    #     return prune_correlated_snps(G_std, n_keep, r2_threshold=min(0.99, r2_threshold + 0.1), rng=rng)
    
    return np.array(selected)


def main():
    # par
    beta = np.loadtxt(snakemake.input.beta)
    V = np.loadtxt(snakemake.input.V)
    W = np.loadtxt(snakemake.input.W)
    C = len(beta)
    L = int(snakemake.params.get("L", snakemake.wildcards.L))
    sig_g = float(snakemake.wildcards.vc.split('_')[1])
    sig_e = float(snakemake.wildcards.vc.split('_')[2])
    mean_nu = float(snakemake.wildcards.vc.split('_')[-1])
    var_nu = (float(snakemake.wildcards.std_nu_scale) * mean_nu) ** 2
    a = np.array([float(x) for x in snakemake.wildcards.a.split('_')])
    ss = int(float(snakemake.wildcards.ss))

    # select genes with enough SNPs
    genes = []
    gene_fs = []
    for gene_f in snakemake.input.genes:
        gene_data = np.load(gene_f, allow_pickle=True).item()
        # exclude genes with variants less than 20
        new_genes = [k for k, v in gene_data.items() if v['nsnps'] >= L]
        genes = genes + new_genes
        gene_fs = gene_fs + [gene_f] * len(new_genes)
    gene_fs = {gene: gene_f for gene, gene_f in zip(genes, gene_fs)}
    print(len(genes))
    # check duplicate gene
    assert len(genes) == len(set(genes))
    # Select evenly spaced genes # TODO: tmp
    n_select = 1000
    indices = np.linspace(0, len(genes)-1, n_select, dtype=int)
    genes = list(np.array(genes)[indices])


    # simulation
    nbatch = len(snakemake.output.data)
    batches = np.array_split(genes, nbatch)
    old_gene_f = ''
    for k, batch in enumerate(batches):
        out_f = snakemake.output.data[k]
        print(f"Processing batch {k+1}/{nbatch}, saving to {out_f}")
        data = {}
        for gene_name in batch:
            gene_f = gene_fs[gene_name]
            if gene_f != old_gene_f:
                print(f"Loading genotype data from {gene_f} for gene {gene_name}")
                old_gene_f = gene_f
                gene_data = np.load(gene_f, allow_pickle=True).item()
            else:
                print(f"Reusing loaded genotype data for gene {gene_name}")
            data[gene_name] = {}
        
            i = genes.index(gene_name)
            rng = np.random.default_rng(snakemake.params.seed + i + ss + L)

            # simulate genotypes
            G = gene_data[gene_name]['G']
            assert G.shape[0] > ss, f"{G.shape[0]} samples in genotype but {ss} needed"
            # random sample individuals
            perm = rng.permutation(G.shape[0])[:ss]
            G = G[perm, :]
            # remove monomorphic SNPs
            stds = np.std(G, axis=0)
            G = G[:, stds > 0]
            assert G.shape[1] >= L, f"{G.shape[1]} SNPs in genotype but {L} needed"

            # data[gene_name]['rawG'] = G

            ## standardize # TODO
            # if np.any(np.std(G, axis=0) == 0):
            #    sys.exit(f'{sum(np.std(G, axis=0) == 0)}')
            # maf = np.mean(G, axis=0) / 2
            # G = (G - 2 * maf) / np.sqrt(2 * maf * (1 - maf))
            G = (G - G.mean(axis=0)) / G.std(axis=0)

            ## save
            data[gene_name]['G'] = G

            # calculate K
            # K = G @ G.T / G.shape[1]
            # data[gene_name]['K'] = K

            # subset causal variants # TODO
            snps = rng.permutation(G.shape[1])[:L]
            G = G[:, snps]
            # TODO: remove correalted SNPs
            # G = G[:, prune_correlated_snps(G[:, rng.permutation(G.shape[1])], L, r2_threshold=0.8, rng=rng)]
            # if G.shape[1] < L:
            #     data.pop(gene_name)
            #     continue
            # else:
            #     G = G[:, :L]

            # TODO: tmp use causal variants to calculate K
            K = G @ G.T / G.shape[1]
            data[gene_name]['K'] = K

            # simulate SNP effect
            ## additive effect
            ### draw from normal distribution of N(0, hom2/L)
            if sig_g == 0:
                add = np.zeros(L)
            else:
                add = rng.normal(0, math.sqrt(sig_g / L), L)
                add = add - np.mean(add)
                add = add * math.sqrt(sig_g / L) / np.std(add)
                if len(add) != L:
                    print(add)
                    print(len(add))
                    sys.exit('Weird')

            ## CT-specific SNP effect
            if np.all(V == np.zeros_like(V)):
                H = np.zeros((L, C))
            else:
                H = rng.multivariate_normal(np.zeros(C), V / L, L)  # of shape SNP x cell type
                H = H - np.mean(H, axis=0)  # NOTE: covariance in Full model is not stded
                H = (H * np.sqrt(np.diag(V)/L)) / np.std(H, axis=0)

            # calculate alpha, shared genetic effect
            alpha_g = G @ add

            # simulate shared noise
            alpha_e = rng.normal(0, math.sqrt(sig_e), ss)

            # simulate cell type proportions
            P = rng.dirichlet(alpha=a, size=ss)
            P = adjust_min_value(P, 0.05)
            assert np.allclose(P.sum(axis=1), np.ones(P.shape[0]))
            data[gene_name]['P'] = P
            pi = np.mean(P, axis=0)
            data[gene_name]['pi'] = pi

            ## estimate S
            ### demean P
            pd = P - pi
            ### covariance
            s = (pd.T @ pd) / ss
            # print(bmatrix(s))
            data[gene_name]['s'] = s

            # calculate ct fixed effect
            ct_main = P @ beta

            # calculate ct-specific genetic
            ct_g = linalg.khatri_rao(G.T, P.T).T @ H.flatten()

            # simulate cell type-specific noise
            gamma_e = rng.multivariate_normal(np.zeros(C), W, ss)
            # calculate cell type-specific noise for OP
            ct_e = linalg.khatri_rao(np.eye(ss), P.T).T @ gamma_e.flatten()

            # draw residual error
            ## draw variance of residual error for each individual from gamma distribution \Gamma(k, theta)
            ## with mean = k * theta, var = k * theta^2, so theta = var / mean, k = mean / theta
            ## since mean = 0.2 and assume var = 0.01, we can get k and theta
            if mean_nu == 0:
                nu = np.zeros(ss)
                data[gene_name]['nu'] = nu
                ctnu = np.zeros((ss, C))
                data[gene_name]['ctnu'] = ctnu
            else:
                theta = var_nu / mean_nu
                k = mean_nu / theta
                ### variance of error for each individual
                nu = rng.gamma(k, scale=theta, size=ss)
                data[gene_name]['nu'] = nu
                #### variance of error for each individual-cell type
                ctnu = nu.reshape(-1, 1) * (1 / P)
                data[gene_name]['ctnu'] = ctnu

            ## draw residual error from normal distribution with variance drawn above
            error = rng.normal(loc=0, scale=np.sqrt(nu))
            ct_error = rng.normal(loc=0, scale=np.sqrt(ctnu))

            # generate overall pseudobulk
            y = ct_main + alpha_g + alpha_e + ct_g + ct_e + error
            Y = np.outer(np.ones(ss), beta) + np.outer(alpha_g + alpha_e, np.ones(C)) + G @ H + gamma_e + ct_error

            # add Extra fixed and random effect
            if 'fixed' in snakemake.wildcards.keys():
                X, b = add_fixed(int(snakemake.wildcards.fixed), ss, rng)
                y = y + X @ b
                Y = Y + (X @ b)[:, np.newaxis]
                data[gene_name]['fixed'] = X
            
            if 'random' in snakemake.wildcards.keys():
                X, b = add_random(int(snakemake.wildcards.random), ss, rng)
                y = y + X @ b
                Y = Y + (X @ b)[:, np.newaxis]
                data[gene_name]['random'] = X

            data[gene_name]['y'] = y
            data[gene_name]['Y'] = Y

        np.save(out_f, data)


if __name__ == '__main__':
    main()
