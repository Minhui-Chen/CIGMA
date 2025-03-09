import os, math, re, sys
from scipy import stats
from scipy import linalg
import numpy as np
import pandas as pd


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


def main():
    batch = snakemake.params.batches[int(snakemake.wildcards.i)]

    # par
    beta = np.loadtxt(snakemake.input.beta)
    V = np.loadtxt(snakemake.input.V)
    W = np.loadtxt(snakemake.input.W)
    C = len(beta)

    sig_g = float(snakemake.wildcards.vc.split('_')[1])
    sig_e = float(snakemake.wildcards.vc.split('_')[2])
    mean_nu = float(snakemake.wildcards.vc.split('_')[-1])
    var_nu = (float(snakemake.wildcards.std_nu_scale) * mean_nu) ** 2
    a = np.array([float(x) for x in snakemake.wildcards.a.split('_')])
    ss = int(float(snakemake.wildcards.ss))


    data = {}
    for i in batch:
        rng = np.random.default_rng(snakemake.params.seed + i)
        data[i] = {}

        # simulate genotypes
        ## draw allele frequency from beta distribution
        frq = rng.beta(snakemake.params.beta[0], snakemake.params.beta[1], snakemake.params.L * 5)
        frq = frq[(frq > snakemake.params.maf) & (frq < (1 - snakemake.params.maf))][:snakemake.params.L]
        G = []
        for frq_ in frq:
            ## draw genotype from binomial distribution based on allele frequency
            G_ = rng.binomial(2, frq_, ss)
            while len(np.unique(G_)) == 1:
                G_ = rng.binomial(2, frq_, ss)
            G.append(G_)
        ## convert SNP x IND to IND x SNP
        G = np.array(G).T
        ## standardize
        # if np.any(np.std(G, axis=0) == 0):
        #    sys.exit(f'{sum(np.std(G, axis=0) == 0)}')
        G = stats.zscore(G)
        ## save
        data[i]['G'] = G

        # calculate K
        K = G @ G.T / G.shape[1]
        data[i]['K'] = K

        # simulate SNP effect
        ## additive effect
        ### draw from normal distribution of N(0, hom2/L)
        if sig_g == 0:
            add = np.zeros(snakemake.params.L)
        else:
            add = rng.normal(0, math.sqrt(sig_g / snakemake.params.L), snakemake.params.L)
            add = add - np.mean(add)
            add = add * math.sqrt(sig_g / snakemake.params.L) / np.std(add)
            if len(add) != snakemake.params.L:
                print(add)
                print(len(add))
                sys.exit('Weird')

        ## CT-specific SNP effect
        if np.all(V == np.zeros_like(V)):
            H = np.zeros((snakemake.params.L, C))
        else:
            H = rng.multivariate_normal(np.zeros(C), V / snakemake.params.L, snakemake.params.L)  # of shape SNP x cell type
            H = H - np.mean(H, axis=0)  # NOTE: covariance in Full model is not stded
            H = (H * np.sqrt(np.diag(V)/snakemake.params.L)) / np.std(H, axis=0)

        # calculate alpha, shared genetic effect
        alpha_g = G @ add

        # simulate shared noise
        alpha_e = rng.normal(0, math.sqrt(sig_e), ss)

        # simulate cell type proportions
        P = rng.dirichlet(alpha=a, size=ss)
        P = adjust_min_value(P, 0.05)
        assert np.allclose(P.sum(axis=1), np.ones(P.shape[0]))
        data[i]['P'] = P
        pi = np.mean(P, axis=0)
        data[i]['pi'] = pi

        ## estimate S
        ### demean P
        pd = P - pi
        ### covariance
        s = (pd.T @ pd) / ss
        # print(bmatrix(s))
        data[i]['s'] = s

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
            data[i]['nu'] = nu
            ctnu = np.zeros((ss, C))
            data[i]['ctnu'] = ctnu
        else:
            theta = var_nu / mean_nu
            k = mean_nu / theta
            ### variance of error for each individual
            nu = rng.gamma(k, scale=theta, size=ss)
            data[i]['nu'] = nu
            #### variance of error for each individual-cell type
            ctnu = nu.reshape(-1, 1) * (1 / P)
            data[i]['ctnu'] = ctnu

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
            data[i]['fixed'] = X
        
        if 'random' in snakemake.wildcards.keys():
            X, b = add_random(int(snakemake.wildcards.random), ss, rng)
            y = y + X @ b
            Y = Y + (X @ b)[:, np.newaxis]
            data[i]['random'] = X

        data[i]['y'] = y
        data[i]['Y'] = Y

    np.save(snakemake.output.data, data)


if __name__ == '__main__':
    main()
