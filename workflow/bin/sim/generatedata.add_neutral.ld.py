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
    L = int(snakemake.params.get("L", snakemake.wildcards.L))
    nL = int(snakemake.params.get("nL", snakemake.wildcards.nL))
    tL = L + nL

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

        # simulate genotypes with LD
        ## draw allele frequency from beta distribution
        frq = rng.beta(snakemake.params.beta[0], snakemake.params.beta[1], tL * 5)
        frq = frq[(frq > snakemake.params.maf) & (frq < (1 - snakemake.params.maf))][:tL]
        
        ## create correlation matrix with r2=0.6 between adjacent SNPs
        r = np.sqrt(float(snakemake.wildcards.ld))  # correlation from r^2 = 0.6
        Sigma = (1-r) * np.eye(tL) + r 
        # for j in range(tL):
        #     if j + 1 < tL:
        #         Sigma[j, j + 1] = Sigma[j + 1, j] = r
                
        ## generate correlated normal variables
        Z = rng.multivariate_normal(mean=np.zeros(tL), cov=Sigma, size=ss)
        
        ## convert to genotypes based on allele frequencies
        G = []
        for snp_idx in range(tL):
            z = Z[:, snp_idx]
            frq_ = frq[snp_idx]
            # Use quantiles to convert to genotypes under HWE
            # P(G=0) = (1-frq)^2, P(G=1) = 2*frq*(1-frq), P(G=2) = frq^2
            threshold1 = stats.norm.ppf((1 - frq_)**2)
            threshold2 = stats.norm.ppf((1 - frq_)**2 + 2 * frq_ * (1 - frq_))
            
            G_ = np.zeros(ss, dtype=int)
            G_[z > threshold1] = 1
            G_[z > threshold2] = 2
            
            # Ensure variation exists
            if len(np.unique(G_)) == 1:
                # Fallback to binomial if no variation
                G_ = rng.binomial(2, frq_, ss)
                while len(np.unique(G_)) == 1:
                    G_ = rng.binomial(2, frq_, ss)
            
            G.append(G_)
            
        ## convert SNP x IND to IND x SNP
        G = np.array(G).T
        data[i]['rawG'] = G
        ## standardize
        # if np.any(np.std(G, axis=0) == 0):
        #    sys.exit(f'{sum(np.std(G, axis=0) == 0)}')
        # maf = np.mean(G, axis=0) / 2
        # G = (G - 2 * maf) / np.sqrt(2 * maf * (1 - maf))
        G = (G - G.mean(axis=0)) / G.std(axis=0)

        ## save
        data[i]['G'] = G

        # calculate K
        # K = G @ G.T / G.shape[1]
        # data[i]['K'] = K

        # subset causal variants
        G = G[:, rng.choice(G.shape[1], L, replace=False)]
        K = G @ G.T / G.shape[1]
        data[i]['K'] = K

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
