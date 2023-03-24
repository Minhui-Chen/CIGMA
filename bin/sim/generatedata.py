import os, math, re
from scipy import stats
from scipy import linalg
import numpy as np

def main():

    batch = snakemake.params.batches[int(snakemake.wildcards.i)]
    G_fs = [snakemake.params.G.replace('repX','rep'+str(x)) for x in batch]
    K_fs = [snakemake.params.K.replace('repX','rep'+str(x)) for x in batch]
    P_fs = [snakemake.params.P.replace('repX','rep'+str(x)) for x in batch]
    pi_fs = [snakemake.params.pi.replace('repX','rep'+str(x)) for x in batch]
    s_fs = [snakemake.params.s.replace('repX','rep'+str(x)) for x in batch]
    nu_fs = [snakemake.params.nu.replace('repX','rep'+str(x)) for x in batch]
    ctnu_fs = [snakemake.params.ctnu.replace('repX','rep'+str(x)) for x in batch]
    y_fs = [snakemake.params.y.replace('repX','rep'+str(x)) for x in batch]
    Y_fs = [snakemake.params.Y.replace('repX','rep'+str(x)) for x in batch]

    with open(snakemake.output.G, 'w') as f: f.write('\n'.join(G_fs))
    with open(snakemake.output.K, 'w') as f: f.write('\n'.join(K_fs))
    with open(snakemake.output.P, 'w') as f: f.write('\n'.join(P_fs))
    with open(snakemake.output.pi, 'w') as f: f.write('\n'.join(pi_fs))
    with open(snakemake.output.s, 'w') as f: f.write('\n'.join(s_fs))
    with open(snakemake.output.nu, 'w') as f: f.write('\n'.join(nu_fs))
    with open(snakemake.output.ctnu, 'w') as f: f.write('\n'.join(ctnu_fs))
    with open(snakemake.output.y, 'w') as f: f.write('\n'.join(y_fs))
    with open(snakemake.output.Y, 'w') as f: f.write('\n'.join(Y_fs))

    # par
    beta = np.loadtxt(snakemake.input.beta)
    V = np.loadtxt(snakemake.input.V)
    W = np.loadtxt(snakemake.input.W)
    C = len( beta )

    sig_g = float( snakemake.wildcards.vc.split('_')[1] )
    sig_e = float( snakemake.wildcards.vc.split('_')[2] )
    mean_nu = float( snakemake.wildcards.vc.split('_')[-1] )
    var_nu = float( snakemake.wildcards.var_nu )
    a = np.array( [ float(x) for x in snakemake.wildcards.a.split('_') ] )
    ss = int( float( snakemake.wildcards.ss ) )

    rng = np.random.default_rng()

    for G_f, K_f, P_f, pi_f, s_f, nu_f, ctnu_f, y_f, Y_f in zip(G_fs,
            K_fs, P_fs, pi_fs, s_fs, nu_fs, ctnu_fs, y_fs, Y_fs):
        os.makedirs(os.path.dirname(P_f), exist_ok=True)
        # simulate genotypes
        ## draw allele frequency from beta distribution
        frq = rng.beta(snakemake.params.beta[0], snakemake.params.beta[1], snakemake.params.L*5)
        frq = frq[(frq > snakemake.params.maf) & (frq < (1-snakemake.params.maf))][:snakemake.params.L]
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
        #if np.any(np.std(G, axis=0) == 0):
        #    sys.exit(f'{sum(np.std(G, axis=0) == 0)}')
        G = stats.zscore(G)
        ## save
        np.savetxt(G_f, G, fmt='%.6e', delimiter='\t')

        # calculate K
        K = G @ G.T / G.shape[1]
        np.savetxt(K_f, K, fmt='%.6e', delimiter='\t')

        # simulate SNP effect
        ## additive effect
        ### draw from normal distribution of N(0, hom2/L)
        add = rng.normal(0, math.sqrt(sig_g/snakemake.params.L), snakemake.params.L)
        if len(add) != snakemake.params.L:
            print(add)
            print(len(add))
            sys.exit('Weird')
        ## interaction effect
        H = rng.multivariate_normal(np.zeros(C), V/snakemake.params.L, snakemake.params.L) # of shape SNP x cell type

        # calculate alpha, shared genetic effect
        alpha_g = G @ add

        # simulate shared noise
        alpha_e = rng.normal(0, math.sqrt(sig_e), ss)

        # simulate cell type proportions
        P = rng.dirichlet(alpha=a, size=ss)
        np.savetxt(P_f, P, fmt='%.6e', delimiter='\t')
        pi = np.mean(P, axis=0)
        np.savetxt(pi_f, pi, delimiter='\t')
        ## estimate S
        ### demean P
        pd = P-pi
        ### covariance
        s = (pd.T @ pd)/ss
        #print(bmatrix(s))
        np.savetxt(s_f, s, delimiter='\t')

        # calculate ct fixed effect
        ct_main = P @ beta

        # calculate ct-specific genetic
        ct_g = linalg.khatri_rao(G.T, P.T).T @ H.flatten()

        # simulate cell type-specific noise
        gamma_e = rng.multivariate_normal(np.zeros(C), W, ss)
        # calculate cell type-specific noise effect
        ct_e = linalg.khatri_rao(np.eye(ss), P.T).T @ gamma_e.flatten()

        # draw residual error
        ## draw variance of residual error for each individual from gamma distribution \Gamma(k, theta)
        ## with mean = k * theta, var = k * theta^2, so theta = var / mean, k = mean / theta
        ## since mean = 0.2 and assume var = 0.01, we can get k and theta
        theta = var_nu / mean_nu
        k = mean_nu / theta
        ### variance of error for each individual
        nu = rng.gamma(k, scale=theta, size=ss)
        np.savetxt(nu_f, nu, delimiter='\t')
        #### variance of error for each individual-cell type
        ctnu = nu.reshape(-1,1) * ( 1/P )
        np.savetxt(ctnu_f, ctnu, delimiter='\t')

        ## draw residual error from normal distribution with variance drawn above
        error = rng.normal(loc=0, scale=np.sqrt(nu))
        ct_error = rng.normal(loc=0, scale=np.sqrt(ctnu))

        # generate overall pseudobulk
        y = ct_main + alpha_g + alpha_e + ct_g + ct_e + error
        Y = np.outer(np.ones(ss), beta) + np.outer(alpha_g+alpha_e, np.ones(C)) + G @ H + gamma_e + ct_error

        np.savetxt(y_f, y, delimiter='\t')
        np.savetxt(Y_f, Y, delimiter='\t')

if __name__ == '__main__':
    main()
