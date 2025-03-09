import math
import numpy as np
import pandas as pd
from scipy import stats

def main():
    # par
    print(snakemake.wildcards)

    # cell type proportion
    alpha = np.array([float(x) for x in snakemake.wildcards.a.split('_')])
    C = len( alpha )
    ## calculate moments of dirichlet distribution
    a0 = np.sum(alpha)
    pi = alpha/a0
    print(pi)
    ## save mean proportions
    np.savetxt(snakemake.output.pi, pi, delimiter='\t')
    
    ## variance and covariance of cell type proportions
    #var = alpha*(a0-alpha)/(a0*a0*(a0+1))
    var = stats.dirichlet(alpha).var()
    s_fun = lambda alpha, a0, i, j: -1*alpha[i]*alpha[j]/(a0*a0*(a0+1))
    s = []
    for i in range(len(alpha)):
        s_e = []
        for j in range(len(alpha)):
            if i == j:
                s_e.append(var[i])
            else:
                s_e.append(s_fun(alpha, a0, i, j))
        s.append(s_e)
    s = np.array(s)
    # save cov
    #print(bmatrix(s))
    np.savetxt(snakemake.output.s, s, delimiter='\t')

    # beta
    ## par
    beta = np.array( [float(x) for x in snakemake.wildcards.beta.split('_')] )
    ct_fixed_vc = float(snakemake.wildcards.vc.split('_')[0]), # variance explained by ct fixed effect \beta^T S \beta

    ## calculate beta based on \beta^T S \beta = 0.25
    ## \beta^2 [1 1/2 1/4 1/8] cov [1 1/2 1/4 1/8]^T = 0.25
    ## \beta^2 = 0.25 / ([1 1/2 1/4 1/8] cov [1 1/2 1/4 1/8]^T)
    ##x = np.array([1*(params.ratio**(i)) for i in range(celltype_no)])
    scale = math.sqrt(ct_fixed_vc / (beta @ s @ beta))
    beta = beta * scale

    ## save
    np.savetxt(snakemake.output.beta, beta, delimiter='\t')

    # V
    # co-variance matrix of interaction effect (gamma) between individual and cell type across cell types
    ct_genetic_vc = float(snakemake.wildcards.vc.split('_')[3]) #variance explained by CT-specific genetic effect
    if ct_genetic_vc == 0:
        V = np.zeros((C,C))
    else:
        if snakemake.wildcards.model == 'iid':
            V = np.eye(C)
        else:
            V_diag = [float(x) for x in snakemake.wildcards.V_diag.split('_')]
            V = np.diag(V_diag)
            if snakemake.wildcards.V_tril != '0':
                V_tril = [float(x) for x in snakemake.wildcards.V_tril.split('_')]

                ## calculate V
                V[np.tril_indices(C, k=-1)] = V_tril
                for i in range(1, C):
                    for j in range(i):
                        V[i,j] = V[i,j] * math.sqrt(V[i,i] * V[j,j])
                V = V + np.tril(V, k=-1).T

    ## \tr{V S} + pi^T V pi = 0.25
    if np.all(V == 0):
        scale = 0
    else:
        scale = ct_genetic_vc / (np.trace(V @ s) + pi @ V @ pi)
    V = V * scale
    #print(bmatrix(V))

    # save
    np.savetxt(snakemake.output.V, V, delimiter='\t')


    # W
    ct_noise_vc = float(snakemake.wildcards.vc.split('_')[4]) 
    if ct_noise_vc == 0:
        W = np.zeros((C,C))
    else:
        if snakemake.wildcards.model == 'iid':
            W = np.eye(C)
        else:
            W_diag = [float(x) for x in snakemake.wildcards.W_diag.split('_')]
            W = np.diag(W_diag)
            if snakemake.wildcards.W_tril != '0':
                W_tril = [float(x) for x in snakemake.wildcards.W_tril.split('_')]

                ## calculate W
                W[np.tril_indices(C, k=-1)] = W_tril
                for i in range(1, C):
                    for j in range(i):
                        W[i,j] = W[i,j] * math.sqrt(W[i,i] * W[j,j])
                W = W + np.tril(W, k=-1).T

    ## \tr{W S} + pi^T W pi = 0.1
    if np.all(W == 0):
        scale = 0
    else:
        scale = ct_noise_vc / (np.trace(W @ s) + pi @ W @ pi)
    W = W * scale

    np.savetxt(snakemake.output.W, W, delimiter='\t')



if __name__ == '__main__':
    main()
