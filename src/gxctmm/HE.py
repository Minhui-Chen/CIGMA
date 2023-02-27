import re, os, sys
import numpy as np, pandas as pd
from ctmm import wald 

def HE(cty, Z, X, ctvs, model):
    '''
    model: hom, free, full
    '''
    N, C = cty.shape
    y = cty.flatten()
    # project out cell type main effect
    proj = np.eye(N * C) - X @ np.linalg.inv(X.T @ X) @ X.T
    y_p = proj @ y
    # y' y'^T - M diag(\nu) M
    Y_m = np.outer(y_p, y_p) - proj @ np.diag(ctvs.flatten()) @ proj ########## optimize
    Y = Y_m.flatten('F')

    # T
    T = []
    if model != 'full':
        T.append( (proj @ np.kron( Z, np.ones((C,C)) ) @ proj).flatten('F') ) # sigma_g2
        T.append( (proj @ np.kron( np.eye(N), np.ones((C,C)) ) @ proj).flatten('F') )

    if model in ['free','full']:
        for c in range(C):
            zero_m = np.zeros((C,C))
            zero_m[c,c] = 1
            T.append( (proj @ np.kron( Z, zero_m ) @ proj).flatten('F') ) # V diagonal
        for c in range(C):
            zero_m = np.zeros((C,C))
            zero_m[c,c] = 1
            T.append( (proj @ np.kron( np.eye(N), zero_m ) @ proj).flatten('F') ) # W diagonal
        if model == 'full':
            for i in range(C-1):
                for j in range(i+1,C):
                    zero_m = np.zeros((C,C))
                    zero_m[i,j] = 1
                    T.append( 2 * (proj @ np.kron( Z, zero_m ) @ proj).flatten('F') ) # V off-diagonal
            for i in range(C-1):
                for j in range(i+1,C):
                    zero_m = np.zeros((C,C))
                    zero_m[i,j] = 1
                    T.append( 2 * (proj @ np.kron( np.eye(N), zero_m ) @ proj).flatten('F') ) # W off-diagonal
    T = np.array(T).T
    
    # theta: sigma_g2, sigma_e2 in Hom; sigma_g2, sigma_e2, V diag, W diag in Free;
    # V diag, W diag, V off diag, and W off diag in Full
    theta = np.linalg.inv(T.T @ T) @ T.T @ Y

    return(theta)

def free_HE(cty_f, Z_f, ctnu_f):
    cty = np.loadtxt(cty_f) 
    Z = np.loadtxt(Z_f)
    ctvs = np.loadtxt(ctnu_f)
    N, C = cty.shape # cell type number
    X = np.kron(np.ones((N,1)), np.eye(C))
    X1 = np.kron(np.ones((N-1,1)), np.eye(C))
    n_par = 2 + 2 * C

    theta = HE(cty, Z, X, ctvs, 'free')
    sigma_g2, sigma_e2 = theta[0], theta[1]
    V, W = np.diag(theta[2:(C+2)]), np.diag(theta[(C+2):(C*2+2)])
    jacks = [HE(np.delete(cty,i,axis=0), np.delete(np.delete(Z,i,axis=0),i,axis=1),
        X1, np.delete(ctvs,i,axis=0), 'free') for i in range(N)]
    covar = (len(jacks)-1.0) * np.cov(np.array(jacks).T, bias=True)
    sigma_g2_var = covar[0,0]
    sigma_e2_var = covar[1,1]
    V_var = covar[2:(C+2), 2:(C+2)]
    W_var = covar[(C+2):(C*2+2), (C+2):(C*2+2)]

    he = {'sigma_g2': sigma_g2, 'sigma_e2':sigma_e2, 'V': V, 'W': W}
    p = {   'sigma_g2': wald.wald_test(sigma_g2, 0, sigma_g2_var, N-n_par),
            'sigma_e2': wald.wald_test(sigma_e2, 0, sigma_e2_var, N-n_par),
            'V': wald.mvwald_test(np.diag(V), np.zeros(C), V_var, n=N, P=n_par),
            'W': wald.mvwald_test(np.diag(W), np.zeros(C), W_var, n=N, P=n_par),
            'VW': wald.mvwald_test( np.append(np.diag(V),np.diag(W)), np.zeros(2*C), 
                covar[2:(2*C+2),2:(2*C+2)], n=N, P=n_par),
            }
    return(he, p)

def full_HE(cty_f, Z_f, ctnu_f):
    cty = np.loadtxt(cty_f) 
    Z = np.loadtxt(Z_f)
    ctvs = np.loadtxt(ctnu_f)
    N, C = cty.shape # cell type number
    ntril = (C-1) * C // 2
    X = np.kron(np.ones((N,1)), np.eye(C))

    theta = HE(cty, Z, X, ctvs, 'full')
    V, W = np.diag(theta[:C]), np.diag(theta[C:(C*2)])
    V[np.triu_indices(C,k=1)] = theta[(C*2):(C*2 + ntril)]
    V = V + V.T - np.diag(theta[:C])
    W[np.triu_indices(C,k=1)] = theta[(C*2 + ntril):(C*2 + ntril*2)]
    W = W + W.T - np.diag(theta[C:(C*2)])

    he = {'V': V, 'W': W}
    return( he )

def main():

    batch = snakemake.params.batches
    outs = [re.sub('/rep/', f'/rep{i}/', snakemake.params.out) for i in batch]
    for cty_f, Z_f, ctnu_f, out_f in zip(
            [line.strip() for line in open(snakemake.input.cty)],
            [line.strip() for line in open(snakemake.input.Z)],
            [line.strip() for line in open(snakemake.input.ctnu)], 
            outs):
        print(cty_f, Z_f, ctnu_f)
        os.makedirs(os.path.dirname(out_f), exist_ok=True)

        ## Free
        free_he, free_he_wald = free_HE(cty_f, Z_f, ctnu_f)

        # Full
        full_he = full_HE(cty_f, Z_f, ctnu_f)

        # save
        np.save(out_f,
                {'free':free_he, 'full':full_he, 
                'wald': {'free':free_he_wald} 
                }
            )
        
    with open(snakemake.output.out, 'w') as f:
        f.write('\n'.join(outs))  


if __name__ == '__main__':
    main()
