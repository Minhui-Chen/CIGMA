from typing import Optional, Tuple
import re, os, sys
import numpy as np, pandas as pd

from ctmm import wald, ctp
from . import log, util

def cal_Vy(A: np.ndarray, B: np.ndarray, K: np.ndarray, ctnu: np.ndarray) -> np.ndarray:
    '''
    Compute covariance matrix of vectorized Cell Type-specific Pseudobulk

    Parameters:
        A:  covariance matrix of genetic effect, sigma_g2 * J_C + V
        B:  covariance matrix of environment effect, signma_e2 * J_C + W
        K:  kinship matrix
        ctnu:   cell type-specific noise variance

    Returns:
        covariance matrix of vectorized Cell Type-specific Pseudobulk V(y)
    '''

    N, C = ctnu.shape
    Vy = np.kron(K, A) + np.kron(np.eye(N), B) + np.diag( vs.flatten() )

    return( Vy )

def he_ols(Y: np.ndarray, K: np.ndarray, X: np.ndarray, ctnu: np.ndarray, model: str
        ) -> np.ndarray:
    '''
    Perform OLS in HE

    Parameters:
        Y:  N * C matrix of Cell Type-specific Pseudobulk
        K:  kinship matrix
        X:  design matrix for fixed effects
        ctnu: cell type-specific noise variance
        model:  free / full
    Returns:
        a tuple of 
            #. 
            #.
    '''

    N, C = Y.shape
    y = Y.flatten()
    D = np.diag( nu.flatten() )
    # projection matrix
    proj = np.eye(N * C) - X @ np.linalg.inv(X.T @ X) @ X.T

    # vec(M @ A @ M)^T @ vec(M @ B @ M) = vec(M @ A)^T @ vec((M @ B)^T)
    # when A, B, and M are symmetrix
    # proj @ y @ y^T @ proj - proj @ D @ proj
    t = np.outer( proj @ y, y ) - proj * ctnu.flatten() 
    t = t.flatten()

    # build Q: list of coefficients
    def L_f(C, c1, c2):
        # fun to build L matrix
        L = np.zeros((C,C))
        L[c1,c2] = 1
        return( L )

    if model == 'free':
        Q = [   proj @ np.kron(K, np.ones((C,C))), 
                proj @ np.kron(np.eye(N), np.ones((C,C))) ]
        for c in range(C):
            L = L_f(C, c, c)
            Q.append( proj @ np.kron(K, L) )
        for c in range(C):
            L = L_f(C, c, c)
            Q.append( proj @ np.kron(np.eye(N), L) )
    elif model == 'full':
        Q = []
        for c in range(C):
            L = L_f(C, c, c)
            Q.append( proj @ np.kron(K, L) )
        for c in range(C):
            L = L_f(C, c, c)
            Q.append( proj @ np.kron(np.eye(N), L) )
        for i in range(C-1):
            for j in range(i+1,C):
                L = L_f(C, i, j) + L_f(C, j, i)
                Q.append( proj @ np.kron(K, L) )
        for i in range(C-1):
            for j in range(i+1,C):
                L = L_f(C, i, j) + L_f(C, j, i)
                Q.append( proj @ np.kron(np.eye(N), L) )
    
    QTQ = np.array([m.flatten('F') for m in Q]) @ np.array([m.flatten() for m in Q]).T
    Qt = np.array([m.flatten('F') for m in Q]) @ t
    theta = np.linalg.inv(QTQ) @ Qt

    return(theta)

def _free_he(Y: np.ndarray, K: np.ndarray, ctnu: np.ndarray, P: np.ndarray) -> dict:
    N, C = Y.shape
    X = ctp.get_X({}, N, C)

    theta = he_ols(Y, K, X, ctnu, 'free')
    sigma_g2, sigma_e2 = theta[0], theta[1]
    V, W = np.diag(theta[2:(C+2)]), np.diag(theta[(C+2):(C*2+2)])

    # GLS to get beta
    A = np.ones((C,C)) * sigma_g2 + V
    B = np.ones((C,C)) * sigma_e2 + W
    Vy = cal_Vy( A, B, K, ctnu )
    beta = util.glse( Vy, X, y )
    # calcualte variance of fixed and random effects, and convert to dict
    beta, fixed_vars = util.cal_variance(beta, P, {}, {}, {})[:2]
    ct_overall_g_var, ct_specific_g_var = util.ct_random_var( V, P )
    ct_overall_e_var, ct_specific_e_var = util.ct_random_var( W, P )

    return( {'sigma_g2':sigma_g2, 'sigma_e2':sigma_e2, 'V':V, 'W':W, 'beta':beta, 'fixed_vars':fixed_vars, 
            'ct_overall_g_var':ct_overall_g_var, 'ct_specific_g_var':ct_specific_g_var, 
            'ct_overall_e_var':ct_overall_e_var, 'ct_specific_e_var':ct_specific_e_var} )

def free_HE(Y_f: str, K_f str, ctnu_f: str, P_f: str) -> Tuple[dict, dict]:
    '''
    Fitting Free model with HE

    Parameters:
        Y_f:    file of cell type-specific pseudobulk (no header no index)
        K_f:    file of kinship matrix 
        ctnu_f: file of cell type-specific noise variance (no header no index)
        P_f:    file of cell type proportions

    Returns:
        a tuple of
            #.  dictionary of parameter estimates
            #.  dictionary of p values 
    '''

    log.logger.info('Fitting Free model with HE')

    Y = np.loadtxt(Y_f) 
    K = np.loadtxt(K_f)
    ctnu = np.loadtxt(ctnu_f)
    P = np.loadtxt(P_f)

    N, C = Y.shape 
    X = ctp.get_X({}, N, C)
    n_par = 2 + 2 * C + X.shape[1]

    out = _free_he(Y, K, ctnu, P)
    out['nu'] = ( ctnu * (P ** 2) ).sum(axis=1) 

    # jackknife
    jacks = {'ct_beta':[], 'V':[], 'W':[], 'VW':[]}
    for i in range(N):
        Y_jk, K_jk, ctnu_jk, _, _, P_jk = util.jk_rmInd(i, Y, K, ctnu, P=P)
        out_jk = _free_he( Y_jk, K_jk, ctnu_jk, P_jk )

        jacks['ct_beta'].append( out_jk['beta']['ct_beta'] )
        jacks['V'].append( np.diag(out_jk['V']) )
        jacks['W'].append( np.diag(out_jk['W']) )
        jacks['VW'].append( np.append( (np.diag(out_jk['V']), np.diag(out_jk['W'])) ) )

    var_V = (N-1) * np.cov( np.array(jacks['V']).T, bias=True )
    var_W = (N-1) * np.cov( np.array(jacks['W']).T, bias=True )
    var_VW = (N-1) * np.cov( np.array(jacks['VW']).T, bias=True )
    var_ct_beta = (N-1) * np.cov( np.array(jacks['ct_beta']).T, bias=True )

    p = {   'V': wald.mvwald_test(np.diag(out['V']), np.zeros(C), var_V, n=N, P=n_par),
            'W': wald.mvwald_test(np.diag(out['W']), np.zeros(C), W_var, n=N, P=n_par),
            'VW': wald.mvwald_test( np.append(np.diag(out['V']),np.diag(out['W'])), np.zeros(2*C), 
                covar[2:(2*C+2),2:(2*C+2)], n=N, P=n_par),
            'ct_beta': util.wald_ct_beta( out['beta']['ct_beta'], var_ct_beta, n=N, P=n_par )
            }
    return(out, p)

def full_HE(Y_f: str, K_f: str, ctnu_f: str, P_f: str) -> dict:
    '''
    Fitting Full model with HE

    Parameters:
        Y_f:    file of cell type-specific pseudobulk (no header no index)
        K_f:    file of kinship matrix 
        ctnu_f: file of cell type-specific noise variance (no header no index)
        P_f:    file of cell type proportions

    Returns:
        a dictionary of parameter estimates
    '''

    log.logger.info('Fitting Full model with HE')

    Y = np.loadtxt(Y_f) 
    K = np.loadtxt(K_f)
    ctnu = np.loadtxt(ctnu_f)

    N, C = Y.shape 
    ntril = (C-1) * C // 2
    X = ctp.get_X({}, N, C)

    theta = he_ols(Y, K, X, ctnu, 'full')
    V, W = np.diag(theta[:C]), np.diag(theta[C:(C*2)])
    V[np.triu_indices(C,k=1)] = theta[(C*2):(C*2 + ntril)]
    V = V + V.T - np.diag(theta[:C])
    W[np.triu_indices(C,k=1)] = theta[(C*2 + ntril):(C*2 + ntril*2)]
    W = W + W.T - np.diag(theta[C:(C*2)])

    he = {'V': V, 'W': W}
    return( he )

