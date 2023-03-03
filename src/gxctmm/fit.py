from typing import Optional, Tuple
import re, os, sys
import numpy as np, pandas as pd
from numpy import linalg

from ctmm import wald, ctp
from . import log, util

def cal_Vy(hom_g2: float, hom_e2: float, V: np.ndarray, W: np.ndarray, K: np.ndarray, ctnu: np.ndarray
        ) -> np.ndarray:
    '''
    Compute covariance matrix of vectorized Cell Type-specific Pseudobulk

    Parameters:
        hom_g2:   variance of genetic effect shared across cell types
        hom_e2:   variance of env effect shared across cell types
        V:  covariance matrix of cell type-specific genetic effect
        W:  covariance matrix of cell type-specific environment effect
        K:  kinship matrix
        ctnu:   cell type-specific noise variance

    Returns:
        covariance matrix of vectorized Cell Type-specific Pseudobulk V(y)
    '''

    N, C = ctnu.shape
    A = hom_g2 * np.ones((C,C)) + V
    B = hom_e2 * np.ones((C,C)) + W
    Vy = np.kron(K, A) + np.kron(np.eye(N), B) + np.diag( ctnu.flatten() )

    return( Vy )

def LL(y: np.ndarray, K: np.ndarray, X: np.ndarray, ctnu: np.ndarray, hom_g2: float, hom_e2: float, V: np.ndarray, 
        W: np.ndarray) -> float:
    '''
    Loglikelihood function

    Parameters:
        y:  vectorized cell type-specific pseudobulk, vec(Y^T)
        K:  kinship matrix
        X:  design matrix for fixed effects
        ctnu:   cell type-specific noise variance
        hom_g2: variance of genetic effect shared across cell types
        hom_e2: variance of env effect shared across cell types
        V:  covariance matrix of cell type-specific genetic effect
        W:  covariance matrix of cell type-specific env effect
    Returns:
        loglikelihood
    '''

    N, C = ctnu.shape
    Vy = cal_Vy( hom_g2, hom_e2, V, W, K, ctnu )

    # inverse variance
    w, v = linalg.eigh(Vy)
    if ( np.amax(w)/np.amin(w) ) > 1e8 or np.amin(w) < 0:
        return(1e12)
    
    # calculate B matrix
    m1 = X.T @ v @ np.diag(1/w) @ v.T
    m2 = m1 @ X

    # calculate loglikelihood
    det_Vy = np.sum( np.log(w) )
    det_XVyX = linalg.slogdet(m2)[1]
    yBy = y @ v @ np.diag(1/w) @ v.T @ y - y @ m1.T @ linalg.inv(m2) @ m1 @ y
    L = det_Vy + det_XVyX + yBy
    L = 0.5 * L
    return( L )

def extract( out: object, model: str, Y: np.ndarray, K: np.ndarray, P: np.ndarray, ctnu: np.ndarray, 
        fixed_covars: dict ) -> dict:
    '''
    Extract REML optimization resutls

    Parameters:
        out:    OptimizationResult from optimize.minimize
        model:  cell type-specific gene expression model, hom/free/full
        Y:  matrix of cell type-specific pseudobulk
        K:  kinship matrix
        P:  matrix of cell type proportions
        ctnu:   cell type-specific noise variance
        fixed_covars:   design matrices for extra fixed effects
    Returns:
        a dict of model parameters and statitics
    '''

    N, C = Y.shape
    ngam = C*(C+1) // 2

    if model == 'hom':
        hom_g2, hom_e2 = out['x']
        V = W = None
        ct_overall_g_var = ct_overall_e_var = 0
    elif model == 'free':
        hom_g2 = out['x'][0]
        hom_e2 = out['x'][1]
        V = np.diag( out['x'][2:(2+C)] )
        W = np.diag( out['x'][(2+C):(2+2*C)] )
        ct_overall_g_var, ct_specific_g_var = util.ct_random_var( V, P )
        ct_overall_e_var, ct_specific_e_var = util.ct_random_var( W, P )
    elif model == 'full':
        V = np.zeros((C,C))
        V[np.tril_indices(C)] = out['x'][:ngam]
        V = V + V.T
        W = np.zeros((C,C))
        W[np.tril_indices(C)] = out['x'][ngam:(2*ngam)]
        W = W + W.T
        ct_overall_g_var, ct_specific_g_var = util.ct_random_var( V, P )
        ct_overall_e_var, ct_specific_e_var = util.ct_random_var( W, P )

    # beta
    y = Y.flatten()
    X = ctp.get_X( fixed_covars, N, C )

    Vy = cal_Vy( hom_g2, hom_e2, V, W, K, ctnu )
    beta = util.glse( Vy, X, y )
    beta, fixed_vars = util.cal_variance( beta, P, fixed_covars, {}, {} )[:2]

    return( {   'hom_g2':hom_g2, 'hom_e2':hom_e2, 'V':V, 'W':W, 'beta':beta, 
                'ct_overall_g_var':ct_overall_g_var, 'ct_overall_e_var':ct_overall_e_var, 
                'fixed_vars':fixed_vars } )

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
    # projection matrix
    proj = np.eye(N * C) - X @ linalg.inv(X.T @ X) @ X.T

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
    theta = linalg.inv(QTQ) @ Qt

    return(theta)

def _reml(loglike_fun: callable, par: list, model: str, Y: np.ndarray, K: np.ndarray, P: np.ndarray, 
        ctnu: np.ndarray, fixed_covars: dict, method: str, nrep: int) -> dict:
    '''
    Wrapper for running REML

    Parameters:
        loglike_fun: loglikelihood function
        par:    initial parameters
        model:  cell type-specific gene expression model
        Y:  cell type-specific pseudobulk
        K:  kinship matrix
        P:  cell type proportions
        ctnu:   cell type-specific noise variance
        fixed_covars:   design matrices for Extra fixed effectts
        method: optimization method, e.g. BFGS
        nrep:   number of optimization repeats when initial optimization failed
    Returns:
        a dict of model parameters and statistics
    '''

    N, C = Y.shape
    y = Y.flatten()
    X = ctp.get_X( fixed_covars, N, C )

    args = (y, K, X, ctnu)

    out, opt = util.optim( loglike_fun, par, args, method )
    res = extract( out, model, Y, K, P, ctnu, fixed_covars )

    if util.check_optim(opt, res['hom_g2'], res['hom_e2'], res['ct_overall_g_var'], res['ct_overall_e_var'], 
            res['fixed_vars']):
        out, opt = util.re_optim(out, opt, loglike_fun, par, args, method, nrep)
        res = extract( out, model, Y, K, P, ctnu, fixed_covars )

    res['opt'] = opt
    #res['out'] = out

    return( res )

def hom_REML_loglike(par: list, y: np.ndarray, K: np.ndarray, X: np.ndarray, ctnu: np.ndarray) -> float:
    '''
    Loglikelihood for REML under Hom model
    '''

    N, C = ctnu.shape
    hom_g2, hom_e2 = par
    V = np.zeros((C,C))
    W = np.zeros((C,C))

    l = LL(y, K, X, ctnu, hom_g2, hom_e2, V, W)
    return( l )

def hom_REML(Y: np.ndarray, K: np.ndarray, P: np.ndarray, ctnu: np.ndarray, fixed_covars: dict={}, 
        par: list=None, method: str=None, nrep: int=10) -> dict:
    '''
    Fit Hom model with REML

    Parameters:
        Y:  cell type-specific pseudobulk
        K:  kinship matrix
        P:  cell type proportioons
        ctnu:   cell type-specific noise variance
        fixed_covars:   design matrix for Extra fixed effects
        par:    initial parameters
        method: optimization method, e.g. BFGS
        nrep:   number of optimization repeats when initial optimization failed
    Returns:
        a dict of model parameters and statistics
    '''

    log.logger.info('Fitting Hom model with REML')

    N, C = Y.shape
    y = Y.flatten()
    X = ctp.get_X( fixed_covars, N, C )

    if par is None:
        beta = linalg.inv( X.T @ X ) @ ( X.T @ y )
        hom_g2 = np.var(y - X @ beta) / 2
        par = [hom_g2] * 2

    res = _reml( hom_REML_loglike, par, 'hom', Y, K, P, ctnu, fixed_covars, method, nrep=nrep )

    return({'hom_g2':res['hom_g2'], 'hom_e2':res['hom_e2'], 'beta':res['beta'], 
        'opt':res['opt']})

def free_REML_loglike(par:list, y:np.ndarray, K:np.ndarray, X:np.ndarray, ctnu:np.ndarray) -> float:
    '''
    Loglikelihood function for REML under Free model
    '''
    N, C = ctnu.shape
    hom_g2 = par[0]
    hom_e2 = par[1]
    V = np.diag(par[2:(C+2)])
    W = np.diag(par[(C+2):(2*C+2)])

    l = LL(y, K, X, ctnu, hom_g2, hom_e2, V, W)
    return( l )

def free_REML(Y:np.ndarray, K:np.ndarray, P:np.ndarray, ctnu:np.ndarray, fixed_covars:dict={}, 
        par:list=None, method:str=None, nrep:int=10, jk:bool=True) -> Tuple[dict,dict]:
    '''
    Fit Free model using REML

    Parameters:
        Y:  cell type-specific pseudobulk
        K:  kinship matrix
        P:  cell type proportions
        ctnu:   cell type-specific noise variance
        fixed_covars:   design matrices for Extra fixed effects
        par:    initial parameters
        method: optimization method, e.g. BFGS
        nrep:   number of optimization repeats when initial optimization failed
    Returns:
        a tuple of 
            #.  dict of model parameters and statistics
            #.  dict of p values
    '''

    log.logger.info('Fitting Free model with REML')

    N, C = Y.shape
    y = Y.flatten()
    X = ctp.get_X( fixed_covars, N, C )
    n_par = 2 + 2 * C + X.shape[1]

    if par is None:
        beta = linalg.inv( X.T @ X ) @ ( X.T @ y )
        hom_g2 = np.var(y - X @ beta) / 4
        par = [hom_g2] * (2 + 2*C) 

    res = _reml(free_REML_loglike, par, 'free', Y, K, P, ctnu, fixed_covars, method, nrep=nrep)

    if jk:
        jacks = {'ct_beta':[], 'V':[], 'W':[], 'VW':[]}
        for i in range(N):
            Y_jk, K_jk, ctnu_jk, fixed_covars_jk, _, P_jk = util.jk_rmInd(i, Y, K, ctnu, fixed_covars, {}, P)

            res_jk = _reml(free_REML_loglike, par, 'free', Y_jk, K_jk, P_jk, ctnu_jk, 
                    fixed_covars_jk, method, nrep=nrep)

            jacks['ct_beta'].append( res_jk['beta']['ct_beta'] )
            jacks['V'].append( np.diag(res_jk['V']) )
            jacks['W'].append( np.diag(res_jk['W']) )
            jacks['VW'].append( np.append( np.diag(res_jk['V']), np.diag(res_jk['W']) ) )

        var_V = (N-1) * np.cov( np.array(jacks['V']).T, bias=True )
        var_W = (N-1) * np.cov( np.array(jacks['W']).T, bias=True )
        vsr_VW = (N-1) * np.cov( np.array(jacks['VW']).T, bias=True )
        var_ct_beta = (N-1) * np.cov( np.array(jacks['ct_beta']).T, bias=True )

        p = {
                'V': wald.mvwald_test(np.diag(res['V']), np.zeros(C), var_V, n=N, P=n_par),
                'W': wald.mvwald_test(np.diag(res['W']), np.zeros(C), var_W, n=N, P=n_par),
                'VW': wald.mvwald_test( np.append(np.diag(res['V']),np.diag(res['W'])), np.zeros(2*C), 
                    var_VW, n=N, P=n_par),
                'ct_beta': util.wald_ct_beta(res['beta']['ct_beta'], var_ct_beta, n=N, P=n_par),
                }
    else:
        p = {}

    return( res, p )

def full_REML_loglike(par:list, y:np.ndarray, K:np.ndarray, X:np.ndarray, ctnu:np.ndarray) -> float:
    '''
    Loglikelihood function for REML under Full model
    '''
    N, C = ctnu.shape
    ngam = C*(C+1) // 2
    hom_g2 = 0
    hom_e2 = 0
    V = np.zeros((C,C))
    V[np.tril_indices(C)] = par[:ngam]
    V = V + V.T
    W = np.zeros((C,C))
    W[np.tril_indices(C)] = par[ngam:(2*ngam)]
    W = W + W.T

    l = LL(y, K, X, ctnu, hom_g2, hom_e2, V, W)
    return( l )

def full_REML(Y:np.ndarray, K:np.ndarray, P:np.ndarray, ctnu:np.ndarray, fixed_covars:dict={}, 
        par:list=None, method:str=None, nrep:int=10) -> dict:
    '''
    Fit Full model using REML

    Parameters:
        Y:  cell type-specific pseudobulk
        K:  kinship matrix
        P:  cell type proportions
        ctnu:   cell type-specific noise variance
        fixed_covars:   design matrices for Extra fixed effects
        par:    initial parameters
        method: optimization method, e.g. BFGS
        nrep:   number of optimization repeats when initial optimization failed
    Returns:
        dict of model parameters and statistics
    '''

    log.logger.info('Fitting Full model with REML')

    N, C = Y.shape
    ngam = C*(C+1) // 2
    y = Y.flatten()
    X = ctp.get_X( fixed_covars, N, C )

    if par is None:
        beta = linalg.inv( X.T @ X ) @ ( X.T @ y )
        hom_g2 = np.var(y - X @ beta) / 4
        V = W = np.diag(np.ones(C))[np.tril_indices(C)] * hom_g2
        par = list(V) + list(W)

    res = _reml(full_REML_loglike, par, 'full', Y, K, P, ctnu, fixed_covars, method, nrep=nrep)

    return( res )

def _free_he(Y: np.ndarray, K: np.ndarray, ctnu: np.ndarray, P: np.ndarray) -> dict:
    N, C = Y.shape
    X = ctp.get_X({}, N, C)

    theta = he_ols(Y, K, X, ctnu, 'free')
    hom_g2, hom_e2 = theta[0], theta[1]
    V, W = np.diag(theta[2:(C+2)]), np.diag(theta[(C+2):(C*2+2)])

    # GLS to get beta
    Vy = cal_Vy( hom_g2, hom_e2, V, W, K, ctnu )
    beta = util.glse( Vy, X, Y.flatten() )
    # calcualte variance of fixed and random effects, and convert to dict
    beta, fixed_vars = util.cal_variance(beta, P, {}, {}, {})[:2]
    ct_overall_g_var, ct_specific_g_var = util.ct_random_var( V, P )
    ct_overall_e_var, ct_specific_e_var = util.ct_random_var( W, P )

    return( {'hom_g2':hom_g2, 'hom_e2':hom_e2, 'V':V, 'W':W, 'beta':beta, 'fixed_vars':fixed_vars, 
            'ct_overall_g_var':ct_overall_g_var, 'ct_specific_g_var':ct_specific_g_var, 
            'ct_overall_e_var':ct_overall_e_var, 'ct_specific_e_var':ct_specific_e_var} )

def free_HE(Y: np.ndarray, K: np.ndarray, ctnu: np.ndarray, P: np.ndarray) -> Tuple[dict, dict]:
    '''
    Fitting Free model with HE

    Parameters:
        Y:    cell type-specific pseudobulk (no header no index)
        K:    kinship matrix 
        ctnu: cell type-specific noise variance (no header no index)
        P:    cell type proportions

    Returns:
        a tuple of
            #.  dictionary of parameter estimates
            #.  dictionary of p values 
    '''

    log.logger.info('Fitting Free model with HE')

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
        jacks['VW'].append( np.append( np.diag(out_jk['V']), np.diag(out_jk['W']) ) )

    var_V = (N-1) * np.cov( np.array(jacks['V']).T, bias=True )
    var_W = (N-1) * np.cov( np.array(jacks['W']).T, bias=True )
    var_VW = (N-1) * np.cov( np.array(jacks['VW']).T, bias=True )
    var_ct_beta = (N-1) * np.cov( np.array(jacks['ct_beta']).T, bias=True )

    p = {   'V': wald.mvwald_test(np.diag(out['V']), np.zeros(C), var_V, n=N, P=n_par),
            'W': wald.mvwald_test(np.diag(out['W']), np.zeros(C), var_W, n=N, P=n_par),
            'VW': wald.mvwald_test( np.append(np.diag(out['V']),np.diag(out['W'])), np.zeros(2*C), 
                var_VW, n=N, P=n_par),
            'ct_beta': util.wald_ct_beta( out['beta']['ct_beta'], var_ct_beta, n=N, P=n_par )
            }
    return(out, p)

def full_HE(Y: np.ndarray, K: np.ndarray, ctnu: np.ndarray, P: np.ndarray) -> dict:
    '''
    Fitting Full model with HE

    Parameters:
        Y:    cell type-specific pseudobulk (no header no index)
        K:    kinship matrix 
        ctnu: cell type-specific noise variance (no header no index)
        P:    cell type proportions

    Returns:
        a dictionary of parameter estimates
    '''

    log.logger.info('Fitting Full model with HE')

    N, C = Y.shape 
    ntril = (C-1) * C // 2
    X = ctp.get_X({}, N, C)

    theta = he_ols(Y, K, X, ctnu, 'full')
    V, W = np.diag(theta[:C]), np.diag(theta[C:(C*2)])
    V[np.triu_indices(C,k=1)] = theta[(C*2):(C*2 + ntril)]
    V = V + V.T - np.diag(theta[:C])
    W[np.triu_indices(C,k=1)] = theta[(C*2 + ntril):(C*2 + ntril*2)]
    W = W + W.T - np.diag(theta[C:(C*2)])

    ct_overall_g_var, ct_specific_g_var = util.ct_random_var( V, P )
    ct_overall_e_var, ct_specific_e_var = util.ct_random_var( W, P )

    he = {'V': V, 'W': W, 'ct_overall_g_var':ct_overall_g_var, 'ct_overall_e_var':ct_overall_e_var}
    return( he )

