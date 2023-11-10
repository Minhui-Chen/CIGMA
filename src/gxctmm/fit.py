from typing import Optional, Tuple, Union
import re, os, sys
import numpy as np, pandas as pd
from numpy import linalg
from scipy import linalg as sla
import dask.array as da
from memory_profiler import profile 

from ctmm import wald
from . import log, util

def _ZZT(random_covars: dict) -> dict:
    """
    Compute Z @ Z.T

    Parameters:
        random_covars:  design matrices for Extra random effects

    Returns:
        # Z @ Z.T of design matrices for Extra random effects
    """
    random_covars_ZZT = {}
    for key, Z in random_covars.items():
        random_covars_ZZT[key] = Z @ Z.T
    return random_covars_ZZT

def cal_Vy(hom_g2: float, hom_e2: float, V: np.ndarray, W: np.ndarray, r2:dict,
           K: np.ndarray, ctnu: np.ndarray, random_covars_ZZT:dict, Kt: Optional[np.ndarray] = None, 
           hom_gt2: Optional[float] = None, Vt: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Compute covariance matrix of vectorized Cell Type-specific Pseudobulk

    Parameters:
        hom_g2:   variance of genetic effect shared across cell types
        hom_e2:   variance of env effect shared across cell types
        V:  covariance matrix of cell type-specific genetic effect
        W:  covariance matrix of cell type-specific environment effect
        r2: variance of Extra random effect
        K:  kinship matrix
        ctnu:   cell type-specific noise variance
        random_covars_ZZT:  Z @ Z.T of design matrices for Extra random effects
        Kt: trans kinship matrix
        hom_gt2:  variance of trans genetic effect shared across cell types
        Vt: covariance matrix of cell type-specific trans genetic effect

    Returns:
        covariance matrix of vectorized Cell Type-specific Pseudobulk V(y)
    """

    N, C = ctnu.shape
    A = hom_g2 * np.ones((C,C), dtype='int8') + V
    B = hom_e2 * np.ones((C,C), dtype='int8') + W
    Vy = np.kron(K, A) + np.kron(np.eye(N, dtype='int8'), B) + np.diag( ctnu.flatten() )
    for key in random_covars_ZZT.keys():
        ZZT = random_covars_ZZT[key]
        if isinstance(r2[key], float):
            Vy += np.kron(ZZT, np.ones((C,C))) * r2[key] # shared random effect
        else:
            Vy += np.kron(ZZT, np.diag(r2[key])) # cell type-specific random effect
    
    if Kt is not None:
        At = hom_gt2 * np.ones((C,C), dtype='int8') + Vt
        Vy += np.kron(Kt, At)

    return Vy


def LL(y: np.ndarray, K: np.ndarray, X: np.ndarray, ctnu: np.ndarray, random_covars_ZZT:dict,
       hom_g2: float, hom_e2: float, V: np.ndarray, W: np.ndarray, r2:dict) -> float:
    """
    Loglikelihood function

    Parameters:
        y:  vectorized cell type-specific pseudobulk, vec(Y^T)
        K:  kinship matrix
        X:  design matrix for fixed effects
        ctnu:   cell type-specific noise variance
        random_covars_ZZT:  Z @ Z.T of design matrices for Extra random effects
        hom_g2: variance of genetic effect shared across cell types
        hom_e2: variance of env effect shared across cell types
        V:  covariance matrix of cell type-specific genetic effect
        W:  covariance matrix of cell type-specific env effect
        r2: variance of Extra random effect
    Returns:
        loglikelihood
    """

    N, C = ctnu.shape
    Vy = cal_Vy( hom_g2, hom_e2, V, W, r2, K, ctnu, random_covars_ZZT )

    # inverse variance
    w, v = linalg.eigh(Vy)
    if ( np.amax(w)/np.amin(w) ) > 1e8 or np.amin(w) < 0:
        return 1e12
    
    # calculate B matrix
    m1 = X.T @ v @ np.diag(1/w) @ v.T 
    m2 = m1 @ X

    # calculate loglikelihood
    det_Vy = np.sum( np.log(w) )
    det_XVyX = linalg.slogdet(m2)[1]
    yBy = y @ v @ np.diag(1/w) @ v.T @ y - y @ m1.T @ linalg.inv(m2) @ m1 @ y
    L = det_Vy + det_XVyX + yBy
    L = 0.5 * L
    return L

def _get_r2(r2:list, random_covars:dict, C:int) -> dict:
    """
    Covert list of r2 to dictionary
    """
    if len(r2) == len(random_covars.keys()):
        shared = True
    else:
        shared = False

    r2_d = {}
    for i, key in enumerate(sorted(random_covars.keys())):
        if shared:
            r2_d[key] = r2[i]
        else:
            r2_d[key] = r2[(i * C):((i + 1) * C)]
    return r2_d

def extract( out: object, model: str, Y: np.ndarray, K: np.ndarray, P: np.ndarray,
             ctnu: np.ndarray, fixed_covars: dict, random_covars:dict ) -> dict:
    """
    Extract REML optimization resutls

    Parameters:
        out:    OptimizationResult from optimize.minimize
        model:  cell type-specific gene expression model, hom/free/full
        Y:  matrix of cell type-specific pseudobulk
        K:  kinship matrix
        P:  matrix of cell type proportions
        ctnu:   cell type-specific noise variance
        fixed_covars:   design matrices for extra fixed effects
        random_covars:  design matrices for extra random effects
    Returns:
        a dict of model parameters and statistics
    """

    N, C = Y.shape
    ngam = C*(C+1) // 2

    if model == 'hom':
        hom_g2, hom_e2 = out['x'][:2]
        V = W = np.zeros((C,C))
        ct_overall_g_var = ct_overall_e_var = 0
        r2 = _get_r2(out['x'][2:], random_covars, C)
    elif model == 'free':
        hom_g2 = out['x'][0]
        hom_e2 = out['x'][1]
        V = np.diag( out['x'][2:(2+C)] )
        W = np.diag( out['x'][(2+C):(2+2*C)] )
        ct_overall_g_var, ct_specific_g_var = util.ct_random_var( V, P )
        ct_overall_e_var, ct_specific_e_var = util.ct_random_var( W, P )
        r2 = _get_r2(out['x'][(2+2*C):], random_covars, C)
    elif model == 'freeW':
        hom_g2 = out['x'][0]
        hom_e2 = out['x'][1]
        W = np.diag( out['x'][2:(2+C)] )
        V = np.zeros((C,C))
        ct_overall_g_var, ct_specific_g_var = util.ct_random_var( V, P )
        ct_overall_e_var, ct_specific_e_var = util.ct_random_var( W, P )
        r2 = _get_r2(out['x'][(2+C):], random_covars, C)
    elif model == 'full':
        hom_g2 = 0
        hom_e2 = 0
        V = np.zeros((C,C))
        V[np.tril_indices(C)] = out['x'][:ngam]
        V = V + V.T
        W = np.zeros((C,C))
        W[np.tril_indices(C)] = out['x'][ngam:(2*ngam)]
        W = W + W.T
        ct_overall_g_var, ct_specific_g_var = util.ct_random_var( V, P )
        ct_overall_e_var, ct_specific_e_var = util.ct_random_var( W, P )
        r2 = _get_r2(out['x'][(2*ngam):], random_covars, C)

    # beta
    if isinstance(list(r2.values())[0], float):
        shared = True
    else:
        shared = False
    y = Y.flatten()
    X = util.get_X( fixed_covars, N, C, fixed_shared=shared )

    Vy = cal_Vy( hom_g2, hom_e2, V, W, r2, K, ctnu, _ZZT(random_covars) )
    beta = util.glse( Vy, X, y )
    beta, fixed_vars = util.cal_variance( beta, P, fixed_covars, r2, random_covars )[:2]

    return( {   'hom_g2':hom_g2, 'hom_e2':hom_e2, 'V':V, 'W':W, 'beta':beta, 
                'ct_overall_g_var':ct_overall_g_var, 'ct_overall_e_var':ct_overall_e_var, 
                'fixed_vars':fixed_vars } )

# def _he_Qloop_full(Q: list, X: np.ndarray, t:np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
#     """
#     Compute QTQ and Qt
#     """
#     QTQ = []
#     Qt = []
#     X_p = X @ linalg.inv(X.T @ X)
#     for m1 in Q:
#         m1 = np.kron(m1[0], m1[1])
#         m1 = ( m1 - X_p @ (X.T @ m1) ).flatten('F')
#         Qt.append( m1 @ t )
#         QTQ_m = []
#         log.logger.info('Timer')
#         for m2 in Q:
#             m2 = np.kron(m2[0], m2[1])
#             m2 = ( m2 - X_p @ (X.T @ m2) ).flatten()
#             QTQ_m.append(m1 @ m2)
#         QTQ.append( QTQ_m )
#     return np.array(QTQ), np.array(Qt)

def _pMp(X:np.ndarray, X_inv: np.ndarray, A:np.ndarray, B:np.ndarray) -> np.ndarray:
    """
    Compute proj @ np.kron(A,B) @ proj
    """
    M = np.kron(A, B)
    p_M = M - X @ X_inv @ (X.T @ M)
    M = p_M - p_M @ X @ X_inv @ X.T
    return M


def he_ols(Y: np.ndarray, K: np.ndarray, X: np.ndarray, ctnu: np.ndarray,
           model: str, random_covars: dict={}, Kt:Optional[np.ndarray] = None, 
           random_shared:bool = True, dtype:str = None) -> np.ndarray:
    """
    Perform OLS in HE

    Parameters:
        Y:  N * C matrix of Cell Type-specific Pseudobulk
        K:  kinship matrix
        X:  design matrix for fixed effects
        ctnu: cell type-specific noise variance
        model:  free / full
        random_covars:  design matrices for Extra random effects
        Kt: kinship matrix for trans-eQTLs (currently only for Free model)
        random_shared: whether Extra random effects are shared across cell types
        dtype:  data type e.g. float32 to save memory
    Returns:
        OLS estimates
    """

    if dtype is None:
        dtype = 'float64'
    else:
        Y, K, X, ctnu = Y.astype(dtype), K.astype(dtype), X.astype(dtype), ctnu.astype(dtype)
        if Kt is not None:
            Kt = Kt.astype(dtype)

    N, C = Y.shape
    y = Y.flatten()
    X_inv = linalg.inv(X.T @ X)
    n_random = len(random_covars.keys()) if random_shared else len(random_covars.keys()) * C

    # projection matrix
    proj = np.eye(N * C, dtype='int8') - X @ X_inv @ X.T # X: 1_N \otimes I_C append sex \otimes 1_C 

    # vec(M @ A @ M)^T @ vec(M @ B @ M) = vec(M @ A)^T @ vec((M @ B)^T)
    # when A, B, and M are symmetric
    # proj @ y @ y^T @ proj - proj @ D @ proj
    y_p = proj @ y
    ctnu_p = proj * ctnu.flatten()
    t = np.outer( y_p, y_p ) - ctnu_p + ctnu_p @ X @ X_inv @ X.T

    # build Q: list of coefficients
    if model == 'free':
        if N * C < (10000 * 100000):
            Q = []
            # shared
            Q.append(_pMp(X, X_inv, K, np.ones((C,C), dtype='int8')))
            if Kt is not None:
                Q.append(_pMp(X, X_inv, Kt, np.ones((C,C), dtype='int8')))
            Q.append(_pMp(X, X_inv, np.eye(N, dtype='int8'), np.ones((C,C), dtype='int8')))

            # cell type-specific
            # cis-eQTLs
            for c in range(C):
                L = util.L_f(C, c, c)
                Q.append(_pMp(X, X_inv, K, L))
            
            # trans-eQTLs
            if Kt is not None:
                for c in range(C):
                    L = util.L_f(C, c, c)
                    Q.append(_pMp(X, X_inv, Kt, L))
            
            # environment
            for c in range(C):
                L = util.L_f(C, c, c)
                Q.append(_pMp(X, X_inv, np.eye(N, dtype='int8'), L))
            
            # extra
            for key in sorted(random_covars.keys()):
                Z = random_covars[key]
                ZZT = Z @ Z.T
                if random_shared:
                    Q.append(_pMp(X, X_inv, ZZT, np.ones((C,C), dtype='int8')))
                else:
                    for c in range(C):
                        L = util.L_f(C, c, c)
                        Q.append(_pMp(X, X_inv, ZZT, L))

            QTQ = np.tensordot(Q, Q, axes=([1, 2], [1, 2]))
            QTt = np.tensordot(Q, t, axes=([1, 2], [0, 1]))
        # else:
        #     log.logger.info('Memory saving mode')
        #     Q_shape = (2 * (C + 1) + n_random, (N * C) ** 2)
        #     # use memmap (don't run in parallel)
        #     tmpfn = util.generate_tmpfn()
        #     Q = np.memmap(tmpfn, dtype=dtype, mode="w+", shape=Q_shape)
        #
        #     Q[0] = _pMp(X, X_inv, K, np.ones((C,C), dtype='int8')).flatten('F') # proj @ np.kron(K, J_C) @ proj
        #     Q[1] = _pMp(X, X_inv, np.eye(N, dtype='int8'), np.ones((C,C), dtype='int8')).flatten('F') # I_N \ot J_C
        #
        #     k = 2
        #     for c in range(C):
        #         L = util.L_f(C, c, c)
        #         Q[k] = _pMp(X, X_inv, K, L).flatten('F')
        #         k += 1
        #     for c in range(C):
        #         L = util.L_f(C, c, c)
        #         Q[k] = _pMp(X, X_inv, np.eye(N, dtype='int8'), L).flatten('F')
        #     k += 1
        #
        #     Q = Q.T
        #     Q.flush()
        #
        #     QTQ = Q.T @ Q
        #     QTt = Q.T @ t.flatten()

    elif model == 'full':
        if N * C < 5000000000:
            log.logger.info('Making Q')
            Q = []
            for c in range(C):
                L = util.L_f(C, c, c)
                Q.append( _pMp(X, X_inv, K, L) )
            for c in range(C):
                L = util.L_f(C, c, c)
                Q.append( _pMp(X, X_inv, np.eye(N, dtype='int8'), L) )
            for i in range(C - 1):
                for j in range(i + 1, C):
                    L = util.L_f(C, i, j) + util.L_f(C, j, i)
                    Q.append( _pMp(X, X_inv, K, L) )
            for i in range(C - 1):
                for j in range(i + 1, C):
                    L = util.L_f(C, i, j) + util.L_f(C, j, i)
                    Q.append( _pMp(X, X_inv, np.eye(N, dtype='int8'), L) )

            for key in sorted(random_covars.keys()):
                Z = random_covars[key]
                ZZT = Z @ Z.T
                if random_shared:
                    Q.append( _pMp(X, X_inv, ZZT, np.ones((C, C), dtype='int8')) )
                else:
                    for c in range(C):
                        L = util.L_f(C, c, c)
                        Q.append( _pMp(X, X_inv, ZZT, L) )

            log.logger.info('Calculating Q products')
            QTQ = np.tensordot(Q, Q, axes=([1, 2], [1, 2]))
            QTt = np.tensordot(Q, t, axes=([1, 2], [0, 1]))

        # else:
        #     log.logger.info('Memmory saving mode')
        #     Q_shape = (C * (C + 1) + n_random, (N * C) ** 2)
        #     # use memmap (don't run in parallel)
        #     tmpfn = util.generate_tmpfn()
        #     Q = np.memmap(tmpfn, dtype=dtype, mode="w+", shape=Q_shape)
        #
        #     k = 0
        #     for c in range(C):
        #         L = util.L_f(C, c, c)
        #         Q[k] = _pMp(X, X_inv, K, L).flatten('F')
        #         k += 1
        #     for c in range(C):
        #         L = util.L_f(C, c, c)
        #         Q[k] = _pMp(X, X_inv, np.eye(N, dtype='int8'), L).flatten('F')
        #         k += 1
        #     for i in range(C-1):
        #         for j in range(i+1,C):
        #             L = util.L_f(C, i, j) + util.L_f(C, j, i)
        #             Q[k] = _pMp(X, X_inv, K, L).flatten('F')
        #             k += 1
        #     for i in range(C-1):
        #         for j in range(i+1,C):
        #             L = util.L_f(C, i, j) + util.L_f(C, j, i)
        #             Q[k] = _pMp(X, X_inv, np.eye(N, dtype='int8'), L).flatten('F')
        #             k += 1
        #
        #     for key in sorted(random_covars.keys()):
        #         Z = random_covars[key]
        #         Q[k] = _pMp(X, X_inv, Z @ Z.T, np.ones((C,C), dtype='int8')).flatten('F')
        #         k += 1
        #
        #     Q = Q.T
        #     Q.flush()
        #
        #     QTQ = Q.T @ Q
        #     QTt = Q.T @ t.flatten()

    # theta
    theta = linalg.inv(QTQ) @ QTt

    if isinstance(Q, np.memmap):
        del Q
        os.remove(tmpfn)

    return theta


def _reml(model:str, par: list, Y: np.ndarray, K: np.ndarray,
          P: np.ndarray, ctnu: np.ndarray, fixed_covars: dict, random_covars: dict,
          shared: bool, method: str, nrep: int) -> dict:
    """
    Wrapper for running REML

    Parameters:
        model:  cell type-specific gene expression model
        par:    initial parameters
        Y:  cell type-specific pseudobulk
        K:  kinship matrix
        P:  cell type proportions
        ctnu:   cell type-specific noise variance
        fixed_covars:   design matrices for Extra fixed effects
        random_covars:   design matrices for Extra random effects
        shared: whether Extra fixed and random effects are shared across cell types
        method: optimization method, e.g. BFGS
        nrep:   number of optimization repeats when initial optimization failed
    Returns:
        a dict of model parameters and statistics
    """

    N, C = Y.shape
    y = Y.flatten()
    X = util.get_X( fixed_covars, N, C, fixed_shared=shared )

    funs = {
        'hom':  hom_REML_loglike,
        'freeW':    freeW_REML_loglike,
        'free': free_REML_loglike,
        'full': full_REML_loglike,
    }
    loglike_fun = funs[model]
    args = (y, K, X, ctnu, _ZZT(random_covars), shared)

    out, opt = util.optim( loglike_fun, par, args, method )
    res = extract( out, model, Y, K, P, ctnu, fixed_covars, random_covars )

    if util.check_optim(opt, res['hom_g2'], res['hom_e2'], res['ct_overall_g_var'], res['ct_overall_e_var'], 
            res['fixed_vars']):
        out, opt = util.re_optim(out, opt, loglike_fun, par, args, method, nrep)
        res = extract( out, model, Y, K, P, ctnu, fixed_covars, random_covars )

    res['opt'] = opt
    #res['out'] = out

    return res

def hom_REML_loglike(par: list, y: np.ndarray, K: np.ndarray, X: np.ndarray,
                     ctnu: np.ndarray, random_covars_ZZT: dict, shared: bool) -> float:
    """
    Loglikelihood for REML under Hom model
    """

    N, C = ctnu.shape
    hom_g2, hom_e2 = par[:2]
    V = np.zeros((C,C))
    W = np.zeros((C,C))
    r2 = {}
    for i, key in enumerate(sorted(random_covars_ZZT.keys())):
        if shared:
            r2[key] = par[2 + i]
        else:
            r2[key] = par[(2 + i * C):(2 + i * (C + 1))]

    l = LL(y, K, X, ctnu, random_covars_ZZT, hom_g2, hom_e2, V, W, r2)
    return l

def hom_REML(Y: np.ndarray, K: np.ndarray, P: np.ndarray, ctnu: np.ndarray,
             fixed_covars: dict={}, random_covars: dict={}, shared: bool=True,
             par: list=None, method: str=None, nrep: int=10) -> dict:
    """
    Fit Hom model with REML

    Parameters:
        Y:  cell type-specific pseudobulk
        K:  kinship matrix
        P:  cell type proportioons
        ctnu:   cell type-specific noise variance
        fixed_covars:   design matrix for Extra fixed effects
        random_covars:   design matrix for Extra random effects
        shared: whether Extra fixed and random effects are shared across cell types
        par:    initial parameters
        method: optimization method, e.g. BFGS
        nrep:   number of optimization repeats when initial optimization failed
    Returns:
        a dict of model parameters and statistics
    """

    log.logger.info('Fitting Hom model with REML')

    N, C = Y.shape
    y = Y.flatten()
    n_random = len(random_covars.keys()) if shared else len(random_covars.keys()) * C
    X = util.get_X( fixed_covars, N, C, fixed_shared=shared )

    if par is None:
        beta = linalg.inv( X.T @ X ) @ ( X.T @ y )
        hom_g2 = np.var(y - X @ beta) / (2 + len(random_covars.keys()))
        par = [hom_g2] * (2 + n_random)

    res = _reml( 'hom', par, Y, K, P, ctnu, fixed_covars,
                 random_covars, shared, method, nrep=nrep )

    return res

def freeW_REML_loglike(par:list, y:np.ndarray, K:np.ndarray, X:np.ndarray, ctnu:np.ndarray,
                       random_covars_ZZT: dict, shared: bool) -> float:
    """
    Loglikelihood function for REML under FreeW model, where env is Free and genetic is Hom
    """
    N, C = ctnu.shape
    hom_g2 = par[0]
    hom_e2 = par[1]
    W = np.diag(par[2:(C+2)])
    V = np.zeros((C,C))
    r2 = {}
    for i, key in enumerate(sorted(random_covars_ZZT.keys())):
        if shared:
            r2[key] = par[C+2+i]
        else:
            r2[key] = par[(C + 2 + i * C):(C + 2 + (i+1) * C)]

    l = LL(y, K, X, ctnu, random_covars_ZZT, hom_g2, hom_e2, V, W, r2)
    return l

def freeW_REML(Y:np.ndarray, K:np.ndarray, P:np.ndarray, ctnu:np.ndarray,
               fixed_covars:dict={}, random_covars:dict={}, shared: bool = True,
               par:list=None, method:str=None, nrep:int=10) -> dict:
    """
    Fit FreeW model using REML, where env is Free and genetic is Hom

    Parameters:
        Y:  cell type-specific pseudobulk
        K:  kinship matrix
        P:  cell type proportions
        ctnu:   cell type-specific noise variance
        fixed_covars:   design matrices for Extra fixed effects
        random_covars:   design matrices for Extra random effects
        shared: whether Extra fixed and random effects are shared across cell types
        par:    initial parameters
        method: optimization method, e.g. BFGS
        nrep:   number of optimization repeats when initial optimization failed
    Returns:
        dict of model parameters and statistics
    """

    log.logger.info('Fitting FreeW model with REML')

    N, C = Y.shape
    y = Y.flatten()
    X = util.get_X( fixed_covars, N, C, fixed_shared=shared )
    n_random = len(random_covars.keys()) if shared else len(random_covars.keys()) * C
    n_par = 2 + 2 * C + X.shape[1] + n_random

    if par is None:
        beta = linalg.inv( X.T @ X ) @ ( X.T @ y )
        hom_g2 = np.var(y - X @ beta) / (3 + len(random_covars.keys()))
        par = [hom_g2] * (2 + 2*C + n_random) # TODO: do I need to stardardize design matrix for extra random effect? Need to fix it later

    res = _reml('freeW', par, Y, K, P, ctnu, fixed_covars,
                random_covars, shared, method, nrep=nrep)

    return res

def free_REML_loglike(par:list, y:np.ndarray, K:np.ndarray, X:np.ndarray, ctnu:np.ndarray,
                      random_covars_ZZT: dict, shared: bool) -> float:
    """
    Loglikelihood function for REML under Free model
    """
    N, C = ctnu.shape
    hom_g2 = par[0]
    hom_e2 = par[1]
    V = np.diag(par[2:(C+2)])
    W = np.diag(par[(C+2):(2*C+2)])
    r2 = {}
    for i, key in enumerate(sorted(random_covars_ZZT.keys())):
        if shared:
            r2[key] = par[2*C+2+i]
        else:
            r2[key] = par[(2*C+2+i*C):(2*C+2+(i+1)*C)]

    l = LL(y, K, X, ctnu, random_covars_ZZT, hom_g2, hom_e2, V, W, r2)
    return l

def free_REML(Y:np.ndarray, K:np.ndarray, P:np.ndarray, ctnu:np.ndarray,
              fixed_covars:dict={}, random_covars:dict={}, shared: bool=True,
              par:list=None, method:str=None, nrep:int=10, jk:bool=True
              ) -> Tuple[dict,dict]:
    '''
    Fit Free model using REML

    Parameters:
        Y:  cell type-specific pseudobulk
        K:  kinship matrix
        P:  cell type proportions
        ctnu:   cell type-specific noise variance
        fixed_covars:   design matrices for Extra fixed effects
        random_covars:   design matrices for Extra random effects
        shared: whether Extra fixed and random effects are shared across cell types
        par:    initial parameters
        method: optimization method, e.g. BFGS
        nrep:   number of optimization repeats when initial optimization failed
        jk: jackknife to estimation dispersion matirx
    Returns:
        a tuple of 
            #.  dict of model parameters and statistics
            #.  dict of p values
    '''

    log.logger.info('Fitting Free model with REML')

    N, C = Y.shape
    y = Y.flatten()
    X = util.get_X( fixed_covars, N, C, fixed_shared=shared )
    n_random = len(random_covars.keys()) if shared else len(random_covars.keys()) * C
    n_par = 2 + 2 * C + X.shape[1] + n_random

    if par is None:
        beta = linalg.inv( X.T @ X ) @ ( X.T @ y )
        hom_g2 = np.var(y - X @ beta) / (4 + len(random_covars.keys()) )
        par = [hom_g2] * (2 + 2*C + n_random)

    res = _reml('free', par, Y, K, P, ctnu, fixed_covars,
                random_covars, shared, method, nrep=nrep)

    if jk:
        jacks = {'ct_beta':[], 'V':[], 'W':[], 'VW':[]}
        for i in range(N):
            Y_jk, K_jk, ctnu_jk, fixed_covars_jk, random_covars_jk, P_jk = util.jk_rmInd(
                i, Y, K, ctnu, fixed_covars, random_covars, P)

            res_jk = _reml('free', par, Y_jk, K_jk, P_jk, ctnu_jk,
                    fixed_covars_jk, random_covars_jk, shared, method, nrep=nrep)

            jacks['ct_beta'].append( res_jk['beta']['ct_beta'] )
            jacks['V'].append( np.diag(res_jk['V']) )
            jacks['W'].append( np.diag(res_jk['W']) )
            jacks['VW'].append( np.append( np.diag(res_jk['V']), np.diag(res_jk['W']) ) )

        var_V = (N-1) * np.cov( np.array(jacks['V']).T, bias=True )
        var_W = (N-1) * np.cov( np.array(jacks['W']).T, bias=True )
        var_VW = (N-1) * np.cov( np.array(jacks['VW']).T, bias=True )
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

    return res, p

def full_REML_loglike(par:list, y:np.ndarray, K:np.ndarray, X:np.ndarray, ctnu:np.ndarray,
                      random_covars_ZZT:dict, shared:bool) -> float:
    """
    Loglikelihood function for REML under Full model
    """
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
    r2 = {}
    for i, key in enumerate(sorted(random_covars_ZZT.keys())):
        if shared:
            r2[key] = par[2*ngam+i]
        else:
            r2[key] = par[(2*ngam+i*C):(2*ngam+(i+1)*C)]

    l = LL(y, K, X, ctnu, random_covars_ZZT, hom_g2, hom_e2, V, W, r2)
    return l

def full_REML(Y:np.ndarray, K:np.ndarray, P:np.ndarray, ctnu:np.ndarray,
              fixed_covars:dict={}, random_covars:dict={}, shared:bool=True,
              par:list=None, method:str=None, nrep:int=10) -> dict:
    """
    Fit Full model using REML

    Parameters:
        Y:  cell type-specific pseudobulk
        K:  kinship matrix
        P:  cell type proportions
        ctnu:   cell type-specific noise variance
        fixed_covars:   design matrices for Extra fixed effects
        random_covars:   design matrices for Extra random effects
        shared: whether Extra fixed and random effects are shared across cell types
        par:    initial parameters
        method: optimization method, e.g. BFGS
        nrep:   number of optimization repeats when initial optimization failed
    Returns:
        dict of model parameters and statistics
    """

    log.logger.info('Fitting Full model with REML')

    N, C = Y.shape
    ngam = C*(C+1) // 2
    y = Y.flatten()
    X = util.get_X( fixed_covars, N, C, fixed_shared=shared )
    n_random = len(random_covars.keys()) if shared else len(random_covars.keys()) * C

    if par is None:
        beta = linalg.inv( X.T @ X ) @ ( X.T @ y )
        hom_g2 = np.var(y - X @ beta) / 2
        V = W = np.diag(np.ones(C))[np.tril_indices(C)] * hom_g2
        par = list(V) + list(W) + [hom_g2] * n_random

    res = _reml('full', par, Y, K, P, ctnu, fixed_covars, random_covars, shared,
                method, nrep=nrep)

    return res

def _free_he(Y: np.ndarray, K: np.ndarray, Kt: Optional[np.ndarray], 
             ctnu: np.ndarray, P: np.ndarray, fixed_covars: dict={}, 
             random_covars: dict={}, fixed_shared: bool=False, random_shared: bool=False, 
             output_beta: bool=True, dtype:str = None) -> dict:

    N, C = Y.shape
    X = util.get_X(fixed_covars, N, C, fixed_shared=fixed_shared)

    theta = he_ols(Y, K, X, ctnu, 'free', random_covars, Kt=Kt, random_shared=random_shared, dtype=dtype)
    if Kt is None:
        hom_g2, hom_e2 = theta[0], theta[1]
        V, W = np.diag(theta[2:(C + 2)]), np.diag(theta[(C + 2):(C * 2 + 2)])
    else:
        hom_g2, hom_gt2, hom_e2 = theta[0], theta[1], theta[2]
        V, Vt, W = np.diag(theta[3:(C + 3)]), np.diag(theta[(C + 3):(C * 2 + 3)]), np.diag(theta[(C * 2 + 3):(C * 3 + 3)])

    # ct specific effect variance
    ct_overall_g_var, ct_specific_g_var = util.ct_random_var(V, P)
    ct_overall_e_var, ct_specific_e_var = util.ct_random_var(W, P)

    out = {'hom_g2':hom_g2, 'hom_e2':hom_e2, 'V':V, 'W':W, 
            'ct_overall_g_var':ct_overall_g_var, 'ct_specific_g_var':ct_specific_g_var, 
            'ct_overall_e_var':ct_overall_e_var, 'ct_specific_e_var':ct_specific_e_var}
    r2 = _get_r2(theta[(C*2+2):], random_covars, C)
    out['r2'] = r2

    # trans effect
    if Kt is not None:
        out['hom_gt2'] = hom_gt2
        out['Vt'] = Vt

        ct_overall_gt_var, ct_specific_gt_var = util.ct_random_var(Vt, P)
        out['ct_overall_gt_var'] = ct_overall_gt_var
        out['ct_specific_gt_var'] = ct_specific_gt_var

    # OLS to get beta
    if output_beta:
        beta = sla.inv(X.T @ X) @ X.T @ Y.flatten()
        beta, fixed_vars, random_vars = util.cal_variance(beta, P, fixed_covars, r2, random_covars)
        
        out['beta'] = beta
        out['fixed_vars'] = fixed_vars
        out['random_vars'] = random_vars

    return out


def free_HE(Y: np.ndarray, K: np.ndarray, ctnu: np.ndarray, P: np.ndarray, fixed_covars: dict={},
        random_covars: dict={}, Kt:Optional[np.ndarray] = None, shared:bool = True, fixed_shared:bool=True, random_shared:bool=True,
        jk:bool = True, dtype:str = None) -> Tuple[dict, dict]:
    """
    Fitting Free model with HE

    Parameters:
        Y:    cell type-specific pseudobulk (no header no index)
        K:    kinship matrix
        ctnu: cell type-specific noise variance (no header no index)
        P:    cell type proportions
        fixed_covars:   design matrices for Extra fixed effects
        random_covars:  design matrices for Extra random effects
        Kt:   kinship matrix for trans effect
        shared: whether Extra fixed and random effects are shared across cell types
        fixed_shared: whether Extra fixed effects are shared across cell types
        random_shared: whether Extra random effects are shared across cell types
        jk: perform jackknife
        dtype:  data type for he_ols, e.g. float32

    Returns:
        a tuple of
            #.  dictionary of parameter estimates
            #.  dictionary of p values
    """

    log.logger.info('Fitting Free model with HE')

    if shared:
        random_shared = True
        fixed_shared = True

    N, C = Y.shape 
    X = util.get_X(fixed_covars, N, C, fixed_shared=shared)
    n_random = len(random_covars.keys()) if shared else len(random_covars.keys()) * C
    n_par = 2 + 2 * C + X.shape[1] + n_random  # TODO: don't include trans


    out = _free_he(Y, K, Kt, ctnu, P, fixed_covars, random_covars=random_covars,
                   fixed_shared=fixed_shared, random_shared=random_shared, 
                   output_beta=False, dtype=dtype)
    out['nu'] = ( ctnu * (P ** 2) ).sum(axis=1)
    log.logger.info(out['nu'].dtype)

    # jackknife
    if jk:
        log.logger.info('Jackknife')
        #jacks = {'ct_beta':[], 'V':[], 'W':[], 'VW':[]}
        jacks = {'hom_g2': [], 'V': [], 'hom_e2': [], 'W': []}
        for i in range(N):
            if i % 10 == 0:
                log.logger.info(i)

            if Kt is None:
                Y_jk, K_jk, ctnu_jk, fixed_covars_jk, random_covars_jk, P_jk = util.jk_rmInd(
                                                i, Y, K, ctnu, fixed_covars, random_covars, P=P)

                out_jk = _free_he(Y_jk, K_jk, None, ctnu_jk, P_jk, fixed_covars_jk,
                               random_covars_jk, fixed_shared, random_shared, output_beta=False, dtype=dtype)

            else:
                Y_jk, K_jk, ctnu_jk, fixed_covars_jk, random_covars_jk, P_jk, Kt_jk = util.jk_rmInd(
                                                i, Y, K, ctnu, fixed_covars, random_covars, P=P, Kt=Kt)

                out_jk = _free_he(Y_jk, K_jk, Kt_jk, ctnu_jk, P_jk, fixed_covars_jk,
                               random_covars_jk, fixed_shared, random_shared, output_beta=False, dtype=dtype)
                
                if 'Vt' not in jacks.keys():
                    jacks['Vt'] = []
                jacks['Vt'].append(np.diag(out_jk['Vt']))

            #jacks['ct_beta'].append(out_jk['beta']['ct_beta'])
            jacks['hom_g2'].append(out_jk['hom_g2'])
            jacks['V'].append(np.diag(out_jk['V']))
            jacks['hom_e2'].append(out_jk['hom_e2'])
            jacks['W'].append(np.diag(out_jk['W']))
            # jacks['VW'].append(np.append( np.diag(out_jk['V']), np.diag(out_jk['W'])))


        var_hom_g2 = (N - 1) * np.var(jacks['hom_g2'])
        var_V = (N - 1) * np.cov(np.array(jacks['V']).T, bias=True)
        var_hom_e2 = (N - 1) * np.var(jacks['hom_e2'])
        var_W = (N - 1) * np.cov(np.array(jacks['W']).T, bias=True)
        # var_VW = (N-1) * np.cov( np.array(jacks['VW']).T, bias=True )
        #var_ct_beta = (N-1) * np.cov( np.array(jacks['ct_beta']).T, bias=True )

        p = {   
                'hom_g2': wald.wald_test(out['hom_g2'], 0, var_hom_g2, N - n_par),
                'V': wald.mvwald_test(np.diag(out['V']), np.zeros(C), var_V, n=N, P=n_par),
                'hom_e2': wald.wald_test(out['hom_e2'], 0, var_hom_e2, N - n_par),
                'W': wald.mvwald_test(np.diag(out['W']), np.zeros(C), var_W, n=N, P=n_par),
                'var_hom_g2': var_hom_g2,
                'var_V': var_V,
                'var_hom_e2': var_hom_e2,
                'var_W': var_W,
                # 'VW': wald.mvwald_test( np.append(np.diag(out['V']),np.diag(out['W'])), np.zeros(2*C), 
                    # var_VW, n=N, P=n_par),
                #'ct_beta': util.wald_ct_beta( out['beta']['ct_beta'], var_ct_beta, n=N, P=n_par )
                }

        if Kt is not None:
            var_Vt = (N-1) * np.cov( np.array(jacks['Vt']).T, bias=True )
            p['Vt'] = wald.mvwald_test(np.diag(out['Vt']), np.zeros(C), var_Vt, n=N, P=n_par)
        
        # OLS to test fixed effect
        y = Y.flatten()
        beta = sla.inv(X.T @ X) @ X.T @ y
        epsilon = y - X @ beta
        s2 = (epsilon.T @ epsilon) / (N * C - X.shape[1]) # df
        var_beta = s2 * sla.inv(X.T @ X)
        p['ct_beta'] = util.wald_ct_beta(beta[:C], var_beta[:C, :C], n=N, P=n_par)  # df
    else:
        p = {}
    return out, p


def full_HE(Y: np.ndarray, K: np.ndarray, ctnu: np.ndarray, P: np.ndarray, fixed_covars: dict={},
        random_covars: dict={}, shared:bool=True, dtype:str=None) -> dict:
    """
    Fitting Full model with HE

    Parameters:
        Y:    cell type-specific pseudobulk (no header no index)
        K:    kinship matrix
        ctnu: cell type-specific noise variance (no header no index)
        P:    cell type proportions
        fixed_covars:   design matrices for Extra fixed effects
        random_covars:  design matrices for Extra random covars
        shared: whether Extra fixed and random effects are shared across cell types
        dtype:  data type for computation in he_ols, e.g. float32

    Returns:
        a dictionary of parameter estimates
    """

    log.logger.info('Fitting Full model with HE')

    N, C = Y.shape 
    ntril = (C-1) * C // 2
    X = util.get_X(fixed_covars, N, C, fixed_shared=shared)

    theta = he_ols(Y, K, X, ctnu, 'full', random_covars=random_covars,
                   random_shared=shared, dtype=dtype)
    log.logger.info(theta.dtype)
    V, W = np.diag(theta[:C]), np.diag(theta[C:(C*2)])
    V[np.triu_indices(C,k=1)] = theta[(C*2):(C*2 + ntril)]
    V = V + V.T - np.diag(theta[:C])
    W[np.triu_indices(C,k=1)] = theta[(C*2 + ntril):(C*2 + ntril*2)]
    W = W + W.T - np.diag(theta[C:(C*2)])

    ct_overall_g_var, ct_specific_g_var = util.ct_random_var( V, P )
    ct_overall_e_var, ct_specific_e_var = util.ct_random_var( W, P )

    he = {'V': V, 'W': W, 'ct_overall_g_var':ct_overall_g_var, 'ct_overall_e_var':ct_overall_e_var}
    return he

