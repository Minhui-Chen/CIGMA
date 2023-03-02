import numpy as np, pandas as pd
from scipy import linalg, optimize
import util, wald

def get_X(fixed_covars, N, C):
    X = np.kron( np.ones((N,1)), np.eye(C) )
    for key in np.sort(list(fixed_covars.keys())):
        m = fixed_covars[key]
        if len( m.shape ) == 1:
            m = m.reshape(-1,1)
        X = np.concatenate( ( X, np.repeat(m, C, axis=0)), axis=1 )
    return( X )

def SIGMA( Z, ctvs, hom_g2=0, hom_e2=0, 
        V=None, W=None ):
    '''
    Variance matrix of y
    '''
    N, C = ctvs.shape
    if V is None:
        V = np.matrix(np.zeros((C,C)))
    if W is None:
        W = np.matrix(np.zeros((C,C)))
    sigma_y2 = np.kron(Z, hom_g2+V)
    sigma_y2 = sigma_y2 + np.kron( np.eye(N), hom_e2+W )
    sigma_y2 = sigma_y2 + np.diag( ctvs.flatten() )

    return sigma_y2

def LL(y, Z, X, ctvs, hom_g2, hom_e2, V, W):
    '''
    Loglikelihood function
    '''
    N, C = ctvs.shape
    sigma_y2 = SIGMA( Z, ctvs, hom_g2, hom_e2, V, W)

    # inverse variance
    w, v = linalg.eigh(sigma_y2)
    if ( np.amax(w)/np.amin(w) ) > 1e8 or np.amin(w) < 0:
        return(1e12)

    # calculate B matrix
    m1 = X.T @ v @ np.diag(1/w) @ v.T
    m2 = m1 @ X

    # calculate loglikelihood
    L = np.sum( np.log(w) )
    L = L + np.linalg.slogdet(m2)[1] # do i need to take care sign of determinant
    yPy = y @ v @ np.diag(1/w) @ v.T @ y - y @ m1.T @ linalg.inv(m2) @ m1 @ y
    L = L + yPy
    L = 0.5 * L
    return( L )

def extract( out, model, cty, Z, P, ctvs, fixed_covars ):
    N, C = cty.shape
    ngam = C*(C+1) // 2

    if model == 'hom':
        hom_g2, hom_e2 = out['x']
        V = W = None
        ct_overall_g2 = ct_overall_e2 = 0
    elif model == 'free':
        hom_g2 = out['x'][0]
        hom_e2 = out['x'][1]
        V = np.diag( out['x'][2:(2+C)] )
        W = np.diag( out['x'][(2+C):(2+2*C)] )
        ct_overall_g2 = util.ct_randomeffect_variance( V, P )
        ct_overall_e2 = util.ct_randomeffect_variance( W, P )
    elif model == 'full':
        V = np.zeros((C,C))
        V[np.tril_indices(C)] = out['x'][:ngam]
        V = V + V.T
        W = np.zeros((C,C))
        W[np.tril_indices(C)] = out['x'][ngam:(2*ngam)]
        W = W + W.T

    # beta
    y = cty.flatten()
    X = get_X( fixed_covars, N, C )

    Vy = SIGMA( Z, ctvs, hom_g2, hom_e2, V, W)
    beta = util.glse( Vy, X, y )
    beta, fixed_vars = util.fixedeffect_vars( beta, P, fixed_covars )

    return( hom_g2, hom_e2, V, W, beta, ct_overall_g2, ct_overall_e2, fixed_vars )

def check_optim(opt, hom_g2, hom_e2, ct_overall_g2, ct_overall_e2, fixed_vars, cut=5):
    if ( ( opt['l'] < -1e10 ) or ( not opt['success'] ) or 
            ( hom_g2 > cut ) or ( hom_e2 > cut ) or
            ( ct_overall_g2 > cut ) or ( ct_overall_e2 > cut ) or
            np.any(np.array(list(fixed_vars.values())) > cut) ):
        return True
    else:
        return False

def reml_f(loglike_fun, par, model, cty, Z, P, ctvs, fixed_covars, method, nrep):
    N, C = cty.shape
    y = cty.flatten()
    X = get_X( fixed_covars, N, C )

    args = (y, Z, X, ctvs)

    out, opt = util.optim( loglike_fun, par, args, method )
    hom_g2, hom_e2, V, W, beta, ct_overall_g2, ct_overall_e2, fixed_vars = extract( 
            out, model, cty, Z, P, ctvs, fixed_covars )

    if check_optim(opt, hom_g2, hom_e2, ct_overall_g2, ct_overall_e2, fixed_vars):
        out, opt = util.re_optim(out, opt, loglike_fun, par, args, method, nrep)
        hom_g2, hom_e2, V, W, beta, ct_overall_g2, ct_overall_e2, fixed_vars = extract( 
                out, model, cty, Z, P, ctvs, fixed_covars )
    return(opt, beta, hom_g2, hom_e2, V, W, ct_overall_g2, ct_overall_e2, fixed_vars)

def hom_REML_loglike(par, y, Z, X, ctvs):
    N, C = ctvs.shape
    hom_g2, hom_e2 = par
    V = np.zeros((C,C))
    W = np.zeros((C,C))

    l = LL(y, Z, X, ctvs, hom_g2, hom_e2, V, W)
    return( l )

def hom_REML(cty, Z, P, ctvs, fixed_covars={}, par=None, method=None, nrep=10):
    N, C = cty.shape
    y = cty.flatten()
    X = get_X( fixed_covars, N, C )

    if par is None:
        beta = linalg.inv( X.T @ X ) @ ( X.T @ y )
        hom_g2 = np.var(y - X @ beta) / 2
        par = [hom_g2] * 2
    
    opt, beta, hom_g2, hom_e2 = reml_f(
            screml_hom_loglike, par, 'hom', cty, Z, P, ctvs, fixed_covars, method, nrep=10)[:4]

    return({'hom_g2':hom_g2, 'hom_e2':hom_e2, 'beta':beta, 'opt':opt})

def free_REML_loglike(par, y, Z, X, ctvs):
    N, C = ctvs.shape
    hom_g2 = par[0]
    hom_e2 = par[1]
    V = np.diag(par[2:(C+2)])
    W = np.diag(par[(C+2):(2*C+2)])

    l = LL(y, Z, X, ctvs, hom_g2, hom_e2, V, W)
    return( l )

def free_REML(cty, Z, P, ctvs, fixed_covars={}, par=None, method=None, nrep=10, jk=True):
    N, C = cty.shape
    y = cty.flatten()
    X = get_X( fixed_covars, N, C )

    if par is None:
        beta = linalg.inv( X.T @ X ) @ ( X.T @ y )
        hom_g2 = np.var(y - X @ beta) / 4
        par = [hom_g2] * (2 + 2*C)

    opt, beta, hom_g2, hom_e2, V, W = reml_f(
            screml_free_loglike, par, 'free', cty, Z, P, ctvs, fixed_covars, method, nrep=10)[:6]
    results = {'hom_g2':hom_g2, 'hom_e2':hom_e2, 'V':V, 'W':W, 'beta':beta, 'opt':opt}
    
    if jk:
        jacks = {'ct_beta':[], 'hom_g2':[], 'hom_e2':[], 'V':[], 'W':[]}
        for i in range(N):
            cty_jk, ctvs_jk, fixed_covars_jk, _, P_jk = util.jk_rmInd(
                    i, cty, ctvs, fixed_covars, {}, P)
            Z_jk = np.delete(np.delete(Z, i, axis=0), i, axis=1)

            _, beta_jk, hom_g2_jk, hom_e2_jk, V_jk, W_jk = reml_f(
                    screml_free_loglike, par, 'free', cty_jk, Z_jk, 
                    P_jk, ctvs_jk, fixed_covars_jk, method, nrep=10)[:6]

            jacks['ct_beta'].append( beta_jk['ct_beta'] )
            jacks['hom_g2'].append( hom_g2_jk )
            jacks['hom_e2'].append( hom_e2_jk )
            jacks['V'].append( np.diag(V_jk) )
            jacks['W'].append( np.diag(W_jk) )

        var_hom_g2 = (N-1) * np.var( jacks['hom_g2'] )
        var_hom_e2 = (N-1) * np.var( jacks['hom_e2'] )
        var_ct_beta = (N-1) * np.cov( np.array(jacks['ct_beta']).T, bias=True )
        var_V = (N-1) * np.cov( np.array(jacks['V']).T, bias=True )
        var_W = (N-1) * np.cov( np.array(jacks['W']).T, bias=True )

        n_par = 2 + 2 * C + X.shape[1]
        p = {
                'hom_g2': wald.wald_test(hom_g2, 0, var_hom_g2, N-n_par),
                'hom_e2': wald.wald_test(hom_e2, 0, var_hom_e2, N-n_par),
                'ct_beta': util.wald_ct_beta(beta['ct_beta'], var_ct_beta, n=N, P=n_par),
                'V': wald.mvwald_test(np.diag(V), np.zeros(C), var_V, n=N, P=n_par)
                'W': wald.mvwald_test(np.diag(W), np.zeros(C), var_W, n=N, P=n_par)
                }
    else:
        p = {}

    return( results, p )


def full_REML_loglike(par, y, Z, X, ctvs):
    N, C = ctvs.shape
    ngam = C*(C+1) // 2
    hom_g2 = 0
    hom_e2 = 0
    V = np.zeros((C,C))
    V[np.tril_indices(C)] = par[:ngam]
    V = V + V.T
    W = np.zeros((C,C))
    W[np.tril_indices(C)] = par[ngam:(2*ngam)]
    W = W + W.T

    l = LL(y, Z, X, ctvs, hom_g2, hom_e2, V, W)
    return( l )

def full_REML(cty, Z, ctvs, fixed_covars={}, par=None, method=None, nrep=10):
    N, C = cty.shape
    ngam = C*(C+1) // 2
    y = cty.flatten()
    X = get_X( fixed_covars, N, C )

    if par is None:
        beta = linalg.inv( X.T @ X ) @ ( X.T @ y )
        hom_g2 = np.var(y - X @ beta) / 4
        V = W = np.diag(np.ones(C))[np.tril_indices(C)] * hom_g2
        par = list(V) + list(W)

    opt, beta, _, _, V, W = reml_f(
            screml_full_loglike, par, 'full', cty, Z, P, ctvs, fixed_covars, method, nrep=10)[:6]
    results = {'V':V, 'W':W, 'beta':beta, 'opt':opt}

    return( results )
