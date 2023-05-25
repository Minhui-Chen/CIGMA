from typing import Tuple, Optional, Union

import os, tempfile, sys, subprocess, re
import numpy as np, pandas as pd
import numpy.typing as npt
import rpy2.robjects as ro
from rpy2.robjects import r, pandas2ri, numpy2ri
from rpy2.robjects.conversion import localconverter
from scipy import stats, linalg, optimize
from numpy.random import default_rng

from ctmm import wald
from . import log

def read_covars(fixed_covars: dict = {}, random_covars: dict = {}, C: Optional[int] = None) -> tuple:
    '''
    Read fixed and random effect design matrices

    Parameters:
        fixed_covars:   files of design matrices for fixed effects,
                        except for cell type-specifc fixed effect, without header
        random_covars:  files of design matrices for random effects,
                        except for cell type-shared and -specifc random effect, without header
        C:  number of cell types
    Returns:
        a tuple of
            #. dict of design matrices for fixed effects
            #. dict of design matrices for random effects
            #. others
    '''
    def read(covars):
        tmp = {}
        for key in covars.keys():
            f = covars[key]
            if isinstance( f, str ):
                tmp[key] = np.loadtxt( f )
            else:
                tmp[key] = f
        return( tmp )

    fixed_covars = read(fixed_covars)
    random_covars = read(random_covars)
    n_fixed, n_random = len( fixed_covars.keys() ), len( random_covars.keys() )
    random_keys = list( np.sort( list(random_covars.keys()) ) )
    Rs = [random_covars[key] for key in random_keys]
    if C:
        random_MMT = [np.repeat( np.repeat(R @ R.T, C, axis=0), C, axis=1 ) for R in Rs]
    else:
        random_MMT = [R @ R.T for R in Rs]

    return fixed_covars, random_covars, n_fixed, n_random, random_keys, Rs, random_MMT

def age_group(age: np.ndarray):
    """
    Separate age groups
    """
    bins = np.arange(25, 91, 5)
    return np.digitize(age, bins)

def optim(fun: callable, par: list, args: tuple, method: str) -> Tuple[object, dict]:
    '''
    Optimization use scipy.optimize.minimize

    Parameters:
        fun:    objective function to minimize (e.g. log-likelihood function)
        par:    initial parameters
        args:   extra arguments passed to objective function
        method: optimization method, e.g. BFGS
    Returns:
        a tuple of
            #. OptimizeResult object from optimize.minimize
            #. dict of optimization parameters and results
    '''
    if method is None:
        method = 'BFGS'
    if method == 'BFGS-Nelder':
        out1 = optimize.minimize( fun, par, args=args, method='BFGS' )
        out = optimize.minimize( fun, out1['x'], args=args, method='Nelder-Mead' )
        opt = {'method1':'BFGS', 'success1':out1['success'], 'status1':out1['status'],
                'message1':out1['message'], 'l1':out1['fun'] * (-1),
                'method':'Nelder-Mead', 'success':out['success'], 'status':out['status'],
                'message':out['message'], 'l':out['fun'] * (-1)}
    else:
        out = optimize.minimize( fun, par, args=args, method=method )
        opt = {'method':method, 'success':out['success'], 'status':out['status'],
                'message':out['message'], 'l':out['fun'] * (-1)}
    return( out, opt )

def check_optim(opt: dict, hom_g2: float, hom_e2: float, ct_overall_g_var: float, ct_overall_e_var: float, 
        fixed_vars: dict, random_vars: dict={}, cut: float=5) -> bool:
    '''
    Check whether optimization converged successfully

    Parameters:
        opt:    dict of optimization results, e.g. log-likelihood
        hom_g2: variance of genetic effect shared across cell types
        hom_e2: variance of env effect shared across cell types
        ct_overall_g_var:  overall variance explained by cell type-specific genetic effect
        ct_overall_e_var:  overall variance explained by cell type-specific env effect
        fixed_var:  dict of variances explained by each fixed effect feature, including cell type-specific fixed effect
        random_var: dict of variances explained by each random effect feature, doesn't include cell type-shared or -specific effect
        cut:    threshold for large variance
    Returns:
        True:   optim failed to converge
        False:  optim successfully to converge
    '''
    if ( (opt['l'] < -1e10) or (not opt['success']) or (hom_g2 > cut) or (hom_e2 > cut) or 
            (ct_overall_g_var > cut) or (ct_overall_e_var > cut) or
            np.any(np.array(list(fixed_vars.values())) > cut) or
            np.any(np.array(list(random_vars.values())) > cut) ):
        return True
    else:
        return False

def re_optim(out: object, opt: dict, fun: callable, par: list, args: tuple, method: str, nrep: int=10) -> Tuple[object, dict]:
    '''
    Rerun optimization

    Parameters:
        out:    OptimizeResult object
        opt:    opmization results, e.g. method used, log-likelihood
        fun:    objective function to minimize
        par:    initial parameters used in the first try of optimization
        args:   extra argument passed to the objective function
        method: optimization method, e.g. BFGS
        nrep:   number of optimization repeats
    Returns:
        a tuple of 
            #. OptimizeResult of the best optimization
            #. results of the best optimization
    '''
    rng = default_rng()
    #print( out['fun'] )
    for i in range(nrep):
        par_ = np.array(par) * rng.gamma(2,1/2,len(par))
        out_, opt_ = optim(fun, par_, args=args, method=method)
        log.logger.info( f"-loglike: {out_['fun']}" )
        if (not out['success']) and out_['success']:
            out, opt = out_, opt_
        elif (out['success'] == out_['success']) and (out['fun'] > out_['fun']):
            out, opt = out_, opt_
    #print( out['fun'] )
    return( out, opt )

def dict2Rlist( X: dict ) -> object:
    '''
    Transform a python dictionary to R list

    Parameters:
        X:  python dictionary
    Returns:
        R list
    '''
    if len( X.keys() ) == 0:
        return( r('NULL') )
    else:
        keys = np.sort( list(X.keys()) )
        rlist = ro.ListVector.from_length( len(keys) )
        for i in range( len(keys) ):
            if isinstance(X[keys[i]], str):
                if os.path.exists( X[keys[i]] ):
                    rlist[i] = r['as.matrix']( r['read.table'](X[keys[i]]) )
                else:
                    try:
                        rlist[i] = np.array( [X[keys[i]]] )
                    except:
                        numpy2ri.activate()
                        rlist[i] = np.array( [X[keys[i]]] )
                        numpy2ri.deactivate() # deactivate would cause numpy2ri deactivated in calling fun
            elif isinstance(X[keys[i]], pd.DataFrame):
                with localconverter(ro.default_converter + pandas2ri.converter):
                    rlist[i] = r['as.matrix']( X[keys[i]] )
            elif isinstance(X[keys[i]], np.ndarray):
                try:
                    rlist[i] = r['as.matrix']( X[keys[i]] )
                except:
                    numpy2ri.activate()
                    rlist[i] = r['as.matrix']( X[keys[i]] )
                    numpy2ri.deactivate()
            elif isinstance(X[keys[i]], int) or isinstance(X[keys[i]], float):
                try:
                    rlist[i] = np.array( [X[keys[i]]] )
                except:
                    numpy2ri.activate()
                    rlist[i] = np.array( [X[keys[i]]] )
                    numpy2ri.deactivate()
        return( rlist )

def generate_HE_initial(he: dict, ML: bool=False, REML: bool=False) -> list:
    '''
    Convert HE estimates to initial parameter for ML / REML

    Parameters:
        he: estiamtes from HE
        ML: generate initial parameters for ML 
        REML:   generate initial parameters for REML
    Returns:
        initial parameters for ML / REML
    '''
    initials_random_effects = []
    # homogeneous effect
    if 'hom2' in he.keys():
        initials_random_effects.append( he['hom2'] )
    # heterogeneous effect
    if 'V' in he.keys():
        C = he['V'].shape[0]
        # determine model based on V
        if np.any( np.diag(np.diag(he['V'])) != he['V'] ):
            # Full model
            initials_random_effects = initials_random_effects + list(he['V'][np.triu_indices(C)])
        elif len( np.unique( np.diag(he['V']) ) ) == 1:
            # IID model
            initials_random_effects.append( he['V'][0,0] )
        else:
            # Free model
            initials_random_effects = initials_random_effects + list(np.diag(he['V']))
    # other random covariates
    if 'r2' in he.keys():
        for key in np.sort( list(he['r2'].keys()) ):
            initials_random_effects.append( he['r2'][key] )

    if REML is True:
        return( initials_random_effects )

    initials_fixed_effects = list(he['beta']['ct_beta'])
    for key in np.sort( list(he['beta'].keys()) ):
        if key != 'ct_beta':
            initials_fixed_effects = initials_fixed_effects + list( he['beta'][key] )

    if ML is True:
        return( initials_fixed_effects + initials_random_effects )

def glse( sig2s: np.ndarray, X: np.ndarray, y: np.ndarray, inverse: bool=False ) -> np.ndarray:
    '''
    Generalized least square estimates

    Parameters:
        sig2s:  covariance matrix of y, pseudobulk
        X:  desing matrix for fixed effects
        y:  pseudobulk
        inverse:    is sig2s inversed
    Returns:
        GLS of fixed effects
    '''
    if not inverse:
        if len( sig2s.shape ) == 1:
            sig2s_inv = 1/sig2s
            A = X.T * sig2s_inv
        else:
            sig2s_inv = np.linalg.inv( sig2s )
            A = X.T @ sig2s_inv
    else:
        sig2s_inv = sig2s
        A = X.T @ sig2s_inv
    B = A @ X
    beta = np.linalg.inv(B) @ A @ y
    return( beta )

def FixedeffectVariance_( beta: np.ndarray, x: np.ndarray ) -> float:
    '''
    Estimate variance explained by fixed effect

    Parameters:
        beta:   fixed effect sizes
        x:  design matrix of fixed effect
    Returns:
        variance explained by fixed effect
    '''
    #xd = x - np.mean(x, axis=0)
    #s = ( xd.T @ xd ) / x.shape[0]
    s = np.cov( x, rowvar=False )
    if len(s.shape) == 0:
        s = s.reshape(1,1)
    return( beta @ s @ beta )

def FixedeffectVariance( beta: np.ndarray, xs: np.ndarray ) -> list:
    '''
    Estimate variance explained by each feature of fixed effect, e.g. cell type, sex

    Parameters:
        beta:   fixed effect sizes
        xs: design matrices for fixed effects
    Returns:
        variances
    '''
    j = 0
    vars = []
    for i,x in enumerate(xs):
        var = FixedeffectVariance_( beta[j:(j+x.shape[1])], x)
        vars.append(var)
        j = j + x.shape[1]
    return( vars )

def fixedeffect_vars(beta: np.ndarray, P: np.ndarray, fixed_covars_d: dict) -> Tuple[dict, dict]:
    '''
    Estimate variance explained by each feature of fixed effect, e.g. cell type, sex

    Parameters:
        beta:   fixed effect sizes
        P:  cell type proportions
        fixed_covars_d: design matrices for fixed effects
    Returns:
        a tuple of 
            #. dict of fixed effects
            #. dict of variances explained
    '''
    # read covars if needed
    fixed_covars_d = read_covars(fixed_covars_d, {})[0]

    beta_d = assign_beta(beta, P, fixed_covars_d)

    fixed_vars_d = {'ct_beta': FixedeffectVariance_( beta_d['ct_beta'], P) }
    for key in fixed_covars_d.keys():
        fixed_vars_d[key] = FixedeffectVariance_( beta_d[key], fixed_covars_d[key] )

#    fixed_covars_l = [P]
#    for key in np.sort(list(fixed_covars_d.keys())):
#        m_ = fixed_covars_d[key]
#        if isinstance( m_, str ):
#            m_ = np.loadtxt( m_ )
#        if len( m_.shape ) == 1:
#            m_ = m_.reshape(-1,1)
#        fixed_covars_l.append( m_ )
#
#    fixedeffect_vars_l = FixedeffectVariance(beta, fixed_covars_l)
#
#    fixedeffect_vars_d = assign_fixedeffect_vars(fixedeffect_vars_l, fixed_covars_d)

    return(beta_d, fixed_vars_d)

def assign_beta(beta_l: list, P: np.ndarray, fixed_covars_d: dict) -> dict:
    '''
    Convert a list of fixed effect to dict for each feature

    Parameters:
        beta_l: fixed effects
        P:  cell type proportions
        fixed_covars_d: design matrices for fixed effects
    Returns:
        dict of fixed effects
    '''
    beta_d = { 'ct_beta': beta_l[:P.shape[1]] }
    beta_l = beta_l[P.shape[1]:]

    for key in np.sort(list(fixed_covars_d.keys())):
        x = fixed_covars_d[key] 
        if len( x.shape ) == 1:
            x = x.reshape(-1,1)
        beta_d[key] = beta_l[:x.shape[1]]
        beta_l = beta_l[x.shape[1]:]

    return(beta_d)

def assign_fixedeffect_vars(fixedeffect_vars_l: list, fixed_covars_d: dict) -> dict:
    '''
    Assign fixed effect variance to each feature

    Parameters:
        fixedeffect_vars_l: fixed effects variances
        fixed_covars_d: design matrices for fixed effects
    Returns:
        fixed effects variances for each feature
    '''
    fixedeffect_vars_d = {'celltype_main_var': fixedeffect_vars_l[0]}
    if len(fixed_covars_d.keys()) > 0:
        for key, value in zip(np.sort(list(fixed_covars_d.keys())), fixedeffect_vars_l[1:]):
            fixedeffect_vars_d[key] = value
    return(fixedeffect_vars_d)

def _random_var( V: np.ndarray, X: np.ndarray ) -> float:
    '''
    Compute variance of random effect

    Parameters:
        V:  covariance matrix of random effect
        X:  design matrix
    Returns:
        variance explained
    '''
    return( np.trace( V @ (X.T @ X) ) / X.shape[0] )

def RandomeffectVariance( Vs: Union[list, dict], Xs: Union[list, dict] ) -> Union[list, dict]:
    if isinstance( Xs, list ):
        if len( np.array( Vs ).shape ) == 1:
            Vs = [V * np.eye(X.shape[1]) for V, X in zip(Vs, Xs)]

        vars = [_random_var(V,X) for V,X in zip(Vs, Xs)]
    elif isinstance( Xs, dict ):
        vars = {}
        for key in Xs.keys():
            V, X = Vs[key], Xs[key]
            if isinstance(V, float):
                V = V  * np.eye(X.shape[1])
            vars[key] = _random_var(V,X)
    return( vars )

def assign_randomeffect_vars(randomeffect_vars_l: list, r2_l: list, random_covars_d: dict) -> Tuple[dict, dict]:
    '''
    Assign variance of random effects
    '''
    randomeffect_vars_d = {}
    r2_d = {}
    keys = np.sort( list(random_covars_d.keys()) )
    if len(keys) != 0:
        for key, v1, v2 in zip( keys, randomeffect_vars_l, r2_l ):
            randomeffect_vars_d[key] = v1
            r2_d[key] = v2

    return( randomeffect_vars_d, r2_d )

def ct_random_var( V: np.ndarray, P: np.ndarray ) -> Tuple[float, np.ndarray]:
    '''
    Compute overall and specific variance of each cell type
    
    Parameters:
        V:  cell type-specific random effect covariance matrix
        P:  cell type proportions
    Returns:
        A tuple of
            #. overall variance
            #. cell type-specific variance
    '''
    N, C = P.shape
    ct_overall_var = _random_var(V, P)
    ct_specific_var = np.array([V[i,i] * ((P[:,i]**2).mean()) for i in range(C)])

    return( ct_overall_var, ct_specific_var )

def cal_variance(beta: np.ndarray, P: np.ndarray, fixed_covars: dict, r2: Union[list, np.ndarray, dict], 
        random_covars: dict) -> Tuple[dict, dict, dict, dict]:
    '''
    Compute variance explained by fixed effects and random effects

    Parameters:
        beta:   fixed effects
        P:  cell type propotions
        fixed_covars: design matrices for additional fixed effects
        r2: variances of additional random effects
        random_covars:  design matrices for additional random effects

    Returns:
        a tuple of 
            #.  dict of fixed effects
            #.  dict of OP variance explained by fixed effects
            #.  dict of random effect variances
            #.  dict of OP variance explained by random effects
    '''
    # calcualte variance of fixed and random effects, and convert to dict
    beta, fixed_vars = fixedeffect_vars( beta, P, fixed_covars ) # fixed effects are always ordered
    if isinstance(r2, list) or isinstance(r2, np.ndarray):
        r2 = dict(zip( np.sort(list(random_covars.keys())), r2 ))
    random_vars = RandomeffectVariance( r2, random_covars )
    return( beta, fixed_vars, r2, random_vars )

def wald_ct_beta(beta: np.ndarray, beta_var: np.ndarray, n: int, P: int) -> float:
    '''
    Wald test on mean expression differentiation

    Parameters:
        beta:   cell type-specific mean expressions
        beta_var:   covariance matrix of cell type-specific mean
        n:  sample size (for Ftest in Wald test)
        P:  number of estimated parameters (for Ftest in Wald test)
    Returns:
        p value for Wald test on mean expression differentiation
    '''
    C = len(beta)
    T = np.concatenate( ( np.eye(C-1), (-1)*np.ones((C-1,1)) ), axis=1 )
    beta = T @ beta
    beta_var = T @ beta_var @ T.T
    return(wald.mvwald_test(beta, np.zeros(C-1), beta_var, n=n, P=P))

def check_R(R: np.ndarray) -> bool:
    '''
    Check R matrix: has to be matrix of 0 and 1
    in the structure of scipy.linalg.block_diag(np.ones((a,1)), np.ones((b,1)), np.ones((c,1))
    '''
    # infer matrix R
    xs = np.sum(R, axis=0).astype('int')
    R_ = np.ones((xs[0],1))
    for i in range(1,len(xs)):
        R_ = linalg.block_diag(R_, np.ones((xs[i],1)))

    if np.any(R != R_):
        print(R[:5,:])
        print(R_[:5,:])
        return( False )
    else:
        return( True )

def order_by_randomcovariate(R: np.ndarray, Xs: list=[], Ys: dict={}
        ) -> Tuple[np.ndarray, np.ndarray, list, dict]:
    '''
    R is the design matrix of 0 and 1 for a random covriate, which we order along by
    Xs or Ys: a list or dict of matrixs we want to order
    '''
    R_df = pd.DataFrame(R)
    index = R_df.sort_values(by=list(R_df.columns), ascending=False).index
    R = np.take_along_axis(R, np.broadcast_to(index, (R.shape[1], R.shape[0])).T, axis=0)
    if not check_R(R):
        sys.exit('Matrix R is wrong!\n')

    new_Xs = []
    for X in Xs:
        if len(X.shape) > 1:
            X = np.take_along_axis(X, np.broadcast_to(index, (X.shape[1], X.shape[0])).T, axis=0)
        else:
            X = np.take_along_axis(X, index, axis=0)
        new_Xs.append(X)

    new_Ys = {}
    for key in Ys.keys():
        Y = Ys[key]
        if len(Y.shape) > 1:
            Y = np.take_along_axis(Y, np.broadcast_to(index, (Y.shape[1], Y.shape[0])).T, axis=0)
        else:
            Y = np.take_along_axis(Y, index, axis=0)
        new_Ys[key] = Y

    return(index, R, new_Xs, new_Ys)

def jk_rmInd(i: int, Y: np.ndarray, K: np.ndarray, ctnu: np.ndarray, fixed_covars: dict={}, 
        random_covars: dict={}, P: Optional[np.ndarray]=None) -> tuple:
    '''
    Remove one individual from the matrices for jackknife

    Parameters:
        i:  index of individual to remove
        Y:  Cell Type-specific Pseudobulk
        K:  Kinship matrix
        ctnu:   cell type-specific noise variance
        fixed_covars:   design matrix for Extra fixed effect
        random_covars:  design matrix for Extra random effect
        P:  cell type proportions
    Returns:
        a tuple of matrices after removing ith individual
    '''
    Y_ = np.delete(Y, i, axis=0)
    K_ = np.delete(np.delete(K,i,axis=0),i,axis=1)
    ctnu_ = np.delete(ctnu, i, axis=0)
    fixed_covars_ = {}
    for key in fixed_covars.keys():
        fixed_covars_[key] = np.delete(fixed_covars[key], i, axis=0)
    random_covars_ = {}
    for key in random_covars.keys():
        random_covars_[key] = np.delete(random_covars[key], i, axis=0)
    if P is None:
        return(Y_, K_, ctnu_, fixed_covars_, random_covars_)
    else:
        P_ = np.delete(P, i, axis=0)
        return(Y_, K_, ctnu_, fixed_covars_, random_covars_, P_)

def lrt(l: float, l0: float, k: int) -> float:
    '''
    Perfomr Likelihood-ration test (LRT)

    Parameters:
        l, l0:  log likelihood for alternative and null hypothesis models
        k:  number of parameters constrained in null model compared to alternative
    Returns:
        p value for LRT
    '''

    Lambda = 2 * (l-l0)
    p = stats.chi2.sf(Lambda, k)
    return(p)

def generate_tmpfn() -> str:
    tmpf = tempfile.NamedTemporaryFile(delete=False)
    tmpfn = tmpf.name
    tmpf.close()
    log.logger.info(tmpfn)
    print(tmpfn)
    return tmpfn

def subprocess_popen(cmd: list, log_fn: str=None) -> None:
    '''
    Run child process using Subprocess.Popen,
    while capture the stdout, stderr, and the exit code of the child process.

    Parameters
    ----------
    cmd : list
        The command for the child process.
        (e.g. ['python', 'test.py'])
    log_fn : str
        The name of log file to keep stdout & stderr

    Returns
    -------
    proc.returncode : str
                    exit code.
    stdout : str
            standard output.
    stderr : str
            strandard error.
    mix : str
            mix of stdout and stderr

    Notes
    -----

    '''
    cmd = [str(x) for x in cmd]
    proc = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE,
            universal_newlines = True)

    # print cmd. proc.args only work for >= python 3.3
    try:
        print(' '.join(proc.args))
    except:
        pass

    stdout, stderr = proc.communicate()

    if log_fn is None:
        sys.stdout.write(stdout)
        if proc.returncode != 0:
            sys.stderr.write('*'*20+'\n'+stderr+'*'*20+'\n')
            print(proc.returncode)
            raise Exception('child exception!')
    else:
        log = open(log_fn, 'w')
        log.write(stdout)
        log.close()
        stdout = open(log_fn).read()
        sys.stdout.write(stdout)
        if proc.returncode != 0:
            log = open(log_fn, 'w')
            log.write(stderr)
            log.close()
            stderr = open(log_fn).read()
            log = open(log_fn, 'w')
            log.write(stdout+'\n')
            log.write('*'*20+'\n'+stderr+'*'*20+'\n')
            log.close()
            sys.stderr.write('*'*20+'\n'+stderr+'*'*20+'\n')
            print(proc.returncode)
            raise Exception('child exception!')

def extract_vcf(input_f:str, panel_f: str=None, samples: npt.ArrayLike=None, samples_f: str=None, 
        pops: npt.ArrayLike=None, snps_f: str=None, 
        snps: npt.ArrayLike=None, maf_filter: bool=True, maf: str='0.000000001', 
        geno: str='1', output_bfile: str=None, bim_only: bool=False, 
        output_vcf_prefix: str=None, update_bim: bool=True,
        ldprune: bool=False, ld_window_size: str='50', ld_step_size: str='5', 
        ld_r2_threshold: str='0.2', memory: str=None, additional_operations: list=None
        ) -> None:
    '''
    Extract samples in pops from 1000G using PLINK (v 1.9).

    Parameters
    ----------
    input_f:    input filename for vcf or bed/bim. If the input file is not end with .vcf/vcf.gz, assume bfile.
    panel_f:    required when pops argu used. 
                Panel filename with at least columns of sample, pop, and super_pop, as in 1000G.
    samples:    a list of sample IDs.
    samples_f:  Sample file used for plink extracting individuals: without header, two columns, both of which are IID.
    pops:   a list of (super)populations to extract.
    snps_f: a file with one snp per line, without header. Only support output bfile.
    snps:   a list of snps. Only support output bfile.
    maf_filter: whether filter variants on maf
    maf:    maf threshold
    geno:   filters out all variants with missing call rates exceeding the provided value
    output_bfile:   prefix for output bfile bed/bim 
    bim_only:   output bim file only --make-just-bim
    update_bim: update the snpname to chr_pos_a2_a1 in bim file.
    output_vcf_prefix:  Prefix for output vcf filename ({output_vcf_prefix}.vcf.bgz).
    memory: memory for plink
    additional_operations:  additional options to give to plink

    # LD prune
    ldprune:    Plink variant pruning using --indep-pairwise
    ld_window_size: a window size in variant count
    ld_step_size:   a variant count to shift the window at the end of each step
    ld_r2_threshold:    r2 threshold

    Notes
    -----
    Exclude monopolymorphic variants.
    vcf-half-cal treated as missing.
    '''

    operations = ['--geno', geno] 
    if memory is not None:
        operations += ['--memory', memory]

    # input file options: vcf file or bfile
    if re.search('(\.vcf$|\.vcf\.gz$)', input_f):
        operations += ['--vcf', input_f, '--double-id', '--keep-allele-order', '--vcf-half-call', 'missing']
    else:
        operations += ['--bfile', input_f]

    # make sample file for plink
    if samples_f:
        operations += ['--keep', samples_f]
    elif samples is not None:
        samples = pd.DataFrame({'sample': samples})
        samples_f = generate_tmpfn()
        samples[['sample', 'sample']].to_csv(samples_f, sep='\t', index=False, header=False)
        operations += ['--keep', samples_f]
    elif pops is not None:
        panel = pd.read_table( panel_f )
        panel = panel.loc[np.isin(panel['pop'], pops) | np.isin(panel['super_pop'], pops)]
        panel = panel.reset_index(drop=True)
        samples_f = generate_tmpfn()
        panel[['sample', 'sample']].to_csv(samples_f, sep='\t', index=False, header=False)
        operations += ['--keep', samples_f]

    # extract snps after updating snp names in bim file
    if snps_f:
        pass
    elif snps is not None:
        snps_f = generate_tmpfn()
        snps = pd.DataFrame({'snp': snps})
        snps[['snp']].to_csv(snps_f, sep='\t', index=False, header=False)

    # maf filter
    if maf_filter:
        operations += ['--maf', str(maf)]

    if additional_operations is not None:
        operations += additional_operations

    # extract samples to bfile
    if output_bfile:
        if bim_only:
            operations += ['--make-just-bim']
        else:
            operations += ['--make-bed']
        subprocess_popen(['plink', '--out', output_bfile]+operations)
        if update_bim:
            update_bim_snpname(output_bfile+'.bim')
        if snps_f:
            cmd = ['plink', '--bfile', output_bfile, '--keep-allele-order',
                    '--extract', snps_f, '--make-bed', '--out', output_bfile]
            if memory:
                cmd += ['--memory', memory]
            subprocess_popen( cmd )
        if ldprune:
            plink_ldprune(output_bfile, ld_window_size=ld_window_size, ld_step_size=ld_step_size,
                    ld_r2_threshold=ld_r2_threshold, output_bfile=output_bfile, memory=memory)

    # extract samples to vcf
    # double check whether ref allele reserved
    if output_vcf_prefix:
        subprocess_popen(['plink2', '--export', 'vcf', 'bgz', 'id-paste=iid',
                '--out', output_vcf_prefix]+operations)

def plink_ldprune(bfile: str, ld_window_size: str='50', ld_step_size: str='5', ld_r2_threshold: str='0.2',
        output_bfile: str=None, memory: str=None) -> None:
    '''
    Plink variance pruning using --indep-pairwise.

    At each step, pairs of variants in the current window with squared correlation greater than
    the threshold are noted, and variants are greedily pruned from the window
    until no such pairs remain.

    Parameters
    ----------
    bfile:  Prefix for bfile {bfile}.bed
    ld_window_size: a window size in variant count
    ld_step_size:   a variant count to shift the window at the end of each step
    ld_r2_threshold:    r2 threshold
    output_bfile: Prefix for pruned bfile. Default: same as {bfile}
    memory: memory of Plink
    '''
    if not output_bfile:
        output_bfile = bfile+'.ldprune'
    prunefile_fn = generate_tmpfn()
    cmd1 = ['plink', '--bfile', bfile, '--indep-pairwise',
            ld_window_size, ld_step_size, ld_r2_threshold,
            '--out', prunefile_fn]
    cmd2 = ['plink', '--bfile', bfile, '--extract',
            prunefile_fn+'.prune.in', '--memory', memory, '--keep-allele-order',
            '--make-bed', '--out', output_bfile]
    if memory:
        cmd1 += ['--memory', memory]
        cmd2 += '--memory', memory
    basic_fun.subprocess_popen( cmd1 )
    basic_fun.subprocess_popen(cmd2)

def update_bim_snpname(bim_fn: str) -> None:
    '''
    Update bim file snpname to chr_pos_a2_a1 (a2 a1 are ref alt if keep-allele-order from vcf)

    Parameters
    ----------
    bim_fn: filename of bim
    '''
    bim = pd.read_table(bim_fn, header=None, names=['chr', 'snp', 'genetic', 'pos', 'a1', 'a2']) # a2 a1 are ref alt
    bim['snp'] = bim['chr'].astype('str')+'_'+bim['pos'].astype('str')+'_'+bim['a2']+'_'+bim['a1']
    bim.to_csv(bim_fn, sep='\t', index=False, header=False)

def grm(bfile: str, chr: int, start: int, end: int, r: int, rel: str, tool: str='plink') -> int:
    """
    Compute kinship matrix for a genomic region (start-r, end+r)

    Parameters:
        bfile:  prefix for chr/genome bed/bim files
        chr:    chromosome
        start:  start position of gene
        end:    end position of gene
        r:  radius to the gene
        rel:    prefix for relationship matrix file (prefix.rel.bin for plink, prefix.grm.bin for gcta)
        tool:   plink or gcta to compute grm
    Returns:
        number of snps in the regions
    """
    start = max(0, start - r)
    end = end + r

    # check number of SNPs in the region
    bim = pd.read_csv(bfile+'.bim', sep='\s+', names=['chr', 'snp', 'cm', 'bp','a1','a2'])
    nsnp = bim.loc[(bim['chr']==chr) & (bim['bp'] > start) & (bim['bp'] < end)].shape[0]

    # compute kinship matrix
    if nsnp > 0:
        if tool == 'plink':
            cmd = ['plink', '--bfile', bfile,
                   '--chr', chr, '--from-bp', start, '--to-bp', end,
                   '--make-rel', 'bin',
                   '--out', rel]
            subprocess_popen( cmd )
        elif tool == 'gcta':
            tmp = generate_tmpfn()
            cmd = ['plink', '--bfile', bfile,
                   '--chr', chr, '--from-bp', start, '--to-bp', end,
                    '--make-bed', '--out', tmp]
            subprocess_popen( cmd )
            cmd = ['gcta', '--bfile', tmp,
                   '--make-grm', '--out', rel]
            subprocess_popen( cmd )

    return nsnp

def design(inds: npt.ArrayLike, pca: pd.DataFrame=None, PC: int=None, cat: pd.Series=None, 
        con: pd.Series=None) -> dict:
    """
    Construct design matrix

    Parameters:
        inds:   order of individuals
        pca:    dataframe of pcs, with index: individuals (sort not required) and columns (PC1-PCx)
        PC: number to PC to adjust
        cat:    series of category elements e.g. sex: male and female
        con:    series of continuous elements e.g. age

    Returns:
        a design matrix
    """

    # pca
    if pca is not None:
        pcs = [f'PC{i}' for i in range(1, int(PC)+1)]
        pca = pca.loc[inds, pcs]
        return( pca )
    elif cat is not None:
        return( pd.get_dummies( cat, drop_first=True ).loc[inds,:] )
    elif con is not None:
        return( con[inds] )

def L_f(C, c1, c2):
    # to build L matrix of 0 and 1
    L = np.zeros((C,C), dtype='int8')
    L[c1,c2] = 1
    return( L )
