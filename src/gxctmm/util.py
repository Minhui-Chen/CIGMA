from typing import Tuple, Optional, Union, List

import os, tempfile, sys, subprocess, re, gzip, time, shutil
import numpy as np, pandas as pd
import numpy.typing as npt
# import jax
# import jax.numpy as jnp
# import jax.scipy as jsp
# import jaxopt
# import tensorflow as tf
# import tensorflow_probability as tfp
import rpy2.robjects as ro
from rpy2.robjects import r, pandas2ri, numpy2ri
from rpy2.robjects.conversion import localconverter
from scipy import stats, linalg, optimize
from numpy.random import default_rng

from . import log, fit, wald



def merge_dicts(reps, out):
    for key, value in reps[0].items():
        print(key)
        if isinstance(value, dict):
            out[key] = {}
            merge_dicts([rep[key] for rep in reps], out[key])
        else:
            out[key] = [rep[key] for rep in reps]


def design(inds: npt.ArrayLike, pca: Optional[pd.DataFrame] = None, 
           PC: Optional[int] = None, cat: Optional[pd.Series] = None,
           con: Optional[pd.Series] = None, drop_first: bool = True) -> np.ndarray:
    """
    Construct design matrix

    Parameters:
        inds:   order of individuals
        pca:    dataframe of pcs, with index: individuals (sort not required) and columns (PC1-PCx)
        PC: number to PC to adjust
        cat:    series of category elements e.g. sex: male and female
        con:    series of continuous elements e.g. age
        drop_first: drop the first column

    Returns:
        a design matrix
    """

    # pca
    if pca is not None:
        pcs = [f'PC{i}' for i in range(1, int(PC) + 1)]
        return pca.loc[inds, pcs].to_numpy()
    elif cat is not None:
        return pd.get_dummies(cat, drop_first=drop_first, dtype='int').loc[inds, :].to_numpy()
    elif con is not None:
        return con[inds].to_numpy()


def get_X(fixed_covars: dict, N: int, C: int, fixed_shared: bool = True) -> np.ndarray:
    """
    Compute the design matrix X for fixed effects.

    Parameters:
        fixed_covars:   a dict of design matrices for each feature of fixed effect,
                        except for cell type-specific fixed effect
        N:  number of individuals
        C:  number of cell types
        fixed_shared: whether fixed effect is shared across cell types
    Returns:
        Design matrix for fixed effects
    """

    X = np.kron(np.ones((N, 1)), np.eye(C))
    for key in sorted(fixed_covars.keys()):
        m = fixed_covars[key]
        if len(m.shape) == 1:
            m = m.reshape(-1, 1)
        # remove columns that have only one unique value
        cols_to_keep = [i for i in range(m.shape[1]) 
                        if np.unique(m[:, i]).size > 1]
        if len(cols_to_keep) < m.shape[1]:
            log.logger.info(f'Remove {m.shape[1] - len(cols_to_keep)} columns in {key}')
        m = m[:, cols_to_keep]
        # if m is empty, skip it
        if m.size == 0:
            continue
        if fixed_shared:
            X = np.concatenate((X, np.kron(m, np.ones((C, 1)))), axis=1)
        else:
            X = np.concatenate((X, np.kron(m, np.eye(C))), axis=1)
    return X


def read_covars(fixed_covars: dict = {}, random_covars: dict = {}, C: Optional[int] = None) -> tuple:
    """
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
    """

    def read(covars):
        tmp = {}
        for key in covars.keys():
            f = covars[key]
            if isinstance(f, str):
                tmp[key] = np.loadtxt(f)
            else:
                tmp[key] = f
        return tmp

    fixed_covars = read(fixed_covars)
    random_covars = read(random_covars)
    n_fixed, n_random = len(fixed_covars.keys()), len(random_covars.keys())
    random_keys = list(np.sort(list(random_covars.keys())))
    Rs = [random_covars[key] for key in random_keys]
    if C:
        random_MMT = [np.repeat(np.repeat(R @ R.T, C, axis=0), C, axis=1) for R in Rs]
    else:
        random_MMT = [R @ R.T for R in Rs]

    return fixed_covars, random_covars, n_fixed, n_random, random_keys, Rs, random_MMT


def age_group(age: pd.Series, bins: npt.NDArray = np.arange(25, 91, 5)):
    """
    Separate age groups
    """
    new = pd.Series(np.digitize(age, bins), index=age.index)
    if age.name is None:
        return new
    else:
        return new.rename(age.name)


def yazar_covars(inds: list, meta_f: str, geno_pca_f: str, op_pca_f: str,
                 geno_pca_n: int = 6, op_pca_n: int = 1
                 ) -> Tuple[dict, dict]:
    """
    Read Yazar covariates

    Parameters:
        inds:  individuals to include in the design matrix
        meta_f:  file of metadata
        geno_pca_f:  file of genotype PCs
        op_pca_f:  file of other PCs
        geno_pca_n:  number of genotype PCs to use
        op_pca_n:  number of other PCs to use
    Returns:
        a tuple of
            #. dict of genotype PCs
            #. dict of other PCs
    """
    
    # read metadata
    meta = pd.read_table(meta_f, usecols=['individual', 'sex', 'age', 'pool'])
    meta = meta.loc[meta['individual'].isin(inds)]
    meta = meta.drop_duplicates()
    meta = meta.set_index('individual')
    geno_pca = pd.read_table(geno_pca_f, index_col=0).drop('IID', axis=1)
    op_pca = pd.read_table(op_pca_f, index_col=0)

    fixed_covars = {
        'op_pca': construct_design(inds, pca=op_pca, PC=op_pca_n).astype('float32'),
        'sex': construct_design(inds, cat=meta['sex']),
        'age': construct_design(inds, cat=age_group(meta['age']))
    }
    if geno_pca_n != 0:
        fixed_covars['geno_pca'] = construct_design(inds, pca=geno_pca, PC=geno_pca_n).astype('float32')

    random_covars = {
        'batch': construct_design(inds, cat=meta['pool'], drop_first=False)
    }

    return fixed_covars, random_covars


def perez_covars(inds: list, meta_f: str, geno_pca_f: str, op_pca_f: str, 
                 batch: Optional[str] = 'shared', geno_pca_n: int = 7,
                 op_pca_n: int = 10, include_dataset: bool = False
                 ) -> Tuple[dict, dict]:
    """
    Read Perez covariates

    Parameters:
        inds:  individuals to include in the design matrix
        meta_f:  file of metadata
        geno_pca_f:  file of genotype PCs
        op_pca_f:  file of other PCs
        batch:  batch effect to include in the design matrix. None, 'fixed', 'shared', or 'specific'
        geno_pca_n:  number of genotype PCs to use
        op_pca_n:  number of other PCs to use
        include_dataset:  whether to include dataset as a fixed effect
    Returns:
        a tuple of
            #. dict of genotype PCs
            #. dict of other PCs
    """
    
    # read metadata
    meta = pd.read_table(meta_f, usecols=['ind_cov', 'Sex', 'Age', 'batch_cov', 
                                          'pop_cov', 'SLE_status', 
                                          'Processing_Cohort', 'dataset'])
    meta = meta.drop_duplicates()
    meta = meta.set_index('ind_cov')
    geno_pca = pd.read_table(geno_pca_f, index_col=0).drop('IID', axis=1)
    op_pca = pd.read_table(op_pca_f, index_col=0)

    fixed_covars = {
        'op_pca': construct_design(inds, pca=op_pca, PC=op_pca_n).astype('float32'),
        'geno_pca': construct_design(inds, pca=geno_pca, PC=geno_pca_n).astype('float32'),
        'sex': construct_design(inds, cat=meta['Sex']),
        'age': construct_design(inds, cat=age_group(meta['Age'], bins=np.arange(25, 71, 5))),
        'sle': construct_design(inds, cat=meta['SLE_status']),
        'cohort': construct_design(inds, cat=meta['Processing_Cohort']),
        'population': construct_design(inds, cat=meta['pop_cov']),
    }
    if include_dataset:
        # only needed in mega analysis. 
        # Euro controls: CLUES are all in Processing cohort 4 and ImmVar are all not in that cohort
        fixed_covars['dataset'] = construct_design(inds, cat=meta['dataset']) 

    random_covars = {}
    if batch is None or batch is False:
        pass
    elif batch == 'shared':
        random_covars = {
            'batch': construct_design(inds, cat=meta['batch_cov'], drop_first=False)
        }

    return fixed_covars, random_covars


def transform_he_to_reml(free_he: dict) -> list:
    """
    Transform HE results to initial parameters for REML
    """
    par = [free_he['hom_g2']] + np.diag(free_he['V']).tolist()
    par += [free_he['hom_e2']] + np.diag(free_he['W']).tolist()
    
    if 'hom_g2_b' in free_he.keys():
        par += [free_he['hom_g2_b']] + np.diag(free_he['V_b']).tolist()
    if len(free_he['r2'].keys()) > 0:
        if len(free_he['r2'].keys()) > 1:
            raise ValueError("Currently only support one extra random effect")
        
        key = list(free_he['r2'].keys())[0]
        if isinstance(free_he['r2'][key], float):
            par += [free_he['r2'][key]]
        else:
            raise ValueError("Currently only support shared extra random effect")

    # change negative values to small positive
    par = np.array(par)
    par[par < 0] = 0.001
    par = par.tolist()

    return par


def optim(fun: callable, par: list, args: tuple, method: Optional[str],
          ) -> Tuple[dict, dict]:
    """
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
    """
    if method is None:
        method = 'BFGS'

    start = time.time()

    if method == 'BFGS-Nelder':
        start1 = time.time()
        out1 = optimize.minimize(fun, par, args=args, method='BFGS')
        start2 = time.time()
        out = optimize.minimize(fun, out1['x'], args=args, method='Nelder-Mead')
        out['l'] = out['fun'] * (-1)
        opt = {'method1': 'BFGS', 'success1': out1['success'], 'status1': out1['status'],
            'message1': out1['message'], 'l1': out1['fun'] * (-1),
            'method': 'Nelder-Mead', 'success': out['success'], 'status': out['status'],
            'message': out['message'], 'l': out['fun'] * (-1),
            'initial': par,
            'time1': start2 - start1, 'time2': time.time() - start2}
    else:
        out = optimize.minimize(fun, par, args=args, method=method)
        out['l'] = out['fun'] * (-1)
        opt = {'method': method, 'success': out['success'], 'status': out['status'],
            'message': out['message'], 'l': out['fun'] * (-1),
            'initial': par}

    opt['time'] = time.time() - start

    return out, opt


def check_optim(opt: dict, vars: dict, fixed_vars: dict, 
                random_vars: dict = {}, cut: float = 1) -> bool:
    """
    Check whether optimization converged successfully

    Parameters:
        opt:    dict of optimization results, e.g. log-likelihood
        var: dict of variances explained by e.g. 
                hom_g2: variance of genetic effect shared across cell types
                hom_e2: variance of env effect shared across cell types
                ct_overall_g_var:  overall variance explained by cell type-specific genetic effect
                ct_overall_e_var:  overall variance explained by cell type-specific env effect
        fixed_vars:  dict of variances explained by each fixed effect feature, including cell type-specific fixed effect
        random_vars:  dict of variances explained by each random effect feature, doesn't include cell type-shared or -specific effect
        cut:    threshold for large variance
    Returns:
        True:   optim failed to converge
        False:  optim successfully to converge
    """
    if ((opt['l'] < -1e10) or (not opt['success']) or 
            np.any(np.array(list(vars.values())) > cut) or
            np.any(np.array(list(fixed_vars.values())) > cut) or
            np.any(np.array(list(random_vars.values())) > cut)):
        info = f"l:{opt['l']}, {opt['success']}\n"
        if 'message' in opt.keys():
            info += f"message: {opt['message']}\n"
        if 'Wolfe' in opt.keys():
            info += f"Wolfe: {opt['Wolfe']}\n"
        info += f"hom_g2: {vars['hom_g2']}, hom_e2: {vars['hom_e2']} \n"
        info += f"overall cell type-specific genetic effect: {vars['ct_overall_g_var']} \n"
        if 'ct_overall_g2_var' in vars.keys():
            info += f"overall cell type-specific genetic effect 2: {vars['ct_overall_g2_var']} \n"
        info += f"overall cell type-specific env effect: {vars['ct_overall_e_var']} \n"
        info += f"cell type-specific fixed effect: {fixed_vars['ct_beta']}\n"
        for key in random_vars.keys():
            info += f"{key}: {random_vars[key]} \n"
        if 'success1' in opt.keys():
            info += f"{opt['success1']}: {opt['message1']} \n"
        log.logger.info(info)

        return True
    else:
        return False


def re_optim(out: dict, opt: dict, fun: callable, par: list, args: tuple, 
             method: Optional[str], nrep: int = 10
             ) -> Tuple[dict, dict]:
    """
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
    """
    rng = default_rng()
    # print( out['fun'] )
    for i in range(nrep):
        par_ = np.array(par) * rng.gamma(2, 1 / 2, len(par))
        out_, opt_ = optim(fun, par_.tolist(), args=args, method=method)
        log.logger.info(f"loglike: {out_['l']}")
        if (not out['success']) and out_['success']:
            out, opt = out_, opt_
        elif (out['success'] == out_['success']) and (out['l'] < out_['l']):
            out, opt = out_, opt_
    # print( out['fun'] )
    return out, opt


def dict2Rlist(X: dict) -> object:
    """
    Transform a python dictionary to R list

    Parameters:
        X:  python dictionary
    Returns:
        R list
    """
    if len(X.keys()) == 0:
        return r('NULL')
    else:
        keys = np.sort(list(X.keys()))
        rlist = ro.ListVector.from_length(len(keys))
        for i in range(len(keys)):
            if isinstance(X[keys[i]], str):
                if os.path.exists(X[keys[i]]):
                    rlist[i] = r['as.matrix'](r['read.table'](X[keys[i]]))
                else:
                    try:
                        rlist[i] = np.array([X[keys[i]]])
                    except:
                        numpy2ri.activate()  # don't think this is useful
                        rlist[i] = np.array([X[keys[i]]])
                        numpy2ri.deactivate()  # deactivate would cause numpy2ri deactivated in calling fun
            elif isinstance(X[keys[i]], pd.DataFrame):
                with localconverter(ro.default_converter + pandas2ri.converter):
                    rlist[i] = r['as.matrix'](X[keys[i]])
            elif isinstance(X[keys[i]], np.ndarray):
                try:
                    rlist[i] = r['as.matrix'](X[keys[i]])
                except:
                    numpy2ri.activate()
                    rlist[i] = r['as.matrix'](X[keys[i]])
                    numpy2ri.deactivate()
            elif isinstance(X[keys[i]], int) or isinstance(X[keys[i]], float):
                try:
                    rlist[i] = np.array([X[keys[i]]])
                except:
                    numpy2ri.activate()
                    rlist[i] = np.array([X[keys[i]]])
                    numpy2ri.deactivate()
        return rlist


def generate_HE_initial(he: dict, ML: bool = False, REML: bool = False) -> list:
    """
    Convert HE estimates to initial parameter for ML / REML

    Parameters:
        he: estiamtes from HE
        ML: generate initial parameters for ML
        REML:   generate initial parameters for REML
    Returns:
        initial parameters for ML / REML
    """
    initials_random_effects = []
    # homogeneous effect
    if 'hom2' in he.keys():
        initials_random_effects.append(he['hom2'])
    # heterogeneous effect
    if 'V' in he.keys():
        C = he['V'].shape[0]
        # determine model based on V
        if np.any(np.diag(np.diag(he['V'])) != he['V']):
            # Full model
            initials_random_effects = initials_random_effects + list(he['V'][np.triu_indices(C)])
        elif len(np.unique(np.diag(he['V']))) == 1:
            # IID model
            initials_random_effects.append(he['V'][0, 0])
        else:
            # Free model
            initials_random_effects = initials_random_effects + list(np.diag(he['V']))
    # other random covariates
    if 'r2' in he.keys():
        for key in np.sort(list(he['r2'].keys())):
            initials_random_effects.append(he['r2'][key])

    if REML is True:
        return initials_random_effects

    initials_fixed_effects = list(he['beta']['ct_beta'])
    for key in np.sort(list(he['beta'].keys())):
        if key != 'ct_beta':
            initials_fixed_effects = initials_fixed_effects + list(he['beta'][key])

    if ML is True:
        return initials_fixed_effects + initials_random_effects


def glse(sig2s: np.ndarray, X: np.ndarray, y: np.ndarray, inverse: bool = False) -> np.ndarray:
    """
    Generalized least square estimates

    Parameters:
        sig2s:  covariance matrix of y, pseudobulk
        X:  desing matrix for fixed effects
        y:  pseudobulk
        inverse:    is sig2s inversed
    Returns:
        GLS of fixed effects
    """
    if not inverse:
        if len(sig2s.shape) == 1:
            sig2s_inv = 1 / sig2s
            A = X.T * sig2s_inv
        else:
            try:
                sig2s_inv = np.linalg.inv(sig2s)
            except np.linalg.LinAlgError as e:
                print(e)
                w = np.linalg.eigvalsh(sig2s)
                print(f"Eigenvalues of sig2s: {w}")
                sys.exit(1)

            A = X.T @ sig2s_inv
    else:
        sig2s_inv = sig2s
        A = X.T @ sig2s_inv
    B = A @ X
    beta = np.linalg.inv(B) @ A @ y
    return beta


def FixedeffectVariance_(beta: np.ndarray, x: np.ndarray) -> float:
    """
    Estimate variance explained by fixed effect

    Parameters:
        beta:   fixed effect sizes
        x:  design matrix of fixed effect
    Returns:
        variance explained by fixed effect
    """
    # xd = x - np.mean(x, axis=0)
    # s = ( xd.T @ xd ) / x.shape[0]
    s = np.cov(x, rowvar=False)
    if len(s.shape) == 0:
        s = s.reshape(1, 1)
    return beta @ s @ beta


def FixedeffectVariance(beta: np.ndarray, xs: np.ndarray) -> list:
    """
    Estimate variance explained by each feature of fixed effect, e.g. cell type, sex

    Parameters:
        beta:   fixed effect sizes
        xs: design matrices for fixed effects
    Returns:
        variances
    """
    j = 0
    vars = []
    for i, x in enumerate(xs):
        var = FixedeffectVariance_(beta[j:(j + x.shape[1])], x)
        vars.append(var)
        j = j + x.shape[1]
    return vars


def fixedeffect_vars(beta: np.ndarray, P: np.ndarray, fixed_covars: dict) -> Tuple[dict, dict]:
    """
    Estimate variance explained by each feature of fixed effect, e.g. cell type, sex

    Parameters:
        beta:   fixed effect sizes
        P:  cell type proportions
        fixed_covars: design matrices for Extra fixed effects
    Returns:
        a tuple of
            #. dict of fixed effects
            #. dict of variances explained
    """
    # read covars if needed
    fixed_covars = read_covars(fixed_covars, {})[0]

    beta_d = assign_beta(beta, P, fixed_covars)

    fixed_vars = {'ct_beta': FixedeffectVariance_(beta_d['ct_beta'], P)}
    for key in fixed_covars.keys():
        fixed_vars[key] = FixedeffectVariance_( beta_d[key], fixed_covars[key] )

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

    return beta_d, fixed_vars


def assign_beta(beta_l: Union[np.ndarray, list], P: np.ndarray, fixed_covars: dict) -> dict:
    """
    Convert a list of fixed effect to dict for each feature

    Parameters:
        beta_l: fixed effects
        P:  cell type proportions
        fixed_covars: design matrices for Extra fixed effects
    Returns:
        dict of fixed effects
    """
    N, C = P.shape

    beta_d = {'ct_beta': beta_l[:C]}
    beta_l = beta_l[C:]

    n = 0
    for value in fixed_covars.values():
        if len(value.shape) == 1:
            n += 1
        else:
            n += value.shape[1]

    if len(beta_l) == n:
        shared = True
    elif len(beta_l) == (n * C):
        shared = False
    else:
        sys.exit('Wrong dimension')
    # log.logger.info(f'{len(beta_l)}, {n}')

    for key in sorted(fixed_covars.keys()):
        x = fixed_covars[key]
        if len(x.shape) == 1:
            x = x.reshape(-1, 1)
        if shared:
            beta_d[key] = beta_l[:x.shape[1]]
            beta_l = beta_l[x.shape[1]:]
        else:
            beta_d[key] = beta_l[:(x.shape[1] * C)]
            beta_l = beta_l[x.shape[1] * C]

    return beta_d


def assign_fixedeffect_vars(fixedeffect_vars_l: list, fixed_covars_d: dict) -> dict:
    """
    Assign fixed effect variance to each feature

    Parameters:
        fixedeffect_vars_l: fixed effects variances
        fixed_covars_d: design matrices for fixed effects
    Returns:
        fixed effects variances for each feature
    """
    fixedeffect_vars_d = {'celltype_main_var': fixedeffect_vars_l[0]}
    if len(fixed_covars_d.keys()) > 0:
        for key, value in zip(np.sort(list(fixed_covars_d.keys())), fixedeffect_vars_l[1:]):
            fixedeffect_vars_d[key] = value
    return fixedeffect_vars_d


def _random_var(V: np.ndarray, X: np.ndarray) -> float:
    """
    Compute variance of random effect

    Parameters:
        V:  covariance matrix of random effect
        X:  design matrix
    Returns:
        variance explained
    """
    return np.trace(V @ (X.T @ X)) / X.shape[0]


def RandomeffectVariance(Vs: Union[list, dict], Xs: Union[list, dict]) -> Union[list, dict]:
    if isinstance(Xs, list):
        if len(np.array(Vs).shape) == 1:
            Vs = [V * np.eye(X.shape[1]) for V, X in zip(Vs, Xs)]

        vars = [_random_var(V, X) for V, X in zip(Vs, Xs)]
    elif isinstance(Xs, dict):
        vars = {}
        for key in Xs.keys():
            V, X = Vs[key], Xs[key]
            if isinstance(V, float):
                V = V * np.eye(X.shape[1])
            vars[key] = _random_var(V, X)
    return vars


def assign_randomeffect_vars(randomeffect_vars_l: list, r2_l: list, random_covars_d: dict) -> Tuple[dict, dict]:
    """
    Assign variance of random effects
    """
    randomeffect_vars_d = {}
    r2_d = {}
    keys = np.sort(list(random_covars_d.keys()))
    if len(keys) != 0:
        for key, v1, v2 in zip(keys, randomeffect_vars_l, r2_l):
            randomeffect_vars_d[key] = v1
            r2_d[key] = v2

    return randomeffect_vars_d, r2_d


def ct_random_var(V: np.ndarray, P: np.ndarray) -> Tuple[float, np.ndarray]:
    """
    Compute overall and specific variance of each cell type

    Parameters:
        V:  cell type-specific random effect covariance matrix
        P:  cell type proportions
    Returns:
        A tuple of
            #. overall variance
            #. cell type-specific variance
    """
    N, C = P.shape
    ct_overall_var = _random_var(V, P)
    ct_specific_var = np.array([V[i, i] * (P[:, i] ** 2).mean() 
                                for i in range(C)])

    return ct_overall_var, ct_specific_var


def cal_variance(beta: np.ndarray, P: np.ndarray, fixed_covars: dict,
                 r2: dict, random_covars: dict
                 ) -> Tuple[dict, dict, dict]:
    """
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
            #.  dict of OP variance explained by random effects
    """
    # calculate variance of fixed and random effects, and convert to dict
    beta_d, fixed_vars = fixedeffect_vars(beta, P, fixed_covars)  # fixed effects are always ordered
    random_vars = {}
    if len(r2.keys()) > 0:
        if isinstance(list(r2.values())[0], float):
            random_vars = RandomeffectVariance(r2, random_covars)

    return beta_d, fixed_vars, random_vars


def wald_ct_beta(beta: np.ndarray, beta_var: np.ndarray, 
                 n: Optional[int]=None, P: Optional[int]=None) -> float:
    """
    Wald test on mean expression differentiation

    Parameters:
        beta:   cell type-specific mean expressions
        beta_var:   covariance matrix of cell type-specific mean
        n:  sample size (for Ftest in Wald test)
        P:  number of estimated parameters (for Ftest in Wald test)
    Returns:
        p value for Wald test on mean expression differentiation
    """
    C = len(beta)
    T = np.concatenate((np.eye(C - 1), (-1) * np.ones((C - 1, 1))), axis=1)
    beta = T @ beta
    beta_var = T @ beta_var @ T.T
    if n:
        p = wald.mvwald_test(beta, np.zeros(C - 1), beta_var, n=n, P=P)
    else:
        p = wald.mvwald_test(beta, np.zeros(C - 1), beta_var, Ftest=False)

    return p


def check_R(R: np.ndarray) -> bool:
    """
    Check R matrix: has to be matrix of 0 and 1
    in the structure of scipy.linalg.block_diag(np.ones((a,1)), np.ones((b,1)), np.ones((c,1))
    """
    # infer matrix R
    xs = np.sum(R, axis=0).astype('int')
    R_ = np.ones((xs[0], 1))
    for i in range(1, len(xs)):
        R_ = linalg.block_diag(R_, np.ones((xs[i], 1)))

    if np.any(R != R_):
        print(R[:5, :])
        print(R_[:5, :])
        return False
    else:
        return True


def order_by_randomcovariate(R: np.ndarray, Xs: list = [], Ys: dict = {}
                             ) -> Tuple[np.ndarray, np.ndarray, list, dict]:
    """
    R is the design matrix of 0 and 1 for a random covriate, which we order along by
    Xs or Ys: a list or dict of matrixs we want to order
    """
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

    return index, R, new_Xs, new_Ys


def jk_rmInd(i: int, Y: np.ndarray, K: np.ndarray, ctnu: np.ndarray, 
             fixed_covars: dict = {}, random_covars: dict = {}, 
             P: Optional[np.ndarray] = None, K2: Optional[np.ndarray] = None
             ) -> list:
    """
    Remove one individual from the matrices for jackknife

    Parameters:
        i:  index of individual to remove
        Y:  Cell Type-specific Pseudobulk
        K:  Kinship matrix
        ctnu:   cell type-specific noise variance
        fixed_covars:   design matrix for Extra fixed effect
        random_covars:  design matrix for Extra random effect
        P:  cell type proportions
        K2: second kinship matrix
    Returns:
        a list of matrices after removing ith individual
    """

    Y_ = np.delete(Y, i, axis=0)
    K_ = np.delete(np.delete(K, i, axis=0), i, axis=1)
    ctnu_ = np.delete(ctnu, i, axis=0)
    fixed_covars_ = {}
    for key in fixed_covars.keys():
        fixed_covars_[key] = np.delete(fixed_covars[key], i, axis=0)
    random_covars_ = {}
    for key in random_covars.keys():
        random_covars_[key] = np.delete(random_covars[key], i, axis=0)
        
    out = [Y_, K_, ctnu_, fixed_covars_, random_covars_]

    if P is not None:
        out.append(np.delete(P, i, axis=0))
    
    if K2 is not None:
        out.append(np.delete(np.delete(K2, i, axis=0), i, axis=1))

    return out


def lrt(l: float, l0: float, k: int) -> float:
    """
    Perfomr Likelihood-ration test (LRT)

    Parameters:
        l:  log likelihood for alternative hypothesis models
        l0: log likelihood for null hypothesis models
        k:  number of parameters constrained in null model compared to alternative
    Returns:
        p value for LRT
    """

    Lambda = 2 * (l - l0)
    p = stats.chi2.sf(Lambda, k)
    return p


def generate_tmpfn() -> str:
    tmpf = tempfile.NamedTemporaryFile(delete=False)
    tmpfn = tmpf.name
    tmpf.close()
    log.logger.info(tmpfn)
    return tmpfn


def subprocess_popen(cmd: list, log_fn: Optional[str] = None) -> None:
    """
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

    """
    cmd = [str(x) for x in cmd]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True)

    # print cmd. proc.args only work for >= python 3.3
    try:
        log.logger.info(' '.join(proc.args))
    except:
        pass

    stdout, stderr = proc.communicate()

    if log_fn is None:
        sys.stdout.write(stdout)
        if proc.returncode != 0:
            sys.stderr.write('*' * 20 + '\n' + stderr + '*' * 20 + '\n')
            print(proc.returncode)
            raise Exception('child exception!')
    else:
        ilog = open(log_fn, 'w')
        ilog.write(stdout)
        ilog.close()
        stdout = open(log_fn).read()
        sys.stdout.write(stdout)
        if proc.returncode != 0:
            ilog = open(log_fn, 'w')
            ilog.write(stderr)
            ilog.close()
            stderr = open(log_fn).read()
            ilog = open(log_fn, 'w')
            ilog.write(stdout + '\n')
            ilog.write('*' * 20 + '\n' + stderr + '*' * 20 + '\n')
            ilog.close()
            sys.stderr.write('*' * 20 + '\n' + stderr + '*' * 20 + '\n')
            print(proc.returncode)
            raise Exception('child exception!')


def extract_vcf(input_f: str, panel_f: Optional[str] = None, samples: Optional[npt.ArrayLike] = None, 
                samples_f: Optional[str] = None, pops: Optional[npt.ArrayLike] = None, 
                snps_f: Optional[str] = None, snps: Optional[npt.ArrayLike] = None, 
                maf_filter: bool = True, maf: Union[float, str] = '0.000000001',
                geno: str = '1', output_bfile: Optional[str] = None, bim_only: bool = False,
                ped: bool = False, output_vcf_prefix: Optional[str] = None, 
                update_bim: bool = True, ldprune: bool = False, ld_window_size: str = '50', 
                ld_step_size: str = '5', ld_r2_threshold: str = '0.2', memory: Optional[str] = None, 
                additional_operations: Optional[list] = None
                ) -> None:
    """
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
    """

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
        panel = pd.read_table(panel_f)
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
        elif ped:
            update_bim = False
        else:
            operations += ['--make-bed']
        subprocess_popen(['plink', '--out', output_bfile] + operations)
        if update_bim:
            update_bim_snpname(output_bfile + '.bim')
        if snps_f:
            cmd = ['plink', '--bfile', output_bfile, '--keep-allele-order',
                   '--extract', snps_f, '--make-bed', '--out', output_bfile]
            if memory:
                cmd += ['--memory', memory]
            subprocess_popen(cmd)
        if ldprune:
            plink_ldprune(output_bfile, ld_window_size=ld_window_size, ld_step_size=ld_step_size,
                          ld_r2_threshold=ld_r2_threshold, output_bfile=output_bfile, memory=memory)

    # extract samples to vcf
    # double check whether ref allele reserved
    if output_vcf_prefix:
        subprocess_popen(['plink2', '--export', 'vcf', 'bgz', 'id-paste=iid',
                          '--out', output_vcf_prefix] + operations)


def plink_ldprune(bfile: str, ld_window_size: str = '50', ld_step_size: str = '5', ld_r2_threshold: str = '0.2',
                  output_bfile: str = None, memory: str = None) -> None:
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
        output_bfile = bfile + '.ldprune'
    prunefile_fn = generate_tmpfn()
    cmd1 = ['plink', '--bfile', bfile, '--indep-pairwise',
            ld_window_size, ld_step_size, ld_r2_threshold,
            '--out', prunefile_fn]
    cmd2 = ['plink', '--bfile', bfile, '--extract',
            prunefile_fn + '.prune.in', '--memory', memory, '--keep-allele-order',
            '--make-bed', '--out', output_bfile]
    if memory:
        cmd1 += ['--memory', memory]
        cmd2 += '--memory', memory
    basic_fun.subprocess_popen(cmd1)
    basic_fun.subprocess_popen(cmd2)


def update_bim_snpname(bim_fn: str) -> None:
    '''
    Update bim file snpname to chr_pos_a2_a1 (a2 a1 are ref alt if keep-allele-order from vcf)

    Parameters
    ----------
    bim_fn: filename of bim
    '''
    bim = pd.read_table(bim_fn, header=None, names=['chr', 'snp', 'genetic', 'pos', 'a1', 'a2'])  # a2 a1 are ref alt
    bim['snp'] = bim['chr'].astype('str') + '_' + bim['pos'].astype('str') + '_' + bim['a2'] + '_' + bim['a1']
    bim.to_csv(bim_fn, sep='\t', index=False, header=False)


def grm(bfile: str, rel: Optional[str]=None, chr: Optional[int]=None, start: Optional[Union[int, List[int]]]=None, 
        end: Optional[Union[int, List[int]]]=None, r: int=0, snps: Optional[list]=None, 
        tool: str = 'plink', format: str = 'bin', nsnp_only: Optional[bool]=False
        ) -> int:
    """
    Compute kinship matrix for a genomic region (start-r, end+r)

    Parameters:
        bfile:  prefix for chr/genome bed/bim files
        rel:    prefix for relationship matrix file (prefix.rel.bin for plink, prefix.grm.bin for gcta)
        chr:    chromosome
        start:  start position of gene
        end:    end position of gene
        r:  radius to the gene
        snps:   SNPs to calculate kinship
        tool:   plink or gcta to compute grm
        format: output format for plink
    Returns:
        number of snps in the regions
    """
    if snps is None:
        bim = pd.read_csv(bfile + '.bim', sep='\s+', names=['chr', 'snp', 'cm', 'bp', 'a1', 'a2'])

        if isinstance(start, int) and isinstance(end, int):
            start = max(0, start - r)
            end = end + r

            # check number of SNPs in the region
            snps = bim.loc[(bim['chr'] == chr) & (bim['bp'] >= start) & (bim['bp'] <= end), 'snp'].tolist()
        elif isinstance(start, list) and isinstance(end, list):
            snps = []
            for x, y in zip(start, end):
                x = max(0, x - r)
                y = y + r

                # check number of SNPs in the region
                snps += bim.loc[(bim['chr'] == chr) & (bim['bp'] >= x) & (bim['bp'] <= y), 'snp'].tolist()
        else:
            log.logger.info('Wrong input!')
            print(start, end)
            sys.exit(-1)

    nsnp = len(snps)

    if nsnp_only:
        return nsnp

    tmpdir = tempfile.mkdtemp()
    tmp = os.path.join(tmpdir, 'tmp')

    snp_f = tmp + '.snps'
    with open(snp_f, 'w') as f:
        f.write('\n'.join(snps))

    if rel is None:
        rel = tmp + '.rel'

    # compute kinship matrix
    if nsnp > 0:
        if tool == 'plink':
            if format == 'bin':
                cmd = ['plink', '--bfile', bfile, '--extract', snp_f,
                    '--make-rel', 'bin',
                    '--out', rel]
            elif format == 'gz':
                cmd = ['plink', '--bfile', bfile, '--extract', snp_f,
                    '--make-rel', 'gz',
                    '--out', rel]
            subprocess_popen(cmd)
        elif tool == 'gcta':
            tmp = generate_tmpfn()
            cmd = ['plink', '--bfile', bfile, '--extract', snp_f,
                   '--make-bed', '--out', tmp]
            subprocess_popen(cmd)
            cmd = ['gcta', '--bfile', tmp,
                   '--make-grm', '--out', rel]
            subprocess_popen(cmd)
    
    shutil.rmtree(tmpdir)

    return nsnp


def transform_grm(grm: np.ndarray) -> np.ndarray:
    """
    Transform 1d grm lower triangle to 2d

    grm:    1d grm lower triangle
    """

    # Calculate the size of the square matrix based on the length of the 1D array
    n = int(np.sqrt(len(grm) * 2))

    # Create a symmetric 2D matrix
    grm2 = np.zeros((n, n), dtype=grm.dtype)

    # Fill the lower triangular part, including the diagonal
    row, col = np.tril_indices(n)
    grm2[row, col] = grm

    # Reflect the lower triangular part to the upper triangular part
    grm2[col, row] = grm

    return grm2


def sort_grm(grm: npt.NDArray[np.float64], old_order: list, new_order: list) -> np.ndarray:
    """
    Sort grm accordding to the order of individuals in new_order

    grm:    2d grm
    old_order:  old order of individuals
    new_order:  new order of individuals
    """

    new_grm = pd.DataFrame(grm, index=old_order, columns=old_order)

    return new_grm.loc[new_order, new_order].to_numpy()
    

def grm_matrix2gcta_grm_gz(grm: np.ndarray, grm_prefix: str, 
                        inds: Optional[list]=None, 
                        nsnp: Union[int, List[int]]=10,) -> None:
    """
    Write grm matrix to gcta grm.gz file

    Parameters:
        grm:    2d grm
        grm_prefix:  prefix for grm file
        inds:   order of individuals in input grm matrix;
                if None, names individuals as ind1, ind2, ind3, etc.
        nsnp:   number of snps to write to grm file
    """

    if inds is None:
        inds = [f'ind{i+1}' for i in range(grm.shape[0])]
    np.savetxt(grm_prefix + '.grm.id', np.array([inds, inds]).T, fmt='%s')

    # Get the lower triangle indices
    indices = np.tril_indices(grm.shape[0])

    # Create a 1D array from the lower triangle elements
    lower_triangle = grm[indices]

    # Merge lower triangle elements with indices
    if isinstance(nsnp, int):
        nsnp = np.ones_like(lower_triangle) * nsnp

    new_grm = np.column_stack((indices[0] + 1, indices[1] + 1, nsnp, lower_triangle))
    np.savetxt(grm_prefix + '.grm.gz', new_grm, fmt='%d %d %d %.6f', delimiter='\t')

    # return inds
    return inds


def read_greml(f: str) -> dict:
    """
    Read GREML output file

    Parameters:
        f:  GREML output file

    Returns:
        GREML results
    """
    
    out = {}
    for line in open(f):
        if re.search(r'V\(', line):
            line = line.strip().split()
            out[line[0]] = {'variance': float(line[1]), 'se': float(line[2])}

    return out
            

def construct_design(inds: list, pca: Optional[pd.DataFrame] = None, PC: Optional[int] = None, cat: Optional[pd.Series] = None,
           con: Optional[pd.Series] = None, drop_first: bool = True) -> np.ndarray:
    """
    Construct design matrix

    Parameters:
        inds:   order of individuals
        pca:    dataframe of pcs, with index: individuals (sort not required) and columns (PC1-PCx)
        PC: number to PC to adjust
        cat:    series of category elements e.g. sex: male and female
        con:    series of continuous elements e.g. age
        drop_first: drop the first column

    Returns:
        a design matrix
    """

    # pca
    if pca is not None and PC is not None:
        pcs = [f'PC{i}' for i in range(1, PC + 1)]
        return pca.loc[inds, pcs].to_numpy()
    elif cat is not None:
        return pd.get_dummies(cat, drop_first=drop_first, dtype=float).loc[inds, :].to_numpy()
    elif con is not None:
        return con[inds].to_numpy()
    else:
        return np.zeros((len(inds), 1))


def L_f(C: int, c1: int, c2: int) -> np.ndarray:
    # to build L matrix of 0 and 1
    L = np.zeros((C, C), dtype='int8')
    L[c1, c2] = 1
    return L


def compute_h2_pergene(hom_g2: float, V: Union[float, np.ndarray], hom_e2: float, 
                       W: Union[float, np.ndarray]) -> Tuple[float, float]:
    '''
    Compute h2 for one gene
    '''
    shared_h2 = None
    specific_h2 = None
    
    if isinstance(V, float):
        shared_h2 = hom_g2 / (hom_g2 + V + hom_e2 + W)
        specific_h2 = V / (hom_g2 + V + hom_e2 + W)
    else:
        mean_V = np.diag(V).mean()
        mean_W = np.diag(W).mean()

        shared_h2 = hom_g2 / (hom_g2 + mean_V + hom_e2 + mean_W)
        specific_h2 = mean_V / (hom_g2 + mean_V + hom_e2 + mean_W)

    return shared_h2, specific_h2


def compute_h2(hom_g2: np.ndarray, V: np.ndarray, hom_e2: np.ndarray, 
               W: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute h2 across transcriptome
    """

    mean_V = np.diagonal(V, axis1=1, axis2=2).mean(axis=1)
    mean_W = np.diagonal(W, axis1=1, axis2=2).mean(axis=1)

    # calculate heritability
    sig2s = hom_g2 + mean_V + hom_e2 + mean_W

    # ct specific
    shared_h2s = hom_g2 / sig2s
    specific_h2s = mean_V / sig2s

    return shared_h2s, specific_h2s


def read_ensembl(f: str) -> pd.DataFrame:
    res = {'gene':[], 'ensembl':[], 'hgnc':[], 'chr':[], 'start':[], 'end':[]}
    for line in gzip.open(f, 'rt'):
        if line[0] != '#':
            line = line.strip().split('\t')
            if line[2] in ['gene', 'pseudogene', 'RNA', 'lincRNA_gene', 'snRNA_gene', 'miRNA_gene', 'processed_transcript'] and line[0].isdigit():
                try:
                    chr, start, end, info = int(line[0]), int(line[3]), int(line[4]), line[-1]
                    # get hgnc id
                    if '[Source:HGNC Symbol' in info:
                        hgnc = info.split('[Source:HGNC Symbol')[1]
                        hgnc = hgnc.split(']')[0]
                        hgnc = hgnc.split(':')[1]
                    else:
                        hgnc = 'NA'
                    res['hgnc'].append(hgnc)

                    info = info.split(';')
                    tmp = {}
                    for x in info:
                        x = x.split('=')
                        tmp[x[0]] = x[1]
                    info = tmp
                    res['gene'].append(info['Name'])
                    if 'gene_id' in info.keys():
                        res['ensembl'].append(info['gene_id'])
                    else:
                        res['ensembl'].append('NA')
                    res['chr'].append(chr)
                    res['start'].append(start) 
                    res['end'].append(end)
                except:
                    log.logger.info('\t'.join(line))
                    sys.exit()

    res = pd.DataFrame(res)

    # sanity check duplicated name or id
    duplicated_gene = (res['gene'].count() != res['gene'].nunique())
    duplicated_ensembl = (res['ensembl'].count() != res['ensembl'].nunique())
    duplicated_hgnc = (res['hgnc'].count() != res['hgnc'].nunique())
    if duplicated_gene:
        genes, counts = np.unique(res['gene'], return_counts=True)
        print('Duplicated genes:', genes[counts > 1])
    if duplicated_ensembl:
        ids, counts = np.unique(res['ensembl'], return_counts=True, equal_nan=False)
        print('Duplicated Ensembl ID:', ids[counts > 1])
    if duplicated_hgnc:
        ids, counts = np.unique(res['hgnc'], return_counts=True, equal_nan=False)
        print('Duplicated HGNC ID:', ids[counts > 1])

    return res


def h2_equal_test(formula: List[str], var: Union[np.ndarray, list], cov: np.ndarray, 
                  ct_h2: np.ndarray) -> float:
    '''
    Delta method to test equal h2 across cell types

    Parameters:
        formula:    a formula mapping from variance to h2, e.g. '~(x1)/(x1+x2)'
                    x1, x2 indicate the first and second element in mean
        var:    estimates of variance parameters, e.g. genetic variance, environment variance
        cov:    estiamted covariance matrix of variance parameters

    Returns:
        p value for ct-specific h2
    '''

    # convert string to R list of formula
    formula_fun = r('function(x) as.formula(x)')
    formula = ro.StrVector(formula)
    formula = r['lapply'](formula, formula_fun)

    numpy2ri.activate()

    h2Covmat = r('msm::deltamethod')(formula, var, cov, ses=False)

    numpy2ri.deactivate()

    # test h2
    p = wald_ct_beta(ct_h2, h2Covmat)  # NOTE: using chi-square test

    return p


def _op_partition_mean(beta: np.ndarray, S: np.ndarray, ) -> float:
    """
    Compute OP partition of cell type fixed effect

    Parameters:
        beta:   estimated cell type means
        S:      estimated covariance of cell type proportions

    Returns:
        OP partition of cell type fixed effect
    """
    return beta @ S @ beta


def _op_partition_specific_effect(X: np.ndarray, S: np.ndarray, pi: np.ndarray, ) -> float:
    """
    Compute OP partition of cell type specific genetic variance

    Parameters:
        X:      estimated genetic variance (V) or environmental variance (W)
        S:      estimated covariance of cell type proportions
        pi:     estimated mean cell type proportions

    Returns:
        OP partition of cell type specific genetic variance
    """
    return np.trace(X @ S) + pi @ X @ pi


def op_partition(beta: np.ndarray, V: np.ndarray, W: np.ndarray, S: np.ndarray, 
                 pi: np.ndarray, ) -> dict:
    """
    Comptued variance partition of OP for cell type fixed effect, 
    ct-specific genetic, and ct-specific environmental variance

    Parameters:
        beta:   estimated cell type means
        V:      estimated genetic variance
        W:      estimated environmental variance
        S:      estimated covariance of cell type proportions
        pi:     estimated mean cell type proportions

    Returns:
        a dictionary of variance components
    """
    return {'ct_fixed_effect': _op_partition_mean(beta, S),
            'ct_specific_gen': _op_partition_specific_effect(V, S, pi),
            'ct_specific_env': _op_partition_specific_effect(W, S, pi)}
