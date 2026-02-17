import time, re, os, logging
import numpy as np, pandas as pd

from cigma import fit, util, log 

def main():
    # par
    R = snakemake.params.get('R', True)
    chol = snakemake.params.get('chol', True)

    #
    outs = []
    data = np.load(snakemake.input.data, allow_pickle=True).item()

    for key in data.keys():
        log.logger.info(f'{key}')

        # read
        Y = data[key]['Y']
        K = data[key]['K']
        ctnu = data[key]['ctnu']
        P = data[key]['P']
        N, C = Y.shape
        if 'fixed' in data[key].keys():
            fixed = {'fixed': data[key]['fixed'][:, :-1]}
        else:
            fixed = {}
        if 'random' in data[key].keys():
            random = {'batch': data[key]['random']}  # here change random to batch for REML in R
        else:
            random = {}

        if snakemake.wildcards.model != 'full':
            # fit free
            free, free_p = fit.free_REML(Y, K, P, ctnu, fixed_covars=fixed, 
                                        random_covars=random, 
                                        method='BFGS', 
                                        jk=False, nrep=2, 
                                        chol=chol, R=R)

            outs.append({'gene': key, 'free': free, 'wald': free_p})

        else:
            # fit full
            full = fit.full_REML(Y, K, P, ctnu, fixed_covars=fixed, 
                                 random_covars=random, method='BFGS-Nelder', 
                                 nrep=2)

            outs.append({'gene': key, 'full': full})

    np.save(snakemake.output.out, outs)


if __name__ == '__main__':
    main()
