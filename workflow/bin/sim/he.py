import re, os, sys
import numpy as np, pandas as pd

from gxctmm import log, fit

def main():

    outs = []
    data = np.load(snakemake.input.data, allow_pickle=True).item()
    free_jk = snakemake.params.get('free_jk', True)

    for key in data.keys():
        log.logger.info(f'{key}')

        # read
        Y = data[key]['Y']
        K = data[key]['K']
        ctnu = data[key]['ctnu']
        P = data[key]['P']
        if 'fixed' in data[key].keys():
            fixed = {'fixed': data[key]['fixed'][:, :-1]}
        else:
            fixed = {}
        if 'random' in data[key].keys():
            random = {'random': data[key]['random']}
        else:
            random = {}

        if snakemake.wildcards.model != 'full':
            ## Free
            free_he, free_he_wald = fit.free_HE(Y, K, ctnu, P, fixed_covars=fixed, random_covars=random, jk=free_jk)

            # Full
            full_he = fit.full_HE(Y, K, ctnu, P, fixed_covars=fixed, random_covars=random)

            # save
            outs.append({'gene': key, 'free':free_he, 'full':full_he, 'wald': {'free':free_he_wald}})
            
        else:
            # Full
            full_he = fit.full_HE(Y, K, ctnu, P, fixed_covars=fixed, random_covars=random)

            # save
            outs.append({'gene': key, 'full':full_he})
        
    np.save(snakemake.output.out, outs)


if __name__ == '__main__':
    main()
