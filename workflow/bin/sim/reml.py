import time, re, os, logging
import numpy as np, pandas as pd

from gxctmm import fit, util, log 

def main():
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

        if snakemake.wildcards.model != 'full':
            # fit hom
            #hom = fit.hom_REML(Y, K, P, ctnu, method='BFGS-Nelder')
            # fit freeW
            #freeW = fit.freeW_REML(Y, K, P, ctnu, method='BFGS-Nelder')
            # fit free
            free, free_p = fit.free_REML(Y, K, P, ctnu, method='BFGS-Nelder', jk=False)

            # LRT
            #free_freeW = util.lrt(free['opt']['l'], freeW['opt']['l'], k=C)
            #free_hom = util.lrt(free['opt']['l'], hom['opt']['l'], k=2*C)

            outs.append({'gene': key, 'free': free, 'wald': free_p})

            #out = { 'hom': hom, 'free':freeW, 'free':free, 'wald':free_p,
                    #'lrt': {'free_freeW':free_freeW, 'free_hom':free_hom}   }
        else:
            # fit full
            full = fit.full_REML(Y, K, P, ctnu, method='BFGS-Nelder', nrep=2)

            outs.append({'gene': key, 'full': full})

    np.save(snakemake.output.out, outs)


if __name__ == '__main__':
    main()
