import re, os, sys
import numpy as np, pandas as pd

from gxctmm import log, fit

def main():

    batch = snakemake.params.batches
    outs = [re.sub('/rep/', f'/rep{i}/', snakemake.params.out) for i in batch]
    for Y_f, K_f, ctnu_f, P_f, out_f in zip(
            [line.strip() for line in open(snakemake.input.Y)],
            [line.strip() for line in open(snakemake.input.K)],
            [line.strip() for line in open(snakemake.input.ctnu)], 
            [line.strip() for line in open(snakemake.input.P)],
            outs):

        log.logger.info(f'{Y_f}\n{K_f}\n{ctnu_f}\n{P_f}\n')

        os.makedirs(os.path.dirname(out_f), exist_ok=True)
        
        # read
        Y = np.loadtxt( Y_f )
        K = np.loadtxt( K_f )
        ctnu = np.loadtxt( ctnu_f )
        P = np.loadtxt( P_f )

        if snakemake.wildcards.model != 'full':
            ## Free
            free_he, free_he_wald = fit.free_HE(Y, K, ctnu, P)

            # Full
            full_he = fit.full_HE(Y, K, ctnu, P)

            # save
            np.save(out_f,
                    {
                        'free':free_he, 
                        'full':full_he, 
                        'wald': {'free':free_he_wald} 
                    }
                )
        else:
            # Full
            full_he = fit.full_HE(Y, K, ctnu, P)

            # save
            np.save(out_f,
                    {
                        'full':full_he, 
                    }
                )
        
    with open(snakemake.output.out, 'w') as f:
        f.write('\n'.join(outs))  


if __name__ == '__main__':
    main()
