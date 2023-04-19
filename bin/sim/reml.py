import time, re, os, logging
import numpy as np, pandas as pd

from gxctmm import fit, util, log 

def main():
    #
    out_fs = []
    for i in snakemake.params.batches:
        out_f = re.sub('/rep/', f'/rep{i}/', snakemake.params.out)
        out_fs.append(out_f)

    Y_fs = [line.strip() for line in open(snakemake.input.Y)]
    K_fs = [line.strip() for line in open(snakemake.input.K)]
    P_fs = [line.strip() for line in open(snakemake.input.P)]
    ctnu_fs = [line.strip() for line in open(snakemake.input.ctnu)]

    for i in range(len(Y_fs)):

        log.logger.info(f'{Y_fs[i]}')
        log.logger.info(f'{K_fs[i]}')
        log.logger.info(f'{ctnu_fs[i]}')
        log.logger.info(f'{P_fs[i]}')

        Y = np.loadtxt( Y_fs[i] )
        K = np.loadtxt( K_fs[i] )
        P = np.loadtxt( P_fs[i] )
        ctnu = np.loadtxt( ctnu_fs[i] )

        N, C = Y.shape

        if snakemake.wildcards.model != 'full':
            # fit hom
            hom = fit.hom_REML(Y, K, P, ctnu, method='BFGS-Nelder')
            # fit freeW
            freeW = fit.freeW_REML(Y, K, P, ctnu, method='BFGS-Nelder')
            # fit free
            free, free_p = fit.free_REML(Y, K, P, ctnu, method='BFGS-Nelder', jk=False)

            # LRT
            free_freeW = util.lrt(free['opt']['l'], freeW['opt']['l'], k=C)
            free_hom = util.lrt(free['opt']['l'], hom['opt']['l'], k=2*C)

            out = { 'hom': hom, 'free':freeW, 'free':free, 'wald':free_p,
                    'lrt': {'free_freeW':free_freeW, 'free_hom':free_hom}   }
        else:
            # fit full
            full = fit.full_REML(Y, K, P, ctnu, method='BFGS-Nelder', nrep=2)

            out = { 'full':full }

        os.makedirs( os.path.dirname(out_fs[i]), exist_ok=True )
        np.save( out_fs[i], out )

    with open(snakemake.output.out, 'w') as f:
        f.write( '\n'.join(out_fs) )

if __name__ == '__main__':
    main()
