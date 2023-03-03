import time, re, os
import numpy as np, pandas as pd

from gxctmm import fit

def main():
    out_fs = []
    for i in snakemake.params.batches:
        out_f = re.sub('/rep/', f'/rep{i}/', snakemake.params.out)
        out_fs.append(out_f)

    Y_fs = [line.strip() for line in open(snakemake.input.Y)]
    K_fs = [line.strip() for line in open(snakemake.input.K)]
    P_fs = [line.strip() for line in open(snakemake.input.P)]
    ctnu_fs = [line.strip() for line in open(snakemake.input.ctnu)]

    for i in range(len(Y_fs)):
        Y = np.loadtxt( Y_fs[i] )
        Z = np.loadtxt( Z_fs[i] )
        P = np.loadtxt( P_fs[i] )
        ctnu = np.loadtxt( ctnu_fs[i] )

        free, free_p = fit.free_REML(Y, Z, P, ctnu, jk=False)

        out = { 'free':free, 'wald':free_p} }

        os.makedirs( os.path.dirname(out_fs[i]), exist_ok=True )
        np.save( out_fs[i], out )

    with open(snakemake.output.out, 'w') as f:
        f.write( '\n'.join(out_fs) )

if __name__ == '__main__':
    main()
