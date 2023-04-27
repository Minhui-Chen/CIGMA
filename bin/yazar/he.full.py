import re, os, sys
import numpy as np, pandas as pd
from pandas_plink import read_rel
from memory_profiler import profile 
from gxctmm import log, fit, util

def main():
    # read
    ctps = pd.read_table(snakemake.input.ctp, index_col=(0,1)).astype('float32')
    ctnus = pd.read_table(snakemake.input.ctnu, index_col=(0,1)).astype('float32')
    P = pd.read_table(snakemake.input.P, index_col=0)
    kinship = pd.read_table(snakemake.input.kinship)

    # collect info and transform
    inds = P.index # order of individuals
    cts = P.columns # order of cts
    genes = ctps.columns
    P = P.to_numpy().astype('float32')

    # check ctp, ctnu, P have the same order
    gene = genes[0]
    ctp = ctps[gene].unstack()
    ctnu = ctnus[gene].unstack()
    if not (ctp.index.equals(inds) and ctnu.index.equals(inds) 
            and ctp.columns.equals(cts) and ctnu.columns.equals(cts)):
        sys.exit('Inds or CTs order not matching!\n')

    # collect covariates
    op_pca = pd.read_table(snakemake.input.op_pca, index_col=0)
    geno_pca = pd.read_table(snakemake.input.geno_pca, index_col=0).drop('IID',axis=1)
    meta = pd.read_table(snakemake.input.meta, usecols=['individual', 'sex', 'age'])
    meta = meta.drop_duplicates()
    meta = meta.set_index('individual')
    fixed_covars = {'op_pca': util.design(inds, pca=op_pca, PC=1).to_numpy().astype('float32'),
            'geno_pca': util.design(inds, pca=geno_pca, PC=6).to_numpy().astype('float32'),
            'sex': util.design(inds, cat=meta['sex']).to_numpy(),
            'age': util.design(inds, con=meta['age']).to_numpy().astype('float32')
            }

    # run
    outs = []
    for gene in genes:
        log.logger.info(f'Fitting {gene}')

        out_f = re.sub('/rep/', f'/{gene}/', snakemake.params.out)
        os.makedirs(os.path.dirname(out_f), exist_ok=True)
        
        # extract gene data
        ctp = ctps[gene].unstack().to_numpy().astype('float32')
        ctnu = ctnus[gene].unstack().to_numpy().astype('float32')
        kinship_gene = kinship.loc[kinship['gene']==gene]
        if kinship_gene.shape[0] == 0:
            continue
        if kinship_gene['snps'].iloc[0] <= snakemake.params.snps:
            continue
        else:
            K = read_rel( kinship_gene['K'].iloc[0] )
            # sort K
            K = K.loc[inds, inds].to_numpy().astype('float32')
        
        # Full
        full_he = fit.full_HE(ctp, K, ctnu, P, fixed_covars, dtype='float32')

        # save
        np.save(out_f,
                {
                    'full':full_he, 
                }
            )
        outs.append( out_f )
        
    with open(snakemake.output.out, 'w') as f:
        f.write('\n'.join(outs))  


if __name__ == '__main__':
    main()
