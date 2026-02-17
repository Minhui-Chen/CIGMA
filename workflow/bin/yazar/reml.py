import time, re, os, logging, sys
import numpy as np, pandas as pd

from cigma import fit, util, log 

def main():
    #
    geno_pca_n = int(snakemake.wildcards.get('geno_pca_n', '6'))
    op_pca_n = int(snakemake.wildcards.get('op_pca_n', '1'))
    batch = snakemake.wildcards.get('batch', 'shared')
    fixed = snakemake.wildcards.get('fixed', 'shared')
    batch_shared = True if batch == 'shared' else False
    fixed_shared = True if fixed == 'shared' else False
    log.logger.info(f'Geno PC: {geno_pca_n}')
    log.logger.info(f'OP PC: {geno_pca_n}')
    log.logger.info(f'{batch} batch effect')
    log.logger.info(f'{fixed} fixed effect')

    # read
    ctps = pd.read_table(snakemake.input.ctp, index_col=(0, 1)).astype('float32')
    ctnus = pd.read_table(snakemake.input.ctnu, index_col=(0, 1)).astype('float32')
    P = pd.read_table(snakemake.input.P, index_col=0)
    inds = P.index  # order of individuals
    cts = P.columns  # order of cts
    genes = ctps.columns
    P = P.to_numpy().astype('float32')

    kinship = np.load(snakemake.input.kinship, allow_pickle=True).item()

    # check ctp, ctnu, P have the same order
    gene = genes[0]
    ctp = ctps[gene].unstack()
    ctnu = ctnus[gene].unstack()
    assert (ctp.index.equals(inds) and ctnu.index.equals(inds)
            and ctp.columns.equals(cts) and ctnu.columns.equals(cts))

    # collect covariates
    fixed_covars, random_covars = util.yazar_covars(inds.to_list(), snakemake.input.obs,
                            snakemake.input.geno_pca, snakemake.input.op_pca, 
                            geno_pca_n=geno_pca_n, op_pca_n=op_pca_n)


    outs = []
    # data = np.load(snakemake.input.data, allow_pickle=True).item()

    for gene in genes:
        log.logger.info(f'Fitting {gene}')

        # extract gene data
        ctp = ctps[gene].unstack().to_numpy().astype('float32')
        ctnu = ctnus[gene].unstack().to_numpy().astype('float32')
        gene_idx = np.nonzero(kinship['gene'] == gene)[0]
        # sanity check
        if len(gene_idx) == 0:
            continue
        elif len(gene_idx) > 1:
            sys.exit('Duplicate gene!')
        gene_idx = gene_idx[0]

        # if 'fixed' in data[key].keys():
        #     fixed = {'fixed': data[key]['fixed'][:, :-1]}
        # else:
        #     fixed = {}
        # if 'random' in data[key].keys():
        #     random = {'batch': data[key]['random']}  # here change random to batch for REML in R
        # else:
        #     random = {}

        # prepare Kinship matrix
        if kinship['nsnp'][gene_idx] <= snakemake.params.snps:
            continue
        else:
            K = util.transform_grm(kinship['K'][gene_idx])
            # sort K
            K = util.sort_grm(K, kinship['ids'], inds.to_list()).astype('float32')


        out = {}
        out['gene'] = gene
        ## Free
        free, free_p = fit.free_REML(ctp, K, P, ctnu, fixed_covars=fixed_covars, 
                                    random_covars=random_covars, 
                                    method='BFGS', 
                                    jk=False, nrep=2, 
                                    chol=True, R=True)

        out['free'] = free
        if len(free_p) != 0:
            out['p'] = {'free': free_p}

        outs.append(out)
        break

    np.save(snakemake.output.out, outs)


if __name__ == '__main__':
    main()
