import re, os, sys
import numpy as np, pandas as pd
from cigma import log, fit, util


def main():
    jk = snakemake.params.get('jk', False)
    geno_pca_n = int(snakemake.wildcards.get('geno_pca_n', '6'))
    op_pca_n = int(snakemake.wildcards.get('op_pca_n', '1'))
    log.logger.info(f'Geno PC: {geno_pca_n}')
    if 'replicates' in snakemake.wildcards.keys():
        prop = float(snakemake.wildcards.replicates)
    elif 'replicates' in snakemake.params.keys():
        prop = float(snakemake.params.replicates)
    else:
        prop = 1
    # batch = snakemake.params.get('batch', 'shared')
    # batch_shared = True if batch == 'shared' else False

    # read
    ctps = pd.read_table(snakemake.input.ctp, index_col=(0, 1)).astype('float32')
    ctnus = pd.read_table(snakemake.input.ctnu, index_col=(0, 1)).astype('float32')
    P = pd.read_table(snakemake.input.P, index_col=0)
    kinship = np.load(snakemake.input.kinship, allow_pickle=True).item()
    permuted_inds = np.loadtxt(snakemake.input.ids, dtype=str).tolist()

    # collect info and transform
    inds = P.index  # order of individuals
    cts = P.columns  # order of cts
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
    fixed_covars, random_covars = util.yazar_covars(inds.to_list(), snakemake.input.obs,
                            snakemake.input.geno_pca, snakemake.input.op_pca, 
                            geno_pca_n=geno_pca_n, op_pca_n=op_pca_n)
    

    # run
    outs = []
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

        if kinship['nsnp'][gene_idx] <= snakemake.params.snps:
            continue
        else:
            K = util.transform_grm(kinship['K'][gene_idx])
            # sort K
            K = util.sort_grm(K, kinship['ids'], permuted_inds).astype('float32')


        ## Free
        out = {'gene': gene}
        if snakemake.params.get('free', True):
            log.logger.info('Fitting cis.')
            free_he, free_he_wald = fit.free_HE(ctp, K, ctnu, P, fixed_covars, 
                                                random_covars, jk=jk, prop=prop, 
                                                dtype='float32')
            out['free'] = free_he
            if free_he_wald:
                if 'p' not in out.keys():
                    out['p'] = {}
                out['p']['free'] = free_he_wald

        ## IID
        if snakemake.params.get('iid', True):
            iid_he, iid_he_wald = fit.iid_HE(ctp, K, ctnu, P, fixed_covars, 
                                             random_covars, jk=jk, dtype='float32')
            out['iid'] = iid_he
            if iid_he_wald:
                if 'p' not in out.keys():
                    out['p'] = {}
                out['p']['iid'] = iid_he_wald


        outs.append(out)

    np.save(snakemake.output.out, outs)


if __name__ == '__main__':
    main()
