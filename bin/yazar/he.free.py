import re, os, sys
import numpy as np, pandas as pd
from gxctmm import log, fit, util


def main():
    jk = snakemake.params.get('jk', False)
    target = snakemake.params.get('target')

    # read
    ctps = pd.read_table(snakemake.input.ctp, index_col=(0, 1)).astype('float32')
    ctnus = pd.read_table(snakemake.input.ctnu, index_col=(0, 1)).astype('float32')
    P = pd.read_table(snakemake.input.P, index_col=0)
    kinship = np.load(snakemake.input.kinship, allow_pickle=True).item()

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
    op_pca = pd.read_table(snakemake.input.op_pca, index_col=0)
    geno_pca = pd.read_table(snakemake.input.geno_pca, index_col=0).drop('IID', axis=1)
    meta = pd.read_table(snakemake.input.obs, usecols=['individual', 'sex', 'age', 'pool'])
    meta = meta.drop_duplicates()
    meta = meta.set_index('individual')
    fixed_covars = {
        'op_pca': util.design(inds, pca=op_pca, PC=1).astype('float32'),
        'geno_pca': util.design(inds, pca=geno_pca, PC=6).astype('float32'),
        'sex': util.design(inds, cat=meta['sex']),
        'age': util.design(inds, cat=util.age_group(meta['age']))
    }


    if snakemake.params.get('batch', True):
        random_covars = {
            'batch': util.design(inds, cat=meta['pool'], drop_first=False)
        }
    else:
        random_covars = {}

    # run
    outs = []
    for gene in genes:
        log.logger.info(f'Fitting {gene}')

        # when only test one Target gene
        if target and gene != target:
            continue

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
            K = util.sort_grm(K, kinship['ids'], inds).astype('float32')


        ## Free
        if 'genome' not in snakemake.input.keys():
            free_he, free_he_wald = fit.free_HE(ctp, K, ctnu, P, fixed_covars, random_covars, jk=jk, dtype='float32')
        else:
            genome = []
            for line in open(snakemake.input.genome):
                genome += line.strip().split()
            genome = np.array(genome)
            genome = util.transform_grm(genome)
            genome_ids = [line.strip().split()[0] for line in open(snakemake.input.genome + '.id')]
            genome = util.sort_grm(genome, genome_ids, inds).astype('float32')

            free_he, free_he_wald = fit.free_HE(ctp, K, ctnu, P, fixed_covars, random_covars, Kt=genome, jk=jk, dtype='float32')

        # save
        out = {'gene': gene, 'free': free_he}
        if len(free_he_wald) != 0:
            out['p'] = {'free': free_he_wald}
        outs.append(out)

    np.save(snakemake.output.out, outs)


if __name__ == '__main__':
    main()
