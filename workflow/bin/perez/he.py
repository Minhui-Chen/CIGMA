import re, os, sys
import numpy as np, pandas as pd
from cigma import log, fit, util


def main():
    jk = snakemake.params.get('jk', False)
    prop = float(snakemake.wildcards.get('replicates', 1))

    # read
    ctps = pd.read_table(snakemake.input.ctp, index_col=(0, 1)).astype('float32')
    ctnus = pd.read_table(snakemake.input.ctnu, index_col=(0, 1)).astype('float32')
    P = pd.read_table(snakemake.input.P, index_col=0)
    inds = P.index  # order of individuals
    cts = P.columns  # order of cts
    genes = ctps.columns
    P = P.to_numpy().astype('float32')

    kinship = np.load(snakemake.input.kinship, allow_pickle=True).item()
    if 'promoter' in snakemake.input.keys():
        promoter_kinship = np.load(snakemake.input.promoter, allow_pickle=True).item()
        # sanity check gene order
        if np.any(kinship['gene'] != promoter_kinship['gene']):
            sys.exit('Mismatching gene order')
    if 'genome' in snakemake.input.keys():
        if isinstance(snakemake.input.genome, str):
            genome = []
            for line in open(snakemake.input.genome):
                genome += line.strip().split()
            genome = np.array(genome)
            genome = util.transform_grm(genome)
            genome_ids = [line.strip().split()[0] for line in open(snakemake.input.genome + '.id')]
            genome = util.sort_grm(genome, genome_ids, inds.to_list()).astype('float32')
        elif isinstance(snakemake.input.genome, list):
            gene_location = pd.read_table(snakemake.input.genes)

            genomes = {}
            for chr, genome_f in zip(snakemake.params.chrs, snakemake.input.genome):
                genome = []
                for line in open(genome_f):
                    genome += line.strip().split()
                genome = np.array(genome)
                genome = util.transform_grm(genome)
                genome_ids = [line.strip().split()[0] for line in open(genome_f + '.id')]
                genome = util.sort_grm(genome, genome_ids, inds.to_list()).astype('float32')
                genomes[chr] = genome


    # check ctp, ctnu, P have the same order
    gene = genes[0]
    ctp = ctps[gene].unstack()
    ctnu = ctnus[gene].unstack()
    if not (ctp.index.equals(inds) and ctnu.index.equals(inds)
            and ctp.columns.equals(cts) and ctnu.columns.equals(cts)):
        sys.exit('Inds or CTs order not matching!\n')

    # collect covariates
    geno_pca_n = 5
    if snakemake.wildcards.anc == 'European' and snakemake.wildcards.status == 'Healthy':
        geno_pca_n = 3
    elif snakemake.wildcards.anc == 'European' and snakemake.wildcards.status == 'SLE':
        geno_pca_n = 4
    elif snakemake.wildcards.anc == 'Asian' and snakemake.wildcards.status == 'Healthy':
        geno_pca_n = 4
    elif snakemake.wildcards.anc == 'Asian' and snakemake.wildcards.status == 'SLE':
        geno_pca_n = 4

    include_dataset = False
    if snakemake.wildcards.anc == 'all' or snakemake.wildcards.status == 'all':
        include_dataset = True
    fixed_covars, random_covars = util.perez_covars(inds.to_list(), snakemake.input.obs,
                            snakemake.input.geno_pca, snakemake.input.op_pca, 
                            batch='shared', geno_pca_n=geno_pca_n, op_pca_n=10,
                            include_dataset=include_dataset)

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

        if 'genome' in snakemake.input.keys():
            if isinstance(snakemake.input.genome, list):
                ichr = gene_location.loc[gene_location['feature'] == gene, 
                                         'chr'].to_numpy()
                if len(ichr) > 1:
                    sys.exit('Wrong chrosome')
                ichr = ichr[0]
                genome = genomes[ichr]

        if kinship['nsnp'][gene_idx] <= snakemake.params.snps:
            continue
        else:
            K = util.transform_grm(kinship['K'][gene_idx])
            # sort K
            K = util.sort_grm(K, kinship['ids'], inds.to_list()).astype('float32')


        ## Free
        free_he, free_he_wald = {}, {}
        if ('genome' not in snakemake.input.keys()) and ('promoter' not in snakemake.input.keys()):
            log.logger.info('Fitting cis...')
            free_he, free_he_wald = fit.free_HE(ctp, K, ctnu, P, fixed_covars, 
                                                random_covars, jk=jk, prop=prop, 
                                                dtype='float32')

        elif 'genome' in snakemake.input.keys():
            # cis vs trans
            log.logger.info('Fitting cis + trans...')
            free_he, free_he_wald = fit.free_HE(ctp, K, ctnu, P, fixed_covars, 
                                                random_covars, Kt=genome, jk=jk, 
                                                prop=prop, dtype='float32')

        elif 'promoter' in snakemake.input.keys():
            # promoter vs enhancer
            log.logger.info('Fitting enhancer + promoter...')
            # doesn't require min nsnp
            promoter_K = util.transform_grm(promoter_kinship['K'][gene_idx])
            promoter_K = util.sort_grm(promoter_K, promoter_kinship['ids'], 
                                       inds.to_list()).astype('float32')
            if np.all(K == promoter_K):
                log.logger.info(f'Equal GRM')
                continue
            free_he, free_he_wald = fit.free_HE(ctp, K, ctnu, P, fixed_covars, 
                                                random_covars, Kt=promoter_K, 
                                                jk=jk, prop=prop, dtype='float32')
            free_he['cis_hom_g2_perSNP'] = free_he['cis_hom_g2'] / kinship['nsnp'][gene_idx]
            free_he['cis_V_perSNP'] = free_he['cis_V'] / kinship['nsnp'][gene_idx]
            free_he['trans_hom_g2_perSNP'] = free_he['trans_hom_g2'] / promoter_kinship['nsnp'][gene_idx]
            free_he['trans_V_perSNP'] = free_he['trans_V'] / promoter_kinship['nsnp'][gene_idx]

        # save
        out = {'gene': gene, 'free': free_he}
        if len(free_he_wald) != 0:
            out['p'] = {'free': free_he_wald}

        # iid
        if snakemake.params.get('iid', False):
            if 'genome' in snakemake.input.keys():
                iid_he = fit.iid_HE(ctp, K, ctnu, P, fixed_covars, random_covars, Kt=genome, dtype='float32')
            else:
                iid_he = fit.iid_HE(ctp, K, ctnu, P, fixed_covars, random_covars, dtype='float32')
            out['iid'] = iid_he

        # Full
        if snakemake.params.get('full', False):
            full_he = fit.full_HE(ctp, K, ctnu, P, fixed_covars, random_covars, dtype='float32')
            out['full'] = full_he

        outs.append(out)

    np.save(snakemake.output.out, outs)


if __name__ == '__main__':
    main()
