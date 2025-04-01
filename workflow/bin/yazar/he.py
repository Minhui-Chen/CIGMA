import re, os, sys
import numpy as np, pandas as pd
from cigma import log, fit, util


def main():
    jk = snakemake.params.get('jk', False)
    target = snakemake.params.get('target')
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

    # Jackknife replicates
    if 'replicates' in snakemake.wildcards.keys():
        prop = float(snakemake.wildcards.replicates)
    elif 'replicates' in snakemake.params.keys():
        prop = float(snakemake.params.replicates)
    else:
        prop = 1

    # read
    ctps = pd.read_table(snakemake.input.ctp, index_col=(0, 1)).astype('float32')
    ctnus = pd.read_table(snakemake.input.ctnu, index_col=(0, 1)).astype('float32')
    P = pd.read_table(snakemake.input.P, index_col=0)
    inds = P.index  # order of individuals
    cts = P.columns  # order of cts
    genes = ctps.columns
    P = P.to_numpy().astype('float32')

    kinship = np.load(snakemake.input.kinship, allow_pickle=True).item()
    
    # prepare trans kinship
    if 'genome' in snakemake.input.keys():
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
            
    if 'promoter' in snakemake.input.keys():
        promoter_kinship = np.load(snakemake.input.promoter, allow_pickle=True).item()
        # sanity check gene order
        assert np.all(kinship['gene'] == promoter_kinship['gene'])

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

        # whether fit trans
        if 'genome' in snakemake.input.keys():
            ichr = gene_location.loc[gene_location['feature'] == gene, 'chr'].to_numpy()
            if len(ichr) > 1:
                sys.exit('Wrong chrosome')
            ichr = ichr[0]

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
        if snakemake.params.get('free', True):
            if ('genome' not in snakemake.input.keys()) and ('promoter' not in snakemake.input.keys()):
                # For fitting cis only
                log.logger.info('Fitting cis...')
                free_he, free_he_wald = fit.free_HE(ctp, K, ctnu, P, 
                    fixed_covars, random_covars, fixed_shared=fixed_shared, 
                    random_shared=batch_shared, jk=jk, prop=prop, dtype='float32')
            elif 'genome' in snakemake.input.keys():
                # For jointly fitting cis and trans
                log.logger.info('Jonitly fitting cis and trans')
                free_he, free_he_wald = fit.free_HE(ctp, K, ctnu, P, 
                    fixed_covars, random_covars, Kt=genomes[ichr], 
                    fixed_shared=fixed_shared, random_shared=batch_shared, 
                    jk=jk, prop=prop, dtype='float32')
            elif 'promoter' in snakemake.input.keys():
                # For jointly fitting promoter and enhancer
                log.logger.info('Jonitly fitting promoter and enhancer')
                if promoter_kinship['nsnp'][gene_idx] <= snakemake.params.snps:
                    continue
                else:
                    promoter_K = util.transform_grm(promoter_kinship['K'][gene_idx])
                    promoter_K = util.sort_grm(promoter_K, promoter_kinship['ids'], inds.to_list()).astype('float32')
                    if np.all(K == promoter_K):
                        log.logger.info(f'Equal GRM')
                        continue
                    free_he, free_he_wald = fit.free_HE(ctp, K, ctnu, P, fixed_covars, random_covars, Kt=promoter_K, jk=jk, prop=prop, dtype='float32')
                    free_he['hom_g2_perSNP'] = free_he['hom_g2'] / kinship['nsnp'][gene_idx]
                    free_he['V_perSNP'] = free_he['V'] / kinship['nsnp'][gene_idx]
                    free_he['hom_g2_perSNP_b'] = free_he['hom_g2_b'] / promoter_kinship['nsnp'][gene_idx]
                    free_he['V_perSNP_b'] = free_he['V_b'] / promoter_kinship['nsnp'][gene_idx]
                    free_he['nsnp'] = kinship['nsnp'][gene_idx]
                    free_he['nsnp_b'] = promoter_kinship['nsnp'][gene_idx]
                    if 'nenhancer' in kinship.keys():
                        free_he['nenhancer'] = kinship['nenhancer'][gene_idx]
                        free_he['nenhancer_b'] = promoter_kinship['nenhancer'][gene_idx]

            # save
            out['free'] = free_he
            if len(free_he_wald) != 0:
                out['p'] = {'free': free_he_wald}

        # IID
        if snakemake.params.get('iid', False):
            if 'genome' in snakemake.input.keys():
                iid_he, iid_he_wald = fit.iid_HE(ctp, K, ctnu, P, fixed_covars, 
                                                 random_covars, Kt=genomes[ichr],
                                                 jk=jk, dtype='float32')
            else:
                iid_he, iid_he_wald = fit.iid_HE(ctp, K, ctnu, P, fixed_covars, 
                                                 random_covars, jk=jk, 
                                                 dtype='float32')
            out['iid'] = iid_he
            if len(iid_he_wald) != 0:
                if 'p' not in out.keys():
                    out['p'] = {}
                out['p']['iid'] = iid_he_wald

        # Full
        if snakemake.params.get('full', False):
            full_he = fit.full_HE(ctp, K, ctnu, P, fixed_covars, random_covars, dtype='float32')
            out['full'] = full_he

        outs.append(out)

    np.save(snakemake.output.out, outs)


if __name__ == '__main__':
    main()
