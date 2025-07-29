######################################################################
# Rare variant analysis for Yazar et al. 2022
######################################################################
rule yazar_rare_he_kinship:
    input:
        ctp = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctp.final.gz',
        genes = 'data/Yazar2022Science/gene_location.txt',
        bed = 'analysis/yazar/data/geno/chr{chr}.bed',
    output:
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/kinship.chr{{chr}}.maf{{maf}}.npy',
    params:
        r = int(float(config['yazar']['radius'])),
    resources:
        mem_mb = lambda wildcards: '20G' if wildcards.chr != '1' else '80G',
    shell: 
        '''
        # load plink 1.9
        if [[ $(hostname) == *midway* ]]; then
            module load plink
        else
            module load gcc/11.3.0 atlas/3.10.3 lapack/3.11.0 plink/1.9
        fi

        python3 workflow/bin/yazar/kinship.npy.maf.py {input.genes} {input.ctp} {params.r} {input.bed} {wildcards.chr} \
                        {wildcards.maf} {output.kinship} 
        '''


use rule yazar_he_kinship_split as yazar_rare_he_kinship_split_part1 with:
    input:
        kinship = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/kinship.chr{chr}.maf{{maf}}.npy'
                    for chr in chrs],
        ctp = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch{i}.gz'
                for i in batch1],
        batch = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch.txt',
    output:
        kinship = [temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/kinship.maf{{maf}}.batch{i}.npy')
                    for i in batch1],


use rule yazar_he_kinship_split as yazar_rare_he_kinship_split_part2 with:
    input:
        kinship = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/kinship.chr{chr}.maf{{maf}}.npy'
                    for chr in chrs],
        ctp = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch{i}.gz'
                for i in batch2],
        batch = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch.txt',
    output:
        kinship = [temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/kinship.maf{{maf}}.batch{i}.npy')
                    for i in batch2],


use rule yazar_he_kinship_split as yazar_rare_he_kinship_split_part3 with:
    input:
        kinship = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/kinship.chr{chr}.maf{{maf}}.npy'
                    for chr in chrs],
        ctp = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch{i}.gz'
                for i in batch3],
        batch = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch.txt',
    output:
        kinship = [temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/kinship.maf{{maf}}.batch{i}.npy')
                    for i in batch3],


use rule yazar_he_kinship_split as yazar_rare_he_kinship_split_part4 with:
    input:
        kinship = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/kinship.chr{chr}.maf{{maf}}.npy'
                    for chr in chrs],
        ctp = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch{i}.gz'
                for i in batch4],
        batch = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch.txt',
    output:
        kinship = [temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/kinship.maf{{maf}}.batch{i}.npy')
                    for i in batch4],


use rule yazar_HE_free as yazar_rare_HE_free_jk with:
    input:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/kinship.maf{{maf}}.batch{{i}}.npy',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.final.gz',
        op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.gz',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he.free.maf{{maf}}.batch{{i}}.jk.npy',
    params:
        jk = True,
        snps = 5, # threshold of snp number per gene
        iid = False,
        full = False,
    resources:
        mem_mb = '20G',


use rule yazar_HE_free_merge as yazar_rare_HE_free_jk_merge with:
    input:
        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he.free.maf{{maf}}.batch{i}.jk.npy'
            for i in range(yazar_reml_batches)],
    output:
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/he.free.maf{{maf}}.jk.npy',


# use rule yazar_HE_free as yazar_rare_HE_free with:
#     input:
#         ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
#         ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
#         kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.maf{{maf}}.batch{{i}}.npy',
#         P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.final.gz',
#         op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.gz',
#         geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
#         obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
#     output:
#         out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he.maf{{maf}}.batch{{i}}.npy',


# use rule yazar_HE_free_merge as yazar_rare_HE_free_merge with:
#     input:
#         out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he.maf{{maf}}.batch{i}.npy'
#             for i in range(yazar_he_batches)],
#     output:
#         out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/he.maf{{maf}}.npy',
