##########################################################################
# Yazar simulation: permute cells
##########################################################################
# yazar_sim_he_batches = config['yazar_sim']['he_nbatch']
# yazar_sim_replicates = config['yazar_sim']['replicates']

rule yazar_sim_permute_ct:
    input:
        h5ad = 'analysis/yazar/data/logp_cp10k.h5ad',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
        ctp = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctp.final.gz',
    output:
        X = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/X.npz',
        obs = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/obs.gz',
        var = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/var.gz',
    params:
        ind_col = yazar_ind_col,
        ct_col = yazar_ct_col,
        pool_col = yazar_pool_col,
        seed = 123,
    resources:
        mem_mb = '90G',
    script: '../bin/yazar/sim.transform.py'


use rule yazar_ctp as yazar_sim_permute_ct_ctp with:
    input:
        X = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/X.npz',
        obs = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/obs.gz',
        var = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/var.gz',
    output:
        ctp = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/ctp.gz',
        ctnu = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/ctnu.gz',
        P = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/P.final.gz',
        n = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/n.final.gz',


use rule yazar_std_op as yazar_sim_permute_ct_std_op with:
    input:
        ctp = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/ctp.gz',
        ctnu = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/ctnu.gz',
        P = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/P.final.gz',
    output:
        op = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/op.final.gz',
        nu = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/nu.final.gz',
        ctp = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/ctp.final.gz',
        ctnu = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/ctnu.final.gz',


use rule yazar_op_pca as yazar_sim_permute_ct_op_pca with:
    input:
        op = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/op.final.gz',
    output:
        evec = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/evec.gz',
        eval = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/eval.gz',
        pca = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/pca.gz',
        png = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/op.pca.png',


use rule yazar_HE_split as yazar_sim_permute_ct_HE_split with:
    input:
        ctp = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/ctp.final.gz',
        ctnu = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/ctnu.final.gz',
        genes = 'data/Yazar2022Science/gene_location.txt',
    output:
        ctp = [f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/he/ctp.batch{i}.gz'
                for i in range(yazar_reml_batches)],
        ctnu = [f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/he/ctnu.batch{i}.gz'
                for i in range(yazar_reml_batches)],
        batch = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/he/ctp.batch.txt',
    

use rule yazar_he_kinship_split as yazar_sim_permute_ct_kinship_split with:
    input:
        kinship = [f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/kinship.chr{chr}.npy'
                    for chr in chrs],
        ctp = [f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/he/ctp.batch{i}.gz'
                for i in range(yazar_reml_batches)],
        batch = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/he/ctp.batch.txt',
    output:
        kinship = [temp(f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/he/kinship.batch{i}.npy')
                    for i in range(yazar_reml_batches)],
    

use rule yazar_HE_free_jk as yazar_sim_permute_ct_HE with:
    input:
        ctp = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/he/kinship.batch{{i}}.npy',
        P = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/P.final.gz',
        op_pca = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/pca.gz',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
    output:
        out = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/he/batch{{i}}.npy',
    params:
        jk = True,
        snps = 5, # threshold of snp number per gene
        iid = False,
        free = True,
        full = False,


use rule yazar_HE_free_merge as yazar_sim_permute_ct_HE_merge with:
    input:
        out = [f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/he/batch{i}.npy'
            for i in range(yazar_reml_batches)],
    output:
        out = f'analysis/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_ct/he.npy',


rule yazar_sim_permute_ct_all:
    input:
        sim = expand('analysis/sim/yazar/{params}/permute_ct/he.npy',
                      params=yazar_paramspace.instance_patterns),





























##########################################################################
# Yazar simulation: permute genotypes
##########################################################################
rule yazar_sim_permute_geno_order_inds:
    input:
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.final.gz',
    output:
        inds = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_geno/inds.txt'
    run:
        rng = np.random.default_rng(123)
        P = pd.read_table(input.P, index_col=0)
        inds = P.index.tolist()

        np.savetxt(output.inds, rng.permutation(inds), fmt='%s')


rule yazar_sim_permute_geno_HE:
    input:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/kinship.batch{{i}}.npy',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.final.gz',
        op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.gz',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
        ids = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_geno/inds.txt',
    output:
        out = f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_geno/he/batch{{i}}.npy',
    params:
        jk = True,
        free = True,
        iid = False,
        snps = 5, # threshold of snp number per gene
    resources:
        mem_mb = '20G',
    script: '../bin/yazar/sim.he.free.py' 


use rule yazar_HE_free_merge as yazar_sim_permute_geno_HE_free_jk_merge with:
    input:
        out = [f'staging/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_geno/he/batch{i}.npy'
            for i in range(yazar_reml_batches)],
    output:
        out = f'analysis/sim/yazar/{yazar_paramspace.wildcard_pattern}/permute_geno/he.npy',


rule yazar_sim_permute_geno_all:
    input:
        out = expand('analysis/sim/yazar/{params}/permute_geno/he.npy',
                      params=yazar_paramspace.instance_patterns),








rule yazar_sim_all:
    input:
        geno = expand('analysis/sim/yazar/{params}/permute_geno/he.npy',
                      params=yazar_paramspace.instance_patterns),
        ct = expand('analysis/sim/yazar/{params}/permute_ct/he.npy',
                      params=yazar_paramspace.instance_patterns),
