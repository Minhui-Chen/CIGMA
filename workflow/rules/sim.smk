#########################################################################################
##   Simulation
#########################################################################################
# par
sim_batches = np.array_split(range(config['sim']['replicates']), 
                            config['sim']['batch_no'])

## paramspace
sim_params = pd.read_table("sim.params.txt", dtype="str", comment='#', na_filter=False)
if sim_params.shape[0] != sim_params.drop_duplicates().shape[0]:
    sys.exit('Duplicated parameters!\n')
sim_par_columns = list(sim_params.columns)
sim_par_columns.remove('model')
sim_paramspace = Paramspace(sim_params[sim_par_columns], filename_params="*")



rule sim_celltype_expectedPInSnBETAnVnW:
    output:
        pi = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/PI.txt',
        s = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/S.txt',
        beta = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/celltypebeta.txt',
        V = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/V.txt',
        W = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/W.txt',
    resources:
        mem_mb = '1G',
    script: '../bin/sim/celltype_expectedPInSnBETAnVnW.py'


rule sim_generatedata:
    input:
        beta = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/celltypebeta.txt',
        V = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/V.txt',
        W = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/W.txt',
    output:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/sim.batch{{i}}.npy',
    params:
        batches = sim_batches,
        beta = (0.5, 0.5), # beta distribution for allele frequency
        maf = 0.05,
        seed = 273672,
    resources:
        mem_mb = '1G',  
    script: '../bin/sim/generatedata.add_neutral.py'


rule sim_HE:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/sim.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/he.batch{{i}}.npy',
    resources:
        partition = lambda wildcards: config['partition1'] if int(wildcards.ss) <= 1000 and len(wildcards.a.split('_')) <=4 else config['partition2'],
        mem_mb = lambda wildcards: '6G' if int(wildcards.ss) <= 1000 and len(wildcards.a.split('_')) <=4 else '20G' if len(wildcards.a.split('_')) <=4 else '30G',
    params:
        free_jk = True,
        full = False,
    script: '../bin/sim/he.py'


rule sim_mergeBatches_HE:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/he.batch{i}.npy' 
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/out.he.npy',
    resources:
        partition = config['partition1'],
    script: '../bin/mergeBatches.py'


sim_benchmark_batches = np.array_split(range(config['sim']['replicates']), 
                            config['sim']['replicates'])

use rule sim_generatedata as sim_generatedata_benchmark with:
    output:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/sim.batch{{i}}.benchmark.npy',
    params:
        batches = sim_benchmark_batches,
        beta = (0.5, 0.5), # beta distribution for allele frequency
        maf = 0.05,
        seed = 273672,


rule sim_HE_benchmark:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/sim.batch{{i}}.benchmark.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/he.batch{{i}}.benchmark.npy',
    resources:
        partition = lambda wildcards: config['partition1'] if int(wildcards.ss) <= 1000 and len(wildcards.a.split('_')) <=4 else config['partition2'],
        mem_mb = lambda wildcards: '6G' if int(wildcards.ss) <= 1000 and len(wildcards.a.split('_')) <=4 else '20G' if len(wildcards.a.split('_')) <=4 else '30G',
    params:
        free_jk = False,
        full = False,
    benchmark: f"benchmark/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/he.batch{{i}}.benchmark.txt",
    script: '../bin/sim/he.benchmark.py'


use rule sim_mergeBatches_HE as sim_mergeBatches_HE_benchmark with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/he.batch{i}.benchmark.npy' 
                for i in range(config['sim']['replicates'])],
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/out.he.benchmark.npy',


def sim_agg_he_truebeta_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/celltypebeta.txt', params=subspace.instance_patterns)


def sim_agg_he_trueV_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/V.txt', params=subspace.instance_patterns)


def sim_agg_he_trueW_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/W.txt', params=subspace.instance_patterns)


def sim_agg_he_truePi_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/PI.txt', params=subspace.instance_patterns)


def sim_agg_he_trueS_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/S.txt', params=subspace.instance_patterns)


def sim_agg_he_data_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/L_{{L}}_nL_{{nL}}/sim.npy', params=subspace.instance_patterns)


def sim_agg_he_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/L_{{L}}_nL_{{nL}}/out.he.npy', params=subspace.instance_patterns)


rule sim_agg_he_out:
    input:
        out = sim_agg_he_out_subspace,
    output:
        out = 'analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.he.npy',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[sim_params['model']==wildcards.model]).iloc[:,:],
    run:
        args = np.array(params.subspace[wildcards.arg])
        data = {}
        for arg, out in zip(args, input.out):
            data[arg] = np.load(out, allow_pickle=True).item()
        np.save(output.out, data)





##########################################
# REML
##########################################
rule sim_REML:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/sim.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/reml.batch{{i}}.npy',
    resources:
        time = '100:00:00',
        mem_mb = '5G',
    script: '../bin/sim/reml.py'


use rule sim_mergeBatches_HE as sim_mergeBatches_REML with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/reml.batch{i}.npy' 
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/out.reml.npy',


use rule sim_REML as sim_REML_benchmark with:
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/reml.batch{{i}}.benchmark.npy',
    benchmark: f"benchmark/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/reml.batch{{i}}.benchmark.txt",


use rule sim_mergeBatches_HE as sim_mergeBatches_REML_benchmark with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/reml.batch{i}.benchmark.npy' 
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/out.reml.benchmark.npy',


def sim_agg_reml_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/L_{{L}}_nL_{{nL}}/out.reml.npy', params=subspace.instance_patterns)


rule sim_agg_reml_out:
    input:
        out = sim_agg_reml_out_subspace,
    output:
        out = 'analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.reml.npy',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[sim_params['model']==wildcards.model]).iloc[:,:],
    run:
        args = np.array(params.subspace[wildcards.arg])
        data = {}
        for arg, out in zip(args, input.out):
            data[arg] = np.load(out, allow_pickle=True).item()
        np.save(output.out, data)


use rule sim_REML as sim_REML_noCholesky with:
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/reml.noCholesky.batch{{i}}.npy',
    params:
        chol = False,


use rule sim_mergeBatches_HE as sim_mergeBatches_REML_noCholesky with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/reml.noCholesky.batch{i}.npy' 
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/out.reml.noCholesky.npy',


def sim_agg_reml_noChol_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/L_{{L}}_nL_{{nL}}/out.reml.noCholesky.npy', params=subspace.instance_patterns)


use rule sim_agg_reml_out as sim_agg_reml_noChol_out with:
    input:
        out = sim_agg_reml_noChol_out_subspace,
    output:
        out = 'analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.reml.noCholesky.npy',



############################
# 1.2 real genotype 
############################
rule sim_realgeno_extract_genotypes:
    input:
        location = 'data/Yazar2022Science/gene_location.txt',
        bed = 'analysis/yazar/data/geno/chr{chr}.bed',
    output:
        genes = 'staging/sim/genes.chr{chr}.genotype.npy',
    params:
        flank = 500000,
        maf = 0.05,
    resources:
        partition = 'tier3q',
        mem_mb = '120G',
    shell:
        '''
        # load plink 1.9
        {config[plink_load]}

        python3 workflow/bin/sim/extract_genotypes.py \
            {input.location} \
            {input.bed} \
            {output.genes} \
            {wildcards.chr} \
            {params.flank} \
            {params.maf}
        '''


rule sim_realgeno_generatedata:
    input:
        beta = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/celltypebeta.txt',
        V = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/V.txt',
        W = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/W.txt',
        genes = expand('staging/sim/genes.chr{chr}.genotype3.npy', chr=range(1,23)),
    output:
        data = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/yazarL_{{L}}/sim.batch{i}.npy'
                for i in range(config['sim']['batch_no'])],
    params:
        out = 'analysis/yazar/ind_min_cellnum~10_ct_min_cellnum~10_prop~0.9_geno_pca_n~6_op_pca_n~1_batch~shared_fixed~shared/he.npy',
        seed = 2713672,
    resources:
        partition = config['partition3'],
        mem_mb = '120G',
    script: '../bin/sim/generatedata.real_genotype.py'


use rule sim_HE as sim_realgeno_HE with:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/yazarL_{{L}}/sim.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/yazarL_{{L}}/he.batch{{i}}.npy',
    params:
        free_jk = False,
        full = False,
    resources:
        partition = config['partition1'],
        mem_mb = '6G',


use rule sim_mergeBatches_HE as sim_realgeno_mergeBatches_HE with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/yazarL_{{L}}/he.batch{i}.npy' 
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/yazarL_{{L}}/out.realgeno.he.npy',


def sim_realgeno_agg_he_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/yazarL_{{L}}/out.realgeno.he.npy', params=subspace.instance_patterns)


use rule sim_agg_he_out as sim_realgeno_agg_he_out with:
    input:
        out = sim_realgeno_agg_he_out_subspace,
    output:
        out = 'analysis/sim/{model}/yazarL_{L}/AGG{arg}.he.npy',


rule sim_agg_L:
    input:
        out = expand('analysis/sim/free5/yazarL_{L}/AGGss.he.npy', L=[0.1, 0.2, 1, 10, 20, 100]),
    output:
        out = 'analysis/sim/free5/yazar/AGGss.he.npy',
    params:
        L = [0.1, 0.2, 1, 10, 20, 100],
    run:
        data = {}
        for L, out in zip(params.L, input.out):
            data[L] = np.load(out, allow_pickle=True).item()
        np.save(output.out, data)





############################
# 1.2 HE without JK 
############################
use rule sim_HE as sim_HE_noJK with:
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/he.noJK.batch{{i}}.npy',
    params:
        free_jk = False,
        full = False,


use rule sim_mergeBatches_HE as sim_mergeBatches_HE_noJK with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/he.noJK.batch{i}.npy' 
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/out.he.noJK.npy',


def sim_agg_he_noJK_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/L_{{L}}_nL_{{nL}}/out.he.noJK.npy', params=subspace.instance_patterns)


use rule sim_agg_he_out as sim_agg_he_noJK_out with:
    input:
        out = sim_agg_he_noJK_out_subspace,
    output:
        out = 'analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.he.noJK.npy',







############################
rule sim_agg_parameters:
    input:
        V = sim_agg_he_trueV_subspace,
        W = sim_agg_he_trueW_subspace,
        Pi = sim_agg_he_truePi_subspace,
        S = sim_agg_he_trueS_subspace,
        beta = sim_agg_he_truebeta_subspace,
    output:
        V = 'analysis/sim/{model}/AGG{arg}.true_V.npy',
        W = 'analysis/sim/{model}/AGG{arg}.true_W.npy',
        Pi = 'analysis/sim/{model}/AGG{arg}.true_Pi.npy',
        S = 'analysis/sim/{model}/AGG{arg}.true_S.npy',
        beta = 'analysis/sim/{model}/AGG{arg}.true_Beta.npy',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[sim_params['model']==wildcards.model]).iloc[:,:],
    run:
        args = np.array(params.subspace[wildcards.arg])
        V = {}
        for arg, V_f in zip(args, input.V):
            V[arg] = np.loadtxt(V_f)
        np.save(output.V, V)

        W = {}
        for arg, W_f in zip(args, input.W):
            W[arg] = np.loadtxt(W_f)
        np.save(output.W, W)

        Pi = {}
        for arg, Pi_f in zip(args, input.Pi):
            Pi[arg] = np.loadtxt(Pi_f)
        np.save(output.Pi, Pi)

        S = {}
        for arg, S_f in zip(args, input.S):
            S[arg] = np.loadtxt(S_f)
        np.save(output.S, S)

        beta = {}
        for arg, beta_f in zip(args, input.beta):
            beta[arg] = np.loadtxt(beta_f)
        np.save(output.beta, beta)


def sim_HE_AGGarg_fun(wildcards):
    effective_args = get_effective_args(sim_params.loc[sim_params['model']==wildcards.model])
    return (expand('analysis/sim/{{model}}/L_{{L}}_nL_{{nL}}/AGG{arg}.he.npy', arg=effective_args) 
    + expand('analysis/sim/{{model}}/AGG{arg}.true_V.npy', arg=effective_args))


rule sim_HE_AGGarg:
    input:
        sim_HE_AGGarg_fun,
    output:
        flag = touch('staging/sim/{model}/L_{L}_nL_{nL}/HE.all.flag'),


rule sim_HE_all:
    input:
        flag = expand('staging/sim/{model}/L_{{L}}_nL_{{nL}}/HE.all.flag', model=['free']),


# cell type number
rule sim_celltype_number_all:
    input:
        out1 = 'analysis/sim/free21/L_{L}_nL_{nL}/AGGss.he.npy',
        V1 = 'analysis/sim/free21/AGGss.true_V.npy',
        out2 = 'analysis/sim/free22/L_{L}_nL_{nL}/AGGss.he.npy',
        V2 = 'analysis/sim/free22/AGGss.true_V.npy',
        out3 = 'analysis/sim/hom21/L_{L}_nL_{nL}/AGGss.he.npy',
        out4 = 'analysis/sim/hom22/L_{L}_nL_{nL}/AGGss.he.npy',
    output:
        flag = touch('analysis/sim/L_{L}_nL_{nL}/celltype_number.flag'),





















###################################################################
# when nu is noisy
###################################################################
rule sim_data_add_noise_to_nu:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/sim.batch{{i}}.npy',
    output:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/sim.noise_{{alpha}}.batch{{i}}.npy',
    run:
        data = np.load(input.data, allow_pickle=True).item()
        rng = np.random.default_rng(125)
        for k in data.keys():
            nu = data[k]['ctnu']
            noise = rng.choice([-1, 1], nu.shape) * rng.beta(float(wildcards.alpha), 1, nu.shape)
            data[k]['ctnu'] = nu * (1 + noise)
        np.save(output.data, data)


use rule sim_HE as sim_HE_noisy_nu with:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/sim.noise_{{alpha}}.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/he.noise_{{alpha}}.batch{{i}}.npy',


use rule sim_mergeBatches_HE as sim_mergeBatches_HE_noisy_nu with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/he.noise_{{alpha}}.batch{i}.npy'
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/out.he.noise_{{alpha}}.npy',


rule sim_agg_he_noisy_nu_out_alpha:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/out.he.noise_{alpha}.npy'
                for alpha in config['sim']['noise_alpha']]
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/out.he.noisy_nu.npy',
    params:
        alpha = config['sim']['noise_alpha'],
    run:
        data = {}
        for alpha, alpha_f in zip(params.alpha, input.out):
            data[alpha] = np.load(alpha_f, allow_pickle=True).item()
        np.save(output.out, data)


# reml
use rule sim_REML as sim_REML_noisy_nu with:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/sim.noise_{{alpha}}.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/reml.noise_{{alpha}}.batch{{i}}.npy',


use rule sim_mergeBatches_HE as sim_mergeBatches_REML_noisy_nu with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/reml.noise_{{alpha}}.batch{i}.npy' 
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/out.reml.noise_{{alpha}}.npy',


use rule sim_agg_he_noisy_nu_out_alpha as sim_agg_reml_noisy_nu_out_alpha with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/out.reml.noise_{alpha}.npy'
                for alpha in config['sim']['noise_alpha']]
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/out.reml.noisy_nu.npy',


def sim_agg_he_noisy_nu_out(wildcards):
    subspace1 = get_subspace('ss', sim_params.loc[sim_params['model']=='hom'].head(1))
    subspace2 = get_subspace('ss', sim_params.loc[sim_params['model']=='free'].head(1))
    return (expand('analysis/sim/hom/{params}/L_{{L}}_nL_{{nL}}/out.he.noisy_nu.npy', params=subspace1.instance_patterns) + 
            expand('analysis/sim/free/{params}/L_{{L}}_nL_{{nL}}/out.he.noisy_nu.npy', params=subspace2.instance_patterns))


rule sim_agg_he_noisy_nu_out_subspace:
    input:
        out = sim_agg_he_noisy_nu_out,
    output:
        out = touch('staging/sim/L_{L}_nL_{nL}/he.noisy_nu.flag'),



















###################################################################
# when nu is not considered when fitting model
###################################################################
#####################################
# CIGMA missing nu
#####################################
rule sim_HE_missing_nu:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/sim.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/he.missing_nu.batch{{i}}.npy',
    resources:
        partition = lambda wildcards: config['partition1'] if int(wildcards.ss) <= 1000 and len(wildcards.a.split('_')) <=4 and wildcards.model != 'full' else config['partition2'],
        mem_mb = lambda wildcards: '6G' if int(wildcards.ss) <= 1000 and len(wildcards.a.split('_')) <=4 and wildcards.model != 'full' else '30G',
    params:
        free_jk = False,
    script: '../bin/sim/he.missing_nu.py'


use rule sim_mergeBatches_HE as sim_mergeBatches_HE_missing_nu with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/he.missing_nu.batch{i}.npy'
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/out.he.missing_nu.npy',


def sim_agg_he_out_subspace_missing_nu(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/L_{{L}}_nL_{{nL}}/out.he.missing_nu.npy', 
                    params=subspace.instance_patterns)


rule sim_agg_he_out_missing_nu:
    input:
        out = sim_agg_he_out_subspace_missing_nu,
    output:
        out = 'analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.he.missing_nu.npy',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[sim_params['model']==wildcards.model]).iloc[:,:],
    run:
        args = np.array(params.subspace[wildcards.arg])
        data = {}
        for arg, out in zip(args, input.out):
            data[arg] = np.load(out, allow_pickle=True).item()
        np.save(output.out, data)









#####################################
# GCTA
#####################################
rule sim_gcta_greml_ctp:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/sim.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/greml.batch{{i}}.npy',
    resources:
        mem_mb = lambda wildcards: '2G' if int(wildcards.ss) <= 2000 else '4G',
    shell:
        '''
        module load gcc/11.3.0 gcta/1.94.1
        python3 workflow/bin/sim/gcta_greml_ctp.py \
                {input.data} {output.out} \
        '''


use rule sim_mergeBatches_HE as sim_gcta_greml_ctp_merge with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/greml.batch{i}.npy'
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/greml.npy',


# TODO: temp
use rule sim_gcta_greml_ctp as sim_gcta_greml_ctp_ld with:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/ld{{ld}}/sim.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/ld{{ld}}/greml.batch{{i}}.npy',
    resources:
        mem_mb = lambda wildcards: '2G' if int(wildcards.ss) <= 2000 else '4G',
        partition = config['partition1'],


use rule sim_mergeBatches_HE as sim_gcta_greml_ctp_ld_merge with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/ld{{ld}}/greml.batch{i}.npy'
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/ld{{ld}}/greml.npy',


def sim_agg_greml_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('staging/sim/{{model}}/{params}/L_{{L}}_nL_{{nL}}/greml.npy', params=subspace.instance_patterns)


rule sim_gcta_greml_ctp_agg:
    input:
        out = sim_agg_greml_out_subspace,
    output:
        out = 'analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.greml.npy',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[sim_params['model']==wildcards.model]).iloc[:,:],
    run:
        args = np.array(params.subspace[wildcards.arg])
        data = {}
        for arg, out in zip(args, input.out):
            data[arg] = np.load(out, allow_pickle=True).item()
        np.save(output.out, data)


rule sim_gcta_all:
    input:
        V = 'analysis/sim/{model}/AGG{arg}.true_V.npy',
        gcta = 'analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.greml.npy',
        missing_nu_out = 'analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.he.missing_nu.npy',
        out = 'analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.he.npy',
    output:
        flag = touch('analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.gcta.flag'),


rule sim_gcta_greml_aireml_ctp:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/sim.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/greml.aireml.batch{{i}}.npy',
    resources:
        mem_mb = lambda wildcards: '2G' if int(wildcards.ss) <= 2000 else '4G',
    priority: 1
    shell:
        '''
        module load gcc/11.3.0 gcta/1.94.1
        python3 workflow/bin/sim/gcta_greml_aireml_ctp.py \
                {input.data} {output.out} \
        '''


use rule sim_mergeBatches_HE as sim_gcta_greml_aireml_ctp_merge with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/greml.aireml.batch{i}.npy'
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/greml.aireml.npy',


def sim_agg_greml_aireml_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('staging/sim/{{model}}/{params}/L_{{L}}_nL_{{nL}}/greml.aireml.npy', params=subspace.instance_patterns)


use rule sim_gcta_greml_ctp_agg as sim_gcta_greml_aireml_ctp_agg with:
    input:
        out = sim_agg_greml_aireml_out_subspace,
    output:
        out = 'analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.greml.aireml.npy',





#####################################
# OTD
#####################################
rule otd_demean_ctp:
    input:
        data = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/sim.batch{i}.npy'
                for i in range(config['sim']['batch_no'])],
    output:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/otd_demean.batch{i}.npy'
                for i in range(config['sim']['batch_no'])],
    resources:
        mem_mb = '2G',
    run:
        for input_f, output_f in zip(input.data, output.out):
            data = np.load(input_f, allow_pickle=True).item()
            for k in data.keys():
                y = data[k]['Y']
                data[k]['y'] = y.mean(axis=1)
                y_demean = y - y.mean(axis=1, keepdims=True)
                data[k]['Y'] = y_demean
            np.save(output_f, data)


use rule sim_gcta_greml_ctp as sim_otd_gcta_greml_ctp with:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/otd_demean.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/otd.batch{{i}}.npy',


use rule sim_mergeBatches_HE as sim_otd_gcta_greml_ctp_merge with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/otd.batch{i}.npy'
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/otd.npy',


def sim_agg_otd_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('staging/sim/{{model}}/{params}/L_{{L}}_nL_{{nL}}/otd.npy', params=subspace.instance_patterns)


use rule sim_gcta_greml_ctp_agg as sim_otd_gcta_greml_ctp_agg with:
    input:
        out = sim_agg_otd_out_subspace,
    output:
        out = 'analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.otd.npy',


rule sim_otd_all:
    input:
        V = 'analysis/sim/free3/AGGvc.true_V.npy',
        otd = 'analysis/sim/free3/L_{L}_nL_{nL}/AGGvc.otd.npy',






















#####################################
# BOLT
#####################################
rule sim_bolt_ctp:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/sim.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/bolt.batch{{i}}.npy',
    resources:
        mem_mb = lambda wildcards: '2G' if int(wildcards.ss) <= 2000 else '4G',
    shell:
        '''
        module load gcc/12.1.0 bolt-lmm/2.4.1 plink/1.9
        python3 workflow/bin/sim/bolt_ctp.py \
                {input.data} {output.out} \
        '''


use rule sim_mergeBatches_HE as sim_bolt_ctp_merge with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/bolt.batch{i}.npy'
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/bolt.npy',

# TODO: temp
use rule sim_bolt_ctp as sim_bolt_ctp_ld with:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/ld{{ld}}/sim.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/ld{{ld}}/bolt.batch{{i}}.npy',
    resources:
        mem_mb = lambda wildcards: '2G' if int(wildcards.ss) <= 2000 else '4G',
        partition = config['partition1'],


use rule sim_mergeBatches_HE as sim_bolt_ctp_ld_merge with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/ld{{ld}}/bolt.batch{i}.npy'
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/ld{{ld}}/bolt.npy',


def sim_agg_bolt_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('staging/sim/{{model}}/{params}/L_{{L}}_nL_{{nL}}/bolt.npy', params=subspace.instance_patterns)


use rule sim_gcta_greml_ctp_agg as sim_bolt_ctp_agg with:
    input:
        out = sim_agg_bolt_out_subspace,
    output:
        out = 'analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.bolt.npy',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[sim_params['model']==wildcards.model]).iloc[:,:],


rule sim_bolt_all:
    input:
        V = 'analysis/sim/{model}/AGG{arg}.true_V.npy',
        bolt = 'analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.bolt.npy',











#####################################
# GxEMM
# "R CMD INSTALL GxEMM_1.0.tar.gz" from https://github.com/andywdahl/gxemm/tree/master
#####################################
rule sim_gxemm:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/sim.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/gxemm.batch{{i}}.txt',
    resources:
        mem_mb = '2G',
    script: '../bin/sim/gxemm.R'


rule sim_gxemm_merge:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/gxemm.batch{i}.txt' 
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/gxemm.txt',
    run:
        with open(output.out, 'w') as outfile:
            for i, infile in enumerate(input.out):
                with open(infile) as f:
                    if i != 0:
                        next(f)  # skip header
                    for line in f:
                        outfile.write(line)


def sim_agg_gxemm_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('staging/sim/{{model}}/{params}/L_{{L}}_nL_{{nL}}/gxemm.txt', params=subspace.instance_patterns)


rule sim_gxemm_ctp_agg:
    input:
        out = sim_agg_gxemm_out_subspace,
    output:
        out = 'analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.gxemm.gz',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[sim_params['model']==wildcards.model]).iloc[:,:],
    run:
        args = np.array(params.subspace[wildcards.arg])
        data = []
        for arg, out in zip(args, input.out):
            tmp = pd.read_table(out)
            tmp['arg'] = arg
            data.append(tmp)
        data = pd.concat(data, axis=0)
        data.to_csv(output.out, sep='\t', index=False)


rule sim_cellcount_all:
    input:
        V = 'analysis/sim/{model}/AGG{arg}.true_V.npy',
        gcta = 'analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.greml.npy',
        otd = 'analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.otd.npy',
        bolt = 'analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.bolt.npy',
        missing_nu_out = 'analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.he.missing_nu.npy',
        out = 'analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.he.npy',
        gxemm = 'analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.gxemm.gz',
    output:
        flag = touch('analysis/sim/{model}/L_{L}_nL_{nL}/AGG{arg}.cellcount.flag'),


















#####################################
# Non Gaussian
#####################################
rule sim_nongaussian_generatedata:
    input:
        beta = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/celltypebeta.txt',
        V = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/V.txt',
        W = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/W.txt',
    output:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/{{dist_g}}_g.{{dist_v}}_e.sim.batch{{i}}.npy',
    params:
        batches = sim_batches,
        beta = (0.5, 0.5), # beta distribution for allele frequency
        maf = 0.05,
        seed = 273672,
    resources:
        mem_mb = '1G',
    script: '../bin/sim/generatedata.non_gaussian.py'


use rule sim_HE as sim_nongaussian_HE with:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/{{dist_g}}_g.{{dist_v}}_e.sim.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/{{dist_g}}_g.{{dist_v}}_e.he.batch{{i}}.npy',
    params:
        free_jk = False, # NOTE: no JK here
        full = False


use rule sim_mergeBatches_HE as sim_nongaussian_mergeBatches_HE with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/{{dist_g}}_g.{{dist_v}}_e.he.batch{i}.npy' 
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/{{dist_g}}_g.{{dist_v}}_e.out.he.npy',


def sim_nongaussian_agg_he_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/L_{{L}}_nL_{{nL}}/{{dist_g}}_g.{{dist_v}}_e.out.he.npy', params=subspace.instance_patterns)


use rule sim_agg_he_out as sim_nongaussian_agg_he_out with:
    input:
        out = sim_nongaussian_agg_he_out_subspace,
    output:
        out = 'analysis/sim/{model}/L_{L}_nL_{nL}/{dist_g}_g.{dist_v}_e/AGG{arg}.he.npy',


rule sim_nongaussian_all:
    input:
        out = expand('analysis/sim/free/L_{L}_nL_{nL}/{dist_g}_g.{dist_v}_e/AGGss.he.npy',
                    zip, dist_g=['gaussian', 'spikeslab', 'spikeslab', 't2'], dist_v=['spikeslab', 'gaussian', 'spikeslab', 't2'],
                    L=[config['sim']['L']] * 4, nL=[config['sim']['nL']] * 4),


# reml
use rule sim_REML as sim_nongaussian_REML with:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/{{dist_g}}_g.{{dist_v}}_e.sim.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/{{dist_g}}_g.{{dist_v}}_e.reml.batch{{i}}.npy',


use rule sim_mergeBatches_HE as sim_nongaussian_mergeBatches_REML with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/{{dist_g}}_g.{{dist_v}}_e.reml.batch{i}.npy' 
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/L_{{L}}_nL_{{nL}}/{{dist_g}}_g.{{dist_v}}_e.out.reml.npy',


# def sim_nongaussian_agg_reml_out_subspace(wildcards):
#     subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
#     return expand('analysis/sim/{{model}}/{params}/L_{{L}}_nL_{{nL}}/{{dist_g}}_g.{{dist_v}}_e.out.reml.npy', params=subspace.instance_patterns)


# use rule sim_agg_he_out as sim_nongaussian_agg_reml_out with:
#     input:
#         out = sim_nongaussian_agg_reml_out_subspace,
#     output:
#         out = 'analysis/sim/{model}/L_{L}_nL_{nL}/{dist_g}_g.{dist_v}_e/AGG{arg}.reml.npy',


# rule sim_nongaussian_reml_all:
#     input:
#         out = expand('analysis/sim/free/L_{L}_nL_{nL}/{dist_g}_g.{dist_v}_e/AGGss.reml.npy',
#                     zip, dist_g=['gaussian', 'spikeslab', 'spikeslab', 't2'], dist_v=['spikeslab', 'gaussian', 'spikeslab', 't2'],
#                     L=[config['sim']['L']] * 4, nL=[config['sim']['nL']] * 4),













#####################################
# number of causal variants
#####################################
# use rule sim_generatedata as sim_L5_generatedata with:
#     input:
#         beta = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/celltypebeta.txt',
#         V = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/V.txt',
#         W = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/W.txt',
#     output:
#         data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/sim.L5.batch{{i}}.npy',
#     params:
#         batches = sim_batches,
#         beta = (0.5, 0.5), # beta distribution for allele frequency
#         maf = 0.05,
#         L = 5, # number of causal SNPs
#         seed = 2736731,


# use rule sim_HE as sim_L5_HE with:
#     input:
#         data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/sim.L5.batch{{i}}.npy',
#     output:
#         out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.L5.batch{{i}}.npy',
#     params:
#         free_jk = True,
#         full = False


# use rule sim_mergeBatches_HE as sim_L5_mergeBatches_HE with:
#     input:
#         out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.L5.batch{i}.npy' 
#                 for i in range(config['sim']['batch_no'])],
#     output:
#         out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/out.he.L5.npy',


# use rule sim_generatedata as sim_L15_generatedata with:
#     input:
#         beta = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/celltypebeta.txt',
#         V = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/V.txt',
#         W = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/W.txt',
#     output:
#         data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/sim.L15.batch{{i}}.npy',
#     params:
#         batches = sim_batches,
#         beta = (0.5, 0.5), # beta distribution for allele frequency
#         maf = 0.05,
#         L = 15, # number of causal SNPs
#         seed = 2736723,


# use rule sim_HE as sim_L15_HE with:
#     input:
#         data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/sim.L15.batch{{i}}.npy',
#     output:
#         out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.L15.batch{{i}}.npy',
#     params:
#         free_jk = True,
#         full = False


# use rule sim_mergeBatches_HE as sim_L15_mergeBatches_HE with:
#     input:
#         out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.L15.batch{i}.npy' 
#                 for i in range(config['sim']['batch_no'])],
#     output:
#         out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/out.he.L15.npy',


# use rule sim_generatedata as sim_L1000_generatedata with:
#     input:
#         beta = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/celltypebeta.txt',
#         V = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/V.txt',
#         W = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/W.txt',
#     output:
#         data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/sim.L1000.batch{{i}}.npy',
#     params:
#         batches = sim_batches,
#         beta = (0.5, 0.5), # beta distribution for allele frequency
#         maf = 0.05,
#         L = 1000, # number of causal SNPs
#         seed = 2736713,


# use rule sim_HE as sim_L1000_HE with:
#     input:
#         data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/sim.L1000.batch{{i}}.npy',
#     output:
#         out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.L1000.batch{{i}}.npy',
#     params:
#         free_jk = True,
#         full = False


# use rule sim_mergeBatches_HE as sim_L1000_mergeBatches_HE with:
#     input:
#         out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.L1000.batch{i}.npy' 
#                 for i in range(config['sim']['batch_no'])],
#     output:
#         out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/out.he.L1000.npy',


# rule sim_L10_nL_990_generatedata:
#     input:
#         beta = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/celltypebeta.txt',
#         V = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/V.txt',
#         W = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/W.txt',
#     output:
#         data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/sim.L10.nL990.batch{{i}}.npy',
#     params:
#         batches = sim_batches,
#         beta = (0.5, 0.5), # beta distribution for allele frequency
#         maf = 0.05,
#         L = 10, # number of causal SNPs
#         nL = 990, # number of non-causal SNPs
#         seed = 273333,
#     resources:
#         mem_mb = '1G',
#     script: '../bin/sim/generatedata.add_neutral.py'


# use rule sim_HE as sim_L10_nL990_HE with:
#     input:
#         data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/sim.L10.nL990.batch{{i}}.npy',
#     output:
#         out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.L10.nL990.batch{{i}}.npy',
#     params:
#         free_jk = True,
#         full = False


# use rule sim_mergeBatches_HE as sim_L10_nL990_mergeBatches_HE with:
#     input:
#         out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.L10.nL990.batch{i}.npy' 
#                 for i in range(config['sim']['batch_no'])],
#     output:
#         out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/out.he.L10.nL990.npy',


# rule sim_L_all:
#     input:
#         L5 = expand('analysis/sim/free/{params}/out.he.L5.npy',
#                     params=Paramspace(sim_params.loc[sim_params['model']=='free'].drop('model', axis=1).reset_index(drop=True).iloc[[0]], 
#                                         filename_params="*").instance_patterns),
#         L15 = expand('analysis/sim/free/{params}/out.he.L15.npy',
#                     params=Paramspace(sim_params.loc[sim_params['model']=='free'].drop('model', axis=1).reset_index(drop=True).iloc[[0]], 
#                                         filename_params="*").instance_patterns),
#         # L15 = expand('analysis/sim/free/{params}/out.he.L15.npy',
#         #             params=Paramspace(sim_params.loc[sim_params['model']=='free'].drop('model', axis=1).reset_index(drop=True).iloc[[0]], 
#         #                                 filename_params="*").instance_patterns),
#         L1000 = expand('analysis/sim/free/{params}/out.he.L1000.npy',
#                     params=Paramspace(sim_params.loc[sim_params['model']=='free'].drop('model', axis=1).reset_index(drop=True).iloc[[0]], 
#                                         filename_params="*").instance_patterns),
#         L10_nL990 = expand('analysis/sim/free/{params}/out.he.L10.nL990.npy',
#                     params=Paramspace(sim_params.loc[sim_params['model']=='free'].drop('model', axis=1).reset_index(drop=True).iloc[[0]], 
#                                         filename_params="*").instance_patterns),















































rule sim_all:
    input:
        hom = expand('staging/sim/hom/L_{L}_nL_{nL}/HE.all.flag',
                     L=config['sim']['L'],
                     nL=config['sim']['nL']),
        free = expand('staging/sim/free/L_{L}_nL_{nL}/HE.all.flag',
                     L=config['sim']['L'],
                     nL=config['sim']['nL']),
        full = expand('analysis/sim/full/L_{L}_nL_{nL}/AGGV_tril.he.noJK.npy',
                     L=config['sim']['L'],
                     nL=config['sim']['nL']),
        full_par = 'analysis/sim/full/AGGV_tril.true_V.npy',
        cellcount = expand('analysis/sim/free3/L_{L}_nL_{nL}/AGGvc.cellcount.flag',  # fig 2
                           L=config['sim']['L'],
                           nL=config['sim']['nL']),
        hom_cellcount = expand('analysis/sim/hom3/L_{L}_nL_{nL}/AGGvc.cellcount.flag',
                           L=config['sim']['L'],
                           nL=config['sim']['nL']),
        ct_num = expand('analysis/sim/L_{L}_nL_{nL}/celltype_number.flag',
                     L=config['sim']['L'],
                     nL=config['sim']['nL']),
        noisy_nu = expand('staging/sim/L_{L}_nL_{nL}/he.noisy_nu.flag',
                     L=config['sim']['L'],
                     nL=config['sim']['nL']),
        out1 = expand('analysis/sim/free/{params}/L_{L}_nL_{nL}/out.he.noJK.npy',
                     params=Paramspace(sim_params.loc[sim_params['model']=='free'].drop('model', axis=1).reset_index(drop=True).iloc[[0]], 
                                        filename_params="*").instance_patterns,
                     L='100',
                     nL='0'),
        out2 = expand('analysis/sim/free/{params}/L_{L}_nL_{nL}/out.he.noJK.npy',
                     params=Paramspace(sim_params.loc[sim_params['model']=='free'].drop('model', axis=1).reset_index(drop=True).iloc[[0]], 
                                        filename_params="*").instance_patterns,
                     L='10',
                     nL='90'),
        non_gaussian = expand('analysis/sim/free/L_{L}_nL_{nL}/{dist_g}_g.{dist_v}_e/AGGss.he.npy',
                    zip, dist_g=['gaussian', 'spikeslab', 'spikeslab', 't2'], dist_v=['spikeslab', 'gaussian', 'spikeslab', 't2'],
                    L=[config['sim']['L']] * 4, nL=[config['sim']['nL']] * 4),



