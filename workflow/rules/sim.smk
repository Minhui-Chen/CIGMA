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

sim_plot_order = {
    'hom':{
        'ss':['50', '100', '200', '300', '500', '1000', '1500', '2000'], 
        'a':['0.5_2_2_2', '1_2_2_2', '2_2_2_2', '4_2_2_2']
        },
    'iid':{
        'ss':['50', '100', '300', '500', '1000'], 'a':['0.5_2_2_2', '1_2_2_2', '2_2_2_2', '4_2_2_2'],
        'vc':['0.3_0.2_0.2_0.0333_0.2667_0', '0.3_0.15_0.25_0.05_0.25_0', '0.3_0.1_0.3_0.0666_0.2334_0',
                '0.3_0.15_0.15_0.05_0.15_0.2', '0.3_0.15_0.05_0.05_0.05_0.4'],
        }, 
    'free': {
        'ss':['50', '100', '200', '300','500', '1000', '1500', '2000'], 
        'a':['0.5_2_2_2', '1_2_2_2', '2_4_4_4', '2_2_2_2', '4_2_2_2', '8_2_2_2'], 
        'vc':['0.3_0.1_0.1_0.05_0.15_0.3', '0.3_0.1_0.1_0.1_0.1_0.3', '0.2_0.1_0.1_0.2_0.2_0.2', 
            '0.1_0.01_0.194_0.001_0.05_0.645', '0.1_0.05_0.15_0.005_0.05_0.645'],
        'V_diag':['0.5_1_1_1', '1_1_1_1', '2_1_1_1', '10_1_1_1'],
        },
    'freeW': {
        'ss':['50', '100', '200', '300','500', '1000'],
        },
    'full':{
        'ss':['50', '100', '200', '300', '500', '1000', '2000'], 'a':['0.5_2_2_2', '1_2_2_2', '2_2_2_2', '4_2_2_2'],
        'vc':['0.3_0.1_0.1_0.1_0.1_0.3', '0.2_0.1_0.1_0.2_0.2_0.2', '0.1_0.1_0.1_0.3_0.3_0.1'],
        'V_diag':['1_1_1_1', '8_4_2_1', '27_9_3_1', '64_16_4_1', '64_64_1_1'],
        'V_tril':['0.25_0.25_0_-0.25_0_0', '0.5_0.5_0_-0.5_0_0', '0.75_0.75_0_-0.75_0_0', '0.95_0.95_0.95_-0.95_-0.95_-0.95']
        },
    }


rule sim_celltype_expectedPInSnBETAnVnW:
    output:
        pi = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/PI.txt',
        s = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/S.txt',
        beta = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/celltypebeta.txt',
        V = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/V.txt',
        W = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/W.txt',
    script: '../bin/sim/celltype_expectedPInSnBETAnVnW.py'


rule sim_generatedata:
    input:
        beta = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/celltypebeta.txt',
        V = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/V.txt',
        W = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/W.txt',
    output:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/sim.batch{{i}}.npy',
    params:
        batches = sim_batches,
        beta = (0.5, 0.5), # beta distribution for allele frequency
        maf = 0.05,
        L = 10, # number of causal SNPs
        seed = 273672,
    script: '../bin/sim/generatedata.py'


rule sim_HE:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/sim.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.batch{{i}}.npy',
    resources:
        # mem_mb = lambda wildcards: '10G' if wildcards.model != 'full' else '150G',
        mem_mb = lambda wildcards: '10G' if int(wildcards.ss) <= 1000 and len(wildcards.a.split('_')) <=4 and wildcards.model != 'full' else '60G',
    params:
        free_jk = True,
        full = lambda wildcards: True if int(wildcards.ss) <= 1000 and len(wildcards.a.split('_')) <=4 and wildcards.model != 'full' else False
    script: '../bin/sim/he.py'


rule sim_mergeBatches_HE:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.batch{i}.npy' 
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/out.he.npy',
    script: '../bin/mergeBatches.py'


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
    return expand('analysis/sim/{{model}}/{params}/sim.npy', params=subspace.instance_patterns)


def sim_agg_he_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/out.he.npy', params=subspace.instance_patterns)


rule sim_agg_he_out:
    input:
        out = sim_agg_he_out_subspace,
    output:
        out = 'analysis/sim/{model}/AGG{arg}.he.npy',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[sim_params['model']==wildcards.model]).iloc[:,:],
    run:
        args = np.array(params.subspace[wildcards.arg])
        data = {}
        for arg, out in zip(args, input.out):
            data[arg] = np.load(out, allow_pickle=True).item()
        np.save(output.out, data)




############################
# 1.2 HE without JK 
############################
use rule sim_HE as sim_HE_noJK with:
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.noJK.batch{{i}}.npy',
    params:
        free_jk = False,


use rule sim_mergeBatches_HE as sim_mergeBatches_HE_noJK with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.noJK.batch{i}.npy' 
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/out.he.noJK.npy',


def sim_agg_he_noJK_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/out.he.noJK.npy', params=subspace.instance_patterns)


use rule sim_agg_he_out as sim_agg_he_noJK_out with:
    input:
        out = sim_agg_he_noJK_out_subspace,
    output:
        out = 'analysis/sim/{model}/AGG{arg}.he.noJK.npy',


############################
# 1.2 HE without Full
############################
use rule sim_HE as sim_HE_noFull with:
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.noFull.batch{{i}}.npy',
    params:
        free_jk = True,
        full = False,


use rule sim_mergeBatches_HE as sim_mergeBatches_HE_noFull with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.noFull.batch{i}.npy' 
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/out.he.noFull.npy',


def sim_agg_he_noFull_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/out.he.noFull.npy', params=subspace.instance_patterns)


use rule sim_agg_he_out as sim_agg_he_noFull_out with:
    input:
        out = sim_agg_he_noFull_out_subspace,
    output:
        out = 'analysis/sim/{model}/AGG{arg}.he.noFull.npy',
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
    return (expand('analysis/sim/{{model}}/AGG{arg}.he.npy', arg=effective_args) 
    + expand('analysis/sim/{{model}}/AGG{arg}.true_V.npy', arg=effective_args))


rule sim_HE_AGGarg:
    input:
        sim_HE_AGGarg_fun,
    output:
        flag = touch('staging/sim/{model}/HE.all.flag'),


rule sim_HE_all:
    input:
        flag = expand('staging/sim/{model}/HE.all.flag', model=['free']),


# cell type number
rule sim_celltype_number_all:
    input:
        out1 = 'analysis/sim/free21/AGGss.he.npy',
        out2 = 'analysis/sim/free22/AGGss.he.noFull.npy',
    output:
        flag = touch('analysis/sim/celltype_number.flag'),





















###################################################################
# when nu is noisy
###################################################################
rule sim_data_add_noise_to_nu:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/sim.batch{{i}}.npy',
    output:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/sim.noise_{{alpha}}.batch{{i}}.npy',
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
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/sim.noise_{{alpha}}.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.noise_{{alpha}}.batch{{i}}.npy',


use rule sim_mergeBatches_HE as sim_mergeBatches_HE_noisy_nu with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.noise_{{alpha}}.batch{i}.npy' 
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/out.he.noise_{{alpha}}.npy',


rule sim_agg_he_noisy_nu_out_alpha:
    input:
        out = [f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/out.he.noise_{alpha}.npy'
                for alpha in config['sim']['noise_alpha']]
    output:
        out = touch(f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.noisy_nu.flag',)


def sim_agg_he_noisy_nu_out(wildcards):
    subspace1 = get_subspace('ss', sim_params.loc[sim_params['model']=='hom'].head(1))
    subspace2 = get_subspace('ss', sim_params.loc[sim_params['model']=='free'].head(1))
    return (expand('staging/sim/hom/{params}/he.noisy_nu.flag', params=subspace1.instance_patterns) + 
            expand('staging/sim/free/{params}/he.noisy_nu.flag', params=subspace2.instance_patterns))


rule sim_agg_he_noisy_nu_out_subspace:
    input:
        out = sim_agg_he_noisy_nu_out,
    output:
        out = touch('staging/sim/he.noisy_nu.flag'),



















###################################################################
# when nu is not considered when fitting model
###################################################################
#####################################
# CIGMA missing nu
#####################################
rule sim_HE_missing_nu:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/sim.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.missing_nu.batch{{i}}.npy',
    resources:
        mem_mb = lambda wildcards: '10G' if int(wildcards.ss) <= 1000 and len(wildcards.a.split('_')) <=4 and wildcards.model != 'full' else '60G',
    params:
        free_jk = False,
    script: '../scripts/sim/he.missing_nu.py'


use rule sim_mergeBatches_HE as sim_mergeBatches_HE_missing_nu with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.missing_nu.batch{i}.npy'
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/out.he.missing_nu.npy',


def sim_agg_he_out_subspace_missing_nu(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/out.he.missing_nu.npy', 
                    params=subspace.instance_patterns)


rule sim_agg_he_out_missing_nu:
    input:
        out = sim_agg_he_out_subspace_missing_nu,
    output:
        out = 'analysis/sim/{model}/AGG{arg}.he.missing_nu.npy',
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
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/sim.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/greml.batch{{i}}.npy',
    resources:
        mem_mb = lambda wildcards: '2G' if int(wildcards.ss) <= 2000 else '60G',
    shell:
        '''
        module load gcc/11.3.0 gcta/1.94.1
        python3 workflow/bin/sim/gcta_greml_ctp.py \
                {input.data} {output.out} \
        '''


use rule sim_mergeBatches_HE as sim_gcta_greml_ctp_merge with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/greml.batch{i}.npy'
                for i in range(config['sim']['batch_no'])],
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/greml.npy',


def sim_agg_greml_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('staging/sim/{{model}}/{params}/greml.npy', params=subspace.instance_patterns)


rule sim_gcta_greml_ctp_agg:
    input:
        out = sim_agg_greml_out_subspace,
    output:
        out = 'analysis/sim/{model}/AGG{arg}.greml.npy',
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
        gcta = 'analysis/sim/{model}/AGG{arg}.greml.npy',
        missing_nu_out = 'analysis/sim/{model}/AGG{arg}.he.missing_nu.npy',
        out = 'analysis/sim/{model}/AGG{arg}.he.npy',
    output:
        flag = touch('analysis/sim/{model}/AGG{arg}.gcta.flag'),






rule sim_all:
    input:
        hom = 'staging/sim/hom/HE.all.flag',
        free = 'staging/sim/free/HE.all.flag',
        gcta = 'analysis/sim/free3/AGGvc.gcta.flag',
        ct_num = 'analysis/sim/celltype_number.flag',
