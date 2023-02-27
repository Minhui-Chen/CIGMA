from snakemake.utils import Paramsapce 
import numpy as np, pandas as pd


#########################################################################################
##   Simulation
#########################################################################################
# par
replicates = 100
batch_no = 10
batches = np.array_split(range(replicates), batch_no)

## paramspace
sim_params = pd.read_table("sim.params.txt", dtype="str", comment='#', na_filter=False)
if sim_params.shape[0] != sim_params.drop_duplicates().shape[0]:
    sys.exit('Duplicated parameters!\n')
sim_par_columns = list(sim_params.columns)
sim_par_columns.remove('model')
sim_paramspace = Paramspace(sim_params[sim_par_columns], filename_params="*")

sim_plot_order = {
    'hom':{
        'ss':['2e1', '5e1', '1e2', '2e2', '3e2', '5e2', '1e3'], 'a':['0.5_2_2_2', '1_2_2_2', '2_2_2_2', '4_2_2_2']
        },
    'iid':{
        'ss':['2e1', '5e1', '1e2', '5e2', '1e3'], 'a':['0.5_2_2_2', '1_2_2_2', '2_2_2_2', '4_2_2_2'],
        'vc':['0.3_0.1_0.1_0.1_0.1_0.3', '0.2_0.1_0.1_0.2_0.2_0.2', '0.1_0.1_0.1_0.3_0.3_0.1'],
        }, 
    'free': {
        'ss':['5e0', '2e1', '5e1', '1e2', '2e2', '3e2','5e2', '1e3'], 'a':['0.5_2_2_2', '1_2_2_2', '2_2_2_2', '4_2_2_2'], 
        'vc':['0.35_0.1_0.1_0.05_0.05_0.35', '0.3_0.1_0.1_0.1_0.1_0.3', '0.2_0.1_0.1_0.2_0.2_0.2', '0.1_0.1_0.1_0.3_0.3_0.1'],
        'V_diag':['1_1_1_1', '8_4_2_1', '27_9_3_1', '64_16_4_1'],
        },
    'full':{
        'ss':['2e1', '5e1', '1e2', '3e2', '5e2', '1e3'], 'a':['0.5_2_2_2', '1_2_2_2', '2_2_2_2', '4_2_2_2'],
        'vc':['0.3_0.1_0.1_0.1_0.1_0.3', '0.2_0.1_0.1_0.2_0.2_0.2', '0.1_0.1_0.1_0.3_0.3_0.1'],
        'V_diag':['1_1_1_1', '8_4_2_1', '27_9_3_1', '64_16_4_1', '64_64_1_1'],
        'V_tril':['0.25_0.25_0_-0.25_0_0', '0.5_0.5_0_-0.5_0_0', '0.75_0.75_0_-0.75_0_0', '0.95_0.95_0.95_-0.95_-0.95_-0.95']
        },
    }

rule ctg_HE:
    input:
        cty = f'staging/og/{{model}}/{sim_paramspace.wildcard_pattern}/cty.batch{{i}}.txt',
        Z = f'staging/og/{{model}}/{sim_paramspace.wildcard_pattern}/Z.batch{{i}}.txt',
        ctnu = f'staging/og/{{model}}/{sim_paramspace.wildcard_pattern}/ctnu.batch{{i}}.txt',
    output:
        out = f'staging/ctg/{{model}}/{sim_paramspace.wildcard_pattern}/he.batch{{i}}',
    params:
        out = f'staging/ctg/{{model}}/{sim_paramspace.wildcard_pattern}/rep/he.npy',
        batches = lambda wildcards: og_batches[int(wildcards.i)],
    resources:
        time = '10:00:00',
        mem_per_cpu = '5000',
    priority: 1
    script: 'bin/ctg_HE.py'

rule ctg_REML:
    input:
        cty = f'staging/og/{{model}}/{sim_paramspace.wildcard_pattern}/cty.batch{{i}}.txt',
        Z = f'staging/og/{{model}}/{sim_paramspace.wildcard_pattern}/Z.batch{{i}}.txt',
        P = f'staging/og/{{model}}/{sim_paramspace.wildcard_pattern}/P.batch{{i}}.txt',
        ctnu = f'staging/og/{{model}}/{sim_paramspace.wildcard_pattern}/ctnu.batch{{i}}.txt',
    output:
        out = f'staging/ctg/{{model}}/{sim_paramspace.wildcard_pattern}/reml.batch{{i}}',
    params:
        out = f'staging/ctg/{{model}}/{sim_paramspace.wildcard_pattern}/rep/reml.npy',
        batches = lambda wildcards: og_batches[int(wildcards.i)],
    resources:
        time = '100:00:00',
        #mem_per_cpu = '5000',
    priority: 1
    script: 'bin/ctg_REML.py'

rule ctg_mergeBatches:
    input:
        he = [f'staging/ctg/{{model}}/{sim_paramspace.wildcard_pattern}/he.batch{i}' for i in range(og_batch_no)],
        reml = [f'staging/ctg/{{model}}/{sim_paramspace.wildcard_pattern}/reml.batch{i}' for i in range(og_batch_no)],
    output:
        out = f'analysis/ctg/{{model}}/{sim_paramspace.wildcard_pattern}/out.npy',
    script: 'bin/og_mergeBatches.py'

def ctg_agg_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/ctg/{{model}}/{params}/out.npy', params=subspace.instance_patterns)

rule ctg_HEestimates_subspace_plot:
    input:
        out = ctg_agg_out_subspace,
        V = og_agg_trueV_subspace,
        W = og_agg_trueW_subspace,
    output:
        png = 'results/ctg/{model}/HE.AGG{arg}.png',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[sim_params['model']==wildcards.model]).iloc[:,:],
        sim_plot_order = sim_plot_order,
        mycolors = mycolors,
        pointcolor = pointcolor,
        colorpalette = colorpalette,
    script: 'bin/ctg_HEestimates_subspace_plot.py'

rule ctg_HEwald_subspace_plot:
    input:
        out = ctg_agg_out_subspace, 
    output:
        png = 'results/ctg/{model}/HE.wald.AGG{arg}.png',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[sim_params['model']==wildcards.model]).iloc[:,:],
        sim_plot_order = sim_plot_order,
        mycolors = mycolors,
    script: 'bin/ctg_HEwald_subspace_plot.py'

rule ctg_REMLestimates_subspace_plot:
    input:
        out = ctg_agg_out_subspace,
        beta = og_agg_truebeta_subspace,
        V = og_agg_trueV_subspace,
        W = og_agg_trueW_subspace,
    output:
        png = 'results/ctg/{model}/REML.AGG{arg}.png',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[sim_params['model']==wildcards.model]).iloc[:,:],
        sim_plot_order = sim_plot_order,
        colorpalette = colorpalette,
        pointcolor = pointcolor,
        mycolors = mycolors,
    script: 'bin/ctg_REMLestimates_subspace_plot.py'

rule ctg_REMLwaldNlrt_subspace_plot:
    input:
        out = ctg_agg_out_subspace, 
    output:
        png = 'results/ctg/{model}/REML.waldNlrt.AGG{arg}.png',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[sim_params['model']==wildcards.model]).iloc[:,:],
        sim_plot_order = sim_plot_order,
        mycolors = mycolors,
    script: 'bin/ctg_REMLwaldNlrt_subspace_plot.py'

rule ctg_collect_HEnREML_subspace_plot:
    input:
        he = 'results/ctg/{model}/HE.AGG{arg}.png',
        he_wald = 'results/ctg/{model}/HE.wald.AGG{arg}.png',
        reml = 'results/ctg/{model}/REML.AGG{arg}.png', # reml is too time-consuming, so not included for now
        reml_waldNlrt = 'results/ctg/{model}/REML.waldNlrt.AGG{arg}.png',
    output:
        flag = touch('staging/ctg/{model}/HEnREML.AGG{arg}.flag'),

def ctg_HEnREML_AGGarg_fun(wildcards):
    effective_args = get_effective_args(sim_params.loc[sim_params['model']==wildcards.model])
    return expand('staging/ctg/{{model}}/HEnREML.AGG{arg}.flag', arg=effective_args)

rule ctg_AGGarg:
    input:
        ctg_HEnREML_AGGarg_fun,
    output:
        flag = touch('staging/ctg/{model}/all.flag'),

rule ctg_all:
    input:
        flag = expand('staging/ctg/{model}/all.flag', model=['hom','free']),
