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
        'ss':['2e1', '5e1', '1e2', '2e2', '3e2', '5e2', '1e3'], 
        'a':['0.5_2_2_2', '1_2_2_2', '2_2_2_2', '4_2_2_2']
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

rule sim_celltype_expectedPInSnBETAnVnW:
    output:
        pi = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/PI.txt',
        s = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/S.txt',
        beta = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/celltypebeta.txt',
        V = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/V.txt',
        W = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/W.txt',
    script: 'bin/sim_celltype_expectedPInSnBETAnVnW.py'

rule sim_generatedata_batch:
    input:
        beta = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/celltypebeta.txt',
        V = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/V.txt',
        W = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/W.txt',
    output:
        G = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/G.batch{{i}}.txt',
        Z = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/Z.batch{{i}}.txt',
        P = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/P.batch{{i}}.txt',
        pi = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/estPI.batch{{i}}.txt',
        s = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/estS.batch{{i}}.txt',
        nu = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/nu.batch{{i}}.txt',
        ctnu = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/ctnu.batch{{i}}.txt',
        y = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/y.batch{{i}}.txt',
        cty = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/cty.batch{{i}}.txt',
    params:
        batches = sim_batches,
        beta = (0.5, 0.5), # beta distribution for allele frequency
        maf = 0.05,
        L = 500, # number of causal SNPs
        G = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/repX/G.txt.gz',
        Z = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/repX/Z.txt.gz',
        P = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/repX/P.txt.gz',
        pi = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/repX/estPI.txt.gz',
        s = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/repX/estS.txt.gz',
        nu = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/repX/nu.txt.gz',
        ctnu = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/repX/ctnu.txt.gz',
        y = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/repX/y.txt.gz',
        cty = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/repX/cty.txt.gz',
    script: 'bin/sim_generatedata_batch.py'

rule sim_HE:
    input:
        cty = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/cty.batch{{i}}.txt',
        Z = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/Z.batch{{i}}.txt',
        ctnu = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/ctnu.batch{{i}}.txt',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.batch{{i}}',
    params:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/rep/he.npy',
        batches = lambda wildcards: sim_batches[int(wildcards.i)],
    resources:
        time = '10:00:00',
        mem_per_cpu = '5000',
    priority: 1
    script: 'bin/sim_HE.py'

rule sim_REML:
    input:
        cty = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/cty.batch{{i}}.txt',
        Z = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/Z.batch{{i}}.txt',
        P = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/P.batch{{i}}.txt',
        ctnu = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/ctnu.batch{{i}}.txt',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/reml.batch{{i}}',
    params:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/rep/reml.npy',
        batches = lambda wildcards: sim_batches[int(wildcards.i)],
    resources:
        time = '100:00:00',
        #mem_per_cpu = '5000',
    priority: 1
    script: 'bin/sim_REML.py'

rule sim_mergeBatches:
    input:
        he = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.batch{i}' for i in range(sim_batch_no)],
        reml = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/reml.batch{i}' for i in range(sim_batch_no)],
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/out.npy',
    script: 'bin/sim_mergeBatches.py'

def sim_agg_truebeta_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/celltypebeta.txt', params=subspace.instance_patterns)

def sim_agg_trueV_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/V.txt', params=subspace.instance_patterns)

def sim_agg_trueW_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/W.txt', params=subspace.instance_patterns)

def sim_agg_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/out.npy', params=subspace.instance_patterns)

rule sim_HEestimates_subspace_plot:
    input:
        out = sim_agg_out_subspace,
        V = sim_agg_trueV_subspace,
        W = sim_agg_trueW_subspace,
    output:
        png = 'results/sim/{model}/HE.AGG{arg}.png',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[sim_params['model']==wildcards.model]).iloc[:,:],
        sim_plot_order = sim_plot_order,
        mycolors = mycolors,
        pointcolor = pointcolor,
        colorpalette = colorpalette,
    script: 'bin/sim_HEestimates_subspace_plot.py'

rule sim_HEwald_subspace_plot:
    input:
        out = sim_agg_out_subspace, 
    output:
        png = 'results/sim/{model}/HE.wald.AGG{arg}.png',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[sim_params['model']==wildcards.model]).iloc[:,:],
        sim_plot_order = sim_plot_order,
        mycolors = mycolors,
    script: 'bin/sim_HEwald_subspace_plot.py'

rule sim_REMLestimates_subspace_plot:
    input:
        out = sim_agg_out_subspace,
        beta = sim_agg_truebeta_subspace,
        V = sim_agg_trueV_subspace,
        W = sim_agg_trueW_subspace,
    output:
        png = 'results/sim/{model}/REML.AGG{arg}.png',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[sim_params['model']==wildcards.model]).iloc[:,:],
        sim_plot_order = sim_plot_order,
        colorpalette = colorpalette,
        pointcolor = pointcolor,
        mycolors = mycolors,
    script: 'bin/sim_REMLestimates_subspace_plot.py'

rule sim_REMLwaldNlrt_subspace_plot:
    input:
        out = sim_agg_out_subspace, 
    output:
        png = 'results/sim/{model}/REML.waldNlrt.AGG{arg}.png',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[sim_params['model']==wildcards.model]).iloc[:,:],
        sim_plot_order = sim_plot_order,
        mycolors = mycolors,
    script: 'bin/sim_REMLwaldNlrt_subspace_plot.py'

rule sim_collect_HEnREML_subspace_plot:
    input:
        he = 'results/sim/{model}/HE.AGG{arg}.png',
        he_wald = 'results/sim/{model}/HE.wald.AGG{arg}.png',
        reml = 'results/sim/{model}/REML.AGG{arg}.png', # reml is too time-consuming, so not included for now
        reml_waldNlrt = 'results/sim/{model}/REML.waldNlrt.AGG{arg}.png',
    output:
        flag = touch('staging/sim/{model}/HEnREML.AGG{arg}.flag'),

def sim_HEnREML_AGGarg_fun(wildcards):
    effective_args = get_effective_args(sim_params.loc[sim_params['model']==wildcards.model])
    return expand('staging/sim/{{model}}/HEnREML.AGG{arg}.flag', arg=effective_args)

rule sim_AGGarg:
    input:
        sim_HEnREML_AGGarg_fun,
    output:
        flag = touch('staging/sim/{model}/all.flag'),

rule sim_all:
    input:
        flag = expand('staging/sim/{model}/all.flag', model=['hom','free']),
