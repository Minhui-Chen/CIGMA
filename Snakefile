from snakemake.utils import Paramspace 
import numpy as np, pandas as pd
import seaborn as sns

# make a logs folders to save log files
os.makedirs('./logs/', exist_ok=True)

# par
colorpalette='bright'
mycolors = sns.color_palette(colorpalette)
pointcolor = 'red' # color for expected values  in estimates plots

def get_subspace(arg, model_params):
    ''' model_params df include or not include the column of model, but the first row is the basemodel'''
    sim_args = list(model_params.columns)
    if 'model' in sim_args:
        sim_args.remove('model')
    model_params = model_params[sim_args].reset_index(drop=True)
    #print(model_params)
    basemodel = model_params.iloc[0].to_dict()
    sim_args.remove(arg)
    evaluate = ' & '.join([f"{_} == '{basemodel[_]}'" for _ in sim_args])
    subspace = model_params[model_params.eval(evaluate)]
    return Paramspace(subspace, filename_params="*")

def get_effective_args(model_params):
    ''' model_params df include or not include the column of model, but the first row is the basemodel'''
    #print(model_params)
    sim_args = list(model_params.columns)
    if 'model' in sim_args:
        sim_args.remove('model')
    model_params = model_params[sim_args]
    effective_args = [] # args with multiple parameters. e.g. ss has '1e2', '2e2', '5e2'
    for arg_ in np.array(model_params.columns):
        subspace = get_subspace(arg_, model_params)
        if subspace.shape[0] > 1:
            effective_args.append(arg_)
    return effective_args

# wildcard constraints
wildcard_constraints: i='[0-9]+'
wildcard_constraints: model='\w+'


#########################################################################################
##   Simulation
#########################################################################################
# par
sim_replicates = 100
sim_batch_no = 50
sim_batches = np.array_split(range(sim_replicates), sim_batch_no)

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
        'vc':['0.35_0.1_0.1_0.05_0.05_0.35', '0.3_0.1_0.1_0.1_0.1_0.3', '0.2_0.1_0.1_0.2_0.2_0.2', 
            '0.1_0.1_0.1_0.3_0.3_0.1', '0.3_0.1_0.1_0.05_0.15_0.3', '0.3_0.1_0.1_0.15_0.05_0.3',
            '0.3_0.1_0.1_0.2_0_0.3'],
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
    script: 'bin/sim/celltype_expectedPInSnBETAnVnW.py'

rule sim_generatedata:
    input:
        beta = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/celltypebeta.txt',
        V = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/V.txt',
        W = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/W.txt',
    output:
        G = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/G.batch{{i}}.txt',
        K = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/K.batch{{i}}.txt',
        P = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/P.batch{{i}}.txt',
        pi = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/estPI.batch{{i}}.txt',
        s = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/estS.batch{{i}}.txt',
        nu = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/nu.batch{{i}}.txt',
        ctnu = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/ctnu.batch{{i}}.txt',
        y = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/y.batch{{i}}.txt',
        Y = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/Y.batch{{i}}.txt',
    params:
        batches = sim_batches,
        beta = (0.5, 0.5), # beta distribution for allele frequency
        maf = 0.05,
        L = 10, # number of causal SNPs
        G = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/repX/G.txt.gz',
        K = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/repX/K.txt.gz',
        P = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/repX/P.txt.gz',
        pi = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/repX/estPI.txt.gz',
        s = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/repX/estS.txt.gz',
        nu = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/repX/nu.txt.gz',
        ctnu = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/repX/ctnu.txt.gz',
        y = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/repX/y.txt.gz',
        Y = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/repX/Y.txt.gz',
    script: 'bin/sim/generatedata.py'

rule sim_HE:
    input:
        Y = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/Y.batch{{i}}.txt',
        K = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/K.batch{{i}}.txt',
        ctnu = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/ctnu.batch{{i}}.txt',
        P = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/P.batch{{i}}.txt',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.batch{{i}}',
    params:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/rep/he.npy',
        batches = lambda wildcards: sim_batches[int(wildcards.i)],
    resources:
        time = '10:00:00',
        mem_per_cpu = '5000',
    priority: 1
    script: 'bin/sim/he.py'

rule sim_REML:
    input:
        Y = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/Y.batch{{i}}.txt',
        K = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/K.batch{{i}}.txt',
        ctnu = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/ctnu.batch{{i}}.txt',
        P = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/P.batch{{i}}.txt',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/reml.batch{{i}}',
    params:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/rep/reml.npy',
        batches = lambda wildcards: sim_batches[int(wildcards.i)],
    resources:
        time = '100:00:00',
        #mem_per_cpu = '5000',
    priority: 1
    script: 'bin/sim/reml.py'

rule sim_mergeBatches:
    input:
        he = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.batch{i}' for i in range(sim_batch_no)],
        reml = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/reml.batch{i}' for i in range(sim_batch_no)],
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/out.npy',
    script: 'bin/sim/mergeBatches.py'

########################################################################
# real data preparation
########################################################################
include: 'real.snake'

########################################################################
# Yazar 2022 Science
########################################################################
include: 'Yazar.snake'

rule yazar_impute_ctp:
    input:
        data = 'data/Yazar2022Science/ctp.gz',
    output:
        data = 'analysis/Yazar2022Science/data/ctp.gz',
    params:
        seed = 1234567,
    resources:
        mem_per_cpu = '5000',
        time = '12:00:00',
    run:
        from ctmm import preprocess
        ctp = pd.read_table(input.data, index_col=(0,1))
        ctp = preprocess.softimpute(ctp, seed=params.seed)
        ctp.to_csv(output.data, sep='\t')

use rule yazar_impute_ctp as yazar_impute_ctnu with:
    input:
        data = 'data/Yazar2022Science/ctnu.gz',
    output:
        data = 'analysis/Yazar2022Science/data/ctnu.gz',
    params:
        seed = 7654321,

rule yazar_std_op:
    input:
        ctp = 'analysis/Yazar2022Science/data/ctp.gz',
        ctnu = 'analysis/Yazar2022Science/data/ctnu.gz',
        P = 'data/Yazar2022Science/P.gz',
    output:
        op = 'staging/Yazar2022Science/data/op.gz',
        nu = 'staging/Yazar2022Science/data/nu.gz',
        ctp = 'staging/Yazar2022Science/data/ctp.gz',
        ctnu = 'staging/Yazar2022Science/data/ctnu.gz',
    run:
        from ctmm import preprocess
        ctp = pd.read_table(snakemake.input.ctp, index_col=(0,1))
        ctnu = pd.read_table(snakemake.input.ctnu, index_col=(0,1))
        P = pd.read_table(snakemake.input.P, index_col=0)

        op, nu, ctp, ctnu = preprocess.std(ctp, ctnu, P)

        op.to_csv(snakemake.output.op, sep='\t')
        nu.to_csv(snakemake.output.nu, sep='\t')
        ctp.to_csv(snakemake.output.ctp, sep='\t')
        ctnu.to_csv(snakemake.output.ctnu, sep='\t')
