#########################################################################################
# overall pseudo-bulk genetic
#########################################################################################
# par
og_replicates = 100
og_batch_no = 10
og_batches = np.array_split(range(og_replicates), og_batch_no)
#og_batchsize = 1
#og_batches = [range(i, min(i+og_batchsize, og_replicates)) 
#        for i in range(0, og_replicates, og_batchsize)]

## declare a dataframe to be a paramspace
og_params = pd.read_table("og.params.txt", dtype="str", comment='#', na_filter=False)
if og_params.shape[0] != og_params.drop_duplicates().shape[0]:
    sys.exit('Duplicated parameters!\n')
og_par_columns = list(og_params.columns)
og_par_columns.remove('model')
og_paramspace = Paramspace(og_params[og_par_columns], filename_params="*")

og_plot_order = {
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

rule og_celltype_expectedPInSnBETAnVnW:
    output:
        pi = f'analysis/og/{{model}}/{og_paramspace.wildcard_pattern}/PI.txt',
        s = f'analysis/og/{{model}}/{og_paramspace.wildcard_pattern}/S.txt',
        beta = f'analysis/og/{{model}}/{og_paramspace.wildcard_pattern}/celltypebeta.txt',
        V = f'analysis/og/{{model}}/{og_paramspace.wildcard_pattern}/V.txt',
        W = f'analysis/og/{{model}}/{og_paramspace.wildcard_pattern}/W.txt',
    script: 'bin/og_celltype_expectedPInSnBETAnVnW.py'

rule og_generatedata_batch:
    input:
        beta = f'analysis/og/{{model}}/{og_paramspace.wildcard_pattern}/celltypebeta.txt',
        V = f'analysis/og/{{model}}/{og_paramspace.wildcard_pattern}/V.txt',
        W = f'analysis/og/{{model}}/{og_paramspace.wildcard_pattern}/W.txt',
    output:
        G = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/G.batch{{i}}.txt',
        Z = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/Z.batch{{i}}.txt',
        P = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/P.batch{{i}}.txt',
        pi = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/estPI.batch{{i}}.txt',
        s = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/estS.batch{{i}}.txt',
        nu = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/nu.batch{{i}}.txt',
        ctnu = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/ctnu.batch{{i}}.txt',
        y = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/y.batch{{i}}.txt',
        cty = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/cty.batch{{i}}.txt',
    params:
        batches = og_batches,
        beta = (0.5, 0.5), # beta distribution for allele frequency
        maf = 0.05,
        L = 500, # number of causal SNPs
        G = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/repX/G.txt.gz', 
        Z = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/repX/Z.txt.gz', 
        P = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/repX/P.txt.gz', 
        pi = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/repX/estPI.txt.gz', 
        s = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/repX/estS.txt.gz', 
        nu = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/repX/nu.txt.gz', 
        ctnu = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/repX/ctnu.txt.gz', 
        y = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/repX/y.txt.gz', 
        cty = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/repX/cty.txt.gz', 
    script: 'bin/og_generatedata_batch.py'

rule og_HE:
    input:
        y = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/y.batch{{i}}.txt',
        G = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/G.batch{{i}}.txt',
        Z = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/Z.batch{{i}}.txt',
        P = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/P.batch{{i}}.txt',
        nu = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/nu.batch{{i}}.txt',
    output:
        out = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/he.batch{{i}}',
    params:
        out = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/rep/he.npy',
        batches = lambda wildcards: og_batches[int(wildcards.i)],
    resources:
        time = '10:00:00',
        mem_per_cpu = '5000',
    priority: 1
    script: 'bin/og_HE.py'

rule og_REML:
    input:
        y = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/y.batch{{i}}.txt',
        G = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/G.batch{{i}}.txt',
        Z = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/Z.batch{{i}}.txt',
        P = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/P.batch{{i}}.txt',
        nu = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/nu.batch{{i}}.txt',
    output:
        out = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/reml.batch{{i}}',
    params:
        out = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/rep/reml.npy',
        batches = lambda wildcards: og_batches[int(wildcards.i)],
    resources:
        time = '100:00:00',
        mem_per_cpu = '5000',
    priority: 1
    script: 'bin/og_REML.R'

rule og_mergeBatches:
    input:
        he = [f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/he.batch{i}' for i in range(og_batch_no)],
        reml = [f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/reml.batch{i}' for i in range(og_batch_no)],
    output:
        out = f'analysis/og/{{model}}/{og_paramspace.wildcard_pattern}/out.npy',
    script: 'bin/og_mergeBatches.py'

#localrules: og_rmFiles
#rule og_rmFiles:
#    input:
#        out = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/out.batch{{i}}',
#        # to compare different ways of REML
#        #reml = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/reml.merge.batch{{i}}',
#        G = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/G.batch{{i}}.txt',
#        Z = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/Z.batch{{i}}.txt',
#    output:
#        flag = touch(f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/rmFiles.batch{{i}}.flag'),
#    priority: 1
#    run:
#        # remove G Z files manually after ML and REML
#        for G_f, Z_f in zip([line.strip() for line in open(input.G)], [line.strip() for line in open(input.Z)]):
#            os.remove(G_f)
#            os.remove(Z_f)

#rule og_MLnREMLnHE_aggReplications:
#    input:
#        s = [f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/estS.batch{i}.txt'
#                for i in range(len(og_batches))],
#        nu = [f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/nu.batch{i}.txt' 
#                for i in range(len(og_batches))],
#        pi = [f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/estPI.batch{i}.txt' 
#                for i in range(len(og_batches))],
#        out = [f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/out.batch{i}' 
#                for i in range(len(og_batches))],
#        flag = [f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/rmFiles.batch{i}.flag'
#                for i in range(len(og_batches))],
#    output:
#        s = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/estS.txt',
#        nu = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/nu.txt',
#        pi = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/estPI.txt',
#        out = f'analysis/og/{{model}}/{og_paramspace.wildcard_pattern}/out.npy',
#    script: "bin/ong_waldNlrt_aggReplications.py"

def og_agg_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, og_params.loc[og_params['model']==wildcards.model])
    return expand('analysis/og/{{model}}/{params}/out.npy', params=subspace.instance_patterns)

def og_agg_truebeta_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, og_params.loc[og_params['model']==wildcards.model])
    return expand('analysis/og/{{model}}/{params}/celltypebeta.txt', params=subspace.instance_patterns)

def og_agg_trueV_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, og_params.loc[og_params['model']==wildcards.model])
    return expand('analysis/og/{{model}}/{params}/V.txt', params=subspace.instance_patterns)

def og_agg_trueW_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, og_params.loc[og_params['model']==wildcards.model])
    return expand('analysis/og/{{model}}/{params}/W.txt', params=subspace.instance_patterns)

#rule og_MLestimates_subspace_plot:
#    input:
#        out = og_agg_out_subspace,
#        beta = og_agg_truebeta_subspace,
#        V = og_agg_trueV_subspace,
#        W = og_agg_trueW_subspace,
#    output:
#        png = 'analysis/og/{model}/ML.AGG{arg}.png',
#        flag = touch('staging/og/{model}/ML.AGG{arg}.flag'),
#    params: 
#        subspace = lambda wildcards: get_subspace(wildcards.arg,
#                og_params.loc[og_params['model']==wildcards.model]).iloc[:,:],
#        og_plot_order = og_plot_order,
#        colorpalette = colorpalette,
#        pointcolor = pointcolor,
#        mycolors = mycolors,
#    script: 'bin/og_MLestimates_subspace_plot.py'
#
#rule og_MLwaldNlrt_subspace_plot:
#    input:
#        out = og_agg_out_subspace, 
#    output:
#        waldNlrt = 'analysis/og/{model}/ML.waldNlrt.AGG{arg}.png',
#    params:
#        subspace = lambda wildcards: get_subspace(wildcards.arg,
#                og_params.loc[og_params['model']==wildcards.model]).iloc[:,:],
#        og_plot_order = og_plot_order,
#        mycolors = mycolors,
#    script: "bin/og_MLwaldNlrt_subspace_plot.py"

rule og_HEestimates_subspace_plot:
    input:
        out = og_agg_out_subspace,
        V = og_agg_trueV_subspace,
        W = og_agg_trueW_subspace,
    output:
        png = 'results/og/{model}/HE.AGG{arg}.png',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                og_params.loc[og_params['model']==wildcards.model]).iloc[:,:],
        og_plot_order = og_plot_order,
        mycolors = mycolors,
        pointcolor = pointcolor,
        colorpalette = colorpalette,
    script: 'bin/ctg_HEestimates_subspace_plot.py'

rule og_HEwald_subspace_plot:
    input:
        out = og_agg_out_subspace, 
    output:
        png = 'results/og/{model}/HE.wald.AGG{arg}.png',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                og_params.loc[og_params['model']==wildcards.model]).iloc[:,:],
        og_plot_order = og_plot_order,
        mycolors = mycolors,
    script: 'bin/ctg_HEwald_subspace_plot.py'

rule og_REMLestimates_subspace_plot:
    input:
        out = og_agg_out_subspace,
        beta = og_agg_truebeta_subspace,
        V = og_agg_trueV_subspace,
        W = og_agg_trueW_subspace,
    output:
        png = 'results/og/{model}/REML.AGG{arg}.png',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                og_params.loc[og_params['model']==wildcards.model]).iloc[:,:],
        og_plot_order = og_plot_order,
        colorpalette = colorpalette,
        pointcolor = pointcolor,
        mycolors = mycolors,
    script: 'bin/ctg_REMLestimates_subspace_plot.py'

rule og_REMLwaldNlrt_subspace_plot:
    input:
        out = og_agg_out_subspace, 
    output:
        png = 'results/og/{model}/REML.waldNlrt.AGG{arg}.png',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                og_params.loc[og_params['model']==wildcards.model]).iloc[:,:],
        og_plot_order = og_plot_order,
        mycolors = mycolors,
    script: 'bin/ctg_REMLwaldNlrt_subspace_plot.py'

rule og_collect_HEnREML_subspace_plot:
    input:
        he = 'results/og/{model}/HE.AGG{arg}.png',
        he_wald = 'results/og/{model}/HE.wald.AGG{arg}.png',
        reml = 'results/og/{model}/REML.AGG{arg}.png',
        reml_waldNlrt = 'results/og/{model}/REML.waldNlrt.AGG{arg}.png',
    output:
        flag = touch('staging/og/{model}/HEnREML.AGG{arg}.flag'),

#def og_MLwaldNlrt_AGGarg_fun(wildcards):
#    effective_args = get_effective_args(og_params.loc[og_params['model']==wildcards.model])
#    return expand('analysis/og/{{model}}/ML.waldNlrt.AGG{arg}.png', arg=effective_args)
#
#def og_MLestimates_AGGarg_fun(wildcards):
#    effective_args = get_effective_args(og_params.loc[og_params['model']==wildcards.model])
#    return expand('staging/og/{{model}}/ML.AGG{arg}.flag', arg=effective_args)

#def og_HEestimates_AGGarg_fun(wildcards):
#    effective_args = get_effective_args(og_params.loc[og_params['model']==wildcards.model])
#    return expand('results/og/{{model}}/HE.AGG{arg}.png', arg=effective_args)
#
#def og_HEwald_AGGarg_fun(wildcards):
#    effective_args = get_effective_args(og_params.loc[og_params['model']==wildcards.model])
#    return expand('results/og/{{model}}/HE.wald.AGG{arg}.png', arg=effective_args)
#
#def og_REMLestimates_AGGarg_fun(wildcards):
#    effective_args = get_effective_args(og_params.loc[og_params['model']==wildcards.model])
#    return expand('results/og/{{model}}/REML.AGG{arg}.png', arg=effective_args)
#
#def og_REMLwaldNlrt_AGGarg_fun(wildcards):
#    effective_args = get_effective_args(og_params.loc[og_params['model']==wildcards.model])
#    return expand('results/og/{{model}}/REML.waldNlrt.AGG{arg}.png', arg=effective_args)

def og_HEnREML_AGGarg_fun(wildcards):
    effective_args = get_effective_args(og_params.loc[og_params['model']==wildcards.model])
    return expand('staging/og/{{model}}/HEnREML.AGG{arg}.flag', arg=effective_args)

rule og_AGGarg:
    input:
        og_HEnREML_AGGarg_fun,
    output:
        flag = touch('staging/og/{model}/all.flag'),

rule og_all:
    input:
        flag = expand('staging/og/{model}/all.flag', model=['hom', 'free']),

#########################################################################################
##   CTG
#########################################################################################
rule ctg_HE:
    input:
        cty = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/cty.batch{{i}}.txt',
        Z = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/Z.batch{{i}}.txt',
        ctnu = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/ctnu.batch{{i}}.txt',
    output:
        out = f'staging/ctg/{{model}}/{og_paramspace.wildcard_pattern}/he.batch{{i}}',
    params:
        out = f'staging/ctg/{{model}}/{og_paramspace.wildcard_pattern}/rep/he.npy',
        batches = lambda wildcards: og_batches[int(wildcards.i)],
    resources:
        time = '10:00:00',
        mem_per_cpu = '5000',
    priority: 1
    script: 'bin/ctg_HE.py'

rule ctg_REML:
    input:
        cty = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/cty.batch{{i}}.txt',
        Z = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/Z.batch{{i}}.txt',
        P = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/P.batch{{i}}.txt',
        ctnu = f'staging/og/{{model}}/{og_paramspace.wildcard_pattern}/ctnu.batch{{i}}.txt',
    output:
        out = f'staging/ctg/{{model}}/{og_paramspace.wildcard_pattern}/reml.batch{{i}}',
    params:
        out = f'staging/ctg/{{model}}/{og_paramspace.wildcard_pattern}/rep/reml.npy',
        batches = lambda wildcards: og_batches[int(wildcards.i)],
    resources:
        time = '100:00:00',
        #mem_per_cpu = '5000',
    priority: 1
    script: 'bin/ctg_REML.py'

use rule og_mergeBatches as ctg_mergeBatches with:
    input:
        he = [f'staging/ctg/{{model}}/{og_paramspace.wildcard_pattern}/he.batch{i}' for i in range(og_batch_no)],
        reml = [f'staging/ctg/{{model}}/{og_paramspace.wildcard_pattern}/reml.batch{i}' for i in range(og_batch_no)],
    output:
        out = f'analysis/ctg/{{model}}/{og_paramspace.wildcard_pattern}/out.npy',

def ctg_agg_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, og_params.loc[og_params['model']==wildcards.model])
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
                og_params.loc[og_params['model']==wildcards.model]).iloc[:,:],
        og_plot_order = og_plot_order,
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
                og_params.loc[og_params['model']==wildcards.model]).iloc[:,:],
        og_plot_order = og_plot_order,
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
                og_params.loc[og_params['model']==wildcards.model]).iloc[:,:],
        og_plot_order = og_plot_order,
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
                og_params.loc[og_params['model']==wildcards.model]).iloc[:,:],
        og_plot_order = og_plot_order,
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
    effective_args = get_effective_args(og_params.loc[og_params['model']==wildcards.model])
    return expand('staging/ctg/{{model}}/HEnREML.AGG{arg}.flag', arg=effective_args)

rule ctg_AGGarg:
    input:
        ctg_HEnREML_AGGarg_fun,
    output:
        flag = touch('staging/ctg/{model}/all.flag'),

rule ctg_all:
    input:
        flag = expand('staging/ctg/{model}/all.flag', model=['hom','free']),
