from snakemake.utils import Paramspace 
import os, re, time, shutil 
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
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
sim_replicates = 1000
sim_batch_no = 150
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
        'vc':['0.3_0.1_0.1_0.05_0.15_0.3', '0.3_0.1_0.1_0.1_0.1_0.3', '0.2_0.1_0.1_0.2_0.2_0.2', 
            '0.1_0.1_0.1_0.3_0.3_0.1', '0.3_0.1_0.1_0.15_0.05_0.3',
            '0.3_0.1_0.1_0.2_0_0.3'],
        'V_diag':['1_1_1_1', '8_4_2_1', '27_9_3_1', '64_16_4_1'],
        },
    'freeW': {
        'ss':['5e0', '2e1', '5e1', '1e2', '2e2', '3e2','5e2', '1e3'],
        },
    'full':{
        'ss':['2e1', '5e1', '1e2', '2e2', '5e2', '1e3'], 'a':['0.5_2_2_2', '1_2_2_2', '2_2_2_2', '4_2_2_2'],
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
        time = '100:00:00',
        mem_mb = lambda wildcards: f'{max(1,int(float(wildcards.ss))//150)}G',
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
        time = '360:00:00',
        mem_mb = lambda wildcards: f'{max(1,int(float(wildcards.ss))//150)}G',
    priority: 1
    script: 'bin/sim/reml.py'

#rule sim_mergeBatches:
#    input:
#        he = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.batch{i}' for i in range(sim_batch_no)],
#        reml = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/reml.batch{i}' for i in range(sim_batch_no)],
#    output:
#        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/out.npy',
#    script: 'bin/sim/mergeBatches.py'

rule sim_mergeBatches_HE:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.batch{i}' for i in range(sim_batch_no)],
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/out.he.npy',
    script: 'bin/mergeBatches.py'

rule sim_mergeBatches_REML:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/reml.batch{i}' for i in range(sim_batch_no)],
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/out.reml.npy',
    script: 'bin/mergeBatches.py'

########################################################################
# Yazar 2022 Science
########################################################################
yazar_ind_col = 'individual'
yazar_ct_col = 'cell_label'

# read parameters
yazar_params = pd.read_table('yazar.params.txt', dtype="str", comment='#')
yazar_paramspace = Paramspace(yazar_params, filename_params="*")

rule yazar_extract_meta:
    input:
        h5ad = 'data/Yazar2022Science/OneK1K_cohort_gene_expression_matrix_14_celltypes.h5ad.gz',
    output:
        obs = 'data/Yazar2022Science/obs.txt', # cells
        var = 'data/Yazar2022Science/var.txt', # genes
    run:
        import scanpy as sc

        data = sc.read_h5ad(input.h5ad, backed='r')
        obs = data.obs.reset_index(drop=False, names='cell')
        obs.to_csv(output.obs, sep='\t', index=False)

        var = data.var.reset_index(drop=False, names='feature')
        var.to_csv(output.var, sep='\t', index=False)

rule yazar_exclude_repeatedpool:
    input:
        obs = 'data/Yazar2022Science/obs.txt',
    output:
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
        dup_inds = 'analysis/yazar/duplicated_inds.txt',
    params:
        ind_col = yazar_ind_col, 
    run:
        obs = pd.read_table(input.obs)
        # id repeated pool: the same individual seqed in more than one pool
        data = obs[['pool',params.ind_col]].drop_duplicates()
        inds, counts = np.unique(data[params.ind_col], return_counts=True)
        inds = inds[counts > 1]
        np.savetxt(output.dup_inds, inds, fmt='%s')

        # for each ind find the largest pool
        for ind in inds:
            pools, counts = np.unique(obs.loc[obs[params.ind_col]==ind,'pool'], return_counts=True)
            excluded_pools = pools[counts < np.amax(counts)]
            obs = obs.loc[~((obs[params.ind_col]==ind) & (obs['pool'].isin(excluded_pools)))]

        obs.to_csv(output.obs, sep='\t', index=False)

rule yazar_age:
    input:
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
    output:
        png = 'results/yazar/age.png',
    run:
        obs = pd.read_table(input.obs)
        obs = obs.drop_duplicates(subset='individual')
        ages, counts = np.unique(obs['age'], return_counts=True)
        fig, ax = plt.subplots(dpi=600)
        plt.bar(ages, counts)
        ax.set_xlabel('Age')
        ax.set_ylabel('Number of individuals')
        fig.savefig(output.png)

rule yazar_ctp_extractX:
    input:
        h5ad = 'data/Yazar2022Science/OneK1K_cohort_gene_expression_matrix_14_celltypes.h5ad.gz',
        var = 'data/Yazar2022Science/var.txt',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
    output:
        X = 'staging/data/yazar/X.npz',
        obs = 'staging/data/yazar/obs.gz',
        var = 'staging/data/yazar/var.gz',
    params:
        ind_col = yazar_ind_col, 
        ct_col = yazar_ct_col,
    resources:
        mem_mb = '40G',
    run:
        import scanpy as sc
        from scipy import sparse

        genes = pd.read_table(input.var)
        if 'feature_is_filtered' in genes.columns:
            genes = genes.loc[~genes['feature_is_filtered'], 'feature'].to_numpy()
        else:
            genes = genes['feature'].to_numpy()

        if 'subset_gene' in params.keys():
            # random select genes
            rng = np.random.default_rng(seed=params.seed)
            genes = rng.choice(genes, params.subset_gene, replace=False)

        obs = pd.read_table(input.obs)
        ind_pool = np.unique(obs[params.ind_col].astype('str')+'+'+obs['pool'].astype('str'))

        ann = sc.read_h5ad(input.h5ad, backed='r')
        data = ann[(~ann.obs[params.ind_col].isna()) 
                & (~ann.obs[params.ct_col].isna()) 
                & (ann.obs[params.ind_col].astype('str')+'+'+ann.obs['pool'].astype('str')).isin(ind_pool), genes]
        # natural logarithm of one plus the input array
        X = data.X.log1p()
        sparse.save_npz(output.X, X)

        data.obs.to_csv(output.obs, sep='\t')
        data.var.to_csv(output.var, sep='\t')

rule yazar_ctp:
    input:
        X = 'staging/data/yazar/X.npz',
        obs = 'staging/data/yazar/obs.gz',
        var = 'staging/data/yazar/var.gz',
    output:
        ctp = 'data/Yazar2022Science/ctp.gz',
        ctnu = 'data/Yazar2022Science/ctnu.gz',
        P = 'data/Yazar2022Science/P.gz',
        n = 'data/Yazar2022Science/n.gz',
    params:
        ind_col = yazar_ind_col, 
        ct_col = yazar_ct_col,
    resources:
        mem_mb = '40G',
    run:
        from scipy import sparse
        from ctmm import preprocess

        X = sparse.load_npz(input.X)
        obs = pd.read_table(input.obs, index_col=0)
        var = pd.read_table(input.var, index_col=0)
        ctp, ctnu, P, n = preprocess.pseudobulk(X=X, obs=obs, var=var, ind_cut=100, ct_cut=10,
                ind_col=params.ind_col, ct_col=params.ct_col)

        # save
        ctp.to_csv(output.ctp, sep='\t')
        ctnu.to_csv(output.ctnu, sep='\t')
        P.to_csv(output.P, sep='\t')
        n.to_csv(output.n, sep='\t')

use rule yazar_ctp_extractX as yazar_var_ctnu_extract_genes with:
# don't know why it takes a lot of memory to extract the X matrix.
# so extract X before compting var of ctnu
    output:
        X = 'staging/data/yazar/var_ctnu/genes.npz',
        obs = 'staging/data/yazar/var_ctnu/obs.gz',
        var = 'staging/data/yazar/var_ctnu/var.gz',
    params:
        ind_col = yazar_ind_col, 
        ct_col = yazar_ct_col,
        seed = 123567,
        subset_gene = 1000,

rule yazar_var_ctnu_split:
    input:
        obs = 'staging/data/yazar/var_ctnu/obs.gz',
    output:
        batches = expand('staging/data/yazar/var_ctnu/ind_ct.batch{i}', i=range(10)),
    params:
        ind_col = yazar_ind_col,
        ct_col = yazar_ct_col,
    run:
        obs = pd.read_table(input.obs)
        obs = obs.rename(columns={params.ind_col:'ind', params.ct_col:'ct'})

        # pairs of ind and ct
        ind_ct = obs.loc[(~obs['ind'].isna()) & (~obs['ct'].isna()), ['ind', 'ct']].drop_duplicates()

        # Split the DataFrame into smaller DataFrames
        ind_ct_batches = np.array_split(ind_ct, len(output.batches))
        
        for batch, batch_f in zip(ind_ct_batches, output.batches):
            batch.to_csv(batch_f, sep='\t', index=False)

rule yazar_var_ctnu:
    input:
        X = 'staging/data/yazar/var_ctnu/genes.npz',
        obs = 'staging/data/yazar/var_ctnu/obs.gz',
        var = 'staging/data/yazar/var_ctnu/var.gz',
        batch = 'staging/data/yazar/var_ctnu/ind_ct.batch{i}',
    output:
        var_ctnu = 'staging/data/yazar/var_ctnu/batch{i}.gz',
    params:
        ind_col = yazar_ind_col,
        ct_col = yazar_ct_col,
        seed = 42,
    resources:
        mem_mb = '2G',
    run:
        from scipy import stats, sparse
        from sklearn.utils import resample

        def cal_ctnu(data):
            ctp = np.squeeze( np.asarray(data.mean(axis=0)) )
            ctp2 = np.squeeze( np.asarray(data.power(2).mean(axis=0)) )
            ctnu = (ctp2 - ctp**2) / (data.shape[0]**2)
            return( ctnu )

        def bootstrap(data, rng, n_resamples=10000):
            ctnus = []
            for i in range(n_resamples):
                sample = resample(data, random_state=rng)
                ctnus.append( cal_ctnu(sample) )
            return( ctnus )

        X = sparse.load_npz(input.X)
        obs = pd.read_table(input.obs)
        obs = obs.rename(columns={params.ind_col:'ind', params.ct_col:'ct'})
        genes = pd.read_table(input.var)['feature'].to_numpy()

        # pairs of ind and ct
        ind_ct = pd.read_table(input.batch)

        # bootstrap
        rng = np.random.RandomState( params.seed )
        boots = {'ind':[], 'ct':[], 'var_ctnu':[]}
        for index, row in ind_ct.iterrows():
            print( index, flush=True )
            ind, ct = row['ind'], row['ct']
            data = X[(obs['ind']==ind) & (obs['ct']==ct), :]
            if data.shape[0] < 10:
                continue
            else:
                ctnus = bootstrap(data, rng)
                var_ctnu = np.std(ctnus, axis=0, ddof=1)**2
            boots['ind'].append( ind )
            boots['ct'].append( ct )
            boots['var_ctnu'].append( var_ctnu )

        var_ctnu = pd.DataFrame(data=boots['var_ctnu'], columns=genes) # double check the order of genes is correct
        var_ctnu.insert(loc=0, column='ct', value=boots['ct'])
        var_ctnu.insert(loc=0, column='ind', value=boots['ind'])
        var_ctnu.to_csv(output.var_ctnu, sep='\t', index=False)

rule yazar_var_ctnu_merge:
    input:
        var_ctnu = expand('staging/data/yazar/var_ctnu/batch{i}.gz', i=range(10)),
    output:
        var_ctnu = 'analysis/yazar/var_ctnu.gz',
    shell:
        "zcat {input.var_ctnu}|awk '!(FNR==1 && NR!=1) {{print}}' |gzip -c > {output.var_ctnu}"

rule yazar_rm_rareINDnCT:
    input:
        ctp = 'data/Yazar2022Science/ctp.gz',
        ctnu = 'data/Yazar2022Science/ctnu.gz',
        n = 'data/Yazar2022Science/n.gz',
    output:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctp.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctnu.gz',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
        n = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/n.gz',
    resources:
        mem_mb = '10G',
    script: 'bin/yazar/rm_rareINDnCT.py'

rule yazar_mvn_ctp:
    input:
        data = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctp.gz',
    output:
        data = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',
    resources:
        mem_mb = '10G',
    run:
        from ctmm import preprocess
        data = pd.read_table(input.data, index_col=(0,1)).astype('float32')
        preprocess.mvn(data).to_csv(output.data, sep='\t')

use rule yazar_mvn_ctp as yazar_mvn_ctnu with:
    input:
        data = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctnu.gz',
    output:
        data = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',

rule yazar_std_op:
    input:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
    output:
        op = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/op.mvn.gz',
        nu = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/nu.mvn.gz',
        ctp = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',
        ctnu = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',
    resources:
        mem_mb = '10G',
    run:
        from ctmm import preprocess
        ctp = pd.read_table(input.ctp, index_col=(0,1)).astype('float32')
        ctnu = pd.read_table(input.ctnu, index_col=(0,1)).astype('float32')
        P = pd.read_table(input.P, index_col=0)

        op, nu, ctp, ctnu = preprocess.std(ctp, ctnu, P)

        op.to_csv(output.op, sep='\t')
        nu.to_csv(output.nu, sep='\t')
        ctp.to_csv(output.ctp, sep='\t')
        ctnu.to_csv(output.ctnu, sep='\t')

rule yazar_op_pca:
    input:
        op = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/op.mvn.gz',
    output:
        evec = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/evec.txt',
        eval = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/eval.txt',
        pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.txt',
        png = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.png',
    resources:
        mem_mb = '4G',
    script: 'bin/yazar/pca.py'

rule yazar_exclude_duplicatedSNPs:
    input:
        vcf = 'data/Yazar2022Science/filter_vcf_r08/chr{chr}.dose.filtered.R2_0.8.vcf.gz',
    output:
        bed = 'analysis/yazar/data/geno/chr{chr}.bed',
        dup = 'analysis/yazar/data/geno/chr{chr}.dup',
    shell:
        '''
        module load gcc/11.3.0 atlas/3.10.3 lapack/3.11.0 plink/1.9
        prefix="$(dirname {output.bed})/$(basename {output.bed} .bed)"
        zcat {input.vcf}|grep -v '#'|cut -f 3|sort|uniq -d > {output.dup}
        plink --vcf {input.vcf} --double-id --keep-allele-order \
                --snps-only \
                --exclude {output.dup} \
                --make-bed --out $prefix
        '''

rule yazar_geno_pca:
    input:
        bed = expand('analysis/yazar/data/geno/chr{chr}.bed',
                chr=range(1,23)),
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
    output:
        eigenvec = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        eigenval = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenval',
    params:
        prefix = lambda wildcards, output: os.path.splitext(output.eigenvec)[0],
        tmp_dir = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/geno_tmp',
    shell:
        '''
        if [ -d {params.tmp_dir} ]; then 
            rm -r {params.tmp_dir} 
        fi
        mkdir -p {params.tmp_dir}
        module load gcc/11.3.0 atlas/3.10.3 lapack/3.11.0 plink/1.9
        ind_f="{params.tmp_dir}/inds.txt" 
        zcat {input.P}|tail -n +2|awk '{{print $1,$1}}' > $ind_f
        merge_f="{params.tmp_dir}/merge.list"
        touch $merge_f
        for bed in {input.bed} 
        do 
            prefix="$(dirname $bed)/$(basename $bed .bed)"
            o_prefix="{params.tmp_dir}/$(basename $bed .bed)"
            echo $o_prefix >> $merge_f
            plink --bfile $prefix \
                --maf 0.05 \
                --keep $ind_f \
                --make-bed --out $o_prefix
        done
        # merge
        merged={params.tmp_dir}/merged
        plink --merge-list $merge_f --out $merged
        # ld prune 
        plink --bfile $merged --indep-pairwise 50 5 0.2 --out $merged
        plink --bfile $merged --extract $merged.prune.in --make-bed --out $merged.ld
        # pca
        plink --bfile $merged.ld --pca 50 header tabs --out {params.prefix}
        rm -r {params.tmp_dir}
        '''

rule yazar_geno_pca_plot:
    input:
        eigenval = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenval',
    output:
        png = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenval.png',
    run:
        import matplotlib.pyplot as plt
        vals = np.loadtxt(input.eigenval)
        plt.rcParams.update({'font.size' :12})
        fig, ax = plt.subplots()
        ax.scatter(np.arange(1,11), vals[:10])
        ax.set_xlabel('PC', fontsize=14)
        ax.set_ylabel('Eigenvalue', fontsize=14)
        fig.savefig(output.png)

rule yazar_gene_location:
    input:
        genes = 'data/Yazar2022Science/genes.txt',
        gff = 'data/gencode.v43lift37.annotation.gff3.gz',
    output:
        genes = 'data/Yazar2022Science/gene_loation.txt',
    script: 'bin/yazar/gene_location.py'

rule yazar_he_kinship:
    input:
        genes = 'data/Yazar2022Science/gene_loation.txt',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
        bed = 'analysis/yazar/data/geno/chr{chr}.bed',
    output:
        kinship = temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.chr{{chr}}.txt'),
    params:
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship/gene.rel.bin', 
        r = int(float(5e5)),
    resources:
        mem_mb = '2G',
    shell: 
        '''
        module load gcc/11.3.0 atlas/3.10.3 lapack/3.11.0 plink/1.9
        mkdir -p $(dirname {params.kinship})
        python3 bin/yazar/kinship.py {input.genes} {input.P} {params.r} {input.bed} {wildcards.chr} \
                        {params.kinship} {output.kinship} 
        '''

rule yazar_he_kinship_merge:
    input:
        kinship = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.chr{chr}.txt'
                for chr in range(1,23)],
    output:
        #kinship = temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt'),
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
        save = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
    shell:
        '''
        awk '!(FNR==1 && NR!=1) {{print}}' {input.kinship} > {output.kinship}
        cp {output.kinship} {output.save}
        '''

yazar_he_batches = 2000
rule yazar_HE_split:
    input:
        ctp = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',
        ctnu = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',
    output:
        ctp = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{i}.gz'
                for i in range(yazar_he_batches)],
        ctnu = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{i}.gz'
                for i in range(yazar_he_batches)],
    resources:
        mem_mb = '10G',
    run:
        ctp = pd.read_table(input.ctp, index_col=(0,1)).astype('float32')
        ctnu = pd.read_table(input.ctnu, index_col=(0,1)).astype('float32')
        genes = ctp.columns
        batches = np.array_split(genes, len(output.ctp))
        for batch, ctp_f, ctnu_f in zip(batches, output.ctp, output.ctnu):
            ctp[batch].to_csv(ctp_f, sep='\t')
            ctnu[batch].to_csv(ctnu_f, sep='\t')

rule yazar_HE_full:
    input:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'staging/data/yazar/obs.gz',
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he.full.batch{{i}}',
    params:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/rep/he.full.npy',
        snps = 5, # threshold of snp number per gene
    resources:
        time = '10:00:00',
        mem_mb = '30G',
    script: 'bin/yazar/he.full.py'

rule yazar_HE_full_merge:
    input:
        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he.full.batch{i}'
            for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/he.full.npy',
    script: 'bin/mergeBatches.py'

rule yazar_HE_Full_plot:
    input:
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/he.full.npy',
    output:
        cov = f'results/yazar/{yazar_paramspace.wildcard_pattern}/he.full.cov.png',
        h2 = f'results/yazar/{yazar_paramspace.wildcard_pattern}/he.full.h2.png',
    script: 'bin/yazar/full.plot.py'

rule yazar_HE_free:
    input:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'staging/data/yazar/obs.gz',
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he.free.batch{{i}}',
    params:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/rep/he.free.npy',
        snps = 5, # threshold of snp number per gene
    resources:
        time = '10:00:00',
        mem_mb = '20G',
    script: 'bin/yazar/he.free.py'

rule yazar_HE_free_merge:
    input:
        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he.free.batch{i}'
            for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/he.free.npy',
    script: 'bin/mergeBatches.py'

rule yazar_HE_Free_plot:
    input:
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/he.free.npy',
    output:
        h2 = f'results/yazar/{yazar_paramspace.wildcard_pattern}/he.free.h2.png',
    script: 'bin/yazar/free.plot.py'

rule yazar_HE_clean:
    input:
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/he.full.npy',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
    output:
        touch(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he.done')
    run:
        kinship = pd.read_table(input.kinship)
        for f in kinship['K'].tolist():
            os.remove(f)

rule yazar_all:
    input:
        full = expand('results/yazar/{params}/he.full.h2.png', params=yazar_paramspace.instance_patterns),
        free = expand('results/yazar/{params}/he.free.h2.png', params=yazar_paramspace.instance_patterns),

#######################################################################################
# Perez 2022 Science
#######################################################################################
perez_ind_col = 'ind_cov'
perez_ct_col = 'author_cell_type'
use rule yazar_extract_meta as perez_extract_meta with:
    input:
        h5ad = 'data/Perez2022Science/local.h5ad',
    output:
        obs = 'data/Perez2022Science/obs.txt',
        var = 'data/Perez2022Science/var.txt',

#use rule yazar_ctp as perez_ctp with:
#    input:
#        h5ad = 'data/Perez2022Science/local.h5ad',
#        genes = 'data/Perez2022Science/genes.txt',
#    output:
#        ctp = 'data/Perez2022Science/ctp.gz',
#        ctnu = 'data/Perez2022Science/ctnu.gz',
#        P = 'data/Perez2022Science/P.gz',
#    params:
#        ind_col = perez_ind_col,
#        ct_col = perez_ct_col,

use rule yazar_var_ctnu_extract_genes as perez_var_ctnu_extract_genes with:
    # don't know why it takes a lot of memory to extract the X matrix.
    # so extract X before compting var of ctnu
    input:
        h5ad = 'data/Perez2022Science/local.h5ad',
        var = 'data/Perez2022Science/var.txt',
    output:
        counts = 'staging/data/Perez2022Science/var_ctnu.genes.npz',
        genes = 'staging/data/Perez2022Science/var_ctnu.genes.txt',
    params:
        seed = 123567,
        gene_no = 10,

use rule yazar_var_ctnu as perez_var_ctnu with:
    input:
        obs = 'data/Perez2022Science/obs.txt',
        counts = 'staging/data/Perez2022Science/var_ctnu.genes.npz',
        genes = 'staging/data/Perez2022Science/var_ctnu.genes.txt',
    output:
        var_ctnu = 'analysis/Perez2022Science/var_nu.gz',
    params:
        ind_col = perez_ind_col,
        ct_col = perez_ct_col,

include: 'dev.snake'
