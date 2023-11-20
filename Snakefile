from snakemake.utils import Paramspace 
import os, re, time, shutil, gzip 
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
wildcard_constraints: l='[^/]+'
wildcard_constraints: r='[^/]+'
wildcard_constraints: disease='[^/]+'
wildcard_constraints: model='\w+'
wildcard_constraints: var_nu='[^/]+'


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
        'ss':['50', '100', '200', '300', '500', '1000'], 
        'a':['0.5_2_2_2', '1_2_2_2', '2_2_2_2', '4_2_2_2']
        },
    'iid':{
        'ss':['50', '100', '500', '1000'], 'a':['0.5_2_2_2', '1_2_2_2', '2_2_2_2', '4_2_2_2'],
        'vc':['0.3_0.1_0.1_0.1_0.1_0.3', '0.2_0.1_0.1_0.2_0.2_0.2', '0.1_0.1_0.1_0.3_0.3_0.1'],
        }, 
    'free': {
        'ss':['50', '100', '200', '300','500', '1000'], 'a':['0.5_2_2_2', '1_2_2_2', '2_2_2_2', '4_2_2_2'], 
        'vc':['0.3_0.1_0.1_0.05_0.15_0.3', '0.3_0.1_0.1_0.1_0.1_0.3', '0.2_0.1_0.1_0.2_0.2_0.2', 
            '0.1_0.1_0.1_0.3_0.3_0.1', '0.3_0.1_0.1_0.15_0.05_0.3',
            '0.3_0.1_0.1_0.2_0_0.3'],
        'V_diag':['1_1_1_1', '8_4_2_1', '27_9_3_1', '64_16_4_1'],
        },
    'freeW': {
        'ss':['50', '100', '200', '300','500', '1000'],
        },
    'full':{
        'ss':['50', '100', '200', '500', '1000'], 'a':['0.5_2_2_2', '1_2_2_2', '2_2_2_2', '4_2_2_2'],
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
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/sim.batch{{i}}.npy',
    params:
        batches = sim_batches,
        beta = (0.5, 0.5), # beta distribution for allele frequency
        maf = 0.05,
        L = 10, # number of causal SNPs
        seed = 273672,
    script: 'bin/sim/generatedata.py'


rule sim_HE:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/sim.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.batch{{i}}.npy',
    resources:
        time = '36:00:00',
        mem_mb = '10G',
    priority: 1
    script: 'bin/sim/he.py'


rule sim_REML:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/sim.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/reml.batch{{i}}.npy',
    resources:
        time = '360:00:00',
        mem_mb = '10G',
    priority: 1
    script: 'bin/sim/reml.py'


rule sim_mergeBatches_HE:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.batch{i}.npy' for i in range(sim_batch_no)],
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/out.he.npy',
    script: 'bin/mergeBatches.py'


rule sim_mergeBatches_REML:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/reml.batch{i}.npy' for i in range(sim_batch_no)],
    output:
        out = f'analysis/sim/{{model}}/{sim_paramspace.wildcard_pattern}/out.reml.npy',
    script: 'bin/mergeBatches.py'


########################################################################
# Yazar 2022 Science
########################################################################
yazar_ind_col = 'individual'
yazar_ct_col = 'cell_label'

yazar_ct_order = np.array(['hom', 'CD4 NC', 'CD8 ET', 'NK', 'CD8 NC', 'B IN', 'CD4 ET', 'B Mem', 'Mono C', 'CD8 S100B', 'Mono NC', 'NK R', 'DC', 'CD4 SOX4', 'Plasma'])
yazar_colors = dict(zip(yazar_ct_order[1:], sns.color_palette()))
yazar_colors['hom'] = '0.7'

# read parameters
yazar_params = pd.read_table('yazar.params.txt', dtype="str", comment='#')
yazar_paramspace = Paramspace(yazar_params, filename_params="*")

include: 'yazar.smk'

rule yazar_rm_rareINDnCT:
    # also select gene expressed in all cts
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


rule yazar_filter_genes:
    input:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctp.gz',
    output:
        genes = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/genes.{{cut}}.txt',
    run:
        def count_expressed(x, threshold=float(wildcards.cut)):
            x = x.unstack()
            prop = (x > 0).sum() / x.count()
            if np.all(prop > threshold):
                return True
            else:
                return False

        ctp = pd.read_table(input.ctp, index_col=(0,1))
        selected = ctp.apply(count_expressed)
        np.savetxt(output.genes, ctp.columns[selected].to_numpy(), fmt='%s')


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
        # load plink 1.9
        if [[ $(hostname) == *midway* ]]; then
            module load plink
        else
            module load gcc/11.3.0 atlas/3.10.3 lapack/3.11.0 plink/1.9
        fi

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
    resources:
        mem_mb = '40G',
    shell:
        '''
        # load plink 1.9
        if [[ $(hostname) == *midway* ]]; then
            module load plink
        else
            module load gcc/11.3.0 atlas/3.10.3 lapack/3.11.0 plink/1.9
        fi

        if [ -d {params.tmp_dir} ]; then 
            rm -r {params.tmp_dir} 
        fi
        mkdir -p {params.tmp_dir}
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
        kinship = temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.chr{{chr}}.npy'),
    params:
        r = int(float(5e5)),
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

        python3 bin/yazar/kinship.npy.py {input.genes} {input.P} {params.r} {input.bed} {wildcards.chr} \
                        {output.kinship} 
        '''


# rule yazar_he_kinship_merge:
#     input:
#         kinship = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.chr{chr}.txt'
#                 for chr in range(1,23)],
#     output:
#         #kinship = temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt'),
#         kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
#         save = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
#     shell:
#         '''
#         awk '!(FNR==1 && NR!=1) {{print}}' {input.kinship} > {output.kinship}
#         cp {output.kinship} {output.save}
#         '''

rule yazar_he_kinship_merge:
    input:
        kinship = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.chr{chr}.npy'
                for chr in range(1,23)],
    output:
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.npz',
    resources:
        mem_mb = '150G',
    run:
        all = {}
        for chr, kinship in zip(range(1, 23), input.kinship):
            all[str(chr)] = np.load(kinship, allow_pickle=True).item()
        
        np.savez_compressed(output.kinship, **all)


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


rule yazar_he_kinship_split:
    input:
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.npz',
        ctp = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{i}.gz'
                for i in range(yazar_he_batches)],
    output:
        kinship = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.batch{i}.npy'
                for i in range(yazar_he_batches)],   # TODO: temp?
    resources:
        mem_mb = '150G',
    run:
        data = {}
        for key, value in np.load(input.kinship, allow_pickle=True).items():
            data[key] = value.item()

        ids = data[list(data.keys())[0]]['ids']
        for value in data.values():
            if np.any(ids != value['ids']):
                sys.exit('Wrong order of individuals!')
        
        for ctp_f, kinship_f in zip(input.ctp, output.kinship):
            genes = gzip.open(ctp_f, 'rt').readline().strip().split()[2:]
            kinship = {'ids': ids, 'gene': [], 'K': [], 'nsnp': []}

            for value in data.values():
                idx = np.isin(value['gene'], genes)
                if idx.sum() == 0:
                    continue
                else:
                    kinship['gene'] += np.array(value['gene'])[idx].tolist()
                    kinship['K'] += np.array(value['K'])[idx].tolist()
                    kinship['nsnp'] += np.array(value['nsnp'])[idx].tolist()

            # convert to array
            for key, value in kinship.items():
                kinship[key] = np.array(value)

            np.save(kinship_f, kinship)


rule yazar_HE_full:
    input:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.batch{{i}}.npy',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he.full.batch{{i}}.npy',
    params:
        snps = 5, # threshold of snp number per gene
    resources:
        partition = 'tier3q',
        mem_mb = '120G',
        burden = 99,
    script: 'bin/yazar/he.full.py'


rule yazar_HE_full_merge:
    input:
        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he.full.batch{i}.npy'
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
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.batch{{i}}.npy',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'staging/data/yazar/obs.gz',  # TODO: try not to use staging
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he.free.batch{{i}}.npy',
    params:
        snps = 5, # threshold of snp number per gene
    resources:
        mem_mb = '40G',
        time = '36:00:00',
    script: 'bin/yazar/he.free.py' 


rule yazar_HE_free_merge:
    input:
        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he.free.batch{i}.npy'
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
    params:
        order = yazar_ct_order,
        colors = yazar_colors,
    script: 'bin/yazar/free.plot.py'


rule yazar_HE_Free_VW_plot:
    input:
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/he.free.npy',
    output:
        png = f'results/yazar/{yazar_paramspace.wildcard_pattern}/he.free.VW.png',
    params:
        order = yazar_ct_order,
        colors = yazar_colors,
    script: 'bin/yazar/free.VW.plot.py'


rule yazar_REML_free:
    input:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.batch{{i}}.npy',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'staging/data/yazar/obs.gz',
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml.free.batch{{i}}.npy',
    params:
        snps = 5, # threshold of snp number per gene
    resources:
        time = '600:00:00',
        mem_mb = '10G',
    script: 'bin/yazar/reml.free.py'


use rule yazar_HE_free_merge as yazar_REML_free_merge with:
    input:
        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml.free.batch{i}.npy'
            for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/reml.free.npy',


use rule yazar_HE_Free_plot as yazar_REML_Free_plot with:
    input:
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/reml.free.npy',
    output:
        h2 = f'results/yazar/{yazar_paramspace.wildcard_pattern}/reml.free.h2.png',


# rule yazar_clean:
#     input:
#         # out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/he.full.npy',
#         kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
#     output:
#         touch(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he.done')
#     run:
#         kinship = pd.read_table(input.kinship)
#         for f in kinship['K'].tolist():
#             if os.path.exists(f):
#                 os.remove(f)


rule yazar_all:
    input:
        full = expand('results/yazar/{params}/he.full.h2.png', params=yazar_paramspace.instance_patterns),
        free = expand('results/yazar/{params}/he.free.h2.png', params=yazar_paramspace.instance_patterns),
        #reml_free = expand('results/yazar/{params}/reml.free.h2.png', params=yazar_paramspace.instance_patterns),




##########################################################################
#  1.2: remove missing individuals
##########################################################################
rule yazar_nomissing_rmIND:
    input:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctp.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctnu.gz',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
        n = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/n.gz',
    output:
        ctp = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/ctp.gz',
        ctnu = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/ctnu.gz',
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        n = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/n.gz',
    resources:
        mem_mb = '10G',
    run:
        ctp = pd.read_table(input.ctp, index_col=(0, 1)).sort_index()
        ctnu = pd.read_table(input.ctnu, index_col=(0, 1)).sort_index()
        P = pd.read_table(input.P, index_col=0)
        n = pd.read_table(input.n, index_col=0)

        # select ids
        ids = n.index.to_numpy()[~(n <= int(wildcards.ct_min_cellnum)).any(axis='columns')]

        # 
        ctp.loc[ctp.index.get_level_values('ind').isin(ids)].to_csv(output.ctp, sep='\t')
        ctnu.loc[ctnu.index.get_level_values('ind').isin(ids)].to_csv(output.ctnu, sep='\t')
        P.loc[P.index.isin(ids)].to_csv(output.P, sep='\t')
        n.loc[n.index.isin(ids)].to_csv(output.n, sep='\t')


use rule yazar_std_op as yazar_nomissing_std_op with:
    input:
        ctp = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/ctp.gz',
        ctnu = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/ctnu.gz',
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
    output:
        op = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/op.std.gz',
        nu = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/nu.std.gz',
        ctp = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/ctp.std.gz',
        ctnu = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/ctnu.std.gz',


use rule yazar_op_pca as yazar_nomissing_op_pca with:
    input:
        op = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/op.std.gz',
    output:
        evec = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/evec.txt',
        eval = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/eval.txt',
        pca = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/pca.txt',
        png = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/pca.png',


use rule yazar_he_kinship as yazar_nomissing_he_kinship with:
    input:
        genes = 'data/Yazar2022Science/gene_loation.txt',
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        bed = 'analysis/yazar/data/geno/chr{chr}.bed',
    output:
        kinship = temp(f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/kinship.chr{{chr}}.npy'),


use rule yazar_he_kinship_merge as yazar_nomissing_he_kinship_merge with:
    input:
        kinship = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/kinship.chr{chr}.npy'
                for chr in range(1,23)],
    output:
        kinship = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/kinship.npz',


use rule yazar_HE_split as yazar_nomissing_he_split with:
    input:
        ctp = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/ctp.std.gz',
        ctnu = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/ctnu.std.gz',
    output:
        # ctp = [temp(f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{i}.gz')
        ctp = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{i}.gz'
                for i in range(yazar_he_batches)],
        # ctnu = [temp(f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{i}.gz')
        ctnu = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{i}.gz'
                for i in range(yazar_he_batches)],


use rule yazar_he_kinship_split as yazar_nomissing_he_kinship_split with:
    input:
        kinship = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/kinship.npz',
        ctp = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{i}.gz'
                for i in range(yazar_he_batches)],
    output:
        # kinship = [temp(f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/kinship.batch{i}.npy')
        kinship = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/kinship.batch{i}.npy'
                for i in range(yazar_he_batches)],   # temp?


use rule yazar_HE_full as yazar_nomissing_HE_full with:
    input:
        ctp = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/kinship.batch{{i}}.npy',
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
    output:
        out = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.full.batch{{i}}.npy',


use rule yazar_HE_full_merge as yazar_nomissing_HE_full_merge with:
    input:
        out = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.full.batch{i}.npy'
            for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.full.npy',


use rule yazar_HE_Full_plot as yazar_nomissing_HE_Full_plot with:
    input:
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.full.npy',
    output:
        cov = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.full.cov.png',
        h2 = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.full.h2.png',


use rule yazar_HE_free as yazar_nomissing_HE_free with:
    input:
        ctp = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/kinship.batch{{i}}.npy',
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',  # TODO: PCA including missing individuals
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
    output:
        out = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.free.batch{{i}}.npy',
    params:
        jk = False,
        snps = 5, # threshold of snp number per gene


use rule yazar_HE_free_merge as yazar_nomissing_HE_free_merge with:
    input:
        out = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.free.batch{i}.npy'
            for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.free.npy',


use rule yazar_HE_Free_plot as yazar_nomissing_HE_Free_plot with:
    input:
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.free.npy',
    output:
        h2 = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.free.h2.png',


use rule yazar_HE_Free_VW_plot as yazar_nomissing_HE_Free_VW_plot with:
    input:
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.free.npy',
    output:
        png = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.free.VW.png',
        p_png = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.free.VW.p.png',


rule yazar_nomissing_all:
    input:
        free = expand('results/yazar/nomissing/{params}/he.free.h2.png', params=yazar_paramspace.instance_patterns),
        VW = expand('results/yazar/nomissing/{params}/he.free.VW.png', params=yazar_paramspace.instance_patterns),
        full = expand('results/yazar/nomissing/{params}/he.full.cov.png', params=yazar_paramspace.instance_patterns),


rule yazar_nomissing_topgenes:
    input:
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.free.npy',
    output:
        tophomg = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.free.tophomg.genes.txt',
        topV = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.free.topV.genes.txt',
    run:
        out = np.load(input.out, allow_pickle=True).item()
        genes = out['gene']
        homgs = out['free']['hom_g2']
        Vs = np.diagonal(out['free']['V'], axis1=1, axis2=2)

        # median V
        V_bar = np.median(Vs, axis=1)

        # top 10
        np.savetxt(output.tophomg, genes[np.argsort(homgs)[-10:][::-1]], fmt='%s')
        np.savetxt(output.topV, genes[np.argsort(V_bar)[-10:][::-1]], fmt='%s')
        
        
rule yazar_nomissing_plottopVgene:
    input:
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.free.npy',
        genes = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.free.topV.genes.txt',
        var = 'data/Yazar2022Science/var.txt', # genes
    output:
        png = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.free.topV.gene.png',
    params:
        order = yazar_ct_order,
        colors = yazar_colors,
    script: 'scripts/yazar/plot.gene.py'


rule yazar_nomissing_plottopVgene_tmp:
    input:
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.free.jk.batch51.npy',
        var = 'data/Yazar2022Science/var.txt', # genes
    output:
        png = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.free.topV.gene.tmp.png',
    params:
        order = yazar_ct_order,
        colors = yazar_colors,
    script: 'scripts/yazar/plot.gene.tmp.py'


#######################################################################################
# Perez 2022 Science
#######################################################################################
perez_ind_col = 'ind_cov'
perez_ct_col = 'author_cell_type'
#use rule yazar_extract_meta as perez_extract_meta with:
#    input:
#        h5ad = 'data/Perez2022Science/local.h5ad',
#    output:
#        obs = 'data/Perez2022Science/obs.txt',
#        var = 'data/Perez2022Science/var.txt',
#
##use rule yazar_ctp as perez_ctp with:
##    input:
##        h5ad = 'data/Perez2022Science/local.h5ad',
##        genes = 'data/Perez2022Science/genes.txt',
##    output:
##        ctp = 'data/Perez2022Science/ctp.gz',
##        ctnu = 'data/Perez2022Science/ctnu.gz',
##        P = 'data/Perez2022Science/P.gz',
##    params:
##        ind_col = perez_ind_col,
##        ct_col = perez_ct_col,
#
#use rule yazar_var_ctnu_extract_genes as perez_var_ctnu_extract_genes with:
#    # don't know why it takes a lot of memory to extract the X matrix.
#    # so extract X before compting var of ctnu
#    input:
#        h5ad = 'data/Perez2022Science/local.h5ad',
#        var = 'data/Perez2022Science/var.txt',
#    output:
#        counts = 'staging/data/Perez2022Science/var_ctnu.genes.npz',
#        genes = 'staging/data/Perez2022Science/var_ctnu.genes.txt',
#    params:
#        seed = 123567,
#        gene_no = 10,
#
#use rule yazar_var_ctnu as perez_var_ctnu with:
#    input:
#        obs = 'data/Perez2022Science/obs.txt',
#        counts = 'staging/data/Perez2022Science/var_ctnu.genes.npz',
#        genes = 'staging/data/Perez2022Science/var_ctnu.genes.txt',
#    output:
#        var_ctnu = 'analysis/Perez2022Science/var_nu.gz',
#    params:
#        ind_col = perez_ind_col,
#        ct_col = perez_ct_col,
#
include: 'dev.smk'
