#########################################################################################
##   Simulation
#########################################################################################


reml_ss_cut = 500

def sim_agg_he_truebeta_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/celltypebeta.txt', params=subspace.instance_patterns)

def sim_agg_he_trueV_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/V.txt', params=subspace.instance_patterns)

def sim_agg_he_trueW_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/W.txt', params=subspace.instance_patterns)

def sim_agg_reml_truebeta_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[(sim_params['model']==wildcards.model) 
        & (sim_params['ss'].astype('float') <= reml_ss_cut)])
    return expand('analysis/sim/{{model}}/{params}/celltypebeta.txt', params=subspace.instance_patterns)

def sim_agg_reml_trueV_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[(sim_params['model']==wildcards.model) 
        & (sim_params['ss'].astype('float') <= reml_ss_cut)])
    return expand('analysis/sim/{{model}}/{params}/V.txt', params=subspace.instance_patterns)

def sim_agg_reml_trueW_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[(sim_params['model']==wildcards.model) 
        & (sim_params['ss'].astype('float') <= reml_ss_cut)])
    return expand('analysis/sim/{{model}}/{params}/W.txt', params=subspace.instance_patterns)

def sim_agg_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/out.npy', params=subspace.instance_patterns)

def sim_agg_he_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('analysis/sim/{{model}}/{params}/out.he.npy', params=subspace.instance_patterns)

def sim_agg_reml_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[(sim_params['model']==wildcards.model) 
        & (sim_params['ss'].astype('float') <= reml_ss_cut)])
    return expand('analysis/sim/{{model}}/{params}/out.reml.npy', params=subspace.instance_patterns)

rule sim_HEestimates_subspace_plot:
    input:
        out = sim_agg_he_out_subspace,
        V = sim_agg_he_trueV_subspace,
        W = sim_agg_he_trueW_subspace,
    output:
        png = 'results/sim/{model}/HE.AGG{arg}.png',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[sim_params['model']==wildcards.model]).iloc[:,:],
        sim_plot_order = sim_plot_order,
        method = None,
        mycolors = mycolors,
        pointcolor = pointcolor,
        colorpalette = colorpalette,
    script: 'scripts/sim/estimates_subspace_plot.py'


rule sim_HEwald_subspace_plot:
    input:
        out = sim_agg_he_out_subspace,
    output:
        png = 'results/sim/{model}/HE.wald.AGG{arg}.png',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[sim_params['model']==wildcards.model]).iloc[:,:],
        sim_plot_order = sim_plot_order,
        mycolors = mycolors,
    script: 'scripts/sim/HEwald_subspace_plot.py'


rule sim_REMLestimates_subspace_plot:
    input:
        out = sim_agg_reml_out_subspace,
        beta = sim_agg_reml_truebeta_subspace, 
        V = sim_agg_reml_trueV_subspace,
        W = sim_agg_reml_trueW_subspace,
    output:
        png = 'results/sim/{model}/REML.AGG{arg}.png',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[(sim_params['model']==wildcards.model) 
                    & (sim_params['ss'].astype('float') <= reml_ss_cut)]).iloc[:,:],
        sim_plot_order = sim_plot_order,
        method = None,
        colorpalette = colorpalette,
        pointcolor = pointcolor,
        mycolors = mycolors,
    script: 'scripts/sim/estimates_subspace_plot.py'


rule sim_REMLwaldNlrt_subspace_plot:
    input:
        out = sim_agg_reml_out_subspace,
    output:
        png = 'results/sim/{model}/REML.waldNlrt.AGG{arg}.png',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[(sim_params['model']==wildcards.model) 
                    & (sim_params['ss'].astype('float') <= reml_ss_cut)]).iloc[:,:],
        sim_plot_order = sim_plot_order,
        mycolors = mycolors,
    script: 'scripts/sim/REMLwaldNlrt_subspace_plot.py'


rule sim_collect_HEnREML_subspace_plot:
    input:
        he = 'results/sim/{model}/HE.AGG{arg}.png',
        he_wald = 'results/sim/{model}/HE.wald.AGG{arg}.png',
        reml = 'results/sim/{model}/REML.AGG{arg}.png', 
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
        flag = expand('staging/sim/{model}/all.flag', model=['free','freeW','full']),


######################################################################
# 1.2 HE without JK 
######################################################################
use rule sim_HE as sim_HE_noJK with:
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.noJK.batch{{i}}.npy',
    params:
        free_jk = False,


use rule sim_mergeBatches_HE as sim_mergeBatches_HE_noJK with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/he.noJK.batch{i}.npy'
                for i in range(sim_batch_no)],
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/out.he.noJK.npy',


def sim_agg_he_noJK_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('staging/sim/{{model}}/{params}/out.he.noJK.npy', params=subspace.instance_patterns)


use rule sim_HEestimates_subspace_plot as sim_HE_noJK_estimates_subspace_plot with:
    input:
        out = sim_agg_he_noJK_out_subspace,
        V = sim_agg_he_trueV_subspace,
        W = sim_agg_he_trueW_subspace,
    output:
        png = 'results/sim/{model}/HE.noJK.AGG{arg}.png',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[sim_params['model']==wildcards.model]).iloc[:,:],
        sim_plot_order = sim_plot_order,
        method = None,
        mycolors = mycolors,
        pointcolor = pointcolor,
        colorpalette = colorpalette,


######################################################################
# 1.2 sim with Extra fixed and random
######################################################################
wildcard_constraints: fixed = '\d+'
wildcard_constraints: random = '\d+'


use rule sim_generatedata as sim_generatedata_extra with:
    output:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/fixed_{{fixed}}/random_{{random}}/sim.batch{{i}}.npy',


use rule sim_HE as sim_HE_extra with:
    input:
        data = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/fixed_{{fixed}}/random_{{random}}/sim.batch{{i}}.npy',
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/fixed_{{fixed}}/random_{{random}}/he.batch{{i}}.npy',
    params:
        free_jk = False,


use rule sim_mergeBatches_HE as sim_mergeBatches_HE_extra with:
    input:
        out = [f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/fixed_{{fixed}}/random_{{random}}/he.batch{i}.npy'
                for i in range(sim_batch_no)],
    output:
        out = f'staging/sim/{{model}}/{sim_paramspace.wildcard_pattern}/fixed_{{fixed}}/random_{{random}}/out.he.npy',


def sim_agg_he_extra_out_subspace(wildcards):
    subspace = get_subspace(wildcards.arg, sim_params.loc[sim_params['model']==wildcards.model])
    return expand('staging/sim/{{model}}/{params}/fixed_{{fixed}}/random_{{random}}/out.he.npy', params=subspace.instance_patterns)


use rule sim_HEestimates_subspace_plot as sim_HE_extra_estimates_subspace_plot with:
    input:
        out = sim_agg_he_extra_out_subspace,
        V = sim_agg_he_trueV_subspace,
        W = sim_agg_he_trueW_subspace,
    output:
        png = 'results/sim/{model}/HE.fixed_{fixed}.random_{random}.AGG{arg}.png',
    params:
        subspace = lambda wildcards: get_subspace(wildcards.arg,
                sim_params.loc[sim_params['model']==wildcards.model]).iloc[:,:],
        sim_plot_order = sim_plot_order,
        method = None,
        mycolors = mycolors,
        pointcolor = pointcolor,
        colorpalette = colorpalette,


########################################################################
# Real data
########################################################################


rule gene_annotation_comparison:
    input:
        v24 = 'data/gencode.v24lift37.annotation.gff3.gz',
        v43 = 'data/gencode.v43lift37.annotation.gff3.gz',
        v82 = 'data/Homo_sapiens.GRCh37.82.gff3.gz',
        v85 = 'data/Homo_sapiens.GRCh37.85.gff3.gz',
        genes = 'data/Yazar2022Science/genes.txt',
    output:
        results = 'data/gene.annotation.txt',
    run:
        import gzip, re
        def read_gencode(f):
            res = {}
            for line in gzip.open(f, 'rt'):
                if line[0] != '#':
                    line = line.strip().split()
                    if line[2] == 'gene' and re.search('chr',line[0]) and line[0][3:].isdigit():
                        try:
                            chr, start, end, info = int(line[0][3:]), int(line[3]), int(line[4]), line[-1]
                            info = info.split(';')
                            gene = info[0].split('.')[0][3:]
                            res[gene] = {'chr':chr, 'start':start, 'end':end}
                        except:
                            print(line)
                            sys.exit()
            return( res )

        def read_ensembl(f):
            res = {}
            chrs = []
            for line in gzip.open(f, 'rt'):
                if line[0] != '#':
                    line = line.strip().split('\t')
                    if 'ID=gene' in line[-1] and line[0].isdigit():
                        try:
                            chr, start, end, info = int(line[0]), int(line[3]), int(line[4]), line[-1]
                            info = info.split(';')
                            gene = info[0].split(':')[1]
                            res[gene] = {'chr':chr, 'start':start, 'end':end}
                            chrs.append( chr )
                        except:
                            print(line)
                            sys.exit()
            print( np.unique( chrs, return_counts=True) )
            return( res )

        def compare_genes(v1, v2):
            v1_v2 = set(v1.keys()) & set(v2.keys())
            match = []
            close = []
            for gene in v1_v2:
                v1_gene, v2_gene = v1[gene], v2[gene]
                if v1_gene['chr'] == v2_gene['chr'] and v1_gene['start'] == v2_gene['start'] and v1_gene['end'] == v2_gene['end']:
                    match.append(gene)
                elif v1_gene['chr'] == v2_gene['chr'] and abs(v1_gene['start'] - v2_gene['start']) < 1000 and abs(v1_gene['end'] - v2_gene['end']) < 1000:
                    close.append(gene)
            return( len(v1_v2), len(match), len(close) )

        v24 = read_gencode(input.v24)
        v43 = read_gencode(input.v43)
        v82 = read_ensembl(input.v82)
        v85 = read_ensembl(input.v85)

        genes = pd.read_table(input.genes)['feature'].tolist()
        with open(output.results, 'w') as f:
            f.write(f'{len(v24.keys())}, {len(v43.keys())}, {len(v85.keys())},  and {len(v82.keys())} genes in v24, v43, v85, and v82.\n')
            f.write(f'{len(set(v24.keys()) & set(genes))} v24 genes in data.\n')
            f.write(f'{len(set(v43.keys()) & set(genes))} v43 genes in data.\n')
            f.write(f'{len(set(v82.keys()) & set(genes))} v82 genes in data.\n')
            f.write(f'{len(set(v85.keys()) & set(genes))} v85 genes in data.\n')
            f.write('%i genes overlap between v24 and v43. %i match, and %i close.\n'%(compare_genes(v24,v43)))
            f.write('%i genes overlap between v24 and v82. %i match, and %i close.\n'%(compare_genes(v24,v82)))
            f.write('%i genes overlap between v43 and v82. %i match, and %i close.\n'%(compare_genes(v43,v82)))
            f.write('%i genes overlap between v43 and v85. %i match, and %i close.\n'%(compare_genes(v43,v85)))



########################################################################
# Yazar 2022 Science
########################################################################


# data check
rule yazar_cell_dist:
    input:
        obs = 'data/Yazar2022Science/obs.txt',
    output:
        png = 'results/yazar/cell.dist.png',
    params:
        ind_col = yazar_ind_col,
        ct_col = yazar_ct_col,
    script: 'scripts/yazar/cell_dist.py'


rule yazar_UMIpercell:
    input:
        h5ad = 'data/Yazar2022Science/OneK1K_cohort_gene_expression_matrix_14_celltypes.h5ad.gz',
        var = 'data/Yazar2022Science/var.txt',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
    output:
        umi = 'staging/data/yazar/umi.gz',
    params:
        ind_col = yazar_ind_col,
        ct_col = yazar_ct_col,
    resources:
        mem_mb = '40G',
    run:
        import scanpy as sc
        from scipy import sparse

        var = pd.read_table(input.var, index_col=0)
        if 'feature_is_filtered' in var.columns:
            genes = var.loc[~var['feature_is_filtered']].index.to_numpy()
        else:
            genes = var.index.to_numpy()

        if 'subset_gene' in params.keys():
            # random select genes
            rng = np.random.default_rng(seed=params.seed)
            genes = rng.choice(genes, params.subset_gene, replace=False)

        obs = pd.read_table(input.obs, index_col=0)
        ind_pool = np.unique(obs[params.ind_col].astype('str')+'+'+obs['pool'].astype('str'))

        ann = sc.read_h5ad(input.h5ad, backed='r')
        data = ann[(~ann.obs[params.ind_col].isna())
                & (~ann.obs[params.ct_col].isna())
                & (ann.obs[params.ind_col].astype('str')+'+'+ann.obs['pool'].astype('str')).isin(ind_pool), genes]

        np.savetxt(output.umi, data.X.sum(axis=1).A1)


rule yazar_ctp_dist:
    input:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctp.gz',
    output:
        png = f'results/yazar/{yazar_paramspace.wildcard_pattern}/ctp.dist.png',
    script: 'scripts/yazar/ctp_dist.py'


rule yazar_impute_ctp_s1:
    input:
        data = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctp.gz',
    output:
        data = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctp.b_imputed.gz',
    resources:
        mem_mb = '10G',
        time = '36:00:00',
    run:
        from ctmm import preprocess
        data = pd.read_table(input.data, index_col=(0,1)).astype('float32')
        #data = preprocess.softimpute(data, seed=params.seed)
        data = data.unstack()
        data.to_csv(output.data, sep='\t', index=False, header=False)


rule yazar_impute_ctp_s2:
    # softimpute is too slow for the whold transcriptome
    input:
        data = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctp.b_imputed.gz',
    output:
        data = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctp.a_imputed.gz',
    params:
        seed = 1234567,
    resources:
        mem_mb = '40G',
        time = '36:00:00',
    script: 'bin/yazar/impute.R'


rule yazar_X:
    input:
        h5ad = 'data/Yazar2022Science/OneK1K_cohort_gene_expression_matrix_14_celltypes.h5ad.gz',
    output:
        out = 'data/Yazar2022Science/OneK1K_cohort_gene_expression_matrix_14_celltypes.X.summary',
    run:
        import scanpy as sc
        ann = sc.read_h5ad(input.h5ad, backed='r')
        umi = np.array(ann.X[:, :].sum(axis=1))
        with open(output.out, 'w') as f:
            f.write(f'Total reads per cell: Max {np.amax(umi)}, Min {np.amin(umi)}, Median {np.median(umi)}. \n')
            if np.all(umi == umi.astype('int')):
                f.write('All ints! \n')


##########################################################################
## 1.2 gcta 
##########################################################################


rule yazar_gcta_grm:
    input:
        genes = 'data/Yazar2022Science/gene_loation.txt',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
        bed = 'analysis/yazar/data/geno/chr{chr}.bed',
    output:
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/kinship.chr{{chr}}.txt',
    params:
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/kinship/gene.grm.bin',
        r = int(float(5e5)),
    resources:
        mem_mb = '2G',
    shell:
        '''
        module load gcc/11.3.0 atlas/3.10.3 lapack/3.11.0 plink/1.9
        module load gcc/11.3.0 gcta/1.94.1
        mkdir -p $(dirname {params.kinship})
        python3 bin/yazar/kinship.py \
                {input.genes} {input.P} {params.r} \
                {input.bed} {wildcards.chr} \
                {params.kinship} {output.kinship}
        '''

#use rule yazar_he_kinship_merge as yazar_gcta_grm_merge with:
#    input:
#        kinship = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/kinship.chr{chr}.txt'
#                for chr in range(1,23)],
#    output:
#        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/kinship.txt',
#        save = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/gcta/kinship.txt',

#rule yazar_gcta_HEreg_op:
#    # covar and qcovar is not working in GCTA HE!!!!
#    input:
#        op = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/op.mvn.gz',
#        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/kinship.chr{{chr}}.txt',
#        op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.txt',
#        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
#        meta = 'data/Yazar2022Science/meta.txt',
#    output:
#        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/op.HEreg.chr{{chr}}',
#    params:
#        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/rep/op.HEreg',
#        snps = 5, # threshold of snp number per gene
#    resources:
#        time = '10:00:00',
#        mem_mb = '2G',
#    shell:
#        '''
#        module load gcc/11.3.0 gcta/1.94.1
#        python3 scripts/yazar/gcta_HEreg_op.py \
#                {input.op} {input.kinship} \
#                {input.op_pca} {input.geno_pca} \
#                {input.meta} {params.snps} {params.out} {output.out}
#        '''
#
#rule yazar_gcta_HEreg_op_merge:
#    input:
#        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/op.HEreg.chr{chr}'
#                for chr in range(1,23)],
#    output:
#        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/op.HEreg',
#    run:
#        with open(output.out, 'w') as f:
#            f.write('gene\th2\n')
#            for chr_out in input.out:
#                for gene_out in open(chr_out):
#                    gene = re.search(r"/(ENSG[^/]+)/", gene_out).group(1)
#                    for line in open(gene_out.strip()):
#                        if re.search(r'V\(G\)\/Vp', line):
#                            h2 = line.strip().split()[1]
#                            break
#                    f.write(f'{gene}\t{h2}\n')

kinship_batches = 500
rule yazar_gcta_split_kinship:
    input:
        kinship = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/kinship.chr{chr}.txt'
                for chr in range(1,23)],
    output:
        kinship = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/kinship.batch{i}.txt'
                for i in range(kinship_batches)],
    run:
        kinship = [pd.read_table(f) for f in input.kinship]
        kinship = pd.concat(kinship, axis=0, ignore_index=True)
        indexs = np.array_split(kinship.index, len(output.kinship))
        for index, f in zip(indexs, output.kinship):
            kinship.loc[index].to_csv(f, sep='\t', index=False)

rule yazar_gcta_greml_op:
    input:
        op = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/op.mvn.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/kinship.batch{{i}}.txt',
        op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/op.greml.batch{{i}}',
    params:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/rep/op.hsq',
        snps = 5, # threshold of snp number per gene
    resources:
        time = '200:00:00',
        mem_mb = '2G',
    shell:
        '''
        module load gcc/11.3.0 gcta/1.94.1
        python3 scripts/yazar/gcta_greml_op.py \
                {input.op} {input.kinship} \
                {input.op_pca} {input.geno_pca} \
                {input.obs} {params.snps} {params.out} {output.out}
        '''

rule yazar_gcta_greml_op_merge:
    input:
        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/op.greml.batch{i}'
                for i in range(kinship_batches)],
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/op.greml',
    run:
        with open(output.out, 'w') as f:
            f.write('gene\th2\n')
            for chr_out in input.out:
                for gene_out in open(chr_out):
                    gene = re.search(r"/(ENSG[^/]+)/", gene_out).group(1)
                    for line in open(gene_out.strip()):
                        if re.search(r'V\(G\)\/Vp', line):
                            h2 = line.strip().split()[1]
                            break
                    f.write(f'{gene}\t{h2}\n')

rule yazar_gcta_greml_op_plot:
    input:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/op.greml',
    output:
        png = f'results/yazar/{yazar_paramspace.wildcard_pattern}/gcta/op.greml.png',
    script: 'scripts/yazar/gcta_op_plot.py'

rule yazar_gcta_greml_ctp:
    input:
        ctp = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/kinship.batch{{i}}.txt',
        op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/ctp.greml.batch{{i}}',
    params:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/rep/ctp.hsq',
        snps = 5, # threshold of snp number per gene
    resources:
        time = '200:00:00',
        mem_mb = '2G',
    shell:
        '''
        module load gcc/11.3.0 gcta/1.94.1
        python3 scripts/yazar/gcta_greml_ctp.py \
                {input.ctp} {input.kinship} \
                {input.op_pca} {input.geno_pca} \
                {input.obs} {params.snps} {params.out} {output.out}
        '''

rule yazar_gcta_greml_ctp_merge:
    input:
        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/ctp.greml.batch{i}'
                for i in range(kinship_batches)],
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/ctp.greml',
    run:
        with open(output.out, 'w') as f:
            f.write('ct\tgene\th2\n')
            for chr_out in input.out:
                for gene_out in open(chr_out):
                    gene_ct = re.search(r"/(ENSG[^/]+)/", gene_out).group(1)
                    gene, ct = gene_ct.split('_')
                    for line in open(gene_out.strip()):
                        if re.search(r'V\(G\)\/Vp', line):
                            h2 = line.strip().split()[1]
                            break
                    f.write(f'{ct}\t{gene}\t{h2}\n')

rule yazar_gcta_greml_ctp_plot:
    input:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/ctp.greml',
    output:
        png = f'results/yazar/{yazar_paramspace.wildcard_pattern}/gcta/ctp.greml.png',
    script: 'scripts/yazar/gcta_ctp_plot.py'


rule yazar_gcta_all:
    input:
        op = expand('results/yazar/{params}/gcta/op.greml.png',
                params=yazar_paramspace.instance_patterns),
        ctp = expand('results/yazar/{params}/gcta/ctp.greml.png',
                params=yazar_paramspace.instance_patterns),



##########################################################################
## 1.2: combine cts 
##########################################################################


combine_cts = {
        'cc1':  {
            'B': ['B IN', 'B Mem'],
            'CD4': ['CD4 ET', 'CD4 NC'],
            'CD8': ['CD8 ET', 'CD8 NC', 'CD8 S100B'],
            },

        'cc2':  {
            'main':  ['B IN', 'B Mem', 'CD4 ET', 'CD4 NC', 'CD8 ET', 'CD8 NC', 'CD8 S100B'],
            }
        }

wildcard_constraints: cc='[\w\d]+'


rule yazar_combine_cts:
    input:
        obs = 'staging/data/yazar/obs.gz',
    output:
        obs = 'staging/data/yazar/combine_cts/{cc}/obs.gz',
    params:
        cc = lambda wildcards: combine_cts[wildcards.cc],
    run:
        obs = pd.read_table(input.obs)
        for key, value in params.cc.items():
            obs.loc[obs['cell_label'].isin(value), 'cell_label'] = key
        obs.to_csv(output.obs, sep='\t', index=False)


use rule yazar_ctp as yazar_cc_ctp with:
    input:
        X = 'staging/data/yazar/X.npz',
        obs = 'staging/data/yazar/combine_cts/{cc}/obs.gz',
        var = 'staging/data/yazar/var.gz',
    output:
        ctp = 'staging/yazar/combine_cts/{cc}/ctp.gz',
        ctnu = 'staging/yazar/combine_cts/{cc}/ctnu.gz',
        P = 'staging/yazar/combine_cts/{cc}/P.gz',
        n = 'staging/yazar/combine_cts/{cc}/n.gz',


use rule yazar_rm_rareINDnCT as yazar_cc_rm_rareINDnCT with:
    input:
        ctp = 'staging/yazar/combine_cts/{cc}/ctp.gz',
        ctnu = 'staging/yazar/combine_cts/{cc}/ctnu.gz',
        n = 'staging/yazar/combine_cts/{cc}/n.gz',
    output:
        ctp = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.gz',
        ctnu = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.gz',
        P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        n = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/n.gz',


use rule yazar_mvn_ctp as yazar_cc_mvn_ctp with:
    input:
        data = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.gz',
    output:
        data = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',


use rule yazar_mvn_ctp as yazar_cc_mvn_ctnu with:
    input:
        data = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.gz',
    output:
        data = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',


use rule yazar_std_op as yazar_cc_std_op with:
    input:
        ctp = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',
        ctnu = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',
        P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
    output:
        op = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/op.std.gz',
        nu = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/nu.std.gz',
        ctp = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.std.gz',
        ctnu = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.std.gz',


use rule yazar_op_pca as yazar_cc_op_pca with:
    input:
        op = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/op.std.gz',
    output:
        evec = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/evec.txt',
        eval = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/eval.txt',
        pca = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/pca.txt',
        png = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/pca.png',


use rule yazar_HE_split as yazar_cc_HE_split with:
    input:
        ctp = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.std.gz',
        ctnu = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.std.gz',
    output:
        ctp = [f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{i}.gz'
                for i in range(yazar_he_batches)],
        ctnu = [f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{i}.gz'
                for i in range(yazar_he_batches)],


use rule yazar_HE_free as yazar_cc_HE_free with:
    input:
        ctp = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
        P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'staging/data/yazar/combine_cts/{cc}/obs.gz',   # TODO: try not to use staging
    output:
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/he.free.batch{{i}}.npy',


use rule yazar_HE_free_merge as yazar_cc_HE_free_merge with:
    input:
        out = [f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/he.free.batch{i}.npy'
                for i in range(yazar_he_batches)],
    output:
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/he.free.npy',


use rule yazar_HE_Free_plot as yazar_cc_HE_Free_plot with:
    input:
        P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/he.free.npy',
    output:
        h2 = f'results/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/he.free.h2.png',


use rule yazar_REML_free as yazar_cc_REML_free with:
    input:
        ctp = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',  # TODO: need to use cc kinship
        P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'staging/data/yazar/combine_cts/{cc}/obs.gz',
    output:
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/reml.free.batch{{i}}.npy',
    params:
        snps = 5, # threshold of snp number per gene


use rule yazar_REML_free_merge as yazar_cc_REML_free_merge with:
    input:
        out = [f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/reml.free.batch{i}.npy'
                for i in range(yazar_he_batches)],
    output:
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/reml.free.npy',


use rule yazar_REML_Free_plot as yazar_cc_REML_Free_plot with:
    input:
        P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/reml.free.npy',
    output:
        h2 = f'results/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/reml.free.h2.png',


rule yazar_cc_all:
    input:
        he = expand('results/yazar/combine_cts/{cc}/{params}/he.free.h2.png', 
                cc=list(combine_cts.keys()), 
                params=yazar_paramspace.instance_patterns),
        # reml = expand('results/yazar/combine_cts/{cc}/{params}/reml.free.h2.png', 
        #         cc='cc1', 
        #         params=yazar_paramspace.instance_patterns),




###################################################################
# 1.2.3 no batch effect 
###################################################################
use rule yazar_HE_free as yazar_cc_xbatch_HE_free with:
    input:
        ctp = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
        P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'staging/data/yazar/combine_cts/{cc}/obs.gz',   # TODO: try not to use staging
    output:
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/xbatch.he.free.batch{{i}}.npy',
    params:
        snps = 5, # threshold of snp number per gene
        batch = False,


use rule yazar_HE_free_merge as yazar_cc_xbatch_HE_free_merge with:
    input:
        out = [f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/xbatch.he.free.batch{i}.npy'
                for i in range(yazar_he_batches)],
    output:
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/xbatch.he.free.npy',


use rule yazar_HE_Free_plot as yazar_cc_xbatch_HE_Free_plot with:
    input:
        P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/xbatch.he.free.npy',
    output:
        h2 = f'results/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/xbatch.he.free.h2.png',


##################### 1.2.3 downsample for REML ############################
rule yazar_cc_downREML_free:
    input:
        ctp = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
        P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'staging/data/yazar/combine_cts/{cc}/obs.gz',
    output:
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/downreml.free.batch{{i}}.npy',
    params:
        snps = 5, # threshold of snp number per gene
    resources:
        time = '400:00:00',
        mem_mb = '10G',
        partition = 'tier2q',
    script: 'scripts/yazar/reml.free.py'


use rule yazar_REML_free_merge as yazar_cc_downREML_free_merge with:
    input:
        out = [f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/downreml.free.batch{i}.npy'
                for i in range(yazar_he_batches)],
    output:
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/downreml.free.npy',


use rule yazar_REML_Free_plot as yazar_cc_downREML_Free_plot with:
    input:
        P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/downreml.free.npy',
    output:
        h2 = f'results/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/downreml.free.h2.png',


rule yazar_cc_downsample_all:
    input:
        reml = expand('results/yazar/combine_cts/{cc}/{params}/downreml.free.h2.png', 
                cc='cc1', 
                params=yazar_paramspace.instance_patterns),


# downsample to 500
rule yazar_cc_down2REML_free:
    input:
        ctp = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
        P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'staging/data/yazar/combine_cts/{cc}/obs.gz',
    output:
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/down2reml.free.batch{{i}}.npy',
    params:
        snps = 5, # threshold of snp number per gene
    resources:
        time = '400:00:00',
        mem_mb = '10G',
        partition = 'tier2q',
    script: 'scripts/yazar/reml2.free.py'


use rule yazar_REML_free_merge as yazar_cc_down2REML_free_merge with:
    input:
        out = [f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/down2reml.free.batch{i}.npy'
                for i in range(yazar_he_batches)],
    output:
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/down2reml.free.npy',


use rule yazar_REML_Free_plot as yazar_cc_down2REML_Free_plot with:
    input:
        P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/down2reml.free.npy',
    output:
        h2 = f'results/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/down2reml.free.h2.png',


rule yazar_cc_downsample2_all:
    input:
        reml = expand('results/yazar/combine_cts/{cc}/{params}/down2reml.free.h2.png', 
                cc='cc1', 
                params=yazar_paramspace.instance_patterns),


##################### 1.2.3 gcta ##################################
use rule yazar_gcta_greml_ctp as yazar_cc_gcta_greml_ctp with:
    input:
        ctp = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.std.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/kinship.batch{{i}}.txt',
        op_pca = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'staging/data/yazar/combine_cts/{cc}/obs.gz',
    output:
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/gcta/ctp.greml.batch{{i}}',


use rule yazar_gcta_greml_ctp_merge as yazar_cc_gcta_greml_ctp_merge with:
    input:
        out = [f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/gcta/ctp.greml.batch{i}'
                for i in range(kinship_batches)],
    output:
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/gcta/ctp.greml',


use rule yazar_gcta_greml_ctp_plot as yazar_cc_gcta_greml_ctp_plot with:
    input:
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/gcta/ctp.greml',
    output:
        png = f'results/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/gcta/ctp.greml.png',


rule yazar_cc_gcta_all:
    input:
        png = expand('results/yazar/combine_cts/{cc}/{params}/gcta/ctp.greml.png',
                cc=list(combine_cts.keys()), 
                params=yazar_paramspace.instance_patterns),


##########################################################################
## 1.2 scRNA preprocessing
##########################################################################
pearson_batch_no = 50

rule yazar_pearson_transform:
    input:
        h5ad = 'data/Yazar2022Science/OneK1K_cohort_gene_expression_matrix_14_celltypes.h5ad.gz',
    output:
        h5ad = 'staging/data/yazar/OneK1K_cohort_gene_expression_matrix_14_celltypes.h5ad',
        seurat = 'staging/data/yazar/OneK1K_cohort_gene_expression_matrix_14_celltypes.rds',
        sct = 'staging/data/yazar/OneK1K_cohort_gene_expression_matrix_14_celltypes.sct.h5ad',
    resources:
        partition = 'tier3q',
        mem_mb = '180G',
    shell:
        '''
        cp {input.h5ad} {output.h5ad}
        Rscript bin/yazar/pearson.R {output.h5ad} {output.seurat} {output.sct}
        '''


rule yazar_pearson_ctp_extractX:
    input:
        h5ad = 'staging/data/yazar/OneK1K_cohort_gene_expression_matrix_14_celltypes.sct.h5ad',
        var = 'data/Yazar2022Science/var.txt',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
    output:
        X = temp('staging/data/yazar/pearson/X.{i}.npz'),
        var = temp('staging/data/yazar/pearson/var.{i}.gz'),
    params:
        batch_no = pearson_batch_no,
        ind_col = yazar_ind_col,
        ct_col = yazar_ct_col,
    resources:
        mem_mb = '40G',
    run:
        import scanpy as sc
        from scipy import sparse
        from ctmm import preprocess

        var = pd.read_table(input.var, index_col=0)
        if 'feature_is_filtered' in var.columns:
            genes = var.loc[~var['feature_is_filtered']].index.to_numpy()
        else:
            genes = var.index.to_numpy()

        # intersect genes
        ann = sc.read_h5ad(input.h5ad, backed='r')
        print(ann)
        print(ann.X[:5, :5].toarray())
        print(ann.obs.head())
        print(ann.var.head())

        genes = genes[np.isin(genes, ann.var.index)]

        # batch
        genes = np.array_split(genes, params.batch_no)[int(wildcards.i)]

        obs = pd.read_table(input.obs, index_col=0)
        ind_pool = np.unique(obs[params.ind_col].astype('str')+'+'+obs['pool'].astype('str'))

        data = ann[(~ann.obs[params.ind_col].isna())
                & (~ann.obs[params.ct_col].isna())
                & (ann.obs[params.ind_col].astype('str')+'+'+ann.obs['pool'].astype('str')).isin(ind_pool), genes]
        # normalize and natural logarithm of one plus the input array
        sparse.save_npz(output.X, data.X)

        data.var.rename_axis('feature').to_csv(output.var, sep='\t')


rule yazar_pearson_ctp:
    input:
        X = 'staging/data/yazar/pearson/X.{i}.npz',
        var = 'staging/data/yazar/pearson/var.{i}.gz',
        obs = 'staging/data/yazar/combine_cts/{cc}/obs.gz',
    output:
        ctp = temp('staging/data/yazar/pearson/{cc}/ctp.{i}.gz'),
        ctnu = temp('staging/data/yazar/pearson/{cc}/ctnu.{i}.gz'),
    params:
        ind_col = yazar_ind_col,
        ct_col = yazar_ct_col,
    resources:
        mem_mb = '40G',
    run:
        import scanpy as sc
        from scipy import sparse
        from ctmm import preprocess

        var = pd.read_table(input.var, index_col=0)
        obs = pd.read_table(input.obs, index_col=0)

        # ctp
        X = sparse.load_npz(input.X)
        ctp, ctnu, P, n= preprocess.pseudobulk(X=X, obs=obs, var=var, ind_cut=100, ct_cut=10,
                ind_col=params.ind_col, ct_col=params.ct_col)
        
        # save
        ctp.to_csv(output.ctp, sep='\t')
        ctnu.to_csv(output.ctnu, sep='\t')


rule yazar_pearson_ctp_merge:
    input:
        ctp = expand('staging/data/yazar/pearson/{{cc}}/ctp.{i}.gz', i=np.arange(pearson_batch_no)),
        ctnu = expand('staging/data/yazar/pearson/{{cc}}/ctnu.{i}.gz', i=np.arange(pearson_batch_no)),
    output:
        ctp = 'staging/data/yazar/pearson/{cc}/ctp.gz',
        ctnu = 'staging/data/yazar/pearson/{cc}/ctnu.gz',
    resources:
        mem_mb = '40G',
    run:
        ctp = [pd.read_table(f, index_col=(0,1)) for f in input.ctp]
        ctp_index = ctp[0].index
        for i in range(1, len(ctp)):
            index = ctp[i].index
            if ctp_index.equals(index):
                continue
            else:
                sys.exit('Wrong index in ctp')
        ctp = pd.concat(ctp, axis=1)
        ctp.to_csv(output.ctp, sep='\t')

        ctnu = [pd.read_table(f, index_col=(0,1)) for f in input.ctnu]
        ctnu_index = ctnu[0].index
        for i in range(1, len(ctnu)):
            index = ctnu[i].index
            if ctnu_index.equals(index):
                continue
            else:
                sys.exit('Wrong index in ctnu')
        ctnu = pd.concat(ctnu, axis=1)
        ctnu.to_csv(output.ctnu, sep='\t')


rule yazar_pearson_rm_rareINDnCT:
    # also select gene expressed in all cts
    input:
        ctp = 'staging/data/yazar/pearson/{cc}/ctp.gz',
        ctnu = 'staging/data/yazar/pearson/{cc}/ctnu.gz',
        ref_ctp = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.gz',
    output:
        ctp = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.gz',
        ctnu = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.gz',
    resources:
        mem_mb = '20G',
    run:
        ref_ctp = pd.read_table(input.ref_ctp, index_col=(0, 1))
        index = ref_ctp.index
        genes = ref_ctp.columns.to_numpy()

        ctp = pd.read_table(input.ctp, index_col=(0, 1))
        ctp = ctp.loc[ctp.index.isin(index)][genes]
        ctp.to_csv(output.ctp, sep='\t')

        ctnu = pd.read_table(input.ctnu, index_col=(0, 1))
        ctnu = ctnu.loc[ctnu.index.isin(index)][genes]
        ctnu.to_csv(output.ctnu, sep='\t')


use rule yazar_mvn_ctp as yazar_pearson_mvn_ctp with:
    input:
        data = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.gz',
    output:
        data = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',


use rule yazar_mvn_ctp as yazar_pearson_mvn_ctnu with:
    input:
        data = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.gz',
    output:
        data = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',


use rule yazar_std_op as yazar_pearson_std_op with:
    input:
        ctp = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',
        ctnu = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',
        P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
    output:
        op = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/op.std.gz',
        nu = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/nu.std.gz',
        ctp = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.std.gz',
        ctnu = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.std.gz',


use rule yazar_op_pca as yazar_pearson_op_pca with:
    input:
        op = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/op.std.gz',
    output:
        evec = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/evec.txt',
        eval = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/eval.txt',
        pca = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/pca.txt',
        png = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/pca.png',


use rule yazar_HE_split as yazar_pearson_HE_split with:
    input:
        ctp = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.std.gz',
        ctnu = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.std.gz',
    output:
        ctp = [f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{i}.gz'
                for i in range(yazar_he_batches)],
        ctnu = [f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{i}.gz'
                for i in range(yazar_he_batches)],


use rule yazar_HE_free as yazar_pearson_HE_free with:
    input:
        ctp = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
        P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'staging/data/yazar/combine_cts/{cc}/obs.gz',   # TODO: try not to use staging
    output:
        out = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/he.free.batch{{i}}.npy',
    params:
        snps = 5, # threshold of snp number per gene


use rule yazar_HE_free_merge as yazar_pearson_HE_free_merge with:
    input:
        out = [f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/he.free.batch{i}.npy'
                for i in range(yazar_he_batches)],
    output:
        out = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/he.free.npy',


use rule yazar_HE_Free_plot as yazar_pearson_HE_Free_plot with:
    input:
        P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/he.free.npy',
    output:
        h2 = f'results/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/he.free.h2.png',


use rule yazar_HE_Free_VW_plot as yazar_pearson_HE_Free_VW_plot with:
    input:
        P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'staging/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/he.free.npy',
    output:
        png = f'results/yazar/pearson/{{cc}}/{yazar_paramspace.wildcard_pattern}/he.free.VW.png',


rule yazar_pearson_all:
    input:
        he = expand('results/yazar/pearson/{cc}/{params}/he.free.h2.png', 
                cc=list(combine_cts.keys()), params=yazar_paramspace.instance_patterns),
        W = expand('results/yazar/pearson/{cc}/{params}/he.free.VW.png', 
                cc=list(combine_cts.keys()), params=yazar_paramspace.instance_patterns),






##########################################################################
## 1.2 quantile normalization
##########################################################################
X_batch_no = 50

rule yazar_qn_ctp_extractX:
    input:
        X = 'staging/data/yazar/X.npz',
    output:
        batch = temp('staging/data/yazar/qn/batch{i}.txt'),
        X = temp('staging/data/yazar/qn/X.{i}.gz'),
    params:
        batch_no = X_batch_no,
    resources:
        mem_mb = '40G',
        burden = 100,
    run:
        from scipy import sparse
        from scipy.stats import norm
        from ctmm import preprocess


        def quantnorm(y):
            ranks = pd.Series(y).rank().to_numpy()
            quantiles = norm.ppf( ( ranks + 1 ) / (len(ranks) + 2) )
            if np.any(ranks != pd.Series(quantiles).rank().to_numpy()):
                sys.exit('Bad rank problem')

            return quantiles


        X = preprocess.normalize(sparse.load_npz(input.X), l=1e4).tocsc()
        batch = np.array_split(np.arange(X.shape[1]), params.batch_no)[int(wildcards.i)]
        np.savetxt(output.batch, batch, '%d')

        X = np.apply_along_axis(quantnorm, 0, X[:, batch].toarray())
        np.savetxt(output.X, X)


rule yazar_qn_ctp:
    input:
        batch = 'staging/data/yazar/qn/batch{i}.txt',
        X = 'staging/data/yazar/qn/X.{i}.gz',
        obs = 'staging/data/yazar/combine_cts/{cc}/obs.gz',
        var = 'staging/data/yazar/var.gz',
    output:
        ctp = temp('staging/data/yazar/qn/{cc}/ctp.{i}.gz'),
        ctnu = temp('staging/data/yazar/qn/{cc}/ctnu.{i}.gz'),
        P = temp('staging/data/yazar/qn/{cc}/P.{i}.gz'),
        n = temp('staging/data/yazar/qn/{cc}/n.{i}.gz'),
    params:
        ind_col = yazar_ind_col,
        ct_col = yazar_ct_col,
    resources:
        mem_mb = '40G',
    run:
        from scipy import sparse
        from ctmm import preprocess

        batch = np.loadtxt(input.batch, 'int')
        X = np.loadtxt(input.X)
        obs = pd.read_table(input.obs, index_col=0)
        var = pd.read_table(input.var, index_col=0)
        genes = var.index.to_numpy()[batch]
        var = var.loc[genes]
        ctp, ctnu, P, n = preprocess.pseudobulk(X=X, obs=obs, var=var, ind_cut=100, ct_cut=10,
                ind_col=params.ind_col, ct_col=params.ct_col)

        # save
        ctp.to_csv(output.ctp, sep='\t')
        ctnu.to_csv(output.ctnu, sep='\t')
        P.to_csv(output.P, sep='\t')
        n.to_csv(output.n, sep='\t')


use rule yazar_pearson_ctp_merge as yazar_qn_ctp_merge with:
    input:
        ctp = expand('staging/data/yazar/qn/{{cc}}/ctp.{i}.gz', i=np.arange(X_batch_no)),
        ctnu = expand('staging/data/yazar/qn/{{cc}}/ctnu.{i}.gz', i=np.arange(X_batch_no)),
    output:
        ctp = 'staging/data/yazar/qn/{cc}/ctp.gz',
        ctnu = 'staging/data/yazar/qn/{cc}/ctnu.gz',


use rule yazar_pearson_rm_rareINDnCT as yazar_qn_rm_rareINDnCT with:
    # also select gene expressed in all cts
    input:
        ctp = 'staging/data/yazar/qn/{cc}/ctp.gz',
        ctnu = 'staging/data/yazar/qn/{cc}/ctnu.gz',
        ref_ctp = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.gz',
    output:
        ctp = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.gz',
        ctnu = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.gz',


use rule yazar_mvn_ctp as yazar_qn_mvn_ctp with:
    input:
        data = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.gz',
    output:
        data = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',


use rule yazar_mvn_ctp as yazar_qn_mvn_ctnu with:
    input:
        data = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.gz',
    output:
        data = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',


use rule yazar_std_op as yazar_qn_std_op with:
    input:
        ctp = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',
        ctnu = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',
        P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
    output:
        op = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/op.std.gz',
        nu = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/nu.std.gz',
        ctp = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.std.gz',
        ctnu = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.std.gz',


use rule yazar_op_pca as yazar_qn_op_pca with:
    input:
        op = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/op.std.gz',
    output:
        evec = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/evec.txt',
        eval = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/eval.txt',
        pca = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/pca.txt',
        png = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/pca.png',


use rule yazar_HE_split as yazar_qn_HE_split with:
    input:
        ctp = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.std.gz',
        ctnu = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.std.gz',
    output:
        ctp = [f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{i}.gz'
                for i in range(yazar_he_batches)],
        ctnu = [f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{i}.gz'
                for i in range(yazar_he_batches)],


use rule yazar_HE_free as yazar_qn_HE_free with:
    input:
        ctp = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
        P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'staging/data/yazar/combine_cts/{cc}/obs.gz',   # TODO: try not to use staging
    output:
        out = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/he.free.batch{{i}}.npy',
    params:
        snps = 5, # threshold of snp number per gene


use rule yazar_HE_free_merge as yazar_qn_HE_free_merge with:
    input:
        out = [f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/he.free.batch{i}.npy'
                for i in range(yazar_he_batches)],
    output:
        out = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/he.free.npy',


use rule yazar_HE_Free_plot as yazar_qn_HE_Free_plot with:
    input:
        P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'staging/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/he.free.npy',
    output:
        h2 = f'results/yazar/qn/{{cc}}/{yazar_paramspace.wildcard_pattern}/he.free.h2.png',


rule yazar_qn_all:
    input:
        he = expand('results/yazar/qn/{cc}/{params}/he.free.h2.png', 
                cc=list(combine_cts.keys()), 
                params=yazar_paramspace.instance_patterns),


##########################################################################
## 1.2: test L = 1e4, 1e6
##########################################################################
rule yazar_L_ctp_extractX:
    input:
        h5ad = 'data/Yazar2022Science/OneK1K_cohort_gene_expression_matrix_14_celltypes.h5ad.gz',
        var = 'data/Yazar2022Science/var.txt',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
    output:
        X = 'staging/data/yazar/L{l}/X.npz',
        obs = 'staging/data/yazar/L{l}/obs.gz',
        var = 'staging/data/yazar/L{l}/var.gz',
    params:
        ind_col = yazar_ind_col,
        ct_col = yazar_ct_col,
    resources:
        mem_mb = '40G',
    run:
        import scanpy as sc
        from scipy import sparse
        from ctmm import preprocess

        var = pd.read_table(input.var, index_col=0)
        if 'feature_is_filtered' in var.columns:
            genes = var.loc[~var['feature_is_filtered']].index.to_numpy()
        else:
            genes = var.index.to_numpy()

        if 'subset_gene' in params.keys():
            # random select genes
            rng = np.random.default_rng(seed=params.seed)
            genes = rng.choice(genes, params.subset_gene, replace=False)

        obs = pd.read_table(input.obs, index_col=0)
        ind_pool = np.unique(obs[params.ind_col].astype('str')+'+'+obs['pool'].astype('str'))

        ann = sc.read_h5ad(input.h5ad, backed='r')
        data = ann[(~ann.obs[params.ind_col].isna())
                & (~ann.obs[params.ct_col].isna())
                & (ann.obs[params.ind_col].astype('str')+'+'+ann.obs['pool'].astype('str')).isin(ind_pool), genes]
        # normalize and natural logarithm of one plus the input array
        X = preprocess.normalize(data.X, float(wildcards.l)).log1p()
        sparse.save_npz(output.X, X)

        data.obs.rename_axis('cell').to_csv(output.obs, sep='\t')
        data.var.rename_axis('feature').to_csv(output.var, sep='\t')


use rule yazar_ctp as yazar_L_ctp with:
    input:
        X = 'staging/data/yazar/L{l}/X.npz',
        obs = 'staging/data/yazar/L{l}/obs.gz',
        var = 'staging/data/yazar/L{l}/var.gz',
    output:
        ctp = 'staging/data/yazar/L{l}/ctp.gz',
        ctnu = 'staging/data/yazar/L{l}/ctnu.gz',
        P = 'staging/data/yazar/L{l}/P.gz',
        n = 'staging/data/yazar/L{l}/n.gz',


use rule yazar_rm_rareINDnCT as yazar_L_rm_rareINDnCT with:
    # also select gene expressed in all cts
    input:
        ctp = 'staging/data/yazar/L{l}/ctp.gz',
        ctnu = 'staging/data/yazar/L{l}/ctnu.gz',
        n = 'staging/data/yazar/L{l}/n.gz',
    output:
        ctp = f'staging/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/ctp.gz',
        ctnu = f'staging/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/ctnu.gz',
        P = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        n = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/n.gz',


use rule yazar_mvn_ctp as yazar_L_mvn_ctp with:
    input:
        data = f'staging/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/ctp.gz',
    output:
        data = f'staging/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',


use rule yazar_mvn_ctp as yazar_L_mvn_ctnu with:
    input:
        data = f'staging/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/ctnu.gz',
    output:
        data = f'staging/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',


use rule yazar_std_op as yazar_L_std_op with:
    input:
        ctp = f'staging/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',
        ctnu = f'staging/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',
        P = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/P.gz',
    output:
        op = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/op.mvn.gz',
        nu = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/nu.mvn.gz',
        ctp = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',
        ctnu = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',


use rule yazar_op_pca as yazar_L_op_pca with:
    input:
        op = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/op.mvn.gz',
    output:
        evec = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/evec.txt',
        eval = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/eval.txt',
        pca = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/pca.txt',
        png = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/pca.png',


use rule yazar_HE_split as yazar_L_HE_split with:
    input:
        ctp = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',
        ctnu = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',
    output:
        ctp = [f'staging/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{i}.gz'
                for i in range(yazar_he_batches)],
        ctnu = [f'staging/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{i}.gz'
                for i in range(yazar_he_batches)],


use rule yazar_HE_free as yazar_L_HE_free with:
    input:
        ctp = f'staging/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
        P = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'staging/data/yazar/obs.gz',   # TODO: try not to use staging
    output:
        out = f'staging/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/he.free.batch{{i}}.npy',
    params:
        snps = 5, # threshold of snp number per gene


use rule yazar_HE_free_merge as yazar_L_HE_free_merge with:
    input:
        out = [f'staging/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/he.free.batch{i}.npy'
            for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/he.free.npy',


use rule yazar_HE_Free_plot as yazar_L_HE_Free_plot with:
    input:
        P = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/he.free.npy',
    output:
        h2 = f'results/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/he.free.h2.png',


use rule yazar_HE_Free_VW_plot as yazar_L_HE_Free_VW_plot with:
    input:
        P = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/he.free.npy',
    output:
        png = f'results/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/he.free.VW.png',


rule yazar_L_all:
    input:
        free = expand('results/yazar/L{l}/{params}/he.free.h2.png', 
                        l=['1e4'], params=yazar_paramspace.instance_patterns),
        free_W = expand('results/yazar/L{l}/{params}/he.free.VW.png',
                        l=['1e4'], params=yazar_paramspace.instance_patterns),




##########################################################################
## 1.2.3: test HE without correcting for batch effect
##########################################################################
use rule yazar_HE_free as yazar_L_nobatch_HE_free with:
    input:
        ctp = f'staging/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
        P = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'staging/data/yazar/obs.gz',   # TODO: try not to use staging
    output:
        out = f'staging/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/he.nobatch.free.batch{{i}}.npy',
    params:
        snps = 5, # threshold of snp number per gene
        batch = False,


use rule yazar_HE_free_merge as yazar_L_nobatch_HE_free_merge with:
    input:
        out = [f'staging/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/he.nobatch.free.batch{i}.npy'
            for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/he.nobatch.free.npy',


use rule yazar_HE_Free_plot as yazar_L_nobatch_HE_Free_plot with:
    input:
        P = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/he.nobatch.free.npy',
    output:
        h2 = f'results/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/he.nobatch.free.h2.png',


use rule yazar_HE_Free_VW_plot as yazar_L_nobatch_HE_Free_VW_plot with:
    input:
        P = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/he.nobatch.free.npy',
    output:
        png = f'results/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/he.nobatch.free.VW.png',


rule yazar_L_nobatch_all:
    input:
        free = expand('results/yazar/L{l}/{params}/he.nobatch.free.h2.png', 
                        l=['1e4'], params=yazar_paramspace.instance_patterns),
        free_W = expand('results/yazar/L{l}/{params}/he.nobatch.free.W.png',
                        l=['1e4'], params=yazar_paramspace.instance_patterns),


##########################################################################
## 1.2.3: test L = 1e4: test imputation scheme (impute target cell types with all cell types)
##########################################################################

rule yazar_L_impute_selectcts:
    input:
        ctp = f'staging/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',
        ctnu = f'staging/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',
        n = f'analysis/yazar/L{{l}}/{yazar_paramspace.wildcard_pattern}/n.gz',
    output:
        ctp = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',
        ctnu = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',
        P = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/P.gz',
        n = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/n.gz',
    params:
        cts = ['CD4 NC', 'CD8 ET', 'NK', 'CD8 NC'],
    run:
        ctp = pd.read_table(input.ctp)
        ctp.loc[ctp['ct'].isin(params.cts)].to_csv(output.ctp, sep='\t', index=False)

        ctnu = pd.read_table(input.ctnu)
        ctnu.loc[ctnu['ct'].isin(params.cts)].to_csv(output.ctnu, sep='\t', index=False)

        n = pd.read_table(input.n, index_col=0)
        cts = n.columns.to_numpy()
        cts = cts[np.isin(cts, params.cts)]
        n = n[cts]
        n.to_csv(output.n, sep='\t')
        
        P = n.divide(n.sum(axis=1), axis=0)
        P.to_csv(output.P, sep='\t')


use rule yazar_std_op as yazar_L_impute_std_op with:
    input:
        ctp = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',
        ctnu = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',
        P = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/P.gz',
    output:
        op = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/op.std.gz',
        nu = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/nu.std.gz',
        ctp = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/ctp.std.gz',
        ctnu = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/ctnu.std.gz',


use rule yazar_op_pca as yazar_L_impute_op_pca with:
    input:
        op = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/op.std.gz',
    output:
        evec = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/evec.txt',
        eval = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/eval.txt',
        pca = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/pca.txt',
        png = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/pca.png',


use rule yazar_HE_split as yazar_L_impute_HE_split with:
    input:
        ctp = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/ctp.std.gz',
        ctnu = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/ctnu.std.gz',
    output:
        ctp = [f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{i}.gz'
                for i in range(yazar_he_batches)],
        ctnu = [f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{i}.gz'
                for i in range(yazar_he_batches)],


use rule yazar_HE_free as yazar_L_impute_HE_free with:
    input:
        ctp = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
        P = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'staging/data/yazar/obs.gz',   # TODO: try not to use staging
    output:
        out = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/he.free.batch{{i}}.npy',
    params:
        snps = 5, # threshold of snp number per gene


use rule yazar_HE_free_merge as yazar_L_impute_HE_free_merge with:
    input:
        out = [f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/he.free.batch{i}.npy'
            for i in range(yazar_he_batches)],
    output:
        out = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/he.free.npy',


use rule yazar_HE_Free_plot as yazar_L_impute_HE_Free_plot with:
    input:
        P = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/he.free.npy',
    output:
        h2 = f'results/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/he.free.h2.png',


use rule yazar_HE_Free_VW_plot as yazar_L_impute_HE_Free_VW_plot with:
    input:
        P = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'staging/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/he.free.npy',
    output:
        png = f'results/yazar/L{{l}}/impute/{yazar_paramspace.wildcard_pattern}/he.free.VW.png',


rule yazar_L_impute_all:
    input:
        # free = expand('results/yazar/L{l}/impute/{params}/he.free.h2.png', 
                        # l=['1e4'], params=yazar_paramspace.instance_patterns),
        free = expand('results/yazar/L{l}/impute/ind_min_cellnum~10_ct_min_cellnum~10_prop~0.9_sex~Y_PC~6_experiment~R_disease~Y/he.free.h2.png', 
                        l=['1e4']),
        # free_W = expand('results/yazar/L{l}/impute/{params}/he.free.VW.png',
                        # l=['1e4'], params=yazar_paramspace.instance_patterns),
        free_W = expand('results/yazar/L{l}/impute/ind_min_cellnum~10_ct_min_cellnum~10_prop~0.9_sex~Y_PC~6_experiment~R_disease~Y/he.free.VW.png',
                        l=['1e4']),




##########################################################################
## 1.2: test gene regions
##########################################################################

use rule yazar_he_kinship as yazar_r_he_kinship with:
    input:
        genes = 'data/Yazar2022Science/gene_loation.txt',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
        bed = 'analysis/yazar/data/geno/chr{chr}.bed',
    output:
        kinship = temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/r{{r}}/kinship.chr{{chr}}.npy'),
    params:
        r = lambda wildcards: int(float(wildcards.r)),


use rule yazar_he_kinship_merge as yazar_r_he_kinship_merge with:
    input:
        kinship = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/r{{r}}/kinship.chr{chr}.npy'
                for chr in range(1,23)],
    output:
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/r{{r}}/kinship.npz',


use rule yazar_he_kinship_split as yazar_r_he_kinship_split with:
    input:
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/r{{r}}/kinship.npz',
        ctp = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{i}.gz'
                for i in range(yazar_he_batches)],
    output:
        kinship = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/r{{r}}/kinship.batch{i}.npy'
                for i in range(yazar_he_batches)],


use rule yazar_HE_free as yazar_r_HE_free with:
    input:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/r{{r}}/kinship.batch{{i}}.npy',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'staging/data/yazar/obs.gz',   # TODO: try not to use staging
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/r{{r}}/he.free.batch{{i}}.npy',


use rule yazar_HE_free_merge as yazar_r_HE_free_merge with:
    input:
        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/r{{r}}/he.free.batch{i}.npy'
            for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/r{{r}}/he.free.npy',


use rule yazar_HE_Free_plot as yazar_r_HE_Free_plot with:
    input:
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/r{{r}}/he.free.npy',
    output:
        h2 = f'results/yazar/{yazar_paramspace.wildcard_pattern}/r{{r}}/he.free.h2.png',


use rule yazar_HE_Free_VW_plot as yazar_r_HE_Free_VW_plot with:
    input:
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/r{{r}}/he.free.npy',
    output:
        png = f'results/yazar/{yazar_paramspace.wildcard_pattern}/r{{r}}/he.free.VW.png',


# use rule yazar_clean as yazar_r_clean with:
#     input:
#         out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/r{{r}}/he.free.npy',
#         kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/r{{r}}/kinship.txt',
#     output:
#         touch(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/r{{r}}/he.done')


rule yazar_r_all:
    input:
        free = expand('results/yazar/{params}/r{r}/he.free.h2.png', 
                        params=yazar_paramspace.instance_patterns, r=[50000, 10000]),
        VW = expand('results/yazar/{params}/r{r}/he.free.VW.png',
                        params=yazar_paramspace.instance_patterns, r=[50000, 10000]),


##########################################################################
## 1.2: trans effect
##########################################################################

rule yazar_trans_genome_kinship:
    input:
        bed = expand('analysis/yazar/data/geno/chr{chr}.bed',
                chr=range(1,23)),
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
    output:
        kinship = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/trans/he/genome.rel',
    params:
        prefix = lambda wildcards, output: os.path.splitext(output.kinship)[0],
        tmp_dir = lambda wildcards, output: os.path.splitext(output.kinship)[0] + '_tmp',
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
        # grm
        plink --bfile $merged.ld --make-rel --out {params.prefix}
        rm -r {params.tmp_dir}
        '''


use rule yazar_HE_free as yazar_trans_HE_free with:
    input:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.batch{{i}}.npy',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'staging/data/yazar/obs.gz',   # TODO: try not to use staging
        genome = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/trans/he/genome.rel',  # don't gz
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/trans/he.free.batch{{i}}.npy',
    resources:
        mem_mb = '30G',


use rule yazar_HE_free_merge as yazar_trans_HE_free_merge with:
    input:
        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/trans/he.free.batch{i}.npy'
            for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/trans/he.free.npy',


use rule yazar_HE_Free_plot as yazar_trans_HE_Free_plot with:
    input:
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/trans/he.free.npy',
    output:
        h2 = f'results/yazar/{yazar_paramspace.wildcard_pattern}/trans/he.free.h2.png',


use rule yazar_HE_Free_VW_plot as yazar_trans_HE_Free_VW_plot with:
    input:
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/trans/he.free.npy',
    output:
        png = f'results/yazar/{yazar_paramspace.wildcard_pattern}/trans/he.free.VW.png',


rule yazar_trans_all:
    input:
        free = expand('results/yazar/{params}/trans/he.free.h2.png', 
                        params=yazar_paramspace.instance_patterns),
        VW = expand('results/yazar/{params}/trans/he.free.VW.png',
                        params=yazar_paramspace.instance_patterns),



















































# ##########################################################################
# #  1.2: remove missing individuals   
# #     .4: eQTL: PADI4 rs10788663
# ##########################################################################
# rule yazar_he_eqtl_kinship:
#     input:
#         genes = 'data/Yazar2022Science/gene_loation.txt',
#         P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
#         bed = expand('analysis/yazar/data/geno/chr{chr}.bed', chr=range(1,23)),
#     output:
#         kinship = temp(f'staging/yazar/nomissing/eqtl/{yazar_paramspace.wildcard_pattern}/he/kinship.{{gene}}.npy'),
#     resources:
#         mem_mb = '20G',
#     shell: 
#         '''
#         module load gcc/11.3.0 atlas/3.10.3 lapack/3.11.0 plink/1.9
#         python3 bin/yazar/kinship.npy.py {input.genes} {input.P} {params.r} {input.bed} {wildcards.chr} \
#                         {output.kinship} 
#         '''










##########################################################################
#  1.2: remove missing individuals   
#     .3: add trans effect
##########################################################################

use rule yazar_trans_genome_kinship as yazar_nomissing_trans_genome_kinship with:
    input:
        bed = expand('analysis/yazar/data/geno/chr{chr}.bed', chr=range(1,23)),
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
    output:
        kinship = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/trans/he/genome.rel',


use rule yazar_HE_free as yazar_nomissing_trans_HE_free with:
    input:
        ctp = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/kinship.batch{{i}}.npy',
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
        genome = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/trans/he/genome.rel',  # don't gz
    output:
        out = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/trans/he.free.batch{{i}}.npy',


use rule yazar_HE_free_merge as yazar_nomissing_trans_HE_free_merge with:
    input:
        out = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/trans/he.free.batch{i}.npy'
            for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/trans/he.free.npy',


use rule yazar_HE_Free_plot as yazar_nomissing_trans_HE_Free_plot with:
    input:
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/trans/he.free.npy',
    output:
        h2 = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/trans/he.free.h2.png',


use rule yazar_HE_Free_VW_plot as yazar_nomissing_trans_HE_Free_VW_plot with:
    input:
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/trans/he.free.npy',
    output:
        png = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/trans/he.free.VW.png',


rule yazar_nomissing_trans_all:
    input:
        free = expand('results/yazar/nomissing/{params}/trans/he.free.h2.png', 
                        params=yazar_paramspace.instance_patterns),
        VW = expand('results/yazar/nomissing/{params}/trans/he.free.VW.png',
                        params=yazar_paramspace.instance_patterns),


##########################################################################
#  1.2: remove missing individuals   
#     .3: add trans effect
#       .4: trans only no cis: use trans K replace cis K
##########################################################################

rule yazar_nomissing_transonly_HE_free:
    input:
        ctp = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
        genome = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/trans/he/genome.rel',  # don't gz
    output:
        out = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/transonly/he.free.batch{{i}}.npy',
    resources:
        mem_mb = '30G',
    script: 'scripts/yazar/he.free.transonly.py' 


use rule yazar_HE_free_merge as yazar_nomissing_transonly_HE_free_merge with:
    input:
        out = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/transonly/he.free.batch{i}.npy'
            for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/transonly/he.free.npy',


use rule yazar_HE_Free_plot as yazar_nomissing_transonly_HE_Free_plot with:
    input:
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/transonly/he.free.npy',
    output:
        h2 = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/transonly/he.free.h2.png',


use rule yazar_HE_Free_VW_plot as yazar_nomissing_transonly_HE_Free_VW_plot with:
    input:
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/transonly/he.free.npy',
    output:
        png = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/transonly/he.free.VW.png',


rule yazar_nomissing_transonly_all:
    input:
        free = expand('results/yazar/nomissing/{params}/transonly/he.free.h2.png', 
                        params=yazar_paramspace.instance_patterns),
        VW = expand('results/yazar/nomissing/{params}/transonly/he.free.VW.png',
                        params=yazar_paramspace.instance_patterns),




##########################################################################
#  1.2: remove missing individuals   
#     .3: cis region size
##########################################################################

use rule yazar_he_kinship as yazar_nomissing_r_he_kinship with:
    input:
        genes = 'data/Yazar2022Science/gene_loation.txt',
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        bed = 'analysis/yazar/data/geno/chr{chr}.bed',
    output:
        kinship = temp(f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/r{{r}}/kinship.chr{{chr}}.npy'),
    params:
        r = lambda wildcards: int(float(wildcards.r)),


use rule yazar_he_kinship_merge as yazar_nomissing_r_he_kinship_merge with:
    input:
        kinship = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/r{{r}}/kinship.chr{chr}.npy'
                for chr in range(1,23)],
    output:
        kinship = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/r{{r}}/kinship.npz',


use rule yazar_he_kinship_split as yazar_nomissing_r_he_kinship_split with:
    input:
        kinship = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/r{{r}}/kinship.npz',
        ctp = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{i}.gz'
                for i in range(yazar_he_batches)],
    output:
        kinship = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/r{{r}}/kinship.batch{i}.npy'
                for i in range(yazar_he_batches)],


use rule yazar_HE_free as yazar_nomissing_r_HE_free with:
    input:
        ctp = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/r{{r}}/kinship.batch{{i}}.npy',
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
    output:
        out = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/r{{r}}/he.free.batch{{i}}.npy',


use rule yazar_HE_free_merge as yazar_nomissing_r_HE_free_merge with:
    input:
        out = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/r{{r}}/he.free.batch{i}.npy'
            for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/r{{r}}/he.free.npy',


use rule yazar_HE_Free_plot as yazar_nomissing_r_HE_Free_plot with:
    input:
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/r{{r}}/he.free.npy',
    output:
        h2 = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/r{{r}}/he.free.h2.png',


use rule yazar_HE_Free_VW_plot as yazar_nomissing_r_HE_Free_VW_plot with:
    input:
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/r{{r}}/he.free.npy',
    output:
        png = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/r{{r}}/he.free.VW.png',


rule yazar_nomissing_r_all:
    input:
        free = expand('results/yazar/nomissing/{params}/r{r}/he.free.h2.png', 
                        params=yazar_paramspace.instance_patterns, r=[100000000, 50000, 10000]),
        VW = expand('results/yazar/nomissing/{params}/r{r}/he.free.VW.png',
                        params=yazar_paramspace.instance_patterns, r=[100000000, 50000, 10000]),

##################################################################################################################
#  1.2: remove missing individuals   
#     .3: ctnu = 0
##################################################################################################################
rule yazar_nomissing_ctnu0:
    input:
        ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
    output:
        ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he_ctnu0/ctnu.batch{{i}}.gz',
    run:
        ctnu = pd.read_table(input.ctnu, index_col=(0, 1))
        ctnu.values[:] = 0
        ctnu.to_csv(output.ctnu, sep='\t')


use rule yazar_HE_free as yazar_nomissing_ctnu0_HE_free with:
    input:
        ctp = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he_ctnu0/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/kinship.batch{{i}}.npy',
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
    output:
        out = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.ctnu0.free.batch{{i}}.npy',


use rule yazar_HE_free_merge as yazar_nomissing_ctnu0_HE_free_merge with:
    input:
        out = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.ctnu0.free.batch{i}.npy'
            for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.ctnu0.free.npy',


use rule yazar_HE_Free_plot as yazar_nomissing_ctnu0_HE_Free_plot with:
    input:
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.ctnu0.free.npy',
    output:
        h2 = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.ctnu0.free.h2.png',


use rule yazar_HE_Free_VW_plot as yazar_nomissing_ctnu0_HE_Free_VW_plot with:
    input:
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.ctnu0.free.npy',
    output:
        png = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.ctnu0.free.VW.png',


rule yazar_nomissing_ctnu0_all:
    input:
        free = expand('results/yazar/nomissing/{params}/he.ctnu0.free.h2.png', params=yazar_paramspace.instance_patterns),
        VW = expand('results/yazar/nomissing/{params}/he.ctnu0.free.VW.png', params=yazar_paramspace.instance_patterns),



##################################################################################################################
#  1.2: remove missing individuals   
#  .3: cell type-specific fixed and random effect
##################################################################################################################
rule yazar_nomissing_ctspecific_HE_free:
    input:
        ctp = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/kinship.batch{{i}}.npy',
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        op_pca = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/pca.txt',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
    output:
        out = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.ctspecific.free.batch{{i}}.npy',
    params:
        snps = 5, # threshold of snp number per gene
    resources:
        mem_mb = '30G',
    script: 'scripts/yazar/he.free.ctspecific.py' 


use rule yazar_HE_free_merge as yazar_nomissing_ctspecific_HE_free_merge with:
    input:
        out = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.ctspecific.free.batch{i}.npy'
            for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.ctspecific.free.npy',


use rule yazar_HE_Free_plot as yazar_nomissing_ctspecific_HE_Free_plot with:
    input:
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.ctspecific.free.npy',
    output:
        h2 = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.ctspecific.free.h2.png',


use rule yazar_HE_Free_VW_plot as yazar_nomissing_ctspecific_HE_Free_VW_plot with:
    input:
        P = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/P.gz',
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.ctspecific.free.npy',
    output:
        png = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he.ctspecific.free.VW.png',


rule yazar_nomissing_ctspecific_all:
    input:
        free = expand('results/yazar/nomissing/{params}/he.ctspecific.free.h2.png', params=yazar_paramspace.instance_patterns),
        VW = expand('results/yazar/nomissing/{params}/he.ctspecific.free.VW.png', params=yazar_paramspace.instance_patterns),



##################################################################################################################
#  1.2: remove missing individuals   
#     .3: choose one ct and copy and permute to simulate 
#       a. no cis genetic effect 
#       b. no genetic correlation
##################################################################################################################
target_ct = 'CD4 NC'


# a. no cis genetic effect
rule yazar_nomissing_permuteA_selectCT:
    input:
        ctp = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/ctp.std.gz',
        ctnu = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/ctnu.std.gz',
    output:
        ctp = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/ctp.gz',
        ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/ctnu.gz',
        P = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/P.gz',
    params:
        ct = target_ct,
        seed = 1983,
    run:
        rng = np.random.default_rng(params.seed)

        ctp = pd.read_table(input.ctp)
        ctnu = pd.read_table(input.ctnu)

        # choose ct
        ctp = ctp.loc[ctp['ct'] == params.ct]
        ctnu = ctnu.loc[ctnu['ct'] == params.ct]

        # santiy check
        if np.any(ctp['ind'] != ctnu['ind']):
            sys.exit('Individual order not matching!')

        # copy
        ctpA = ctp.copy()
        ctpA['ct'] = ctpA['ct'] + ':A'
        ctnuA = ctnu.copy()
        ctnuA['ct'] = ctnuA['ct'] + ':A'
        ctpB = ctp.copy()
        ctpB['ct'] = ctpB['ct'] + ':B'
        ctnuB = ctnu.copy()
        ctnuB['ct'] = ctnuB['ct'] + ':B'

        # permute
        # TODO: double check        
        genes = ctp.columns.tolist()[2:]
        N = len(np.unique(ctp['ind']))
        for gene in genes:
            shuffled_indices = rng.permutation(N)
            ctpA[gene] = ctpA[gene].to_numpy()[shuffled_indices]
            ctnuA[gene] = ctnuA[gene].to_numpy()[shuffled_indices]
            shuffled_indices = rng.permutation(N)
            ctpB[gene] = ctpB[gene].to_numpy()[shuffled_indices]
            ctnuB[gene] = ctnuB[gene].to_numpy()[shuffled_indices]
        # A_inds = ctp['ind'].sample(frac=1, random_state=rng.integers(10000)).to_numpy()
        # ctpA['ind'] = A_inds
        # ctnuA['ind'] = A_inds
        # B_inds = ctp['ind'].sample(frac=1, random_state=rng.integers(10000)).to_numpy()
        # ctpB['ind'] = B_inds
        # ctnuB['ind'] = B_inds

        # merge
        ctp = pd.concat([ctpA, ctpB], ignore_index=True)
        ctnu = pd.concat([ctnuA, ctnuB], ignore_index=True)

        # sort
        ctp = ctp.sort_values(by=['ind', 'ct'])
        ctnu = ctnu.sort_values(by=['ind', 'ct'])

        ctp.to_csv(output.ctp, sep='\t', index=False)
        ctnu.to_csv(output.ctnu, sep='\t', index=False)

        # make P of 1 / # cts
        p = 1 / len(np.unique(ctp['ct']))
        inds = np.unique(ctp['ind'])
        P = pd.DataFrame({'ind': inds})
        for ct in np.unique(ctp['ct']):
            P[ct] = p
        P.to_csv(output.P, sep='\t', index=False)


use rule yazar_HE_split as yazar_nomissing_permuteA_HE_split with:
    input:
        ctp = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/ctp.gz',
        ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/ctnu.gz',
    output:
        ctp = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/he/ctp.batch{i}.gz'
                for i in range(yazar_he_batches)],
        ctnu = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/he/ctnu.batch{i}.gz'
                for i in range(yazar_he_batches)],


rule yazar_nomissing_permuteA_HE_free:
    input:
        ctp = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/kinship.batch{{i}}.npy',
        P = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/P.gz',
    output:
        out = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/he.free.batch{{i}}.npy',
    params:
        snps = 5, # threshold of snp number per gene
    resources:
        mem_mb = '30G',
    script: 'scripts/yazar/nomissing.permute.he.free.py' 


use rule yazar_HE_free_merge as yazar_nomissing_permuteA_HE_free_merge with:
    input:
        out = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/he.free.batch{i}.npy'
            for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/he.free.npy',


# use rule yazar_HE_Free_plot as yazar_nomissing_permuteA_HE_Free_plot with:
#     input:
#         P = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/P.gz',
#         out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/he.free.npy',
#     output:
#         h2 = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/he.free.h2.png',


use rule yazar_HE_Free_VW_plot as yazar_nomissing_permuteA_HE_Free_VW_plot with:
    input:
        P = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/P.gz',
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/he.free.npy',
    output:
        png = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/he.free.VW.png',


rule yazar_nomissing_permuteA_he_free_VW_topvsbot_plot:
    input:
        P = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/P.gz',
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/he.free.npy',
        ctp = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/ctp.gz',
    output:
        png = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/he.free.VW.topvsbot.png',
    params:
        order = yazar_ct_order,
        colors = yazar_colors,
    script: 'scripts/yazar/free.VW.plot.topvsbot.py'


# # a2. no cis genetic effect (break ctp-ctnu connection)
# rule yazar_nomissing_permuteA2_selectCT:
#     input:
#         ctp = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/ctp.std.gz',
#         ctnu = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/ctnu.std.gz',
#     output:
#         ctp = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis2/ctp.gz',
#         ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis2/ctnu.gz',
#         P = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis2/P.gz',
#     params:
#         ct = target_ct,
#         seed = 1989,
#     run:
#         rng = np.random.default_rng(params.seed)

#         ctp = pd.read_table(input.ctp)
#         ctnu = pd.read_table(input.ctnu)

#         # choose ct
#         ctp = ctp.loc[ctp['ct'] == params.ct]
#         ctnu = ctnu.loc[ctnu['ct'] == params.ct]

#         # santiy check
#         if np.any(ctp['ind'] != ctnu['ind']):
#             sys.exit('Individual order not matching!')

#         # copy
#         ctpA = ctp.copy()
#         ctpA['ct'] = ctpA['ct'] + ':A'
#         ctnuA = ctnu.copy()
#         ctnuA['ct'] = ctnuA['ct'] + ':A'
#         ctpB = ctp.copy()
#         ctpB['ct'] = ctpB['ct'] + ':B'
#         ctnuB = ctnu.copy()
#         ctnuB['ct'] = ctnuB['ct'] + ':B'

#         # permute
#         A_inds = ctp['ind'].sample(frac=1, random_state=rng.integers(10000)).to_numpy()
#         ctpA['ind'] = A_inds
#         A_inds = ctnu['ind'].sample(frac=1, random_state=rng.integers(10000)).to_numpy()
#         ctnuA['ind'] = A_inds
#         B_inds = ctp['ind'].sample(frac=1, random_state=rng.integers(10000)).to_numpy()
#         ctpB['ind'] = B_inds
#         B_inds = ctnu['ind'].sample(frac=1, random_state=rng.integers(10000)).to_numpy()
#         ctnuB['ind'] = B_inds

#         # merge
#         ctp = pd.concat([ctpA, ctpB], ignore_index=True)
#         ctnu = pd.concat([ctnuA, ctnuB], ignore_index=True)

#         # sort
#         ctp = ctp.sort_values(by=['ind', 'ct'])
#         ctnu = ctnu.sort_values(by=['ind', 'ct'])

#         ctp.to_csv(output.ctp, sep='\t', index=False)
#         ctnu.to_csv(output.ctnu, sep='\t', index=False)

#         # make P of 1 / # cts
#         p = 1 / len(np.unique(ctp['ct']))
#         inds = np.unique(ctp['ind'])
#         P = pd.DataFrame({'ind': inds})
#         for ct in np.unique(ctp['ct']):
#             P[ct] = p
#         P.to_csv(output.P, sep='\t', index=False)


# use rule yazar_HE_split as yazar_nomissing_permuteA2_HE_split with:
#     input:
#         ctp = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis2/ctp.gz',
#         ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis2/ctnu.gz',
#     output:
#         ctp = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis2/he/ctp.batch{i}.gz'
#                 for i in range(yazar_he_batches)],
#         ctnu = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis2/he/ctnu.batch{i}.gz'
#                 for i in range(yazar_he_batches)],


# use rule yazar_nomissing_permuteA_HE_free as yazar_nomissing_permuteA2_HE_free with:
#     input:
#         ctp = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis2/he/ctp.batch{{i}}.gz',
#         ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis2/he/ctnu.batch{{i}}.gz',
#         kinship = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/kinship.batch{{i}}.npy',
#         P = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis2/P.gz',
#     output:
#         out = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis2/he.free.batch{{i}}.npy',


# use rule yazar_HE_free_merge as yazar_nomissing_permuteA2_HE_free_merge with:
#     input:
#         out = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis2/he.free.batch{i}.npy'
#             for i in range(yazar_he_batches)],
#     output:
#         out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis2/he.free.npy',


# use rule yazar_HE_Free_VW_plot as yazar_nomissing_permuteA2_HE_Free_VW_plot with:
#     input:
#         P = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis2/P.gz',
#         out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis2/he.free.npy',
#     output:
#         png = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis2/he.free.VW.png',


# b. no genetic correlation
rule yazar_nomissing_permuteB_selectCT:
    input:
        ctp = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/ctp.std.gz',
        ctnu = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/ctnu.std.gz',
    output:
        ctp = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocor/ctp.gz',
        ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocor/ctnu.gz',
        P = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocor/P.gz',
    params:
        ct = target_ct,
        seed = 1983,
    run:
        rng = np.random.default_rng(params.seed)

        ctp = pd.read_table(input.ctp)
        ctnu = pd.read_table(input.ctnu)

        # choose ct
        ctp = ctp.loc[ctp['ct'] == params.ct]
        ctnu = ctnu.loc[ctnu['ct'] == params.ct]

        # santiy check
        if np.any(ctp['ind'] != ctnu['ind']):
            sys.exit('Individual order not matching!')

        # copy
        ctpA = ctp.copy()
        ctpA['ct'] = ctpA['ct'] + ':A'
        ctnuA = ctnu.copy()
        ctnuA['ct'] = ctnuA['ct'] + ':A'
        ctpB = ctp.copy()
        ctpB['ct'] = ctpB['ct'] + ':B'
        ctnuB = ctnu.copy()
        ctnuB['ct'] = ctnuB['ct'] + ':B'

        # permute
        B_inds = ctp['ind'].sample(frac=1, random_state=rng.integers(10000)).to_numpy()
        ctpB['ind'] = B_inds
        ctnuB['ind'] = B_inds

        # merge
        ctp = pd.concat([ctpA, ctpB], ignore_index=True)
        ctnu = pd.concat([ctnuA, ctnuB], ignore_index=True)

        # sort
        ctp = ctp.sort_values(by=['ind', 'ct'])
        ctnu = ctnu.sort_values(by=['ind', 'ct'])

        ctp.to_csv(output.ctp, sep='\t', index=False)
        ctnu.to_csv(output.ctnu, sep='\t', index=False)

        # make P of 1 / # cts
        p = 1 / len(np.unique(ctp['ct']))
        inds = np.unique(ctp['ind'])
        P = pd.DataFrame({'ind': inds})
        for ct in np.unique(ctp['ct']):
            P[ct] = p
        P.to_csv(output.P, sep='\t', index=False)
 

use rule yazar_HE_split as yazar_nomissing_permuteB_HE_split with:
    input:
        ctp = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocor/ctp.gz',
        ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocor/ctnu.gz',
    output:
        ctp = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocor/he/ctp.batch{i}.gz'
                for i in range(yazar_he_batches)],
        ctnu = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocor/he/ctnu.batch{i}.gz'
                for i in range(yazar_he_batches)],


use rule yazar_nomissing_permuteA_HE_free as yazar_nomissing_permuteB_HE_free with:
    input:
        ctp = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocor/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocor/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/kinship.batch{{i}}.npy',
        P = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocor/P.gz',
    output:
        out = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocor/he.free.batch{{i}}.npy',


use rule yazar_HE_free_merge as yazar_nomissing_permuteB_HE_free_merge with:
    input:
        out = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocor/he.free.batch{i}.npy'
            for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocor/he.free.npy',


# use rule yazar_HE_Free_plot as yazar_nomissing_permuteB_HE_Free_plot with:
#     input:
#         P = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocor/P.gz',
#         out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocor/he.free.npy',
#     output:
#         h2 = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocor/he.free.h2.png',


use rule yazar_HE_Free_VW_plot as yazar_nomissing_permuteB_HE_Free_VW_plot with:
    input:
        P = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocor/P.gz',
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocor/he.free.npy',
    output:
        png = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocor/he.free.VW.png',


# 0. hom genetic
rule yazar_nomissing_permute0_selectCT:
    input:
        ctp = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/ctp.std.gz',
        ctnu = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/ctnu.std.gz',
    output:
        ctp = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/hom/ctp.gz',
        ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/hom/ctnu.gz',
        P = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/hom/P.gz',
    params:
        ct = target_ct,
        seed = 1983,
    run:
        rng = np.random.default_rng(params.seed)

        ctp = pd.read_table(input.ctp)
        ctnu = pd.read_table(input.ctnu)

        # choose ct
        ctp = ctp.loc[ctp['ct'] == params.ct]
        ctnu = ctnu.loc[ctnu['ct'] == params.ct]

        # santiy check
        if np.any(ctp['ind'] != ctnu['ind']):
            sys.exit('Individual order not matching!')

        # copy
        ctpA = ctp.copy()
        ctpA['ct'] = ctpA['ct'] + ':A'
        ctnuA = ctnu.copy()
        ctnuA['ct'] = ctnuA['ct'] + ':A'
        ctpB = ctp.copy()
        ctpB['ct'] = ctpB['ct'] + ':B'
        ctnuB = ctnu.copy()
        ctnuB['ct'] = ctnuB['ct'] + ':B'

        # merge
        ctp = pd.concat([ctpA, ctpB], ignore_index=True)
        ctnu = pd.concat([ctnuA, ctnuB], ignore_index=True)

        # sort
        ctp = ctp.sort_values(by=['ind', 'ct'])
        ctnu = ctnu.sort_values(by=['ind', 'ct'])

        ctp.to_csv(output.ctp, sep='\t', index=False)
        ctnu.to_csv(output.ctnu, sep='\t', index=False)

        # make P of 1 / # cts
        p = 1 / len(np.unique(ctp['ct']))
        inds = np.unique(ctp['ind'])
        P = pd.DataFrame({'ind': inds})
        for ct in np.unique(ctp['ct']):
            P[ct] = p
        P.to_csv(output.P, sep='\t', index=False)
 

use rule yazar_HE_split as yazar_nomissing_permute0_HE_split with:
    input:
        ctp = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/hom/ctp.gz',
        ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/hom/ctnu.gz',
    output:
        ctp = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/hom/he/ctp.batch{i}.gz'
                for i in range(yazar_he_batches)],
        ctnu = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/hom/he/ctnu.batch{i}.gz'
                for i in range(yazar_he_batches)],


use rule yazar_nomissing_permuteA_HE_free as yazar_nomissing_permute0_HE_free with:
    input:
        ctp = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/hom/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/hom/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/kinship.batch{{i}}.npy',
        P = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/hom/P.gz',
    output:
        out = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/hom/he.free.batch{{i}}.npy',


use rule yazar_HE_free_merge as yazar_nomissing_permute0_HE_free_merge with:
    input:
        out = [f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/hom/he.free.batch{i}.npy'
            for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/hom/he.free.npy',


use rule yazar_HE_Free_VW_plot as yazar_nomissing_permute0_HE_Free_VW_plot with:
    input:
        P = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/hom/P.gz',
        out = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/hom/he.free.npy',
    output:
        png = f'results/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/hom/he.free.VW.png',


rule yazar_nomissing_permute_all:
    input:
        hom_VW = expand('results/yazar/nomissing/{params}/permute/hom/he.free.VW.png', params=yazar_paramspace.instance_patterns),
        nocis_VW = expand('results/yazar/nomissing/{params}/permute/nocis/he.free.VW.png', params=yazar_paramspace.instance_patterns),
        nocor_VW = expand('results/yazar/nomissing/{params}/permute/nocor/he.free.VW.png', params=yazar_paramspace.instance_patterns),


##################################################################################################################
#  1.2: remove missing individuals   
#     .3: choose one ct and copy and permute to simulate 
#       .4: using only individuals with >0 expression
#         a. no cis genetic effect 
##################################################################################################################
# # a. no cis genetic effect
# rule yazar_nomissing_nozero_permuteA_HE_free:
#     input:
#         raw_ctp = f'analysis/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/ctp.gz',
#         ctp = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/he/ctp.batch{{i}}.gz',
#         ctnu = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/he/ctnu.batch{{i}}.gz',
#         kinship = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/he/kinship.batch{{i}}.npy',
#         P = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/P.gz',
#     output:
#         out = f'staging/yazar/nomissing/nozero/{yazar_paramspace.wildcard_pattern}/permute/nocis/he.free.batch{{i}}.npy',
#     params:
#         snps = 5, # threshold of snp number per gene
#     resources:
#         mem_mb = '30G',
#     script: 'scripts/yazar/nomissing.nozero.permute.he.free.py' 


# # use rule yazar_HE_free_merge as yazar_nomissing_nozero_permuteA_HE_free_merge with:
# rule yazar_nomissing_nozero_permuteA_HE_free_merge:
#     input:
#         out = [f'staging/yazar/nomissing/nozero/{yazar_paramspace.wildcard_pattern}/permute/nocis/he.free.batch{i}.npy'
#             for i in range(yazar_he_batches)],
#     output:
#         out = f'analysis/yazar/nomissing/nozero/{yazar_paramspace.wildcard_pattern}/permute/nocis/he.free.npy',
#     script: 'scripts/mergeBatches.py'


# use rule yazar_HE_Free_VW_plot as yazar_nomissing_nozero_permuteA_HE_Free_VW_plot with:
#     input:
#         P = f'staging/yazar/nomissing/{yazar_paramspace.wildcard_pattern}/permute/nocis/P.gz',
#         out = f'analysis/yazar/nomissing/nozero/{yazar_paramspace.wildcard_pattern}/permute/nocis/he.free.npy',
#     output:
#         png = f'results/yazar/nomissing/nozero/{yazar_paramspace.wildcard_pattern}/permute/nocis/he.free.VW.png',













































# ##########################################################################
# ## 1.2: raw counts
# ##########################################################################


# rule yazar_ctp_extractX_raw:
#     input:
#         h5ad = 'data/Yazar2022Science/OneK1K_cohort_gene_expression_matrix_14_celltypes.h5ad.gz',
#         var = 'data/Yazar2022Science/var.txt',
#         obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
#     output:
#         X = 'staging/data/yazar/raw/X.npz',
#         obs = 'staging/data/yazar/raw/obs.gz',
#         var = 'staging/data/yazar/raw/var.gz',
#     params:
#         ind_col = yazar_ind_col,
#         ct_col = yazar_ct_col,
#     resources:
#         mem_mb = '40G',
#     run:
#         import scanpy as sc
#         from scipy import sparse

#         genes = pd.read_table(input.var)
#         if 'feature_is_filtered' in genes.columns:
#             genes = genes.loc[~genes['feature_is_filtered'], 'feature'].to_numpy()
#         else:
#             genes = genes['feature'].to_numpy()

#         if 'subset_gene' in params.keys():
#             # random select genes
#             rng = np.random.default_rng(seed=params.seed)
#             genes = rng.choice(genes, params.subset_gene, replace=False)

#         obs = pd.read_table(input.obs)
#         ind_pool = np.unique(obs[params.ind_col].astype('str')+'+'+obs['pool'].astype('str'))

#         ann = sc.read_h5ad(input.h5ad, backed='r')
#         data = ann[(~ann.obs[params.ind_col].isna())
#                 & (~ann.obs[params.ct_col].isna())
#                 & (ann.obs[params.ind_col].astype('str')+'+'+ann.obs['pool'].astype('str')).isin(ind_pool), genes]
#         # raw without normalize and logarithm
#         X = data.X
#         sparse.save_npz(output.X, X)

#         data.obs.rename_axis('cell').to_csv(output.obs, sep='\t')
#         data.var.rename_axis('feature').to_csv(output.var, sep='\t')


# use rule yazar_ctp as yazar_raw_ctp with:
#     input:
#         X = 'staging/data/yazar/raw/X.npz',
#         obs = 'staging/data/yazar/raw/obs.gz',
#         var = 'staging/data/yazar/raw/var.gz',
#     output:
#         ctp = 'data/Yazar2022Science/raw/ctp.gz',
#         ctnu = 'data/Yazar2022Science/raw/ctnu.gz',
#         P = 'data/Yazar2022Science/raw/P.gz',
#         n = 'data/Yazar2022Science/raw/n.gz',


# use rule yazar_rm_rareINDnCT as yazar_raw_rm_rareINDnCT with:
#     # also select gene expressed in all cts
#     input:
#         ctp = 'data/Yazar2022Science/raw/ctp.gz',
#         ctnu = 'data/Yazar2022Science/raw/ctnu.gz',
#         n = 'data/Yazar2022Science/raw/n.gz',
#     output:
#         ctp = f'staging/yazar/raw/{yazar_paramspace.wildcard_pattern}/ctp.gz',
#         ctnu = f'staging/yazar/raw/{yazar_paramspace.wildcard_pattern}/ctnu.gz',
#         P = f'analysis/yazar/raw/{yazar_paramspace.wildcard_pattern}/P.gz',
#         n = f'analysis/yazar/raw/{yazar_paramspace.wildcard_pattern}/n.gz',


# use rule yazar_filter_genes as yazar_raw_filter_genes with:
#     input:
#         ctp = f'staging/yazar/raw/{yazar_paramspace.wildcard_pattern}/ctp.gz',
#     output:
#         genes = f'staging/yazar/raw/{yazar_paramspace.wildcard_pattern}/genes.{{cut}}.txt',


# use rule yazar_mvn_ctp as yazar_raw_mvn_ctp with:
#     input:
#         data = f'staging/yazar/raw/{yazar_paramspace.wildcard_pattern}/ctp.gz',
#     output:
#         data = f'staging/yazar/raw/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',


# use rule yazar_mvn_ctp as yazar_raw_mvn_ctnu with:
#     input:
#         data = f'staging/yazar/raw/{yazar_paramspace.wildcard_pattern}/ctnu.gz',
#     output:
#         data = f'staging/yazar/raw/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',


# use rule yazar_std_op as yazar_raw_std_op with:
#     input:
#         ctp = f'staging/yazar/raw/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',
#         ctnu = f'staging/yazar/raw/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',
#         P = f'analysis/yazar/raw/{yazar_paramspace.wildcard_pattern}/P.gz',
#     output:
#         op = f'analysis/yazar/raw/{yazar_paramspace.wildcard_pattern}/op.std.gz',
#         nu = f'analysis/yazar/raw/{yazar_paramspace.wildcard_pattern}/nu.std.gz',
#         ctp = f'analysis/yazar/raw/{yazar_paramspace.wildcard_pattern}/ctp.std.gz',
#         ctnu = f'analysis/yazar/raw/{yazar_paramspace.wildcard_pattern}/ctnu.std.gz',


# use rule yazar_op_pca as yazar_raw_op_pca with:
#     input:
#         op = f'analysis/yazar/raw/{yazar_paramspace.wildcard_pattern}/op.std.gz',
#     output:
#         evec = f'analysis/yazar/raw/{yazar_paramspace.wildcard_pattern}/evec.txt',
#         eval = f'analysis/yazar/raw/{yazar_paramspace.wildcard_pattern}/eval.txt',
#         pca = f'analysis/yazar/raw/{yazar_paramspace.wildcard_pattern}/pca.txt',
#         png = f'analysis/yazar/raw/{yazar_paramspace.wildcard_pattern}/pca.png',


# use rule yazar_he_kinship as yazar_raw_he_kinship with:
#     input:
#         genes = 'data/Yazar2022Science/gene_loation.txt',
#         P = f'analysis/yazar/raw/{yazar_paramspace.wildcard_pattern}/P.gz',
#         bed = 'analysis/yazar/data/geno/chr{chr}.bed',
#     output:
#         kinship = temp(f'staging/yazar/raw/{yazar_paramspace.wildcard_pattern}/he/kinship.chr{{chr}}.txt'),


# use rule yazar_he_kinship_merge as yazar_raw_he_kinship_merge with:
#     input:
#         kinship = [f'staging/yazar/raw/{yazar_paramspace.wildcard_pattern}/he/kinship.chr{chr}.txt'
#                 for chr in range(1,23)],
#     output:
#         kinship = f'staging/yazar/raw/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
#         save = f'analysis/yazar/raw/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',


# use rule yazar_HE_split as yazar_raw_HE_split with:
#     input:
#         ctp = f'analysis/yazar/raw/{yazar_paramspace.wildcard_pattern}/ctp.std.gz',
#         ctnu = f'analysis/yazar/raw/{yazar_paramspace.wildcard_pattern}/ctnu.std.gz',
#     output:
#         ctp = [f'staging/yazar/raw/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{i}.gz'
#                 for i in range(yazar_he_batches)],
#         ctnu = [f'staging/yazar/raw/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{i}.gz'
#                 for i in range(yazar_he_batches)],


# use rule yazar_HE_free as yazar_raw_HE_free with:
#     input:
#         ctp = f'staging/yazar/raw/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
#         ctnu = f'staging/yazar/raw/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
#         kinship = f'staging/yazar/raw/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
#         P = f'analysis/yazar/raw/{yazar_paramspace.wildcard_pattern}/P.gz',
#         op_pca = f'analysis/yazar/raw/{yazar_paramspace.wildcard_pattern}/pca.txt',
#         geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
#         obs = 'staging/data/yazar/raw/obs.gz',
#     output:
#         out = f'staging/yazar/raw/{yazar_paramspace.wildcard_pattern}/he.free.batch{{i}}.npy',


# use rule yazar_HE_free_merge as yazar_raw_HE_free_merge with:
#     input:
#         out = [f'staging/yazar/raw/{yazar_paramspace.wildcard_pattern}/he.free.batch{i}.npy'
#             for i in range(yazar_he_batches)],
#     output:
#         out = f'analysis/yazar/raw/{yazar_paramspace.wildcard_pattern}/he.free.npy',


# use rule yazar_HE_Free_plot as yazar_raw_HE_Free_plot with:
#     input:
#         P = f'analysis/yazar/raw/{yazar_paramspace.wildcard_pattern}/P.gz',
#         out = f'analysis/yazar/raw/{yazar_paramspace.wildcard_pattern}/he.free.npy',
#     output:
#         h2 = f'results/yazar/raw/{yazar_paramspace.wildcard_pattern}/he.free.h2.png',


# rule yazar_raw_all:
#     input:
#         free = expand('results/yazar/raw/{params}/he.free.h2.png', params=yazar_paramspace.instance_patterns),




# ##################### 1.2 subset cts ###############################
# subset_cts = {
#         #'B': ['B IN', 'B Mem'],
#         'CD4': ['CD4 ET', 'CD4 NC'],
#         'CD8': ['CD8 ET', 'CD8 NC', 'CD8 S100B'],
#         }

# wildcard_constraints: sc='[\w\d]+'
# rule yazar_subset_cts:
#     # TODO: ctp is stded on OP of all cts
#     input:
#         ctp = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',
#         ctnu = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',
#         P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
#     output:
#         ctp = f'staging/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',
#         ctnu = f'staging/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',
#         P = f'staging/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
#     params:
#         subset_cts = lambda wildcards: subset_cts[wildcards.sc],
#     run:
#         ctp = pd.read_table(input.ctp)
#         ctp = ctp.loc[ctp['ct'].isin(params.subset_cts)]
#         ctp.to_csv(output.ctp, sep='\t', index=False)
#         ctnu = pd.read_table(input.ctnu)
#         ctnu.loc[ctnu['ct'].isin(params.subset_cts)].to_csv(output.ctnu, sep='\t', index=False)
#         P = pd.read_table(input.P, index_col=0)
#         cts = np.unique(params.subset_cts)
#         cts = cts[np.isin(cts, P.columns)]
#         P = P.loc[np.unique(ctp['ind']), cts]
#         P = P.div(P.sum(axis=1), axis=0)
#         P.to_csv(output.P, sep='\t')

# use rule yazar_HE_split as yazar_subset_cts_HE_split with:
#     input:
#         ctp = f'staging/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/ctp.mvn.gz',
#         ctnu = f'staging/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/ctnu.mvn.gz',
#     output:
#         ctp = [f'staging/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{i}.gz'
#                 for i in range(yazar_he_batches)],
#         ctnu = [f'staging/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{i}.gz'
#                 for i in range(yazar_he_batches)],

# use rule yazar_HE_free as yazar_subset_cts_HE_free with:
#     input:
#         ctp = f'staging/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
#         ctnu = f'staging/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
#         kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
#         P = f'staging/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
#         op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.txt',
#         geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
#         obs = 'staging/data/yazar/obs.gz',
#     output:
#         out = f'staging/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/he.free.batch{{i}}.npy',
#     params:
#         snps = 5, # threshold of snp number per gene

# use rule yazar_HE_free_merge as yazar_subset_cts_HE_free_merge with:
#     input:
#         out = [f'staging/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/he.free.batch{i}.npy'
#                 for i in range(yazar_he_batches)],
#     output:
#         out = f'analysis/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/he.free.npy',


# use rule yazar_HE_Free_plot as yazar_subset_cts_HE_Free_plot with:
#     input:
#         P = f'staging/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
#         out = f'analysis/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/he.free.npy',
#     output:
#         h2 = f'results/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/he.free.h2.png',


# rule yazar_subset_cts_HE_Free_subsetGenes:
#     input:
#         out = f'analysis/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/he.free.npy',
#         genes = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/genes.{{cut}}.txt',
#     output:
#         out = f'staging/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/he.{{cut}}.free.txt',
#     run:
#         from gxctmm import util
#         out = np.load(input.out, allow_pickle=True).item()
#         hom_g2 = out['free']['hom_g2']
#         V = out['free']['V']
#         hom_e2 = out['free']['hom_e2']
#         W = out['free']['W']
#         h2 = util.compute_h2(hom_g2, V, hom_e2, W)
#         genes = np.loadtxt(input.genes, dtype='str')
#         selected = np.isin(out['gene'], genes)
#         h2 = h2[selected]
#         mean = "\t".join(np.mean(h2, axis=0).astype('str'))
#         median = "\t".join(np.median(h2, axis=0).astype('str'))
#         with open(output.out, 'w') as f:
#             f.write(f'mean:{mean}\nmedian:{median}\n')


# use rule yazar_REML_free as yazar_subset_cts_REML_free with:
#     input:
#         ctp = f'staging/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
#         ctnu = f'staging/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
#         kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
#         P = f'staging/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
#         op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.txt',
#         geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
#         obs = 'staging/data/yazar/obs.gz',
#     output:
#         out = f'staging/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/reml.free.batch{{i}}.npy',
#     params:
#         snps = 5, # threshold of snp number per gene


# use rule yazar_HE_free_merge as yazar_subset_cts_REML_free_merge with:
#     input:
#         out = [f'staging/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/reml.free.batch{i}.npy'
#                 for i in range(yazar_he_batches)],
#     output:
#         out = f'analysis/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/reml.free.npy',


# use rule yazar_HE_Free_plot as yazar_subset_cts_REML_Free_plot with:
#     input:
#         P = f'staging/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
#         out = f'analysis/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/reml.free.npy',
#     output:
#         h2 = f'results/yazar/subset_cts/{{sc}}/{yazar_paramspace.wildcard_pattern}/reml.free.h2.png',


# rule yazar_subset_cts_all:
#     input:
#         he = expand('results/yazar/subset_cts/{sc}/{params}/he.free.h2.png',
#                 sc=list(subset_cts.keys()),
#                 params=yazar_paramspace.instance_patterns),


# ##################### 1.2 get rid of ctnu ##########################
# rule yazar_ccNOctnu_HE_free:
#     input:
#         ctp = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
#         ctnu = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
#         kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
#         P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
#         op_pca = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/pca.txt',
#         geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
#         obs = 'staging/data/yazar/combine_cts/{cc}/obs.gz',
#     output:
#         out = f'staging/yazar/combine_cts/{{cc}}/noctnu/{yazar_paramspace.wildcard_pattern}/he.free.batch{{i}}.npy',
#     params:
#         snps = 5, # threshold of snp number per gene
#     script: 'scripts/yazar/he.free.ccNOctnu.py'


# use rule yazar_HE_free_merge as yazar_ccNOctnu_HE_free_merge with:
#     input:
#         out = [f'staging/yazar/combine_cts/{{cc}}/noctnu/{yazar_paramspace.wildcard_pattern}/he.free.batch{i}.npy'
#                 for i in range(yazar_he_batches)],
#     output:
#         out = f'staging/yazar/combine_cts/{{cc}}/noctnu/{yazar_paramspace.wildcard_pattern}/he.free.npy',


# use rule yazar_HE_Free_plot as yazar_ccNOctnu_HE_Free_plot with:
#     input:
#         P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
#         out = f'staging/yazar/combine_cts/{{cc}}/noctnu/{yazar_paramspace.wildcard_pattern}/he.free.npy',
#     output:
#         h2 = f'results/yazar/combine_cts/{{cc}}/noctnu/{yazar_paramspace.wildcard_pattern}/he.free.h2.png',


# rule yazar_ccNOctnu_all:
#     input:
#         h2 = expand('results/yazar/combine_cts/{cc}/noctnu/{params}/he.free.h2.png', 
#                 cc=list(combine_cts.keys()), 
#                 params=yazar_paramspace.instance_patterns),


# ########### 1.2 correct for P cell type proportions (not implemented) ################
# #rule yazar_HE_free_adjP:
# #    input:
# #        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
# #        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
# #        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.txt',
# #        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.gz',
# #        op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.txt',
# #        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
# #        obs = 'staging/data/yazar/obs.gz',
# #    output:
# #        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he.free.adjP.batch{{i}}',
# #    params:
# #        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/rep/he.free.adjP.npy',
# #        snps = 5, # threshold of snp number per gene
# #    resources:
# #        time = '10:00:00',
# #        mem_mb = '20G',
# #    script: 'scripts/yazar/he.free.adjP.py'
# #
# #use rule yazar_HE_free_merge as yazar_HE_free_merge_adjP with:
# #    input:
# #        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he.free.adjP.batch{i}'
# #            for i in range(yazar_he_batches)],
# #    output:
# #        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/he.free.adjP.npy',


