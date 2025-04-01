########################################################################
# Yazar 2022 Science
########################################################################
yazar_ind_col = config['yazar']['yazar_ind_col']
yazar_ct_col = config['yazar']['yazar_ct_col']
yazar_pool_col = config['yazar']['yazar_pool_col']
yazar_he_batches = config['yazar']['he_nbatch']
yazar_reml_batches = config['yazar']['reml_nbatch']

yazar_ct_order = np.array(['hom', 'CD4 NC', 'CD8 ET', 'NK', 'CD8 NC', 'B IN', 'CD4 ET', 'B Mem', 'Mono C', 'CD8 S100B', 'Mono NC', 'NK R', 'DC', 'CD4 SOX4', 'Plasma'])
yazar_colors = dict(zip(yazar_ct_order[1:], sns.color_palette(config['colorpalette'])))
yazar_colors['hom'] = '0.7'
# paried light colors
colors = [colorsys.rgb_to_hsv(*color) for color in sns.color_palette(config['colorpalette'])]
colors = [(color[0], min(1, max(0, color[1] - 0.3)), color[2]) for color in colors]
colors = [colorsys.hsv_to_rgb(*color) for color in colors]
yazar_light_colors = dict(zip(yazar_ct_order[1:], colors))
yazar_light_colors['hom'] = '0.9'

# read parameters
yazar_params = pd.read_table('yazar.params.txt', dtype="str", comment='#')
yazar_paramspace = Paramspace(yazar_params, filename_params="*")



rule get_biomart_v75:
    # only keep autosomes
    output:
        gene_info = 'analysis/data/Homo_sapiens.GRCh37.75.txt',
    script: '../bin/utils/get_biomart_v75.R'


# data preprocessing
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
        for column in var.columns:
            var[column] = var[column].str.replace(' ', '')
        var.to_csv(output.var, sep='\t', index=False)


rule yazar_exclude_repeatedpool:
    input:
        obs = 'data/Yazar2022Science/obs.txt',
    output:
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
        dup_inds = 'analysis/yazar/duplicated_inds.txt',
    params:
        ind_col = yazar_ind_col,
        pool_col = yazar_pool_col,
    run:
        obs = pd.read_table(input.obs)
        # id repeated pool: the same individual seqed in more than one pool
        data = obs[[params.pool_col, params.ind_col]].drop_duplicates()
        inds, counts = np.unique(data[params.ind_col], return_counts=True)
        inds = inds[counts > 1]
        np.savetxt(output.dup_inds, inds, fmt='%s')

        # for each ind find the largest pool
        for ind in inds:
            pools, counts = np.unique(obs.loc[obs[params.ind_col]==ind, params.pool_col], return_counts=True)
            excluded_pools = pools[counts < np.amax(counts)]
            obs = obs.loc[~((obs[params.ind_col]==ind) & (obs[params.pool_col].isin(excluded_pools)))]

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


rule yazar_ctp_transform:
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
        pool_col = yazar_pool_col,
    resources:
        mem_mb = '120G',
    script: '../bin/yazar/transform.py'


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
        mem_mb = '120G',
    run:
        from scipy import sparse
        from cigma import preprocess

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


rule yazar_P_plot:
    input:
        P = 'data/Yazar2022Science/P.gz',
    output:
        png = 'results/yazar/P.png',
    run:
        P = pd.read_table(input.P, index_col=0)
        P = P.drop(['Erythrocytes', 'Platelets'], axis=1)
        P = P.div(P.sum(axis=1), axis=0)
        P = P[P.mean().sort_values(ascending=False).index]
        P.columns = P.columns.str.replace(' ', '_')

        plt.rcParams.update({'font.size' : 6})
        fig, ax = plt.subplots(figsize=(8,4), dpi=600)
        sns.violinplot(data=P, scale='width', cut=0)
        ax.axhline(y=0, color='0.9', ls='--', zorder=0)
        ax.set_xlabel('Cell type', fontsize=10)
        ax.set_ylabel('Cell type proportion', fontsize=10)
        plt.tight_layout()
        fig.savefig(output.png)


use rule yazar_ctp_transform as yazar_ctp_transform_var_ctnu with:
    output:
        X = 'staging/data/yazar/var_ctnu/X.npz',
        obs = 'staging/data/yazar/var_ctnu/obs.gz',
        var = 'staging/data/yazar/var_ctnu/var.gz',


var_ctnu_batches = 500
rule yazar_var_ctnu_split:
    input:
        obs = 'staging/data/yazar/var_ctnu/obs.gz',
    output:
        batches = expand('staging/data/yazar/var_ctnu/ind_ct.batch{i}', i=range(var_ctnu_batches)),
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
        X = 'staging/data/yazar/var_ctnu/X.npz',
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
        mem_mb = '50G',
    run:
        from scipy import stats, sparse
        from sklearn.utils import resample

        def cal_ctnu(data):
            ctp = np.squeeze( np.asarray(data.mean(axis=0)) )
            ctp2 = np.squeeze( np.asarray(data.power(2).mean(axis=0)) )
            ctnu = (ctp2 - ctp**2) / data.shape[0]
            return ctnu 

        def bootstrap(data, rng, n_resamples=10000):
            ctnus = []
            for i in range(n_resamples):
                sample = resample(data, random_state=rng)
                ctnus.append( cal_ctnu(sample) )
            return ctnus 

        X = sparse.load_npz(input.X).tocsr()
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
        print(np.array(boots['var_ctnu']).shape)
        var_ctnu = pd.DataFrame(data=boots['var_ctnu'], columns=genes)
        var_ctnu.insert(loc=0, column='ct', value=boots['ct'])
        var_ctnu.insert(loc=0, column='ind', value=boots['ind'])
        var_ctnu.to_csv(output.var_ctnu, sep='\t', index=False)


rule yazar_var_ctnu_merge:
    input:
        var_ctnu = expand('staging/data/yazar/var_ctnu/batch{i}.gz', i=range(var_ctnu_batches)),
    output:
        var_ctnu = 'analysis/yazar/var_ctnu.gz',
    run:
        var_ctnu = [pd.read_table(f) for f in input.var_ctnu]
        pd.concat(var_ctnu, ignore_index=True).to_csv(output.var_ctnu,
                sep='\t', index=False)


rule yazar_varNU_dist:
    input:
        var_nu = 'analysis/yazar/var_ctnu.gz',
        nu = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctnu.rmid.gz',
    output:
        cv = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/cv.gz',
        png = f'results/yazar/{yazar_paramspace.wildcard_pattern}/var_nu.png',
    script: '../bin/yazar/varNU_dist.py'


rule yazar_gene_location:
    input:
        genes = 'data/Yazar2022Science/var.txt',
        gene_info = 'analysis/data/Homo_sapiens.GRCh37.75.txt',
        # gff = 'data/Homo_sapiens.GRCh37.82.gff3.gz',
    output:
        location = 'data/Yazar2022Science/gene_location.txt',
        meta = 'data/Yazar2022Science/gene_meta.txt',
    script: '../bin/yazar/gene_location.py'


# data check
rule yazar_cell_dist:
    input:
        obs = 'data/Yazar2022Science/obs.txt',
    output:
        png = 'results/yazar/cell.dist.png',
    params:
        ind_col = yazar_ind_col,
        ct_col = yazar_ct_col,
        colors = yazar_colors,
        analyzed_ct_num = 7,  # number of main cell types included for analysis
    script: '../bin/yazar/cell_dist.py'


rule yazar_rm_rareINDnCT:
    # also select gene expressed in all cts
    input:
        ctp = 'data/Yazar2022Science/ctp.gz',
        ctnu = 'data/Yazar2022Science/ctnu.gz',
        n = 'data/Yazar2022Science/n.gz',
    output:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctp.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctnu.gz',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/mvn/P.final.gz',
        n = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/mvn/n.final.gz',
    resources:
        mem_mb = '20G',
    script: '../bin/yazar/rm_rareINDnCT.py'


rule yazar_rm_missingIND:
    input:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctp.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctnu.gz',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/mvn/P.final.gz',
        n = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/mvn/n.final.gz',
    output:
        ctp = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctp.rmid.gz',
        ctnu = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctnu.rmid.gz',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.final.gz',
        n = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/n.final.gz',
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


rule yazar_std_op:
    input:
        ctp = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctp.rmid.gz',
        ctnu = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctnu.rmid.gz',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.final.gz',
    output:
        op = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/op.std.gz',
        nu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/nu.std.gz',
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctp.std.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctnu.std.gz',
    resources:
        mem_mb = '10G',
    run:
        from cigma import preprocess
        ctp = pd.read_table(input.ctp, index_col=(0,1)).astype('float32')
        ctnu = pd.read_table(input.ctnu, index_col=(0,1)).astype('float32')
        P = pd.read_table(input.P, index_col=0)

        op, nu, ctp, ctnu = preprocess.std(ctp, ctnu, P)

        op.to_csv(output.op, sep='\t')
        nu.to_csv(output.nu, sep='\t')
        ctp.to_csv(output.ctp, sep='\t')
        ctnu.to_csv(output.ctnu, sep='\t')


rule yazar_exclude_sexchr:
    input:
        genes = 'data/Yazar2022Science/gene_location.txt',
        op = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/op.std.gz',
        nu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/nu.std.gz',
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctp.std.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ctnu.std.gz',
    output:
        op = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/op.final.gz',
        nu = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/nu.final.gz',
        ctp = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctp.final.gz',
        ctnu = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctnu.final.gz',
    resources:
        mem_mb = '10G',
    run:
        genes = pd.read_table(input.genes)
        genes = np.unique(genes['feature'])

        op = pd.read_table(input.op, index_col=0)
        op[genes[np.isin(genes, op.columns)]].to_csv(output.op, sep='\t')
        nu = pd.read_table(input.nu, index_col=0)
        nu[genes[np.isin(genes, nu.columns)]].to_csv(output.nu, sep='\t')

        ctp = pd.read_table(input.ctp, index_col=(0, 1))
        ctp[genes[np.isin(genes, ctp.columns)]].to_csv(output.ctp, sep='\t')
        ctnu = pd.read_table(input.ctnu, index_col=(0, 1))
        ctnu[genes[np.isin(genes, ctnu.columns)]].to_csv(output.ctnu, sep='\t')


rule yazar_op_pca:
    input:
        op = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/op.final.gz',
    output:
        evec = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/evec.gz',
        eval = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/eval.gz',
        pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.gz',
        png = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/op.pca.png',
    resources:
        mem_mb = '4G',
    script: '../bin/yazar/pca.py'


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
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.final.gz',
    output:
        eigenvec = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        eigenval = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenval',
    params:
        prefix = lambda wildcards, output: os.path.splitext(output.eigenvec)[0],
        tmp_dir = lambda wildcards, output: os.path.splitext(output.eigenvec)[0] + '_tmp',
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
        ax.set_title('OneK1K')
        fig.savefig(output.png)


rule yazar_he_kinship:
    input:
        ctp = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctp.final.gz',
        genes = 'data/Yazar2022Science/gene_location.txt',
        bed = 'analysis/yazar/data/geno/chr{chr}.bed',
    output:
        kinship = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/kinship.chr{{chr}}.npy',
        # kinship = temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/kinship.chr{{chr}}.npy'),
    params:
        r = int(float(config['yazar']['radius'])),
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

        python3 workflow/bin/yazar/kinship.npy.py {input.genes} {input.ctp} {params.r} {input.bed} {wildcards.chr} \
                        {output.kinship} 
        '''


rule yazar_HE_split:
    input:
        ctp = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctp.final.gz',
        ctnu = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctnu.final.gz',
        genes = 'data/Yazar2022Science/gene_location.txt',
    output:
        ctp = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{i}.gz'
                for i in range(yazar_he_batches)],
        ctnu = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{i}.gz'
                for i in range(yazar_he_batches)],
        batch = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctp.batch.txt',
    resources:
        mem_mb = '10G',
    run:
        meta = pd.read_table(input.genes, usecols=['feature', 'chr'])
        meta = meta.drop_duplicates()
        ctp = pd.read_table(input.ctp, index_col=(0,1))
        ctnu = pd.read_table(input.ctnu, index_col=(0,1))

        meta = meta.loc[meta['feature'].isin(ctp.columns)]
        ngene = meta['feature'].nunique()
        # sanity check
        if ctp.shape[1] != ngene or meta.shape[0] != ngene:
            sys.exit('Missing genes')

        chrs, counts = np.unique(meta['chr'], return_counts=True)
        nbatch = len(output.ctp)
        ngene_per_batch = ngene / nbatch
        chr_nbatchs = []
        for chr in chrs:
            chr_genes = meta.loc[meta['chr'] == chr, 'feature']
            chr_nbatch = int(len(chr_genes) // ngene_per_batch)
            chr_nbatchs.append(chr_nbatch)
        chr_nbatchs[0] = nbatch - sum(chr_nbatchs[1:])

        k = 0
        note = open(output.batch, 'w')
        for i, chr in enumerate(chrs):
            chr_genes = meta.loc[meta['chr'] == chr, 'feature']
            chr_nbatch = chr_nbatchs[i]
            if chr_nbatch == 0:
                continue
            chr_batches = np.array_split(chr_genes, chr_nbatch)
                
            chr_ctps = output.ctp[k:(k + chr_nbatch)]
            chr_ctnus = output.ctnu[k:(k + chr_nbatch)]
            for batch, ctp_f, ctnu_f in zip(chr_batches, chr_ctps, chr_ctnus):
                if len(batch) == 0:
                    sys.exit('Empty batch')
                ctp[batch].to_csv(ctp_f, sep='\t')
                ctnu[batch].to_csv(ctnu_f, sep='\t')
                note.write(f'{chr}\t{ctp_f}\n')

            k += chr_nbatch
        note.close()


rule yazar_he_kinship_split:
    input:
        kinship = [f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/kinship.chr{chr}.npy'
                    for chr in chrs],
        ctp = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{i}.gz'
                for i in range(yazar_he_batches)],
        batch = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctp.batch.txt',
    output:
        kinship = [temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.batch{i}.npy')
                    for i in range(yazar_he_batches)],
    params:
        chrs = chrs,
    resources:
        mem_mb = '20G',
    script: '../bin/yazar/kinship.split.py'


rule yazar_HE_free:
    input:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.batch{{i}}.npy',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.final.gz',
        op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.gz',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he.batch{{i}}.npy',
    params:
        snps = 5, # threshold of snp number per gene
        iid = True,
        full = True,
    resources:
        mem_mb = '80G',
    script: '../bin/yazar/he.py' 


rule yazar_HE_free_merge:
    input:
        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he.batch{i}.npy'
            for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/he.npy',
    script: '../bin/mergeBatches.py'


rule yazar_HE_togz:
    input:
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/he.npy',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.final.gz',
    output:
        out = f'results/yazar/{yazar_paramspace.wildcard_pattern}/he.gz',
    script: '../bin/yazar/he.togz.py'




#################  HE Free with JK  #####################
use rule yazar_HE_split as yazar_reml_split with:
    output:
        ctp = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch{i}.gz'
                for i in range(yazar_reml_batches)],
        ctnu = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctnu.batch{i}.gz'
                for i in range(yazar_reml_batches)],
        batch = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch.txt',


batch1, batch2, batch3, batch4 = np.array_split(np.arange(yazar_reml_batches), 4)
use rule yazar_he_kinship_split as yazar_reml_kinship_split_part1 with:
    input:
        kinship = [f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/kinship.chr{chr}.npy'
                    for chr in chrs],
        ctp = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch{i}.gz'
                for i in batch1],
        batch = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch.txt',
    output:
        # kinship = [temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/kinship.batch{i}.npy')
        kinship = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/kinship.batch{i}.npy'
                    for i in batch1],


use rule yazar_he_kinship_split as yazar_reml_kinship_split_part2 with:
    input:
        kinship = [f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/kinship.chr{chr}.npy'
                    for chr in chrs],
        ctp = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch{i}.gz'
                for i in batch2],
        batch = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch.txt',
    output:
        # kinship = [temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/kinship.batch{i}.npy')
        kinship = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/kinship.batch{i}.npy'
                    for i in batch2],


use rule yazar_he_kinship_split as yazar_reml_kinship_split_part3 with:
    input:
        kinship = [f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/kinship.chr{chr}.npy'
                    for chr in chrs],
        ctp = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch{i}.gz'
                for i in batch3],
        batch = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch.txt',
    output:
        # kinship = [temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/kinship.batch{i}.npy')
        kinship = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/kinship.batch{i}.npy'
                    for i in batch3],


use rule yazar_he_kinship_split as yazar_reml_kinship_split_part4 with:
    input:
        kinship = [f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/kinship.chr{chr}.npy'
                    for chr in chrs],
        ctp = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch{i}.gz'
                for i in batch4],
        batch = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch.txt',
    output:
        # kinship = [temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/kinship.batch{i}.npy')
        kinship = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/kinship.batch{i}.npy'
                    for i in batch4],


use rule yazar_HE_free as yazar_HE_free_jk with:
    input:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/kinship.batch{{i}}.npy',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.final.gz',
        op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.gz',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he.free.batch{{i}}.jk.npy',
    params:
        jk = True,
        snps = 5, # threshold of snp number per gene
        iid = False,
        full = False,
    resources:
        mem_mb = '20G',


use rule yazar_HE_free_merge as yazar_HE_free_jk_merge with:
    input:
        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he.free.batch{i}.jk.npy'
            for i in range(yazar_reml_batches)],
    output:
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/he.free.jk.npy',


use rule yazar_HE_togz as yazar_HE_free_jk_togz with:
    input:
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/he.free.jk.npy',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.final.gz',
    output:
        out = f'results/yazar/{yazar_paramspace.wildcard_pattern}/he.free.jk.gz',










##########################################################################
# correlation with gene annotations
##########################################################################
rule yazar_eds:
    input:
        var = 'data/Yazar2022Science/var.txt', 
        eds = 'data/Wang_Goldstein.tableS1.txt',
        blood_connected = 'data/gene_annots/saha_et_al/twns/processed/WholeBlood.degree.txt',
        combined_connected = 'data/gene_annots/saha_et_al/twns/processed/genes_by_combined_rank.txt',
        gnomad = 'data/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz',
    output:
        eds = 'analysis/yazar/eds.txt',
    run:
        var = pd.read_table(input.var, sep='\s+')
        var = var.rename(columns={'feature': 'gene_id'})
        eds = pd.read_table(input.eds)
        eds = eds.rename(columns={'GeneSymbol': 'gene_id'})
        blood_connected = pd.read_table(input.blood_connected)
        combined_connected = pd.read_table(input.combined_connected)

        print(var.shape, eds.shape)
        print(len(np.unique(var['gene_id'])))
        print(len(np.unique(eds['gene_id'])))

        # NOTE: drop dozens of duplicated genes
        eds = eds.drop_duplicates(subset=['gene_id'], keep=False, ignore_index=True)

        # drop pLI in eds, instead using pLI from gnomad
        eds = eds.drop('pLI', axis=1)

        # read gene length from gnomad
        gnomad = pd.read_table(input.gnomad, usecols=['gene_id', 'gene_length', 'pLI', 'oe_lof_upper'])
        gnomad = gnomad.rename(columns={'oe_lof_upper': 'LOEUF'})
        eds = eds.merge(gnomad, how='outer')

        # add connectedness rank
        print(blood_connected.shape, combined_connected.shape)
        
        blood_connected = blood_connected.merge(var, left_on='gene', right_on='GeneSymbol')
        blood_connected = blood_connected[['gene_id']]
        blood_connected['blood_connected_rank'] = np.arange(blood_connected.shape[0]) + 1

        combined_connected = combined_connected.merge(var, left_on='gene', right_on='GeneSymbol')
        combined_connected = combined_connected[['gene_id']]
        combined_connected['combined_connected_rank'] = np.arange(combined_connected.shape[0]) + 1  # lower rank -> more connection
        
        print(blood_connected.shape, combined_connected.shape)
        eds = eds.merge(blood_connected, how='outer')
        eds = eds.merge(combined_connected, how='outer')
        
        # keep genes in the dataset
        eds = eds.loc[eds['gene_id'].isin(var['gene_id'])]
        eds.to_csv(output.eds, sep='\t', index=False)


# rule yazar_eds_all:
#     input:
#         h2 = expand('results/yazar/{params}/matrix.cor.png', params=yazar_paramspace.instance_patterns),




#######################################################################################
# LDSC
#######################################################################################
gwass = config['yazar']['gwas']
neg_traits = ['Height', 'CAD', 'SCZ'] 
pos_traits = ['UC', 'RA', 'PBC', 'MS', 'Crohns', 'Celiac', 'Lupus']
traits = neg_traits + pos_traits

# gcta
rule yazar_gcta_grm:
    input:
        genes = 'data/Yazar2022Science/gene_location.txt',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.final.gz',
        bed = 'analysis/yazar/data/geno/chr{chr}.bed',
    output:
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/kinship.chr{{chr}}.txt',
    params:
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/kinship/gene.grm.bin',  # file size: each gene has ~1.7 MB bin and ~1.7 MB N.bin
        r = int(float(config['yazar']['radius'])),
    resources:
        mem_mb = '20G',
    shell:
        '''
        # load plink 1.9 and gcta
        module load gcc/11.3.0 atlas/3.10.3 lapack/3.11.0 plink/1.9 gcta/1.94.1

        mkdir -p $(dirname {params.kinship})
        python3 workflow/bin/yazar/kinship.py \
                {input.genes} {input.P} {params.r} \
                {input.bed} {wildcards.chr} \
                {params.kinship} {output.kinship}
        '''

rule yazar_gcta_split_kinship:
    input:
        kinship = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/kinship.chr{chr}.txt'
                for chr in range(1,23)],
    output:
        kinship = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/kinship.batch{i}.txt'
                for i in range(yazar_he_batches)],
    run:
        kinship = [pd.read_table(f) for f in input.kinship]
        kinship = pd.concat(kinship, axis=0, ignore_index=True)
        indexs = np.array_split(kinship.index, len(output.kinship))
        for index, f in zip(indexs, output.kinship):
            kinship.loc[index].to_csv(f, sep='\t', index=False)


rule yazar_gcta_greml_op:
    input:
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.final.gz',
        op = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/op.final.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/kinship.batch{{i}}.txt',
        op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.gz',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/op.greml.batch{{i}}',
    params:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/rep/op.hsq',
        snps = 5, # threshold of snp number per gene
    resources:
        mem_mb = '2G',
    shell:
        '''
        module load gcc/11.3.0 gcta/1.94.1
        python3 workflow/bin/yazar/gcta_greml_op.py \
                {input.P} {input.op} {input.kinship} \
                {input.op_pca} {wildcards.op_pca_n} \
                {input.geno_pca} {wildcards.geno_pca_n} \
                {input.obs} {params.snps} {params.out} {output.out}
        '''


rule yazar_gcta_greml_op_merge:
    input:
        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/gcta/op.greml.batch{i}'
                for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/gcta/op.greml',
    run:
        with open(output.out, 'w') as f:
            f.write('gene\th2\tp\n')
            for chr_out in input.out:
                for gene_out in open(chr_out):
                    gene = re.search(r"/(ENSG[^/]+)/", gene_out).group(1)
                    h2, v_g, v_e, pval = None, None, None, None
                    for line in open(gene_out.strip()):
                        line = line.strip().split()
                        if len(line) == 0:
                            continue
                        elif line[0] == 'V(G)/Vp':
                            h2 = line[1]
                        elif line[0] == 'V(G1)':
                            v_g = float(line[1])
                        elif line[0] == 'V(e)':
                            v_e = float(line[1])
                        elif line[0] == 'Pval':
                            pval = float(line[1])
                    if h2 is None:
                        h2 = v_g / (v_g + v_e)
                    f.write(f'{gene}\t{h2}\t{pval}\n')


rule yazar_ldsc_make_geneset:
    input:
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/he.free.jk.npy',
        location = 'data/Yazar2022Science/gene_location.txt',
        gcta = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/gcta/op.greml',
    output:
        all = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/all.genes.txt',
        random = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/random.genes.txt',
        shared = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/shared.genes.txt',
        var = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/var.genes.txt',
        mean = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/mean.genes.txt',
        gcta = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/gcta.genes.txt',
    params:
        seed = 123,
        chr = config['ldsc']['mhc']['chr'],
        start = config['ldsc']['mhc']['start'],
        end = config['ldsc']['mhc']['end'],
    run:
        out = np.load(input.out, allow_pickle=True).item()
        gene_df = pd.read_table(input.location)

        # remove hla genes: https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37
        gene_df = gene_df.loc[~((gene_df['chr'] == params.chr) & (gene_df['end'] > params.start) & (gene_df['start'] < params.end))]
        filter = np.isin(out['gene'], gene_df['feature'])

        # filter genes
        genes = out['gene'][filter]
        p = out['p']['free']['V'][filter]
        shared_p = out['p']['free']['hom_g2'][filter]
        var_beta = np.var(out['free']['ct_beta'], axis=1)[filter]

        # control
        np.savetxt(output.all, genes, fmt='%s')

        # number of genes to select
        ngene = int(wildcards.ngene)
        if ngene == 0:
            ngene = (p < (0.05 / len(out['gene']))).sum()

        # find top genes
        var_genes = genes[np.argsort(p)][:ngene]
        np.savetxt(output.var, var_genes, fmt='%s')

        # shared
        shared_genes = genes[np.argsort(shared_p)][:ngene]
        np.savetxt(output.shared, shared_genes, fmt='%s')

        # beta
        mean_genes = genes[np.argsort(var_beta)[::-1][:ngene]]
        np.savetxt(output.mean, mean_genes, fmt='%s')

        # random
        rng = np.random.default_rng(seed=params.seed)
        np.savetxt(output.random, rng.choice(genes, size=ngene, replace=False), 
                fmt='%s')

        # gcta
        gcta = pd.read_table(input.gcta, usecols=['gene', 'p'])
        print(gcta.shape)
        gcta = gcta.loc[gcta['gene'].isin(genes)]
        print(gcta.shape)
        gcta = gcta.sort_values('p')
        if gcta.loc[gcta['p'] == 0].shape[0] > ngene:
            gcta = gcta.loc[gcta['p'] == 0].sample(n=ngene, replace=False, random_state=rng)
        else:
            gcta = gcta.iloc[:ngene]
        np.savetxt(output.gcta, gcta['gene'], fmt='%s')


rule yazar_ldsc_make_genecoord:
    input:
        location = 'data/Yazar2022Science/gene_location.txt',
    output:
        gene_coord = 'staging/data/ldsc/genecoord.txt',
    run:
        data = pd.read_table(input.location, usecols=['feature', 'chr', 'start', 'end'])
        data = data.rename(columns={'feature': 'GENE', 'chr': 'CHR', 'start': 'START', 'end': 'END'})
        data[['GENE', 'CHR', 'START', 'END']].to_csv(output.gene_coord, sep='\t', index=False)


rule yazar_ldsc_rmMHC:
    input:
        bim = 'data/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.{chr}.bim',
    output:
        bim = 'analysis/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.noMHC.{chr}.bim',
        bed = 'analysis/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.noMHC.{chr}.bed',
    params:
        chr = config['ldsc']['mhc']['chr'],
        start = config['ldsc']['mhc']['start'],
        end = config['ldsc']['mhc']['end'],
    shell:
        '''
        # load plink 1.9
        module load gcc/11.3.0 atlas/3.10.3 lapack/3.11.0 plink/1.9

        input_prefix="$(dirname {input.bim})/$(basename {input.bim} .bim)"
        output_prefix="$(dirname {output.bim})/$(basename {output.bim} .bim)"

        if [ {wildcards.chr} -eq {params.chr} ]; then
            awk -v start={params.start} -v end={params.end} '($4 < start || $4 > end) {{print $2}}' {input.bim} > {output.bed}.tmp.snps

            plink --bfile $input_prefix  --extract {output.bed}.tmp.snps \
                    --make-bed --out $output_prefix
            rm {output.bed}.tmp.snps
        else
            cp $input_prefix.bed $output_prefix.bed
            cp $input_prefix.bim $output_prefix.bim
            cp $input_prefix.fam $output_prefix.fam
        fi
        '''


rule yazar_ldsc_make_annot:
    input:
        all = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/all.genes.txt',
        random = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/random.genes.txt',
        shared = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/shared.genes.txt',
        var = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/var.genes.txt',
        mean = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/mean.genes.txt',
        gcta = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/gcta.genes.txt',
        gene_coord = 'staging/data/ldsc/genecoord.txt',
        bim = 'analysis/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.noMHC.{chr}.bim',
    output:
        tmp_all = temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/all.{{chr}}.annot.tmp.gz'),
        tmp_random = temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/random.{{chr}}.annot.tmp.gz'),
        tmp_shared = temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/shared.{{chr}}.annot.tmp.gz'),
        tmp_var = temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/var.{{chr}}.annot.tmp.gz'),
        tmp_mean = temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/mean.{{chr}}.annot.tmp.gz'),
        tmp_gcta = temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/gcta.{{chr}}.annot.tmp.gz'),
        bim = temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/{{chr}}.bim'),
        all = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/all.{{chr}}.annot.gz',
        random = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/random.{{chr}}.annot.gz',
        shared = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/shared.{{chr}}.annot.gz',
        var = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/var.{{chr}}.annot.gz',
        mean = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/mean.{{chr}}.annot.gz',
        gcta = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/gcta.{{chr}}.annot.gz',
    conda: '../../ldsc/environment.yml'
    shell:
        '''
        python ldsc/make_annot.py \
            --gene-set-file {input.all} \
            --gene-coord-file {input.gene_coord} \
            --windowsize {wildcards.window} \
            --bimfile {input.bim} \
            --annot-file {output.tmp_all}

        python ldsc/make_annot.py \
            --gene-set-file {input.random} \
            --gene-coord-file {input.gene_coord} \
            --windowsize {wildcards.window} \
            --bimfile {input.bim} \
            --annot-file {output.tmp_random}

        python ldsc/make_annot.py \
            --gene-set-file {input.shared} \
            --gene-coord-file {input.gene_coord} \
            --windowsize {wildcards.window} \
            --bimfile {input.bim} \
            --annot-file {output.tmp_shared}

        python ldsc/make_annot.py \
            --gene-set-file {input.var} \
            --gene-coord-file {input.gene_coord} \
            --windowsize {wildcards.window} \
            --bimfile {input.bim} \
            --annot-file {output.tmp_var}

        python ldsc/make_annot.py \
            --gene-set-file {input.mean} \
            --gene-coord-file {input.gene_coord} \
            --windowsize {wildcards.window} \
            --bimfile {input.bim} \
            --annot-file {output.tmp_mean}

        python ldsc/make_annot.py \
            --gene-set-file {input.gcta} \
            --gene-coord-file {input.gene_coord} \
            --windowsize {wildcards.window} \
            --bimfile {input.bim} \
            --annot-file {output.tmp_gcta}

        # add snp info
        echo -e 'CHR\\tBP\\tSNP\\tCM' > {output.bim}
        awk 'BEGIN {{OFS = "\\t"}} {{print $1, $4, $2, $3}}' {input.bim} >> {output.bim}
        paste {output.bim} <(zcat {output.tmp_all}) | gzip -c > {output.all}
        paste {output.bim} <(zcat {output.tmp_random}) | gzip -c > {output.random}
        paste {output.bim} <(zcat {output.tmp_shared}) | gzip -c > {output.shared}
        paste {output.bim} <(zcat {output.tmp_var}) | gzip -c > {output.var}
        paste {output.bim} <(zcat {output.tmp_mean}) | gzip -c > {output.mean}
        paste {output.bim} <(zcat {output.tmp_gcta}) | gzip -c > {output.gcta}
        '''


rule yazar_ldsc_compldscore:
    input:
        all = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/all.{{chr}}.annot.gz',
        random = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/random.{{chr}}.annot.gz',
        shared = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/shared.{{chr}}.annot.gz',
        var = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/var.{{chr}}.annot.gz',
        mean = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/mean.{{chr}}.annot.gz',
        gcta = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/gcta.{{chr}}.annot.gz',
        bim = 'analysis/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.noMHC.{chr}.bim',
        hapmap3 = 'data/ldsc/hm3_no_MHC.list.txt',
    output:
        all = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/all.{{chr}}.l2.ldscore.gz',
        random = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/random.{{chr}}.l2.ldscore.gz',
        shared = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/shared.{{chr}}.l2.ldscore.gz',
        var = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/var.{{chr}}.l2.ldscore.gz',
        mean = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/mean.{{chr}}.l2.ldscore.gz',
        gcta = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/gcta.{{chr}}.l2.ldscore.gz',
    params:
        bfile = lambda wildcards, input: os.path.splitext(input.bim)[0],
        all = lambda wildcards, output: re.sub('\.l2.ldscore.gz$', '', output.all),
        random = lambda wildcards, output: re.sub('\.l2.ldscore.gz$', '', output.random),
        shared = lambda wildcards, output: re.sub('\.l2.ldscore.gz$', '', output.shared),
        var = lambda wildcards, output: re.sub('\.l2.ldscore.gz$', '', output.var),
        mean = lambda wildcards, output: re.sub('\.l2.ldscore.gz$', '', output.mean),
        gcta = lambda wildcards, output: re.sub('\.l2.ldscore.gz$', '', output.gcta),
    conda: '../../ldsc/environment.yml'
    shell:
        '''
        python ldsc/ldsc.py \
            --l2 \
            --bfile {params.bfile} \
            --ld-wind-cm 1 \
            --annot {input.all} \
            --out {params.all} \
            --print-snps {input.hapmap3}

        python ldsc/ldsc.py \
            --l2 \
            --bfile {params.bfile} \
            --ld-wind-cm 1 \
            --annot {input.random} \
            --out {params.random} \
            --print-snps {input.hapmap3}

        python ldsc/ldsc.py \
            --l2 \
            --bfile {params.bfile} \
            --ld-wind-cm 1 \
            --annot {input.shared} \
            --out {params.shared} \
            --print-snps {input.hapmap3}

        python ldsc/ldsc.py \
            --l2 \
            --bfile {params.bfile} \
            --ld-wind-cm 1 \
            --annot {input.var} \
            --out {params.var} \
            --print-snps {input.hapmap3}

        python ldsc/ldsc.py \
            --l2 \
            --bfile {params.bfile} \
            --ld-wind-cm 1 \
            --annot {input.mean} \
            --out {params.mean} \
            --print-snps {input.hapmap3}

        python ldsc/ldsc.py \
            --l2 \
            --bfile {params.bfile} \
            --ld-wind-cm 1 \
            --annot {input.gcta} \
            --out {params.gcta} \
            --print-snps {input.hapmap3}
        '''


rule yazar_ldsc_gwas_addP:
    input:
        gwas = lambda wildcards: gwass[wildcards.gwas],
    output:
        gwas = 'staging/data/ldsc/{gwas}.gz',
    run:
        from scipy.stats import norm

        gwas = pd.read_table(input.gwas)
        if 'p' in gwas.columns or 'P' in gwas.columns:
            pass 
        else:
            gwas['P'] = gwas['Z'].apply(lambda z: norm.sf(abs(z)) * 2)

        # for gwas with negative Z
        if np.any(gwas['Z'] < 0):
            gwas.loc[gwas['Z'] < 0, ['A1', 'A2']] = gwas.loc[gwas['Z'] < 0, ['A2', 'A1']].values
            gwas['Z'] = np.abs(gwas['Z'])

        gwas.to_csv(output.gwas, sep='\t', index=False)


rule yazar_ldsc_format_gwas:
    input:
        gwas = 'staging/data/ldsc/{gwas}.gz',
        hapmap3 = 'data/ldsc/broad/GINGER/ginger_vc_c1_year3/lab_11.13.19/w_hm3.snplist',  # snp list with A1 and A2
    output:
        gwas = 'staging/data/ldsc/{gwas}.sumstats.gz',
    params:
        prefix = lambda wildcards, output: re.sub('\.sumstats.gz$', '', output.gwas),
    resources:
        mem_mb = '10G',
    conda: '../../ldsc/environment.yml'
    shell:
        # NOTE: assmuning A1 is trait-increasing
        '''
        python ldsc/munge_sumstats.py \
            --sumstats {input.gwas} \
            --merge-alleles {input.hapmap3} \
            --chunksize 500000 \
            --a1-inc \
            --out {params.prefix}
        '''


rule yazar_ldsc_seg_make_ldcts:
    input:
        shared = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/shared.{chr}.l2.ldscore.gz'
                    for chr in chrs],
        var = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/var.{chr}.l2.ldscore.gz'
                    for chr in chrs],
        mean = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/mean.{chr}.l2.ldscore.gz'
                    for chr in chrs],
        gcta = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/gcta.{chr}.l2.ldscore.gz'
                    for chr in chrs],
        control = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/all.{chr}.l2.ldscore.gz'
                    for chr in chrs],
    output: 
        ldcts = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/he.ldcts',
    run:
        with open(output.ldcts, 'w') as f:
            f.write('shared\t{},{}\n'.format( re.sub('1.l2.ldscore.gz', '', input.shared[0]), 
                                            re.sub('1.l2.ldscore.gz', '', input.control[0])))
            f.write('var\t{},{}\n'.format( re.sub('1.l2.ldscore.gz', '', input.var[0]), 
                                            re.sub('1.l2.ldscore.gz', '', input.control[0])))
            f.write('mean\t{},{}\n'.format( re.sub('1.l2.ldscore.gz', '', input.mean[0]), 
                                            re.sub('1.l2.ldscore.gz', '', input.control[0])))
            f.write('gcta\t{},{}\n'.format( re.sub('1.l2.ldscore.gz', '', input.gcta[0]), 
                                            re.sub('1.l2.ldscore.gz', '', input.control[0])))


rule yazar_ldsc_seg_regression:
    input:
        gwas = 'staging/data/ldsc/{gwas}.sumstats.gz',
        ldcts = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/he.ldcts',
        ref_ld = expand('data/ldsc/baseline_v1.2/baseline.{chr}.l2.ldscore.gz', chr=chrs),  # readme_baseline_versions: use baselineLD v2.2 for estimating heritability enrichment; baseline v1.2 for identifying critical tissues/cell-types via P-value of tau  
        weight = expand('data/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.{chr}.l2.ldscore.gz', chr=chrs),
    output:
        out = f'results/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/he.{{gwas}}.cell_type_results.txt',
    params:
        ref_ld = lambda wildcards, input: re.sub('1.l2.ldscore.gz', '', input.ref_ld[0]),
        weight = lambda wildcards, input: re.sub('1.l2.ldscore.gz', '', input.weight[0]),
        out = lambda wildcards, output: re.sub('\.cell_type_results.txt', '', output.out),
    resources:
        mem_mb = '10G',
    conda: '../../ldsc/environment.yml'
    shell:
        '''
        ldsc/ldsc.py \
            --h2-ct {input.gwas} \
            --ref-ld-chr {params.ref_ld} \
            --out {params.out} \
            --ref-ld-chr-cts {input.ldcts} \
            --w-ld-chr {params.weight}
        '''


rule yazar_ldsc_seg_regression_summary:
    input:
        stacked = [f'results/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/he.{gwas}.cell_type_results.txt'
                    for gwas in traits],
    output:
        png = f'results/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/he.h2_cts.png',
    params:
        neg_traits = neg_traits,
        pos_traits = pos_traits,
        traits = traits,
    script: '../bin/yazar/ldsc.h2_cts.py'


rule yazar_ldsc_controlmean_seg_make_ldcts:
    input:
        shared = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/shared.{chr}.l2.ldscore.gz'
                    for chr in chrs],
        var = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/var.{chr}.l2.ldscore.gz'
                    for chr in chrs],
        mean = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/mean.{chr}.l2.ldscore.gz'
                    for chr in chrs],
        gcta = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/gcta.{chr}.l2.ldscore.gz'
                    for chr in chrs],
        control = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/all.{chr}.l2.ldscore.gz'
                    for chr in chrs],
    output: 
        ldcts = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/he.controlmean.ldcts',
    run:
        with open(output.ldcts, 'w') as f:
            f.write('shared\t{},{},{}\n'.format( re.sub('1.l2.ldscore.gz', '', input.shared[0]), 
                                            re.sub('1.l2.ldscore.gz', '', input.mean[0]),
                                            re.sub('1.l2.ldscore.gz', '', input.control[0])))
            f.write('var\t{},{},{}\n'.format( re.sub('1.l2.ldscore.gz', '', input.var[0]), 
                                            re.sub('1.l2.ldscore.gz', '', input.mean[0]),
                                            re.sub('1.l2.ldscore.gz', '', input.control[0])))
            f.write('gcta\t{},{},{}\n'.format( re.sub('1.l2.ldscore.gz', '', input.gcta[0]), 
                                            re.sub('1.l2.ldscore.gz', '', input.mean[0]),
                                            re.sub('1.l2.ldscore.gz', '', input.control[0])))


use rule yazar_ldsc_seg_regression as yazar_ldsc_controlmean_seg_regression with:
    input:
        gwas = 'staging/data/ldsc/{gwas}.sumstats.gz',
        ldcts = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/he.controlmean.ldcts',
        ref_ld = expand('data/ldsc/baseline_v1.2/baseline.{chr}.l2.ldscore.gz', chr=chrs),  # readme_baseline_versions: use baselineLD v2.2 for estimating heritability enrichment; baseline v1.2 for identifying critical tissues/cell-types via P-value of tau  
        weight = expand('data/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.{chr}.l2.ldscore.gz', chr=chrs),
    output:
        out = f'results/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/he.{{gwas}}.controlmean.cell_type_results.txt',


use rule yazar_ldsc_seg_regression_summary as yazar_ldsc_controlmean_seg_regression_summary with:
    input:
        stacked = [f'results/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/he.{gwas}.controlmean.cell_type_results.txt'
                    for gwas in traits],
    output:
        png = f'results/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/he.h2_cts.controlmean.png',


rule yazar_ldsc_all:
    input:
        # h2_cts = expand('results/yazar/ind_min_cellnum~10_ct_min_cellnum~10_prop~0.9_geno_pca_n~6_op_pca_n~1_batch~shared_fixed~shared/ldsc/top_{ngene}/window_{window}/he.h2_cts.png',
        h2_cts = expand('results/yazar/{params}/ldsc/top_{ngene}/window_{window}/he.h2_cts.png',
                        params=yazar_paramspace.instance_patterns, 
                        ngene = [0, 100, 200, 300],
                        window=[300000, 500000, 700000]),
        h2_cts_controlmean = expand('results/yazar/{params}/ldsc/top_{ngene}/window_{window}/he.h2_cts.controlmean.png',
                        params=yazar_paramspace.instance_patterns, 
                        ngene = [0, 100, 200, 300],
                        window=[300000, 500000, 700000]),






































                






















































































##########################################################################
## 1.2: trans-genetic effect
##########################################################################

rule yazar_trans_genome_kinship:
    input:
        bed = expand('analysis/yazar/data/geno/chr{chr}.bed',
                chr=range(1,23)),
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.final.gz',
    output:
        kinship = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/trans/he/genome.rel',
    params:
        prefix = lambda wildcards, output: os.path.splitext(output.kinship)[0],
        tmp_dir = lambda wildcards, output: os.path.splitext(output.kinship)[0] + '_tmp',
    resources:
        mem_mb = '20G',
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
        # grm
        plink --bfile $merged.ld --make-rel --out {params.prefix}
        rm -r {params.tmp_dir}
        '''


use rule yazar_trans_genome_kinship as yazar_trans_rm_target_chr_kinship with:
    input:
        bed = lambda wildcards: expand('analysis/yazar/data/geno/chr{chr}.bed',
                chr=[x for x  in chrs if int(x) != int(wildcards.chr)]),
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.final.gz',
    output:
        kinship = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/trans/he/chr{{chr}}.rel',


use rule yazar_HE_free as yazar_trans_HE_free with:
    input:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.batch{{i}}.npy',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.final.gz',
        op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.gz',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
        genome = [f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/trans/he/chr{chr}.rel'
                    for chr in chrs],
        genes = 'data/Yazar2022Science/gene_location.txt',
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/trans/he.batch{{i}}.npy',
    params:
        chrs = chrs,
        snps = 5, # threshold of snp number per gene
        iid = True,
        full = False,
    resources:
        mem_mb = lambda wildcards:'40G' if wildcards.batch != 'shared' else '20G',


use rule yazar_HE_free_merge as yazar_trans_HE_free_merge with:
    input:
        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/trans/he.batch{i}.npy'
            for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/trans/he.npy',


use rule yazar_HE_togz as yazar_trans_HE_togz with:
    input:
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/trans/he.npy',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.final.gz',
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/trans/he.gz',


use rule yazar_trans_HE_free as yazar_trans_HE_free_jk with:
    input:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/reml/kinship.batch{{i}}.npy',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.final.gz',
        op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.gz',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
        genome = [f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/trans/he/chr{chr}.rel'
                    for chr in chrs],
        genes = 'data/Yazar2022Science/gene_location.txt',
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/trans/he.jk.batch{{i}}.npy',
    params:
        chrs = chrs,
        snps = 5, # threshold of snp number per gene
        replicates = 0.3,
        iid = False,
        full = False,
        jk = True,
    resources:
        mem_mb = '20G',
        time = '200:00:00',


use rule yazar_HE_free_merge as yazar_trans_HE_free_jk_merge with:
    input:
        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/trans/he.jk.batch{i}.npy'
            for i in range(yazar_reml_batches)],
    output:
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/trans/he.jk.npy',



rule yazar_trans_all:
    input:
        jk = expand('analysis/yazar/{params}/trans/he.jk.npy',
                        params=yazar_paramspace.instance_patterns),











###########################################################
# trans-genetic effect restricting to main cell types
###########################################################
rule yazar_extract_mainct:
    input:
        ctp = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctp.final.gz',
        ctnu = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctnu.final.gz',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.final.gz',
    output:
        op = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/mainct/op.final.gz',
        ctp = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/mainct/ctp.final.gz',
        ctnu = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/mainct/ctnu.final.gz',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/mainct/P.final.gz',
    run:
        P = pd.read_table(input.P, index_col=0)
        cts = P.columns.values
        # find main cell types
        pi = P.mean(axis=0)
        main_cts = pi.sort_values(ascending=False).index[:4].to_numpy()
        print(main_cts)

        # P
        P = P.loc[:, P.columns.isin(main_cts)]
        print(P.head())
        P = P.div(P.sum(axis=1), axis=0)
        P.to_csv(output.P, sep='\t')

        # ctp
        ctp = pd.read_table(input.ctp, index_col=(0, 1))
        ctp = ctp.loc[ctp.index.get_level_values('ct').isin(main_cts)]
        ctp.to_csv(output.ctp, sep='\t')

        # ctnu
        ctnu = pd.read_table(input.ctnu)
        ctnu = ctnu.loc[ctnu['ct'].isin(main_cts)]
        ctnu.to_csv(output.ctnu, sep='\t', index=False)

        # op
        ctp = ctp.unstack()
        # check ct order
        assert ctp[ctp.columns.levels[0][0]].columns.equals(P.columns)
        # ctp * P
        for gene in ctp.columns.levels[0]:
            ctp[gene] = ctp[gene].mul(P)
        ctp = ctp.stack()
        # sum ctp to op
        op = ctp.groupby(level=0).sum()
        op.to_csv(output.op, sep='\t')


use rule yazar_op_pca as yazar_op_pca_mainct with:
    input:
        op = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/mainct/op.final.gz',
    output:
        evec = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/mainct/evec.gz',
        eval = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/mainct/eval.gz',
        pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/mainct/pca.gz',
        png = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/mainct/op.pca.png',


use rule yazar_HE_split as yazar_HE_split_mainct with:
    input:
        ctp = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/mainct/ctp.final.gz',
        ctnu = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/mainct/ctnu.final.gz',
        genes = 'data/Yazar2022Science/gene_location.txt',
    output:
        ctp = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/mainct/he/ctp.batch{i}.gz'
                for i in range(yazar_he_batches)],
        ctnu = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/mainct/he/ctnu.batch{i}.gz'
                for i in range(yazar_he_batches)],
        batch = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/mainct/he/ctp.batch.txt',


use rule yazar_he_kinship_split as yazar_he_kinship_split_mainct with:
    input:
        kinship = [f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/kinship.chr{chr}.npy'
                    for chr in chrs],
        ctp = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/mainct/he/ctp.batch{i}.gz'
                for i in range(yazar_he_batches)],
        batch = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/mainct/he/ctp.batch.txt',
    output:
        kinship = [temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/mainct/he/kinship.batch{i}.npy')
                    for i in range(yazar_he_batches)],


use rule yazar_HE_free as yazar_trans_HE_free_mainct with:
    input:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/mainct/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/mainct/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/mainct/he/kinship.batch{{i}}.npy',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/mainct/P.final.gz',
        op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/mainct/pca.gz',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
        genome = [f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/trans/he/chr{chr}.rel'
                    for chr in chrs],
        genes = 'data/Yazar2022Science/gene_location.txt',
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/trans/mainct/he.batch{{i}}.npy',
    params:
        chrs = chrs,
        snps = 5, # threshold of snp number per gene
        iid = False,
        full = False,
    resources:
        mem_mb = '20G',


use rule yazar_HE_free_merge as yazar_trans_HE_free_mainct_merge with:
    input:
        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/trans/mainct/he.batch{i}.npy'
            for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/trans/mainct/he.npy',








































##########################################################################
## 1.2 transformation
##########################################################################
use rule yazar_ctp_transform as yazar_transform_ctp_transform with:
    output:
        X = 'staging/data/yazar/transform_{transform}/X.npz',
        obs = 'staging/data/yazar/transform_{transform}/obs.gz',
        var = 'staging/data/yazar/transform_{transform}/var.gz',
    params:
        ind_col = yazar_ind_col,
        ct_col = yazar_ct_col,
        pool_col = yazar_pool_col,
    resources:
        mem_mb = lambda wildcards: '300G' if wildcards.transform == 'pearson' else '90G',


use rule yazar_ctp as yazar_transform_ctp with:
    input:
        X = 'staging/data/yazar/transform_{transform}/X.npz',
        obs = 'staging/data/yazar/transform_{transform}/obs.gz',
        var = 'staging/data/yazar/transform_{transform}/var.gz',
    output:
        ctp = 'staging/data/yazar/transform_{transform}/ctp.gz',
        ctnu = 'staging/data/yazar/transform_{transform}/ctnu.gz',
        P = 'staging/data/yazar/transform_{transform}/P.gz',
        n = 'staging/data/yazar/transform_{transform}/n.gz',


use rule yazar_rm_rareINDnCT as yazar_transform_rm_rareINDnCT with:
    # also select gene expressed in all cts
    input:
        ctp = 'staging/data/yazar/transform_{transform}/ctp.gz',
        ctnu = 'staging/data/yazar/transform_{transform}/ctnu.gz',
        n = 'staging/data/yazar/transform_{transform}/n.gz',
    output:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/ctp.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/ctnu.gz',
        P = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/P.gz',
        n = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/n.gz',


use rule yazar_rm_missingIND as yazar_transform_rm_missingIND with:
    input:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/ctp.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/ctnu.gz',
        P = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/P.gz',
        n = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/n.gz',
    output:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/ctp.rmid.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/ctnu.rmid.gz',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/P.final.gz',
        n = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/n.final.gz',


use rule yazar_std_op as yazar_transform_std_op with:
    input:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/ctp.rmid.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/ctnu.rmid.gz',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/P.final.gz',
    output:
        op = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/op.std.gz',
        nu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/nu.std.gz',
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/ctp.std.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/ctnu.std.gz',


use rule yazar_exclude_sexchr as yazar_transform_exclude_sexchr with:
    input:
        genes = 'data/Yazar2022Science/gene_location.txt',
        op = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/op.std.gz',
        nu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/nu.std.gz',
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/ctp.std.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/ctnu.std.gz',
    output:
        op = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/op.final.gz',
        nu = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/nu.final.gz',
        ctp = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/ctp.final.gz',
        ctnu = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/ctnu.final.gz',


use rule yazar_op_pca as yazar_transform_op_pca with:
    input:
        op = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/op.final.gz',
    output:
        evec = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/evec.txt',
        eval = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/eval.txt',
        pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/pca.gz',
        png = f'results/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/pca.png',


# use rule yazar_he_kinship as yazar_transform_he_kinship with:
#     input:
#         ctp = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/ctp.final.gz',
#         genes = 'data/Yazar2022Science/gene_location.txt',
#         bed = 'analysis/yazar/data/geno/chr{chr}.bed',
#     output:
#         kinship = temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/he/kinship.chr{{chr}}.npy'),


use rule yazar_HE_split as yazar_transform_he_split with:
    input:
        ctp = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/ctp.final.gz',
        ctnu = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/ctnu.final.gz',
        genes = 'data/Yazar2022Science/gene_location.txt',
    output:
        ctp = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/he/ctp.batch{i}.gz'
                for i in range(yazar_he_batches)],
        ctnu = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/he/ctnu.batch{i}.gz'
                for i in range(yazar_he_batches)],
        batch = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/he/ctp.batch.txt',


# use rule yazar_he_kinship_split as yazar_transform_he_kinship_split with:
#     input:
#         kinship = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/he/kinship.chr{chr}.npy'
#                     for chr in chrs],
#         ctp = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/he/ctp.batch{i}.gz'
#                 for i in range(yazar_he_batches)],
#         batch = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/he/ctp.batch.txt',
#     output:
#         kinship = [temp(f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/he/kinship.batch{i}.npy')
#                 for i in range(yazar_he_batches)],


use rule yazar_HE_free as yazar_transform_trans_HE_free with:
    # NOTE: assume transformation doesn't change individuals
    input:
        ctp = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/he/kinship.batch{{i}}.npy',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/P.final.gz',
        op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/pca.gz',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
        genome = [f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/trans/he/chr{chr}.rel'
                    for chr in chrs],
        genes = 'data/Yazar2022Science/gene_location.txt',
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/trans/he.batch{{i}}.npy',
    params:
        chrs = chrs,
        snps = 5, # threshold of snp number per gene
        iid = True,
        full = False,


use rule yazar_HE_free_merge as yazar_transform_trans_HE_free_merge with:
    input:
        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/trans/he.batch{i}.npy'
            for i in range(yazar_he_batches)],
    output:
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/transform_{{transform}}/trans/he.npy',



rule yazar_transform_all:
    input:
        out = expand('analysis/yazar/{params}/transform_{transform}/trans/he.npy',
                        params=yazar_paramspace.instance_patterns,
                        transform=['logp_cpm', 'pearson']),









rule yazar_all:
    input:
        # cis = expand('analysis/yazar/{params}/he.npy',
        #                 params=yazar_paramspace.instance_patterns),
        trans = expand('analysis/yazar/{params}/trans/he.npy',
                        params=yazar_paramspace.instance_patterns),


