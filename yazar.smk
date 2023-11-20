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
        mem_mb = '90G',
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
        X = preprocess.normalize(data.X, 1e4).log1p()
        sparse.save_npz(output.X, X)

        data.obs.rename_axis('cell').to_csv(output.obs, sep='\t')
        data.var.rename_axis('feature').to_csv(output.var, sep='\t')


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
        mem_mb = '80G',
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


var_ctnu_batches = 500
rule yazar_var_ctnu_split:
    input:
        obs = 'staging/data/yazar/obs.gz',
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
        X = 'staging/data/yazar/X.npz',
        obs = 'staging/data/yazar/obs.gz',
        var = 'staging/data/yazar/var.gz',
        batch = 'staging/data/yazar/var_ctnu/ind_ct.batch{i}',
    output:
        var_ctnu = 'staging/data/yazar/var_ctnu/batch{i}.gz',
    params:
        ind_col = yazar_ind_col,
        ct_col = yazar_ct_col,
        seed = 42,
    resources:
        mem_mb = '2G',
        partition = 'tier3q',
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
        var_ctnu = expand('staging/data/yazar/var_ctnu/batch{i}.gz', i=range(var_ctnu_batches)),
    output:
        var_ctnu = 'analysis/yazar/var_ctnu.gz',
    run:
        var_ctnu = [pd.read_table(f) for f in input.var_ctnu]
        pd.concat(var_ctnu, ignore_index=True).to_csv(output.var_ctnu,
                sep='\t', index=False)


