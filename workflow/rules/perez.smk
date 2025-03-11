#######################################################################################
# Perez 2022 Science
#######################################################################################
perez_ind_col = 'ind_cov'
perez_ct_col = 'cg_cov'
perez_pool_col = 'batch_cov'


perez_ct_order = np.array(['hom', 'T4', 'cM', 'T8', 'B', 'NK', 'ncM', 'cDC', 'Prolif', 'pDC', 'PB', 'Progen'])
perez_colors = dict(zip(perez_ct_order[1:], sns.color_palette(config['colorpalette'])))
perez_colors['hom'] = '0.7'
# paried light colors
colors = [colorsys.rgb_to_hsv(*color) for color in sns.color_palette(config['colorpalette'])]
colors = [(color[0], min(1, max(0, color[1] - 0.3)), color[2]) for color in colors]
colors = [colorsys.hsv_to_rgb(*color) for color in colors]
perez_light_colors = dict(zip(perez_ct_order[1:], colors))
perez_light_colors['hom'] = '0.9'

# read parameters
perez_params = pd.read_table('perez.params.txt', dtype="str", comment='#')
perez_paramspace = Paramspace(perez_params, filename_params="*")

perez_he_batches = config['perez']['he_nbatch']
# perez_reml_batches = config['perez']['reml_nbatch']

rule perez_merge_vcfs:
    input:
        vcf1 = 'data/Perez2022Science/clues/all_clues.processed.vcf.gz',
        vcf2 = 'data/Perez2022Science/clues/immvar.processed.vcf.gz',
    output:
        vcf = temp('staging/perez/data/geno/perez.vcf.gz'),
        bed = expand('analysis/perez/data/geno/chr{chr}.bed', chr=chrs),
    params:
        dir = lambda wildcards, output: os.path.dirname(output.vcf),
        prefix = lambda wildcards, output: os.path.splitext(output.bed[0])[0][:-1],
    shell:
        """
        module load gcc/11.3.0 htslib/1.17 bcftools/1.17
        module load gcc/11.3.0 atlas/3.10.3 lapack/3.11.0 plink/1.9

        mkdir -p {params.dir}_tmp
        bcftools isec -p {params.dir}_tmp -O z {input.vcf1} {input.vcf2}
        bcftools merge {params.dir}_tmp/0002.vcf.gz {params.dir}_tmp/0003.vcf.gz -o {output.vcf}

        for chr in {{1..22}}
        do
            plink --vcf {output.vcf} --double-id \
                    --keep-allele-order --chr $chr \
                    --make-bed --out {params.prefix}$chr
        done

        rm -rf {params.dir}_tmp
        """


rule perez_check_numSNPs:
    input:
        bed = expand('analysis/perez/data/geno/chr{chr}.bed',
                chr=range(1,23)),
    output:
        txt = 'analysis/perez/data/geno/num_snps.txt',
    params:
        tmp_dir = 'staging/perez/data/geno/tmp',
    shell:
        """
        module load gcc/11.3.0 atlas/3.10.3 lapack/3.11.0 plink/1.9

        mkdir -p {params.tmp_dir}
        for bed in {input.bed} 
        do 
            prefix="$(dirname $bed)/$(basename $bed .bed)"
            o_prefix="{params.tmp_dir}/$(basename $bed .bed)"
            plink --bfile $prefix \
                --maf 0.1 \
                --make-bed --out $o_prefix
        done
        wc -l {params.tmp_dir}/*bim > {output.txt}
        rm -r {params.tmp_dir}
        """

rule perez_udpate_anndata:
    # extract raw counts
    # change gene index to ensembl
    input:
        h5ad = 'data/Perez2022Science/GSE174188_CLUES1_adjusted.h5ad',
    output:
        h5ad = 'data/Perez2022Science/adjusted.all.h5ad',
        obs = 'data/Perez2022Science/adjusted.all.obs.gz',
        var = 'data/Perez2022Science/adjusted.all.var.gz',
    run:
        import anndata as ad
        import scanpy as sc

        data = sc.read_h5ad(input.h5ad, backed='r')
        adata = ad.AnnData(data.raw.X)
        adata.obs = data.obs
        var = data.raw.var.reset_index(names='GeneSymbol')
        var = var.rename(columns={'gene_ids': 'feature'})
        var = var.set_index('feature')
        adata.var = var

        adata.write_h5ad(output.h5ad, compression='gzip', compression_opts=9)

        obs = adata.obs.reset_index(drop=False, names='cell')
        obs.to_csv(output.obs, sep='\t', index=False)

        var = adata.var.reset_index(drop=False, names='feature')
        for column in var.columns:
            var[column] = var[column].str.replace(' ', '')
        var.to_csv(output.var, sep='\t', index=False)


rule perez_gene_location_allgenes:
    input:
        genes = 'data/Perez2022Science/adjusted.all.var.gz',
        gene_info = 'analysis/data/Homo_sapiens.GRCh37.75.txt',
    output:
        location = 'analysis/perez/data/gene_location.txt',
        meta = 'analysis/perez/data/gene_meta.txt',
    script: '../bin/perez/gene_location.py'


use rule yazar_cell_dist as perez_cell_dist with:
    input:
        obs = 'data/Perez2022Science/adjusted.all.obs.gz',
    output:
        png = 'results/perez/cell.dist.png',
    params:
        ind_col = perez_ind_col,
        ct_col = perez_ct_col,
        colors = perez_colors,
        analyzed_ct_num = 7,  # number of main cell types included for analysis


rule perez_age:
    input:
        obs = 'data/Perez2022Science/adjusted.all.obs.gz',
    output:
        png = 'results/perez/age.png',
    params:
        ind_col = perez_ind_col,
    run:
        obs = pd.read_table(input.obs)
        obs = obs.drop_duplicates(subset=params.ind_col)
        ages, counts = np.unique(obs['Age'], return_counts=True)
        fig, ax = plt.subplots(dpi=600)
        plt.bar(ages, counts)
        ax.set_xlabel('Age')
        ax.set_ylabel('Number of individuals')
        fig.savefig(output.png)


rule perez_filter_samples:
    input:
        h5ad = 'data/Perez2022Science/adjusted.all.h5ad',
        bed = 'analysis/perez/data/geno/chr1.bed',
    output:
        h5ad = 'analysis/perez/data/{anc}.{status}.h5ad',
        inds = 'analysis/perez/data/{anc}.{status}.ids',
        log = 'analysis/perez/data/{anc}.{status}.log',
    resources:
        mem_mb = '30G',
    run:
        import scanpy as sc

        adata = sc.read_h5ad(input.h5ad, backed='r')

        log = open(output.log, 'w')

        # only keep Asian Case and European Control + Case
        cells = adata.obs['pop_cov'].isin(['Asian', 'European'])
        asian_control = (adata.obs['pop_cov'] == 'Asian') & (adata.obs['SLE_status'] == 'Healthy')
        cells = cells & (~asian_control)
        nind = adata[cells, :].obs['ind_cov'].nunique()
        log.write(f'Number of individuals: {nind}\n')

        # only keep individuals with genotype
        fam_f = input.bed.replace('.bed', '.fam')
        ind_ids = [line.split()[0] for line in open(fam_f)]
        cells = cells & adata.obs['ind_cov'].isin(ind_ids)
        nind = adata[cells, :].obs['ind_cov'].nunique()
        log.write(f'Number of individuals with genotype: {nind}\n')

        # keep individuals within pop and with status
        if wildcards.anc == 'all':
            pass
        else:
            cells = cells & (adata.obs['pop_cov'] == wildcards.anc)
        if wildcards.status == 'all':
            pass
        else:
            cells = cells & (adata.obs['SLE_status'] == wildcards.status)

        adata = adata[cells, :]
        # save filtered adata
        adata.write_h5ad(output.h5ad, compression='gzip', compression_opts=9)
        nind = adata.obs['ind_cov'].nunique()
        log.write(f'Number of individuals after filtering: {nind}\n')

        # save individual ids
        inds = adata.obs['ind_cov'].unique().tolist()
        with open(output.inds, 'w') as f:
            for ind in inds:
                f.write(f'{ind}\t{ind}\n')


rule perez_filter_plink:
    input:
        bed = 'analysis/perez/data/geno/chr{chr}.bed',
        inds = 'analysis/perez/data/{anc}.{status}.ids',
    output:
        bed = 'analysis/perez/data/geno/{anc}.{status}/chr{chr}.bed',
    shell:
        """
        module load gcc/11.3.0 atlas/3.10.3 lapack/3.11.0 plink/1.9

        prefix="$(dirname {input.bed})/$(basename {input.bed} .bed)"
        o_prefix="$(dirname {output.bed})/$(basename {output.bed} .bed)"
        plink --bfile $prefix \
            --keep {input.inds} \
            --make-bed --out $o_prefix
        """


use rule yazar_extract_meta as perez_extract_meta with:
   input:
        h5ad = 'analysis/perez/data/{anc}.{status}.h5ad',
   output:
        obs = 'analysis/perez/data/{anc}.{status}.obs.gz',
        var = 'analysis/perez/data/{anc}.{status}.var.gz',


rule perez_add_dataset_info:
    input:
        obs = 'analysis/perez/data/{anc}.{status}.obs.gz',
        vcf1 = 'data/Perez2022Science/clues/all_clues.processed.vcf.gz',
        vcf2 = 'data/Perez2022Science/clues/immvar.processed.vcf.gz',
    output:
        obs = temp('analysis/perez/data/{anc}.{status}.obs.dataset.gz'),
    run:
        obs = pd.read_table(input.obs)
        # collect individual ids from vcf
        for line in gzip.open(input.vcf1, 'rt'):
            if line.startswith('#CHROM'):
                clues_ids = line.strip().split()[9:]
                break 
        for line in gzip.open(input.vcf2, 'rt'):
            if line.startswith('#CHROM'):
                immvar_ids = line.strip().split()[9:]
                break
        # add dataset info to obs
        obs['dataset'] = 'NA'
        obs.loc[obs['ind_cov'].isin(clues_ids), 'dataset'] = 'CLUES'
        obs.loc[obs['ind_cov'].isin(immvar_ids), 'dataset'] = 'IMMVAR'
        assert obs.loc[obs['dataset'] == 'NA'].shape[0] == 0
        obs.to_csv(output.obs, sep='\t', index=False)


use rule yazar_exclude_repeatedpool as perez_exclude_repeatedpool with:
    input:
        obs = 'analysis/perez/data/{anc}.{status}.obs.dataset.gz',
    output:
        # obs = 'staging/data/perez/{anc}.{status}.exclude_repeatedpool.obs.gz',
        obs = 'analysis/perez/data/{anc}.{status}.exclude_repeatedpool.obs.gz',
        dup_inds = 'analysis/perez/data/{anc}.{status}.duplicated_inds.txt',
    params:
        ind_col = perez_ind_col,
        pool_col = perez_pool_col,


use rule yazar_ctp_transform as perez_ctp_transform with:
    input:
        h5ad = 'analysis/perez/data/{anc}.{status}.h5ad',
        var = 'analysis/perez/data/{anc}.{status}.var.gz',
        obs = 'analysis/perez/data/{anc}.{status}.exclude_repeatedpool.obs.gz',
    output:
        X = 'staging/data/perez/{anc}.{status}.X.npz',
        obs = 'staging/data/perez/{anc}.{status}.obs.gz',
        var = 'staging/data/perez/{anc}.{status}.var.gz',
    params:
        ind_col = perez_ind_col,
        ct_col = perez_ct_col,
        pool_col = perez_pool_col,


use rule perez_gene_location_allgenes as perez_gene_location with:
    input:
        genes = 'staging/data/perez/{anc}.{status}.var.gz',
        gene_info = 'analysis/data/Homo_sapiens.GRCh37.75.txt',
    output:
        location = 'analysis/perez/data/{anc}.{status}.gene_location.txt',
        meta = 'analysis/perez/data/{anc}.{status}.gene_meta.txt',


use rule yazar_ctp as perez_ctp with:
   input:
        X = 'staging/data/perez/{anc}.{status}.X.npz',
        obs = 'staging/data/perez/{anc}.{status}.obs.gz',
        var = 'staging/data/perez/{anc}.{status}.var.gz',
   output:
        ctp = 'analysis/perez/data/{anc}.{status}.ctp.gz',
        ctnu = 'analysis/perez/data/{anc}.{status}.ctnu.gz',
        P = 'analysis/perez/data/{anc}.{status}.P.gz',
        n = 'analysis/perez/data/{anc}.{status}.n.gz',
   params:
       ind_col = perez_ind_col,
       ct_col = perez_ct_col,


use rule yazar_rm_rareINDnCT as perez_rm_rareINDnCT with:
    # also select gene expressed in all cts
    input:
        ctp = 'analysis/perez/data/{anc}.{status}.ctp.gz',
        ctnu = 'analysis/perez/data/{anc}.{status}.ctnu.gz',
        n = 'analysis/perez/data/{anc}.{status}.n.gz',
    output:
        ctp = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/ctp.gz',
        ctnu = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/ctnu.gz',
        P = f'analysis/perez/{perez_paramspace.wildcard_pattern}/mvn/{{anc}}.{{status}}/P.final.gz',
        n = f'analysis/perez/{perez_paramspace.wildcard_pattern}/mvn/{{anc}}.{{status}}/n.final.gz',


use rule yazar_rm_missingIND as perez_rm_missingIND with:
    input:
        ctp = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/ctp.gz',
        ctnu = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/ctnu.gz',
        P = f'analysis/perez/{perez_paramspace.wildcard_pattern}/mvn/{{anc}}.{{status}}/P.final.gz',
        n = f'analysis/perez/{perez_paramspace.wildcard_pattern}/mvn/{{anc}}.{{status}}/n.final.gz',
    output:
        ctp = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/ctp.rmid.gz',
        ctnu = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/ctnu.rmid.gz',
        P = f'analysis/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/P.final.gz',
        n = f'analysis/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/n.final.gz',


use rule yazar_std_op as perez_std_op with:
    input:
        ctp = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/ctp.rmid.gz',
        ctnu = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/ctnu.rmid.gz',
        P = f'analysis/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/P.final.gz',
    output:
        op = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/op.std.gz',
        nu = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/nu.std.gz',
        ctp = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/ctp.std.gz',
        ctnu = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/ctnu.std.gz',


use rule yazar_exclude_sexchr as perez_exclude_sexchr with:
    input:
        genes = 'analysis/perez/data/{anc}.{status}.gene_location.txt',
        op = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/op.std.gz',
        nu = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/nu.std.gz',
        ctp = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/ctp.std.gz',
        ctnu = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/ctnu.std.gz',
    output:
        op = f'analysis/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/op.final.gz',
        nu = f'analysis/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/nu.final.gz',
        ctp = f'analysis/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/ctp.final.gz',
        ctnu = f'analysis/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/ctnu.final.gz',


use rule yazar_op_pca as perez_op_pca with:
    input:
        op = f'analysis/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/op.final.gz',
    output:
        evec = f'analysis/perez/{yazar_paramspace.wildcard_pattern}/{{anc}}.{{status}}/evec.gz',
        eval = f'analysis/perez/{yazar_paramspace.wildcard_pattern}/{{anc}}.{{status}}/eval.gz',
        pca = f'analysis/perez/{yazar_paramspace.wildcard_pattern}/{{anc}}.{{status}}/pca.gz',
        png = f'results/perez/{yazar_paramspace.wildcard_pattern}/{{anc}}.{{status}}/pca.png',


use rule yazar_geno_pca as perez_geno_pca with:
    input:
        bed = expand('analysis/perez/data/geno/{{anc}}.{{status}}/chr{chr}.bed',
                chr=range(1,23)),
        P = f'analysis/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/P.final.gz',
    output:
        eigenvec = f'analysis/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/geno.eigenvec',
        eigenval = f'analysis/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/geno.eigenval',


use rule yazar_geno_pca_plot as perez_geno_pca_plot with:
    input:
        eigenval = f'analysis/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/geno.eigenval',
    output:
        png = f'results/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/geno.eigenval.png',


use rule yazar_he_kinship as perez_he_kinship with:
    input:
        ctp = f'analysis/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/ctp.final.gz',
        genes = 'analysis/perez/data/{anc}.{status}.gene_location.txt',
        bed = 'analysis/perez/data/geno/{anc}.{status}/chr{chr}.bed',
    output:
        kinship = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/he/kinship.chr{{chr}}.npy',


use rule yazar_HE_split as perez_HE_split with:
    input:
        ctp = f'analysis/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/ctp.final.gz',
        ctnu = f'analysis/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/ctnu.final.gz',
        genes = 'analysis/perez/data/{anc}.{status}.gene_location.txt',
    output:
        ctp = [f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/he/ctp.batch{i}.gz'
                for i in range(perez_he_batches)],
        ctnu = [f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/he/ctnu.batch{i}.gz'
                for i in range(perez_he_batches)],
        batch = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/he/ctp.batch.txt',


use rule yazar_he_kinship_split as perez_he_kinship_split with:
    input:
        kinship = [f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/he/kinship.chr{chr}.npy'
                    for chr in chrs],
        ctp = [f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/he/ctp.batch{i}.gz'
                for i in range(perez_he_batches)],
        batch = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/he/ctp.batch.txt',
    output:
        kinship = [temp(f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/he/kinship.batch{i}.npy')
                    for i in range(perez_he_batches)],


rule perez_HE_free:
    input:
        ctp = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/he/ctp.batch{{i}}.gz',
        ctnu = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/he/ctnu.batch{{i}}.gz',
        kinship = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/he/kinship.batch{{i}}.npy',
        P = f'analysis/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/P.final.gz',
        op_pca = f'analysis/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/pca.gz',
        geno_pca = f'analysis/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/geno.eigenvec',
        obs = 'analysis/perez/data/{anc}.{status}.exclude_repeatedpool.obs.gz',
    output:
        out = f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/he.batch{{i}}.npy',
    params:
        snps = 5, # threshold of snp number per gene
        iid = False,
        full = False,
        jk = True,
    script: '../bin/perez/he.py' 


use rule yazar_HE_free_merge as perez_HE_free_merge with:
    input:
        out = [f'staging/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/he.batch{i}.npy'
            for i in range(perez_he_batches)],
    output:
        out = f'analysis/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/he.npy',


use rule yazar_HE_togz as perez_HE_free_togz with:
    input:
        out = f'analysis/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/he.npy',
        P = f'analysis/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/P.final.gz',
    output:
        out = f'results/perez/{perez_paramspace.wildcard_pattern}/{{anc}}.{{status}}/he.gz',


rule perez_HE_Free_meta:
    input:
        P = [f'analysis/perez/{perez_paramspace.wildcard_pattern}/{anc}.{status}/P.final.gz'
                    for anc, status in zip(['European', 'European', 'Asian'], ['Healthy', 'SLE', 'SLE'])],
        out = [f'analysis/perez/{perez_paramspace.wildcard_pattern}/{anc}.{status}/he.npy' 
                    for anc, status in zip(['European', 'European', 'Asian'], ['Healthy', 'SLE', 'SLE'])],
    output:
        out = f'analysis/perez/{perez_paramspace.wildcard_pattern}/he.free.meta.npy', 
    script: '../bin/perez/he.free.meta.py'


rule perez_HE_Free_meta_togz:
    input:
        out = f'analysis/perez/{perez_paramspace.wildcard_pattern}/he.free.meta.npy',
        P = f'analysis/perez/{perez_paramspace.wildcard_pattern}/European.SLE/P.final.gz',
    output:
        out = f'results/perez/{perez_paramspace.wildcard_pattern}/he.free.meta.gz',
    script: '../bin/perez/he.togz.py'


use rule yazar_eds as perez_eds with:
    input:
        var = 'data/Perez2022Science/adjusted.all.var.gz',
        eds = 'data/Wang_Goldstein.tableS1.txt',
        blood_connected = 'data/gene_annots/saha_et_al/twns/processed/WholeBlood.degree.txt',
        combined_connected = 'data/gene_annots/saha_et_al/twns/processed/genes_by_combined_rank.txt',
        gnomad = 'data/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz',
    output:
        eds = 'analysis/perez/eds.txt',










rule perez_all:
    input:
        all_all_pca = expand('results/perez/{params}/all.all/geno.eigenval.png',
                        params=perez_paramspace.instance_patterns, ),
        eur_case_pca = expand('results/perez/{params}/European.SLE/geno.eigenval.png',
                        params=perez_paramspace.instance_patterns, ),
        eur_ctrl_pca = expand('results/perez/{params}/European.Healthy/geno.eigenval.png',
                        params=perez_paramspace.instance_patterns, ),
        asia_case_pca = expand('results/perez/{params}/Asian.SLE/geno.eigenval.png',
                        params=perez_paramspace.instance_patterns, ),
        cis_eds = expand('results/perez/{params}/matrix.cor.main.png',
                    params=perez_paramspace.instance_patterns), 
