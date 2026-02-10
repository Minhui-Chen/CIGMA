##########################################################################
# coding heritability
# variance effect https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html
##########################################################################
rule variant_annotation_download:
    output:
        vanno = 'data/benneale/variants.tsv.gz',
        vanno2 = 'analysis/data/variants.tsv.gz',
    params:
        link = 'https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz',
    resources: mem_mb = '1G'
    shell:
        '''
        mkdir -p $(dirname {output.vanno})
        wget -O - {params.link} | gunzip -c | gzip -c > {output.vanno}
        zcat {output.vanno} | awk '$2 != "X"' | gzip -c > {output.vanno2}
        '''

rule yazar_ldsc_coding_make_annot:
    input:
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/he.free.jk.npy',
        op = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/op.rmid.gz',
        location = 'data/Yazar2022Science/gene_location.txt',
        vanno = 'analysis/data/variants.tsv.gz',
        bim = 'data/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.{chr}.bim',
    output:
        var_annot = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/window_{{window}}/var.{{chr}}.annot.gz',
        mean_annot = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/window_{{window}}/mean.{{chr}}.annot.gz',
    params: 
        coding_variants = ['splice_donor_5th_base_variant', 'missense_variant', 'splice_region_variant', 
                            'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 
                            'start_lost', 'stop_lost', 'frameshift_variant',   
                            'protein_altering_variant'],
    resources: 
        mem_mb = '8G',
    script: '../bin/yazar/ldsc_coding_make_annot.py'


rule yazar_ldsc_coding_compldscore:
    input:
        var = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/window_{{window}}/var.{{chr}}.annot.gz',
        mean = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/window_{{window}}/mean.{{chr}}.annot.gz',
        bim = 'data/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.{chr}.bim',
        hapmap3 = 'data/ldsc/hm3_no_MHC.list.txt',
    output:
        var = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/window_{{window}}/var.{{chr}}.l2.ldscore.gz',
        mean = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/window_{{window}}/mean.{{chr}}.l2.ldscore.gz',
    params:
        bfile = lambda wildcards, input: os.path.splitext(input.bim)[0],
        var = lambda wildcards, output: re.sub('\.l2.ldscore.gz$', '', output.var),
        mean = lambda wildcards, output: re.sub('\.l2.ldscore.gz$', '', output.mean),
    conda: '../../ldsc/environment.yml'
    # resources: 
        # partition = 'tier2q',
    shell:
        '''
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
        '''


rule yazar_ldsc_coding_seg_regression:
    input:
        gwas = 'staging/data/ldsc/{gwas}.sumstats.gz',
        var = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/window_{{window}}/var.{chr}.l2.ldscore.gz'
                    for chr in chrs],
        mean = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/window_{{window}}/mean.{chr}.l2.ldscore.gz'
                    for chr in chrs],
        ref_ld = expand('data/ldsc/baseline_v1.2/baseline.{chr}.l2.ldscore.gz', chr=chrs),  # readme_baseline_versions: use baselineLD v2.2 for estimating heritability enrichment; baseline v1.2 for identifying critical tissues/cell-types via P-value of tau  
        weight = expand('data/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.{chr}.l2.ldscore.gz', chr=chrs),
        frq = expand('data/ldsc/1000G_Phase3_frq/1000G.EUR.QC.{chr}.frq', chr=chrs),
    output:
        var = f'results/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/window_{{window}}/var.{{gwas}}.results',
        var_log = f'results/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/window_{{window}}/var.{{gwas}}.log',
        mean = f'results/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/window_{{window}}/mean.{{gwas}}.results',
        mean_log = f'results/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/window_{{window}}/mean.{{gwas}}.log',
    params:
        var1 = lambda wildcards, input: re.sub('1.l2.ldscore.gz', '', input.var[0]),
        mean1 = lambda wildcards, input: re.sub('1.l2.ldscore.gz', '', input.mean[0]),
        ref_ld = lambda wildcards, input: re.sub('1.l2.ldscore.gz', '', input.ref_ld[0]),
        weight = lambda wildcards, input: re.sub('1.l2.ldscore.gz', '', input.weight[0]),
        frq = lambda wildcards, input: re.sub('1.frq$', '', input.frq[0]),
        var = lambda wildcards, output: re.sub('\.results', '', output.var),
        mean = lambda wildcards, output: re.sub('\.results', '', output.mean),
    conda: '../../ldsc/environment.yml'
    # resources: 
        # partition = 'tier2q',
    shell:
        '''
        ldsc/ldsc.py \
            --h2 {input.gwas} \
            --ref-ld-chr {params.var1},{params.ref_ld} \
            --w-ld-chr {params.weight} \
            --overlap-annot \
            --frqfile-chr {params.frq} \
            --print-coefficients \
            --out {params.var} 
        
        ldsc/ldsc.py \
            --h2 {input.gwas} \
            --ref-ld-chr {params.mean1},{params.ref_ld} \
            --w-ld-chr {params.weight} \
            --overlap-annot \
            --frqfile-chr {params.frq} \
            --print-coefficients \
            --out {params.mean}
        '''


rule yazar_ldsc_coding_seg_regression_agg1:
    input:
        var = expand('results/yazar/{params}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/window_{{window}}/var.{gwas}.results',
                        params=yazar_paramspace.instance_patterns,
                        gwas=traits),
        var_log = expand('results/yazar/{params}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/window_{{window}}/var.{gwas}.log',
                        params=yazar_paramspace.instance_patterns,
                        gwas=traits),
        mean = expand('results/yazar/{params}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/window_{{window}}/mean.{gwas}.results',
                        params=yazar_paramspace.instance_patterns,
                        gwas=traits),
        mean_log = expand('results/yazar/{params}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/window_{{window}}/mean.{gwas}.log',
                        params=yazar_paramspace.instance_patterns,
                        gwas=traits),
    output:
        var = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/window_{{window}}/var.results',
        mean = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/window_{{window}}/mean.results',
    params:
        traits = traits,
    resources: 
        mem_mb = '1G',
    run:
        dfs_var = []
        for f1, f2, gwas in zip(input.var, input.var_log, params.traits):
            df = pd.read_table(f1)
            df = df.assign(gwas=gwas)
            with open(f2, 'r') as logf:
                for line in logf:
                    if line.startswith('Total Observed scale h2'):
                        h2 = float(line.strip().split()[-2])
                        df = df.assign(h2=h2)
            dfs_var.append(df)
        df_var = pd.concat(dfs_var)
        df_var.to_csv(output.var, sep='\t', index=False)


        dfs_mean = []
        for f1, f2, gwas in zip(input.mean, input.mean_log, params.traits):
            df = pd.read_table(f1)
            df = df.assign(gwas=gwas)
            with open(f2, 'r') as logf:
                for line in logf:
                    if line.startswith('Total Observed scale h2'):
                        h2 = float(line.strip().split()[-2])
                        df = df.assign(h2=h2)
            dfs_mean.append(df)
        df_mean = pd.concat(dfs_mean)
        df_mean.to_csv(output.mean, sep='\t', index=False)


rule yazar_ldsc_coding_seg_regression_agg2:
    input:
        var = expand('staging/yazar/{params}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/window_{window}/var.results',
                        params=yazar_paramspace.instance_patterns, window=[500000]),
        mean = expand('staging/yazar/{params}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/window_{window}/mean.results',
                        params=yazar_paramspace.instance_patterns, window=[500000]),
    output:
        var = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/var.results',
        mean = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_coding/genefeature_{{feature}}/nbin_{{n}}/mean.results',
    params:
        window = [500000],
    resources: 
        # partition = 'tier2q',
        mem_mb = '1G',
    run:
        dfs_var = [pd.read_table(f).assign(window=window) for f, window in zip(input.var, params.window)]
        df_var = pd.concat(dfs_var)
        df_var.to_csv(output.var, sep='\t', index=False)

        dfs_mean = [pd.read_table(f).assign(window=window) for f, window in zip(input.mean, params.window)]
        df_mean = pd.concat(dfs_mean)
        df_mean.to_csv(output.mean, sep='\t', index=False)


rule yazar_ldsc_coding_seg_regression_agg3:
    input:
        var = expand('staging/yazar/{params}/ldsc_coding/genefeature_{{feature}}/nbin_{n}/var.results',
                        params=yazar_paramspace.instance_patterns, n=[5]),
        mean = expand('staging/yazar/{params}/ldsc_coding/genefeature_{{feature}}/nbin_{n}/mean.results',
                        params=yazar_paramspace.instance_patterns, n=[5]),
    output:
        var = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_coding/genefeature_{{feature}}/var.results',
        mean = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_coding/genefeature_{{feature}}/mean.results',
    params:
        n = [5],
    resources: 
        # partition = 'tier2q',
        mem_mb = '1G',
    run:
        dfs_var = [pd.read_table(f).assign(nbin=n) for f, n in zip(input.var, params.n)]
        df_var = pd.concat(dfs_var)
        df_var.to_csv(output.var, sep='\t', index=False)

        dfs_mean = [pd.read_table(f).assign(nbin=n) for f, n in zip(input.mean, params.n)]
        df_mean = pd.concat(dfs_mean)
        df_mean.to_csv(output.mean, sep='\t', index=False)


tmp_gene_features = ['specificity', 'hom_g2', 'v', 'var_beta']
rule yazar_ldsc_coding_seg_regression_agg4:
    input:
        var = expand('staging/yazar/{params}/ldsc_coding/genefeature_{feature}/var.results',
                      params=yazar_paramspace.instance_patterns, feature=tmp_gene_features),
        mean = expand('staging/yazar/{params}/ldsc_coding/genefeature_{feature}/mean.results',
                        params=yazar_paramspace.instance_patterns, feature=tmp_gene_features),
    output:
        var = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_coding/var.results',
        mean = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_coding/mean.results',
    params:
        feature = tmp_gene_features,
    resources: 
        # partition = 'tier2q',
        mem_mb = '1G',
    run:
        dfs_var = [pd.read_table(f).assign(feature=n) for f, n in zip(input.var, params.feature)]
        df_var = pd.concat(dfs_var)
        df_var.to_csv(output.var, sep='\t', index=False)

        dfs_mean = [pd.read_table(f).assign(feature=n) for f, n in zip(input.mean, params.feature)]
        df_mean = pd.concat(dfs_mean)
        df_mean.to_csv(output.mean, sep='\t', index=False)


rule yazar_ldsc_coding_all:
    input:
        var = expand('analysis/yazar/{params}/ldsc_coding/var.results',
                        params=yazar_paramspace.instance_patterns),










































##########################################################################
# More traits for LDSC (used in Cuomo 2025, from Margoliash)
##########################################################################
biomarkers = ['alanine_aminotransferase', 'albumin', 'alkaline_phosphatase', 'apolipoprotein_a', 'apolipoprotein_b', 'aspartate_aminotransferase',
               'c_reactive_protein', 'calcium', 'cholesterol', 'creatinine', 'cystatin_c', 'gamma_glutamyltransferase', 'glucose', 
               'glycated_haemoglobin', 'hdl_cholesterol', 'igf_1', 'ldl_cholesterol_direct', 'phosphate', 'shbg', 'total_bilirubin', 
               'total_protein', 'triglycerides', 'urate', 'urea', 'vitamin_d']
blood_cell_counts = ['eosinophil_count', 'eosinophil_percent', 'haematocrit', 'haemoglobin_concentration', 'lymphocyte_count', 'lymphocyte_percent', 
                'mean_corpuscular_haemoglobin', 'mean_corpuscular_haemoglobin_concentration', 'mean_corpuscular_volume', 'mean_platelet_volume',
                'mean_sphered_cell_volume', 'neutrophil_count', 'neutrophil_percent', 'platelet_count', 'platelet_crit', 'platelet_distribution_width',
                'red_blood_cell_count', 'red_blood_cell_distribution_width', 'white_blood_cell_count']
moregwas = biomarkers + blood_cell_counts

rule yazar_ldscmore_downloadgwas:
    output:
        gwas = 'data/ldsc/gwas/{gwas}.gz',
    params:
        link = lambda wildcards: re.sub('phenotype',wildcards.gwas, 
        'https://margoliash-et-al-2023.s3.amazonaws.com/associations/white_british_phenotype_snp_gwas_results.tab.gz'),
    shell:
        'wget {params.link} -O {output.gwas}'


rule yazar_ldscmore_gwas_addP:
    input:
        gwas = 'data/ldsc/gwas/{gwas}.gz',
    output:
        gwas = 'staging/data/ldscmore/{gwas}.gz',
    resources:
        partition = 'tier3q',
        mem_mb = '60G',
    run:
        gwas = pd.read_table(input.gwas, comment=None)

        # remove error rows
        gwas = gwas[gwas['ERRCODE'] == '.']

        # rename columns: #CHROM to CHROM
        gwas.rename(columns={'#CHROM':'CHROM', 'ID':'SNP', 'OBS_CT':'N', 'T_STAT':'Z'}, inplace=True)

        # add A2
        gwas['A2'] = gwas['REF']
        gwas.loc[gwas['REF'] == gwas['A1'], 'A2'] = gwas.loc[gwas['REF'] == gwas['A1'], 'ALT']

        # make A1 trait-incresing
        mask = gwas['BETA'] < 0
        tmp = gwas.loc[mask, 'A1'].copy()
        gwas.loc[mask, 'A1'] = gwas.loc[mask, 'A2']
        gwas.loc[mask, 'A2'] = tmp
        gwas.loc[mask, 'BETA'] = -gwas.loc[mask, 'BETA']
        gwas.loc[mask, 'Z'] = -gwas.loc[mask, 'Z']
        if np.all(gwas['Z'] >= 0):
            pass
        else:
            print(gwas[~(gwas['Z'] >= 0)])
            sys.exit()

        gwas.to_csv(output.gwas, sep='\t', index=False)


use rule yazar_ldsc_format_gwas as yazar_ldscmore_format_gwas with:
    input:
        gwas = 'staging/data/ldscmore/{gwas}.gz',
        hapmap3 = 'data/ldsc/broad/GINGER/ginger_vc_c1_year3/lab_11.13.19/w_hm3.snplist',  # snp list with A1 and A2
    output:
        gwas = 'staging/data/ldscmore/{gwas}.sumstats.gz',


use rule yazar_ldsc_seg_regression as yazar_ldscmore_seg_regression with:
    input:
        gwas = 'staging/data/ldscmore/{gwas}.sumstats.gz',
        ldcts = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/he.ldcts',
        ref_ld = expand('data/ldsc/baseline_v1.2/baseline.{chr}.l2.ldscore.gz', chr=chrs),  # readme_baseline_versions: use baselineLD v2.2 for estimating heritability enrichment; baseline v1.2 for identifying critical tissues/cell-types via P-value of tau  
        weight = expand('data/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.{chr}.l2.ldscore.gz', chr=chrs),
    output:
        out = f'results/yazar/{yazar_paramspace.wildcard_pattern}/ldscmore/top_{{ngene}}/window_{{window}}/he.{{gwas}}.cell_type_results.txt',


rule yazar_ldscmore_partition_stacked_summary:
    input:
        stacked = [f'results/yazar/{yazar_paramspace.wildcard_pattern}/ldscmore/top_{{ngene}}/window_{{window}}/he.{gwas}.cell_type_results.txt'
                    for gwas in moregwas],
    output:
        png = f'results/yazar/{yazar_paramspace.wildcard_pattern}/ldscmore/top_{{ngene}}/window_{{window}}/he.h2_cts.png',
    params:
        traits1 = biomarkers,
        traits2 = blood_cell_counts,
    script: '../scripts/yazar/ldsc.moregwas.h2_cts.py'



rule yazar_ldscmore_all:
    input:
        h2_cts = expand('results/yazar/{params}/ldscmore/top_{ngene}/window_{window}/he.h2_cts.png',
                params=yazar_paramspace.instance_patterns, 
                ngene = [0, 100, 200, 300],
                window=[300000, 500000, 700000]),
















































##########################################################################
# LDSC continuous annotation 
# TODO
##########################################################################
rule yazar_ldsc_make_continuous_annot:
    input:
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/he.free.jk.npy',
        location = 'data/Yazar2022Science/gene_location.txt',
        gcta = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/gcta/op.greml',
        bim = 'data/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.{chr}.bim',
    output:
        # shared = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/shared.{{chr}}.annot.gz',
        var = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_continuous/window_{{window}}/var.{{chr}}.annot.gz',
        # mean = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/mean.{{chr}}.annot.gz',
        # gcta = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/gcta.{{chr}}.annot.gz',
    params:
        chr = config['ldsc']['mhc']['chr'],
        start = config['ldsc']['mhc']['start'],
        end = config['ldsc']['mhc']['end'],
    script: '../bin/yazar/ldsc_make_continuous_annot.py'


rule yazar_ldsc_continuous_compldscore:
    input:
        # all = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/all.{{chr}}.annot.gz',
        # random = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/random.{{chr}}.annot.gz',
        # shared = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/shared.{{chr}}.annot.gz',
        var = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_continuous/window_{{window}}/var.{{chr}}.annot.gz',
        # mean = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/mean.{{chr}}.annot.gz',
        # gcta = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/gcta.{{chr}}.annot.gz',
        bim = 'data/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.{chr}.bim',
        hapmap3 = 'data/ldsc/hm3_no_MHC.list.txt',
    output:
        # all = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/all.{{chr}}.l2.ldscore.gz',
        # random = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/random.{{chr}}.l2.ldscore.gz',
        # shared = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/shared.{{chr}}.l2.ldscore.gz',
        var = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_continuous/window_{{window}}/var.{{chr}}.l2.ldscore.gz',
        # mean = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/mean.{{chr}}.l2.ldscore.gz',
        # gcta = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc/top_{{ngene}}/window_{{window}}/gcta.{{chr}}.l2.ldscore.gz',
    params:
        bfile = lambda wildcards, input: os.path.splitext(input.bim)[0],
        # all = lambda wildcards, output: re.sub('\.l2.ldscore.gz$', '', output.all),
        # random = lambda wildcards, output: re.sub('\.l2.ldscore.gz$', '', output.random),
        # shared = lambda wildcards, output: re.sub('\.l2.ldscore.gz$', '', output.shared),
        var = lambda wildcards, output: re.sub('\.l2.ldscore.gz$', '', output.var),
        # mean = lambda wildcards, output: re.sub('\.l2.ldscore.gz$', '', output.mean),
        # gcta = lambda wildcards, output: re.sub('\.l2.ldscore.gz$', '', output.gcta),
    conda: '../../ldsc/environment.yml'
    shell:
        '''
        python ldsc/ldsc.py \
            --l2 \
            --bfile {params.bfile} \
            --ld-wind-cm 1 \
            --annot {input.var} \
            --out {params.var} \
            --print-snps {input.hapmap3}
        '''


rule yazar_ldsc_continuous_partition:
    input:
        gwas = 'staging/data/ldsc/{gwas}.sumstats.gz',
        ref_ld = expand('data/ldsc/1000G_Phase3_baselineLD_v2.2/baseline.{chr}.l2.ldscore.gz', chr=chrs),  # readme_baseline_versions: use baselineLD v2.2 for estimating heritability enrichment; baseline v1.2 for identifying critical tissues/cell-types via P-value of tau  
        weight = expand('data/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.{chr}.l2.ldscore.gz', chr=chrs),
        frq = expand('data/ldsc/1000G_Phase3_frq/1000G.EUR.QC.{chr}.frq', chr=chrs),
        var = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_continuous/window_{{window}}/var.{chr}.l2.ldscore.gz'
                for chr in chrs],
    output:
        var = f'results/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_continuous/window_{{window}}/he.var.{{gwas}}.results',
    params:
        ref_ld = lambda wildcards, input: re.sub('1.l2.ldscore.gz', '', input.ref_ld[0]),
        weight = lambda wildcards, input: re.sub('1.l2.ldscore.gz', '', input.weight[0]),
        frq = lambda wildcards, input: re.sub('1.frq', '', input.frq[0]),
        var = lambda wildcards, input: re.sub('1.l2.ldscore.gz', '', input.var[0]),
        var_out = lambda wildcards, output: re.sub('\.results', '', output.var),
    conda: '../../ldsc/environment.yml'
    shell:
        '''
        ldsc/ldsc.py \
            --h2 {input.gwas} \
            --ref-ld-chr {params.var},{params.ref_ld} \
            --w-ld-chr {params.weight} \
            --overlap-annot \
            --frqfile-chr {params.frq} \
            --print-coefficients \
            --out {params.var_out} \
        '''


rule yazar_ldsc_continuous_partition_agg:
    input:
        var = [f'results/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_continuous/window_{{window}}/he.var.{gwas}.results'
                    for gwas in traits],
    output:
        ldsc = f'results/yazar/{yazar_paramspace.wildcard_pattern}/ldsc_continuous/window_{{window}}/ldsc.gz',
    params:
        traits = traits,
    run:
        data = []
        for gwas, f in zip(params.traits, input.var):
            # TODO check file results file
            df = pd.read_table(f).assign(trait=gwas)
            data.append(df)
        data = pd.concat(data, ignore_index=True)
        data.to_csv(output.ldsc, sep='\t', index=False)


rule yazar_ldsc_continuous_all:
    input:
        ldsc = expand('results/yazar/{params}/ldsc_continuous/window_{window}/ldsc.gz',
                        params=yazar_paramspace.instance_patterns, 
                        window=[500000]),








###########################################################
# For MIRA
###########################################################
##########################################################################
## 1.2: combine cts 
##########################################################################


combine_cts = {
        'cc1':  {
            'B': ['B IN', 'B Mem'],
            'T': ['CD4 ET', 'CD4 NC', 'CD8 ET', 'CD8 NC'],
            },
        }

wildcard_constraints: cc='[\w\d]+'


rule yazar_combine_cts:
    input:
        obs = 'staging/data/yazar/obs.gz',
    output:
        obs = 'analysis/data/yazar/combine_cts/{cc}/obs.gz',
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
        obs = 'analysis/data/yazar/combine_cts/{cc}/obs.gz',
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


use rule yazar_rm_missingIND as yazar_cc_rm_missingIND with:
    input:
        ctp = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.gz',
        ctnu = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.gz',
        P = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.gz',
        n = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/n.gz',
    output:
        ctp = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.rmid.gz',
        ctnu = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.rmid.gz',
        P = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.final.gz',
        n = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/n.final.gz',


use rule yazar_std_op as yazar_cc_std_op with:
    input:
        ctp = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.rmid.gz',
        ctnu = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.rmid.gz',
        P = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.final.gz',
    output:
        op = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/op.std.gz',
        nu = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/nu.std.gz',
        ctp = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.std.gz',
        ctnu = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.std.gz',


use rule yazar_exclude_sexchr as yazar_cc_exclude_sexchr with:
    input:
        genes = 'data/Yazar2022Science/gene_location.txt',
        op = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/op.std.gz',
        nu = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/nu.std.gz',
        ctp = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.std.gz',
        ctnu = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.std.gz',
    output:
        op = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/op.final.gz',
        nu = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/nu.final.gz',
        ctp = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.final.gz',
        ctnu = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.final.gz',


use rule yazar_op_pca as yazar_cc_op_pca with:
    input:
        op = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/op.final.gz',
    output:
        evec = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/evec.gz',
        eval = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/eval.gz',
        pca = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/pca.gz',
        png = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/op.pca.png',


use rule yazar_geno_pca as yazar_cc_geno_pca with:
    input:
        bed = expand('analysis/yazar/data/geno/chr{chr}.bed',
                chr=range(1,23)),
        P = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.final.gz',
    output:
        eigenvec = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        eigenval = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/geno.eigenval',


use rule yazar_he_kinship as yazar_cc_he_kinship with:
    input:
        ctp = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.final.gz',
        genes = 'data/Yazar2022Science/gene_location.txt',
        bed = 'analysis/yazar/data/geno/chr{chr}.bed',
    output:
        kinship = temp(f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/kinship.chr{{chr}}.npy'),


network = ['STAT5B', 'STAT5A', 'GATA3', 'IRF4', 'KMT2A', 'JAK3', 'FOXP1', 'PTEN', 'ETS1', 'RELA', 'MBD2', 'YY1', 'MYB', 
           'IRF2', 'HIVEP2', 'ZNF217', 'ATXN7L3', 'FOXK1', 'KLF2', 'TNFAIP3', 'CBFB']
rule yazar_kinship_network:
    input:
        genes = 'data/Yazar2022Science/gene_location.txt',
        bed = expand('analysis/yazar/data/geno/chr{chr}.bed', chr=chrs),
        ctp = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.final.gz',
    output:
        kinship = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/network.kinship.npy',
    resources:
        partition = 'tier3q',
        mem_mb = '80G',
    params:
        bed = lambda wildcards, input: ','.join(input.bed),
        chrs = ','.join([str(chr) for chr in chrs]),
        genes = ','.join(network),
        r = int(float(config['yazar']['radius'])),
    shell:
        '''
        {config[plink_load]}

        python3 workflow/scripts/yazar/kinship.py \
            {input.genes} \
            {params.genes} \
            {params.r} \
            {params.bed} \
            {params.chrs} \
            {input.ctp} \
            {output.kinship}
        '''


rule yazar_n2g:
    input:
        ctp = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.final.gz',
        ctnu = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.final.gz',
        kinship = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/network.kinship.npy',
        P = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.final.gz',
        op_pca = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/pca.gz',
        geno_pca = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
        genes = 'data/Yazar2022Science/gene_location.txt',
    output:
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/n2g.npy',
    params:
        gene = 'CTLA4',
    resources:
        partition = 'tier2q',
        mem_mb = '16G',
    script: '../scripts/yazar/network.py'


use rule yazar_kinship_network as yazar_kinship_gene with:
    input:
        genes = 'data/Yazar2022Science/gene_location.txt',
        bed = expand('analysis/yazar/data/geno/chr{chr}.bed', chr=chrs),
        ctp = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.final.gz',
    output:
        kinship = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/gene.kinship.npy',
    resources:
        partition = 'tier3q',
        mem_mb = '80G',
    params:
        bed = lambda wildcards, input: ','.join(input.bed),
        chrs = ','.join([str(chr) for chr in chrs]),
        genes = 'IL2RA',
        r = int(float(config['yazar']['radius'])),


use rule yazar_n2g as yazar_g2n with:
    input:
        ctp = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.final.gz',
        ctnu = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.final.gz',
        kinship = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/gene.kinship.npy',
        P = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.final.gz',
        op_pca = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/pca.gz',
        geno_pca = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
        genes = 'data/Yazar2022Science/gene_location.txt',
    output:
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/regulator.{{gene}}.npy',



rule yazar_g2n_all:
    input:
        n2g = expand('analysis/yazar/combine_cts/{cc}/{params}/n2g.npy',
                      cc=['cc1'], params=yazar_paramspace.instance_patterns),
        g2n = expand('staging/yazar/combine_cts/{cc}/{params}/regulator.{gene}.npy',
                      cc=['cc1'], params=yazar_paramspace.instance_patterns,
                      gene=['ZNF217', 'KLF2', 'TNFAIP3', 'GATA3']),


rule sim_n2g:
    input:
        V = 'analysis/sim/free4/AGGvc.true_V.npy',
        out = 'analysis/sim/free4/L_20_nL_80/AGGvc.he.npy',
        out_pergene = 'analysis/sim/free41/L_20_nL_80/AGGvc.he.npy',


################ trans-PCO gene modules ##########################
rule yazar_kinship_module:
    input:
        modules = 'ngx/transpco_tableS1.csv',
        genes = 'data/Yazar2022Science/gene_location.txt',
        bed = expand('analysis/yazar/data/geno/chr{chr}.bed', chr=chrs),
        ctp = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.final.gz',
    output:
        kinship = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/{{gene}}.module{{i}}.kinship.npy',
    resources:
        partition = 'tier3q',
        mem_mb = '80G',
    params:
        bed = lambda wildcards, input: ','.join(input.bed),
        chrs = ','.join([str(chr) for chr in chrs]),
        genes = lambda wildcards, input: pd.read_csv(input.modules).loc[
            pd.read_csv(input.modules)['gene_module'] == int(wildcards.i), 'genes_in_module'
            ].values[0].replace(';', ','),
        r = int(float(config['yazar']['radius'])),
    shell:
        '''
        {config[plink_load]}

        python3 workflow/scripts/yazar/kinship.rmchr.py \
            {input.genes} \
            {params.genes} \
            {params.r} \
            {params.bed} \
            {params.chrs} \
            {input.ctp} \
            {wildcards.gene} \
            {output.kinship}
        '''


use rule yazar_n2g as yazar_module_n2g with:
    input:
        ctp = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctp.final.gz',
        ctnu = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/ctnu.final.gz',
        kinship = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/{{gene}}.module{{i}}.kinship.npy',
        P = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/P.final.gz',
        op_pca = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/pca.gz',
        geno_pca = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
        genes = 'data/Yazar2022Science/gene_location.txt',
    output:
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/n2g.{{gene}}.module{{i}}.npy',
    params:
        jk = True,


rule yazar_module_n2g_add_module:
    input:
        out = [f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/n2g.{{gene}}.module{i}.npy' for i in range(1,167)],
    output:
        out = [f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/n2g.{{gene}}.module{i}.tmp.npy' for i in range(1,167)],
    params:
        module = range(1,167),
    run:
        for in_f, out_f, module_idx in zip(input.out, output.out, params.module):
            n2g = np.load(in_f, allow_pickle=True).item()
            n2g['module'] = module_idx
            np.save(out_f, [n2g])


use rule yazar_HE_free_merge as yazar_module_n2g_merge_modules with:
    input:
        out = [f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/n2g.{{gene}}.module{i}.tmp.npy' for i in range(1,167)],
    output:
        out = f'staging/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/n2g.{{gene}}.npy',


rule yazar_n2g_merge_genes:
    input:
        out = expand('staging/yazar/combine_cts/{{cc}}/{params}/n2g.{gene}.npy',
                      params=yazar_paramspace.instance_patterns,
                      gene=['CTLA4', 'IL2RA']),
    output:
        out = f'analysis/yazar/combine_cts/{{cc}}/{yazar_paramspace.wildcard_pattern}/n2g.allgenes.npy',
    params:
        gene = ['CTLA4', 'IL2RA'],
    run:
        n2g_all = {}
        for f, gene in zip(input.out, params.gene):
            n2g = np.load(f, allow_pickle=True).item()
            n2g_all[gene] = n2g
        np.save(output.out, n2g_all)



use rule yazar_kinship_module as yazar_allcts_kinship_module with:
    input:
        modules = 'ngx/transpco_tableS1.csv',
        genes = 'data/Yazar2022Science/gene_location.txt',
        bed = expand('analysis/yazar/data/geno/chr{chr}.bed', chr=chrs),
        ctp = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctp.final.gz',
    output:
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/{{gene}}.module{{i}}.kinship.npy',


use rule yazar_n2g as yazar_allcts_module_n2g with:
    input:
        ctp = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctp.final.gz',
        ctnu = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/ctnu.final.gz',
        kinship = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/{{gene}}.module{{i}}.kinship.npy',
        P = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/P.final.gz',
        op_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/pca.gz',
        geno_pca = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/geno.eigenvec',
        obs = 'analysis/yazar/exclude_repeatedpool.obs.txt',
        genes = 'data/Yazar2022Science/gene_location.txt',
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/n2g.{{gene}}.module{{i}}.npy',
    params:
        jk = True,


use rule yazar_module_n2g_add_module as yazar_allcts_module_n2g_add_module with:
    input:
        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/n2g.{{gene}}.module{i}.npy' for i in range(1,167)],
    output:
        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/n2g.{{gene}}.module{i}.tmp.npy' for i in range(1,167)],


use rule yazar_HE_free_merge as yazar_allcts_module_n2g_merge_modules with:
    input:
        out = [f'staging/yazar/{yazar_paramspace.wildcard_pattern}/n2g.{{gene}}.module{i}.tmp.npy' for i in range(1,167)],
    output:
        out = f'staging/yazar/{yazar_paramspace.wildcard_pattern}/n2g.{{gene}}.npy',


use rule yazar_n2g_merge_genes as yazar_allcts_n2g_merge_genes with:
    input:
        out = expand('staging/yazar/{params}/n2g.{gene}.npy',
                      params=yazar_paramspace.instance_patterns,
                      gene=['CTLA4', 'IL2RA']),
    output:
        out = f'analysis/yazar/{yazar_paramspace.wildcard_pattern}/n2g.allgenes.npy',


rule yazar_module_n2g_all:
    input:
        n2g = expand('analysis/yazar/combine_cts/{cc}/{params}/n2g.allgenes.npy',
                      cc=['cc1'], params=yazar_paramspace.instance_patterns,
                      ),
        n2g_allcts = expand('analysis/yazar/{params}/n2g.allgenes.npy',
                      params=yazar_paramspace.instance_patterns,
                      ),
