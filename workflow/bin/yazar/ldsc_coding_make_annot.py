import numpy as np
import pandas as pd
from cigma import util
from pybedtools import BedTool


def main():
    # par
    out = np.load(snakemake.input.out, allow_pickle=True).item()
    op = pd.read_table(snakemake.input.op, index_col=0)
    location = pd.read_csv(snakemake.input.location, sep='\t')
    vanno = pd.read_csv(snakemake.input.vanno, sep='\t', usecols=['chr', 'pos', 'ref', 'alt', 'rsid', 'consequence'])
    bim = pd.read_table(snakemake.input.bim, names=['chr', 'rsid', 'cm', 'pos', 'a1', 'a2'])
    annot = bim[['chr', 'rsid', 'cm', 'pos']].copy()
    nbin = int(snakemake.wildcards['n'])
    window = int(float(snakemake.wildcards['window']))
    chr = int(snakemake.wildcards['chr'])
    feature = snakemake.wildcards.feature
    coding_variants = snakemake.params.coding_variants

    # data process
    data = util.read_out(out, ['gene', 'free:hom_g2', 'free:v'])
    data = data.rename(columns={'free:hom_g2': 'hom_g2', 'free:v': 'v'})
    data['var_beta'] = out['free']['ct_beta'].var(axis=1)
    data['g'] = data['hom_g2'] + data['v']
    data = data.loc[data['g'] > 0] # remove negative genetic variance
    print(data.shape[0], 'genes with positive genetic variance')
    data['specificity'] = data['v'] / data['g']
    data = data.merge(location[['feature', 'chr', 'start', 'end']], left_on='gene', right_on='feature').drop(columns=['feature'])
    op = op.mean(axis=0)
    op.name = 'mean_expr'
    op.index.name = 'gene'
    op = op.reset_index()
    data = data.merge(op, on='gene')

    ## vanno process
    vanno = vanno.loc[(vanno['chr'] == chr) & (vanno['consequence'].isin(coding_variants))]
    print(f'number of coding variants on chr{chr}: {vanno.shape[0]}')
    vanno = vanno.merge(bim[['chr', 'rsid', 'pos', 'a1', 'a2']])
    print(f'number of coding variants on chr{chr}: {vanno.shape[0]}')
    vanno = vanno.loc[((vanno['ref'] == vanno['a1']) & (vanno['alt'] == vanno['a2'])) | 
                      ((vanno['ref'] == vanno['a2']) & (vanno['alt'] == vanno['a1']))]
    print(f'number of coding variants on chr{chr}: {vanno.shape[0]}')

    # divide into bins
    data[f'{feature}_bin'] = pd.qcut(data[feature], nbin, labels=False)
    print(data.groupby(f'{feature}_bin')[feature].median())
    data['mean_expr_bin'] = pd.qcut(data['mean_expr'], nbin, labels=False)
    print(data.groupby('mean_expr_bin')['mean_expr'].median())

    # add window
    data['raw_start'] = data['start']
    data['raw_end'] = data['end']
    data['start'] = np.maximum(1, data['start'] - window)
    data['end'] = data['end'] + window

    # make annot
    iter_bim = [['chr'+str(x1), x2 - 1, x2] for (x1, x2) in np.array(bim[['chr', 'pos']])]
    bimbed = BedTool(iter_bim)

    ## gene bin
    for var in [f'{feature}_bin', 'mean_expr_bin']:
        for i in range(nbin):
            print(f'making annot for {var} bin {i}')
            data_subset = data.loc[data[var] == i]
            iter_df = [['chr'+(str(x1).lstrip('chr')), x2 - 1, x3] for (x1,x2,x3) in np.array(data_subset[['chr', 'start', 'end']])]
            bed_for_annot = BedTool(iter_df).sort().merge()

            annotbed = bimbed.intersect(bed_for_annot)
            bp = [x.start + 1 for x in annotbed]
            df_int = pd.DataFrame({'pos': bp, f'{var}_{i}':1})
            annot = pd.merge(annot, df_int, how='left', on='pos')
            annot.fillna(0, inplace=True)
            annot[[f'{var}_{i}']] = annot[[f'{var}_{i}']].astype(int)

    ## coding variants
    print('making annot for coding variants')
    ### sanity check rsids are unique and info matching
    tmp = annot.merge(vanno[['rsid', 'pos']], on='rsid')
    assert all(np.sort(tmp['rsid'].values) == np.sort(vanno['rsid'].values))
    assert tmp['pos_x'].equals(tmp['pos_y'])

    annot = annot.rename(columns={'chr':'CHR', 'pos':'BP', 'rsid':'SNP', 'cm':'CM'})
    for var, output_f in zip([f'{feature}_bin', 'mean_expr_bin'], [snakemake.output.var_annot, snakemake.output.mean_annot]):
        annot_tmp = annot.copy()
        # drop the first bin to avoid collinearity
        annot_tmp = annot_tmp.drop(columns=[f'mean_expr_bin_0'])

        # coding variants in mean expression bins
        for i in range(1, nbin):
            data_subset = data.loc[data['mean_expr_bin'] == i]
            iter_df = [['chr'+(str(x1).lstrip('chr')), x2 - 1, x3] for (x1,x2,x3) in np.array(data_subset[['chr', 'raw_start', 'raw_end']])]
            bed_for_annot = BedTool(iter_df).sort().merge()
            annotbed = bimbed.intersect(bed_for_annot)
            bp = [x.start + 1 for x in annotbed]

            annot_tmp[f'coding_mean_expr_bin_{i}'] = 0
            annot_tmp.loc[(annot_tmp['SNP'].isin(vanno['rsid'])) & (annot_tmp['BP'].isin(bp)), f'coding_mean_expr_bin_{i}'] = 1

        # coding variants in other bins
        if var != 'mean_expr_bin':
            for i in range(nbin):
                data_subset = data.loc[data[var] == i]
                iter_df = [['chr'+(str(x1).lstrip('chr')), x2 - 1, x3] for (x1,x2,x3) in np.array(data_subset[['chr', 'raw_start', 'raw_end']])]
                bed_for_annot = BedTool(iter_df).sort().merge()
                annotbed = bimbed.intersect(bed_for_annot)
                bp = [x.start + 1 for x in annotbed]

                annot_tmp[f'coding_{var}_{i}'] = 0
                annot_tmp.loc[(annot_tmp['SNP'].isin(vanno['rsid'])) & (annot_tmp['BP'].isin(bp)), f'coding_{var}_{i}'] = 1

        # save
        annot_tmp.to_csv(output_f, sep='\t', index=False)


if __name__ == '__main__':
    main()