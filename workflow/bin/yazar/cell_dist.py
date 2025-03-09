import numpy as np, pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt

def main():
    # set seaborn style
    #sns.set_theme()

    # data
    meta = pd.read_table(snakemake.input.obs, usecols=['cell', snakemake.params.ind_col, snakemake.params.ct_col])
    meta = meta.rename( columns={snakemake.params.ind_col:'ind', snakemake.params.ct_col:'ct'} )

    # count for each individiual
    ind_grouped = meta.groupby('ind')
    ind_sizes = ind_grouped['cell'].count()
    ## inds in order of decreasing number of cells
    inds = ind_sizes.sort_values(ascending=False).index.to_numpy()

    # count for cell types
    ct_grouped = meta.groupby('ct')
    ct_sizes = ct_grouped['cell'].count()
    ## cts in order of decreasing number of cells
    cts = ct_sizes.sort_values(ascending=False).index.to_numpy()
    ## exclude Plasma Erythrocytes
    cts = cts[~np.isin(cts,['Platelets','Erythrocytes'])]

    # count for ind-ct pairs
    grouped = meta.groupby(['ind','ct'])
    group_size = grouped['cell'].count()
    group_size = group_size.unstack(fill_value=0).stack()

    # collect data for plot
    x = range( 1, len(inds)+1 )
    y = []
    for ct in cts:
        y_tmp = group_size.loc[group_size.index.get_level_values('ct')==ct]
        y_tmp = y_tmp[inds]
        print( y_tmp )
        y.append( y_tmp.tolist() )

    # plot
    colors = sns.color_palette()[:len(cts)]
    plt.rcParams.update({'font.size': 14})
    fig, axes = plt.subplots(nrows=2, figsize=(8,10), dpi=600)

    axes[0].stackplot(x, y, colors=colors, labels=cts)
    axes[0].legend(loc='upper right', ncols=2)
    axes[0].set_xlabel( 'Individual' )
    axes[0].set_ylabel( 'Number of cells' )

    colors = snakemake.params.colors
    for ct in cts[int(snakemake.params.analyzed_ct_num):]:
        colors[ct] = colors['hom']

    data = pd.DataFrame(np.log10(np.array(y).T+1), columns=cts)
    sns.violinplot(data=data, cut=0, density_norm='width', palette=colors, ax=axes[1])
    plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right', fontsize=10)
    axes[1].axhline(y=np.log10(10+1), color='0.9', ls='--', zorder=0)
    axes[1].set_xlabel( 'Cell type', fontsize=16 )
    axes[1].set_ylabel( '$log_{10}$ (Number of cells + 1)', fontsize=16 )

    fig.tight_layout()
    plt.savefig( snakemake.output.png )

if __name__ == '__main__':
    main()
