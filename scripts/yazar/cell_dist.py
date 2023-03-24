import numpy as np, pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main():
    # set seaborn style
    sns.set_theme()

    # data
    meta = pd.read_table(snakemake.input.meta, usecols=['cell','donor_id','predicted.celltype.l2'])
    meta = meta.rename( columns={'donor_id':'ind', 'predicted.celltype.l2':'ct'} )

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

    # count for ind-ct pairs
    grouped = meta.groupby(['ind','ct'])
    group_size = grouped['cell'].count()

    # collect data for plot
    x = range( 1, len(inds)+1 )
    y = []
    for ct in cts:
        y_tmp = grouped.loc[grouped.index.get_level_values('ct')==ct]
        y_tmp = y_tmp[inds].tolist()
        y.append( y_ tmp )

    # plot
    fig, ax = plt.subplots()

    plt.stackplot(x, y, labels=cts)
    plt.legend('upper right')
    plt.axhline(y=100, color='0.9', ls='--', zorder=0)
    plt.set_xlabel( 'Individual' )
    plt.set_ylabel( 'Number of cells' )

    plt.savefig( snakemake.output.png )

if __name__ == '__main__':
    main()
