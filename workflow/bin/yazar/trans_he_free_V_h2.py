import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from gxctmm import util


def main():
    #
    P = pd.read_table(snakemake.input.P, index_col=0)
    out = np.load(snakemake.input.out, allow_pickle=True).item()

    CTs = P.columns.tolist()
    order = snakemake.params.order
    order = order[np.isin(order, CTs)]
    missing_cts = [ct for ct in CTs if ct not in order]
    order = ['hom'] + order.tolist() + missing_cts
    # CTs = [re.sub(' ', '_', ct) for ct in CTs]
    colors1 = snakemake.params.colors1
    colors2 = snakemake.params.colors2
    # print(colors)
    # for ct in missing_cts:
    #     colors[0][ct] = '0.5'
    #     colors[1][ct] = '0.6'

    cs1 = [colors2[ct] for ct in order]
    cs2 = [colors1[ct] for ct in order]

    #
    cis_V = np.diagonal(out['free']['cis_V'], axis1=1, axis2=2)
    cis_V = pd.DataFrame(cis_V, columns=CTs)
    cis_V['hom'] = out['free']['cis_hom_g2']
    cis_V = cis_V.melt(var_name='CT', value_name='variance')
    cis_V['Regulation'] = 'cis'
    trans_V = np.diagonal(out['free']['trans_V'], axis1=1, axis2=2)
    trans_V = trans_V.clip(-10, 10)
    trans_V = pd.DataFrame(trans_V, columns=CTs)
    trans_V['hom'] = out['free']['trans_hom_g2']
    trans_V = trans_V.melt(var_name='CT', value_name='variance')
    trans_V['Regulation'] = 'trans'
    V = pd.concat([cis_V, trans_V], axis=0)


    cis_h2s = util.compute_h2(out['free']['cis_hom_g2'], out['free']['cis_V'], 
                            out['free']['hom_e2'] + out['free']['trans_hom_g2'], 
                            out['free']['W'] + out['free']['trans_V'])
    cis_h2s = pd.DataFrame(cis_h2s, columns=CTs)
    cis_h2s = cis_h2s.clip(-5, 5)
    cis_h2s = cis_h2s.melt(var_name='CT', value_name='h2')
    cis_h2s['Regulation'] = 'cis'
    trans_h2s = util.compute_h2(out['free']['trans_hom_g2'], out['free']['trans_V'], 
                            out['free']['hom_e2'] + out['free']['cis_hom_g2'], 
                            out['free']['W'] + out['free']['cis_V'])
    trans_h2s = trans_h2s.clip(-50, 50)
    trans_h2s = pd.DataFrame(trans_h2s, columns=CTs)
    trans_h2s = trans_h2s.melt(var_name='CT', value_name='h2')
    trans_h2s['Regulation'] = 'trans'
    h2s = pd.concat([cis_h2s, trans_h2s], axis=0)
    
    # plot
    fig, axes = plt.subplots(nrows=2, figsize=(12, 8))

    ax = sns.barplot(data=V, x='CT', y='variance', hue='Regulation', order=order, 
                hue_order=['cis', 'trans'], palette=[cs1[0], cs2[0]],
                estimator='mean', errorbar='se', ax=axes[0])
    for bars, colors in zip(ax.containers, (cs1, cs2)):
        for bar, color in zip(bars, colors):
            bar.set_facecolor(color)

    ax = sns.barplot(data=h2s, x='CT', y='h2', hue='Regulation', order=order[1:], 
                hue_order=['cis', 'trans'], palette=[cs1[0], cs2[0]],
                estimator='mean', errorbar='se', ax=axes[1])  # TODO: using mean
    for bars, colors in zip(ax.containers, (cs1[1:], cs2[1:])):
        for bar, color in zip(bars, colors):
            bar.set_facecolor(color)

    fig.tight_layout()
    fig.savefig(snakemake.output.png)


if __name__ == '__main__':
    main()