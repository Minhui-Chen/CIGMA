import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from gxctmm import util, plot


def cor(V):
    std = np.sqrt(np.diag(V))
    return (V / np.outer(std, std))


def main():
    # 
    P = pd.read_table(snakemake.input.P, index_col=0)
    cts = P.columns
    out = np.load(snakemake.input.out, allow_pickle=True).item()
    order = snakemake.params.order
    order = order[np.isin(order, cts)].tolist()
    # order = order[np.isin(order, cts.to_list() + ['hom'])].tolist()
    colors = snakemake.params.colors
    colors = {key: colors[key] for key in order}

    #
    hom_g2 = out['free']['cis_hom_g2']
    hom_e2 = out['free']['hom_e2']

    # calculate heritability
    if 'trans_V' in out['free'].keys():
        h2s = util.compute_h2(hom_g2, out['free']['cis_V'], hom_e2 + out['free']['trans_hom_g2'], out['free']['W'] + out['free']['trans_V'])
        trans_h2s = util.compute_h2(out['free']['trans_hom_g2'], out['free']['trans_V'], hom_e2 + hom_g2, out['free']['W'] + out['free']['cis_V'])
    else: 
        h2s = util.compute_h2(hom_g2, out['free']['cis_V'], hom_e2, out['free']['W'])
    # mean_h2s = np.mean(h2s, axis=0)  # mean across genes

    # get median
    # h2_median = np.median(h2s, axis=0)
    # e2_median = 1 - h2_median
    # V_median = np.median(Vs, axis=0)
    # W_median = np.median(Ws, axis=0)

    # plot for h2
    if 'trans_V' in out['free'].keys():
        fig, axes = plt.subplots(figsize=(10, 8), nrows=2)
        ax = axes[0]
    else:
        fig, ax = plt.subplots(figsize=(10, 4))

    data = pd.DataFrame(h2s, columns=cts)
    # data = pd.DataFrame(h2s, columns=['hom'] + cts.to_list())
    if ((data > 1) | (data < -1)).sum().sum() / (data.shape[0] * data.shape[1]) > 0.2:
        clipped = data.clip(lower=-5, upper=5)
    elif ((data > 0.5) | (data < -0.5)).sum().sum() / (data.shape[0] * data.shape[1]) > 0.2:
        clipped = data.clip(lower=-1, upper=1)
    elif ((data > 0.2) | (data < -0.1)).sum().sum() / (data.shape[0] * data.shape[1]) > 0.2:
        clipped = data.clip(lower=-.5, upper=.5)
    else:
        clipped = data.clip(lower=-0.1, upper=0.2)
    sns.violinplot(data=clipped, cut=0, order=order, palette=colors, ax=ax)
    ax.axhline(y=0, color='0.9', ls='--', zorder=0)
    ax.set_ylabel('cis h2', fontsize=14)

    # add mean h2
    y_loc = -0.15
    medians = data.median(axis=0)[order]
    ax.text(-0.05, y_loc, "median:", ha='center', va='center', fontsize=10, transform=ax.transAxes)
    for xtick, median in zip(ax.get_xticks(), medians):
        x = ax.transLimits.transform((xtick, median))[0]
        ax.text(x, y_loc, f"{median:.3f}", ha='center', va='center', fontsize=10, transform=ax.transAxes)
    # means = data.mean(axis=0)[order]
    # ax.text(-0.05, y_loc, "mean:", ha='center', va='center', fontsize=10, transform=ax.transAxes)
    # for xtick, mean in zip(ax.get_xticks(), means):
    #     x = ax.transLimits.transform((xtick, mean))[0]
    #     ax.text(x, y_loc, f"{mean:.3f}", ha='center', va='center', fontsize=10, transform=ax.transAxes)

    if 'trans_V' in out['free'].keys():
        ax = axes[1]

        data = pd.DataFrame(trans_h2s, columns=cts)
        # data = pd.DataFrame(trans_h2s, columns=['hom'] + cts.to_list())
        if ((data > 5) | (data < -5)).sum().sum() / (data.shape[0] * data.shape[1]) > 0.2:
            clipped = data.clip(lower=-10, upper=10)
        elif ((data > 0.2) | (data < -0.1)).sum().sum() / (data.shape[0] * data.shape[1]) > 0.2:
            clipped = data.clip(lower=-5, upper=5)
        else:
            clipped = data.clip(lower=-0.1, upper=0.2)
        sns.violinplot(data=clipped, cut=0, order=order, palette=colors, ax=ax)
        ax.axhline(y=0, color='0.9', ls='--', zorder=0)
        ax.set_ylabel('trans h2', fontsize=14)

        # add mean h2
        y_loc = -0.15
        medians = data.median(axis=0)[order]
        ax.text(-0.05, y_loc, "median:", ha='center', va='center', fontsize=10, transform=ax.transAxes)
        for xtick, median in zip(ax.get_xticks(), medians):
            x = ax.transLimits.transform((xtick, median))[0]
            ax.text(x, y_loc, f"{median:.3f}", ha='center', va='center', fontsize=10, transform=ax.transAxes)
            
    # 
    fig.tight_layout(h_pad=5)
    fig.savefig(snakemake.output.h2)


if __name__ == '__main__':
    main()
