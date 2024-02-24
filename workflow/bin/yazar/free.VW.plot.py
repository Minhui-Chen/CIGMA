import os, math, re, sys
import numpy as np, pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from ctmm import draw


def main():
    # read
    P = pd.read_table(snakemake.input.P, index_col=0)
    ctp = np.load(snakemake.input.out, allow_pickle=True).item()

    #
    V = ctp['free']['cis_V']
    n = V.shape[0]
    C = V.shape[1]
    CTs = P.columns.tolist()
    order = snakemake.params.order
    order = order[np.isin(order, CTs)]
    missing_cts = [ct for ct in CTs if ct not in order]
    order = ['hom'] + order.tolist() + missing_cts
    # CTs = [re.sub(' ', '_', ct) for ct in CTs]
    colors = snakemake.params.colors
    colors['hom'] = '0.8'
    for ct in missing_cts:
        colors[ct] = '0.5'

    gen = pd.DataFrame()
    env = pd.DataFrame()
    # Free model
    ## gen
    gen['hom'] = ctp['free']['cis_hom_g2']
    V = np.diagonal(ctp['free']['cis_V'], axis1=1, axis2=2)
    V = pd.DataFrame(V, columns=CTs)
    gen = pd.concat((gen, V), axis=1)

    ## env
    env['hom'] = ctp['free']['hom_e2']
    W = np.diagonal(ctp['free']['W'], axis1=1, axis2=2)
    W = pd.DataFrame(W, columns=CTs)
    env = pd.concat((env, W), axis=1)

    ## trans
    if 'trans_V' in ctp['free'].keys():
        trans = pd.DataFrame()

        trans['hom'] = ctp['free']['trans_hom_g2']
        Vt = np.diagonal(ctp['free']['trans_V'], axis1=1, axis2=2)
        Vt = pd.DataFrame(Vt, columns=CTs)
        trans = pd.concat((trans, Vt), axis=1)

    # for ct in CTs:
        # env[f'{ct}+hom'] = env[ct] + env['hom']

    # plot
    mpl.rcParams.update({'font.size': 14})
    if 'trans_V' not in ctp['free'].keys():
        fig, axes = plt.subplots(nrows=2, figsize=(10, 10), dpi=600)
    else:
        fig, axes = plt.subplots(nrows=3, figsize=(10, 16), dpi=600)

    # cis gen
    if ((gen > 1) | (gen < -1)).sum().sum() / (gen.shape[0] * gen.shape[1]) > .2:
        gen_clip = gen.clip(lower=-5, upper=5)
    elif ((gen > .15) | (gen < -.15)).sum().sum() / (gen.shape[0] * gen.shape[1]) > .2:
        gen_clip = gen.clip(lower=-1, upper=1)
    else:
        gen_clip = gen.clip(lower=-.15, upper=.15)

    sns.violinplot(data=gen_clip, order=order, palette=colors, cut=0, scale='width', ax=axes[0])

    # Add median values as text annotations
    # TODO: try mean    
    y_loc = -0.12
    # gen_medians = gen.median(axis=0)[order]
    # axes[0].text(-0.05, y_loc, "median:", ha='center', va='center', fontsize=10, transform=axes[0].transAxes)
    # for xtick, median in zip(axes[0].get_xticks(), gen_medians):
    #     x = axes[0].transLimits.transform((xtick, median))[0]
    #     axes[0].text(x, y_loc, f"{median:.4f}", ha='center', va='center', fontsize=10, transform=axes[0].transAxes)
    gen_means = gen.mean(axis=0)[order]
    axes[0].text(-0.05, y_loc, "mean:", ha='center', va='center', transform=axes[0].transAxes)
    for xtick, mean in zip(axes[0].get_xticks(), gen_means):
        x = axes[0].transLimits.transform((xtick, mean))[0]
        axes[0].text(x, y_loc, f"{mean:.3f}", ha='center', va='center', transform=axes[0].transAxes)

    # env
    if ((env > 1) | (env < -1)).sum().sum() / (env.shape[0] * env.shape[1]) > .2:
        env_clip = env.clip(lower=-5, upper=5)
    else:
        env_clip = env.clip(lower=-1, upper=1)
    sns.violinplot(data=env_clip, order=order, palette=colors, cut=0, scale='width', ax=axes[1])

    # Add median values as text annotations
    # change to mean    
    # env_medians = env.median(axis=0)[order]
    # axes[1].text(-0.05, y_loc, "median:", ha='center', va='center', fontsize=10, transform=axes[1].transAxes)
    # for xtick, median in zip(axes[1].get_xticks(), env_medians):
    #     x = axes[1].transLimits.transform((xtick, median))[0]
    #     axes[1].text(x, y_loc, f"{median:.2f}", ha='center', va='center', fontsize=10, transform=axes[1].transAxes)
    env_means = env.mean(axis=0)[order]
    axes[1].text(-0.05, y_loc, "mean:", ha='center', va='center', transform=axes[1].transAxes)
    for xtick, mean in zip(axes[1].get_xticks(), env_means):
        x = axes[1].transLimits.transform((xtick, mean))[0]
        axes[1].text(x, y_loc, f"{mean:.3f}", ha='center', va='center', transform=axes[1].transAxes)

    axes[0].axhline(0, ls='--', color='0.8', zorder=0)
    axes[1].axhline(0, ls='--', color='0.8', zorder=0)
    axes[1].set_xlabel('')
    axes[0].set_ylabel('Cis-genetic variance (V)', fontsize=16)
    axes[1].set_ylabel('Environment variance (W)', fontsize=16)

    # trans
    if 'trans_V' in ctp['free'].keys():
        trans_clip = trans.clip(lower=-5, upper=5)
        sns.violinplot(data=trans_clip, order=order, palette=colors, cut=0, scale='width', ax=axes[2])

        # Add median values as text annotations
        # change to mean    
        # trans_medians = trans.median(axis=0)[order]
        # axes[2].text(-0.05, y_loc, "median:", ha='center', va='center', fontsize=10, transform=axes[2].transAxes)
        # for xtick, median in zip(axes[2].get_xticks(), trans_medians):
        #     x = axes[2].transLimits.transform((xtick, median))[0]
        #     axes[2].text(x, y_loc, f"{median:.2f}", ha='center', va='center', fontsize=10, transform=axes[2].transAxes)

        trans_means = trans.mean(axis=0)[order]
        axes[2].text(-0.05, y_loc, "mean:", ha='center', va='center', transform=axes[2].transAxes)
        for xtick, mean in zip(axes[2].get_xticks(), trans_means):
            x = axes[2].transLimits.transform((xtick, mean))[0]
            axes[2].text(x, y_loc, f"{mean:.3f}", ha='center', va='center', transform=axes[2].transAxes)

        axes[2].axhline(0, ls='--', color='0.8', zorder=0)
        axes[2].set_xlabel('')
        axes[2].set_ylabel('Trans-genetic variance (trans-V)', fontsize=16)



    # # change ticks in Free model
    # fig.canvas.draw_idle()
    # # plt.sca(ax)
    # locs, labels = plt.xticks()
    # print(locs)
    # for i in range(len(labels)):
    #     label = labels[i].get_text()
    #     print(label)
    #     if '-' in label:
    #         ct = label.split('-')[1]
    #         ct = re.sub('_', '\\_', ct)
    #         labels[i] = r'$W_{%s}$'%(ct)
    #     elif label == 'hom2':
    #         labels[i] = r'$\sigma_e^2$'
    # plt.xticks(locs, labels)
    fig.tight_layout(h_pad=2.5)
    fig.savefig(snakemake.output.png)


    # p value plot
    # if 'cis_V' not in ctp['p'].keys():
    #     with open(snakemake.output.p_png, 'w') as f:
    #         pass
    # else:
    #     fig, axes = plt.subplots(ncols=3, figsize=(18, 4))

    #     draw.scatter(-1 * np.log10(ctp['p']['free']['ct_beta']), -1 * np.log10(ctp['p']['free']['cis_V']), 
    #                  s=5, heatscatter=True, linregress=False, xlab='-log10 p(mean differentiation)', 
    #                  ylab='-log10 p(genetic variance differentiation)', ax=axes[0])

    #     draw.scatter(-1 * np.log10(ctp['p']['free']['ct_beta']), -1 * np.log10(ctp['p']['free']['W']), 
    #                  s=5, heatscatter=True, linregress=False, xlab='-log10 p(mean differentiation)', 
    #                  ylab='-log10 p(noise variance differentiation)', ax=axes[1])

    #     draw.scatter(-1 * np.log10(ctp['p']['free']['cis_V']), -1 * np.log10(ctp['p']['free']['W']), 
    #                  s=5, heatscatter=True, linregress=False, xlab='-log10 p(genetic variance differentiation)', 
    #                  ylab='-log10 p(noise variance differentiation)', ax=axes[2])

    #     fig.savefig(snakemake.output.p_png)
        

if __name__ == '__main__':
    main()
