import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from gxctmm import util, plot
from ctmm import draw


# def violinplot(data, ax, property, feature):
#     '''
#     feature: h2 or V
#     '''

#     # clip h2
#     data[feature] = data[feature].clip(-1, 1)

#     sns.violinplot(data=data, x=property + '_bin', y=feature, cut=0)
#     # sns.boxplot(data=data, x=property + '_bin', y='h2')
#     plt.xticks(fontsize=7, rotation=10)
#     ax.set_ylabel(f'Mean {feature} across cell types')

#     # add mean h2
#     data_grouped = data.groupby(property + '_bin')
#     medians = data_grouped[feature].median()
#     print(medians)  # TODO: check 
#     y_loc = -0.15
#     ax.text(-0.05, y_loc, "median:", ha='center', va='center', fontsize=10, transform=ax.transAxes)
#     for xtick, median in zip(ax.get_xticks(), medians):
#         x = ax.transLimits.transform((xtick, median))[0]
#         ax.text(x, y_loc, f"{median:.3f}", ha='center', va='center', fontsize=10, transform=ax.transAxes)


def scatter_bin(data, property, num_bins, png): 
    # divide EDS into bins
    data[property + '_bin'], bins = pd.qcut(data[property], q=num_bins, retbins=True)
    # print(bins)

    h2_lower = -2
    h2_upper = 2
    V_upper = 1
    V_prop_lower = -10
    V_prop_upper = 10

    data['h2'] = data['h2'].clip(h2_lower, h2_upper)
    data['V'] = data['V'].clip(upper=V_upper)
    data['V_prop'] = data['V_prop'].clip(V_prop_lower, V_prop_upper)
    grouped = data.groupby(property + '_bin', observed=True)

    # print(grouped['h2'].std() / grouped.size().apply(lambda x: x**0.5))

    # plot
    fig, axes = plt.subplots(ncols=2, nrows=6, figsize=(18, 30))
    i = -1
    fs = 10

    # h2
    i += 1
    draw.scatter(data[property], data['h2'], ax=axes[i, 0], 
                 xlab=property, ylab='Mean h2 across cell types', 
                 linregress_label=True, coeff_determination=False, heatscatter=True, s=6)

    # lineplot
    sns.pointplot(data=data, x=property + '_bin', y='h2', estimator='mean', 
                  errorbar='se', label='mean', ax=axes[i, 1])
    sns.pointplot(data=data, x=property + '_bin', y='h2', estimator='median', 
                  errorbar=None, label='median', color='green', ax=axes[i, 1])
    axes[i, 1].set_ylabel('Mean h2 across cell types')
    axes[i, 1].tick_params(axis='x', labelsize=fs, rotation=30)
    axes[i, 1].legend()

    # add meta regression
    line, p = plot.meta_regression(grouped, 'h2', 'mean')
    line2, p2 = plot.meta_regression(grouped, 'h2', 'median')
    axes[i, 1].text(0, 1.02, f'mean: {line}, p(slope) = {p:.3e}; median: {line2}, p(slope) = {p2:.3e}', 
                    fontsize=fs, transform=axes[i, 1].transAxes)


    # hom_g2 + V
    i += 1
    draw.scatter(data[property], data['g'], ax=axes[i, 0], xlab=property, 
                 ylab='Mean genetic variance (hom_g + V) across cell types', 
                 linregress_label=True, coeff_determination=False, heatscatter=True, s=6)

    # lineplot
    sns.pointplot(data=data, x=property + '_bin', y='g', estimator='mean', errorbar='se', 
                  label='mean', ax=axes[i, 1])
    sns.pointplot(data=data, x=property + '_bin', y='g', estimator='median', errorbar=None, 
                  label='median', color='green', ax=axes[i, 1])
    axes[i, 1].set_ylabel('Mean genetic variance (hom_g + V) across cell types')
    axes[i, 1].tick_params(axis='x', labelsize=fs, rotation=30)
    axes[i, 1].legend()

    # add meta regression
    line, p = plot.meta_regression(grouped, 'g', 'mean')
    line2, p2 = plot.meta_regression(grouped, 'g', 'median')
    axes[i, 1].text(0, 1.02, f'mean: {line}, p(slope) = {p:.3e}; median: {line2}, p(slope) = {p2:.3e}', 
                    fontsize=fs, transform=axes[i, 1].transAxes)

    # hom_g2
    i += 1
    draw.scatter(data[property], data['hom_g2'], ax=axes[i, 0], xlab=property, 
                 ylab='Shared genetic variance (hom_g2)', 
                 linregress_label=True, coeff_determination=False, heatscatter=True, s=6)

    # lineplot
    sns.pointplot(data=data, x=property + '_bin', y='hom_g2', estimator='mean', errorbar='se', 
                  label='mean', ax=axes[i, 1])
    sns.pointplot(data=data, x=property + '_bin', y='hom_g2', estimator='median', errorbar=None, 
                  label='median', color='green', ax=axes[i, 1])
    axes[i, 1].set_ylabel('Shared genetic variance (hom_g2)')
    axes[i, 1].tick_params(axis='x', labelsize=fs, rotation=30)
    axes[i, 1].legend()

    # add meta regression
    line, p = plot.meta_regression(grouped, 'hom_g2', 'mean')
    line2, p2 = plot.meta_regression(grouped, 'hom_g2', 'median')
    axes[i, 1].text(0, 1.02, f'mean: {line}, p(slope) = {p:.3e}; median: {line2}, p(slope) = {p2:.3e}', 
                    fontsize=fs, transform=axes[i, 1].transAxes)


    # V
    i += 1
    draw.scatter(data[property], data['V'], ax=axes[i, 0], xlab=property, 
                 ylab='Mean ct-specific genetic variance (V) across cell types', 
                 linregress_label=True, coeff_determination=False, heatscatter=True, s=6)

    # lineplot
    sns.pointplot(data=data, x=property + '_bin', y='V', estimator='mean', errorbar='se', 
                  label='mean', ax=axes[i, 1])
    sns.pointplot(data=data, x=property + '_bin', y='V', estimator='median', errorbar=None, 
                  label='median', color='green', ax=axes[i, 1])
    axes[i, 1].set_ylabel('Mean ct-specific genetic variance (V) across cell types')
    axes[i, 1].tick_params(axis='x', labelsize=fs, rotation=30)
    axes[i, 1].legend()

    # add meta regression
    line, p = plot.meta_regression(grouped, 'V', 'mean')
    line2, p2 = plot.meta_regression(grouped, 'V', 'median')
    axes[i, 1].text(0, 1.02, f'mean: {line}, p(slope) = {p:.3e}; median: {line2}, p(slope) = {p2:.3e}', 
                    fontsize=fs, transform=axes[i, 1].transAxes)


    # V / (hom_g2 + V)
    i += 1
    draw.scatter(data[property], data['V_prop'], ax=axes[i, 0], xlab=property, 
                 ylab='Proportion of ct-specific genetic variance (V / [hom_g2 + V])', 
                 linregress_label=True, coeff_determination=False, heatscatter=True, s=6)

    # lineplot
    sns.pointplot(data=data, x=property + '_bin', y='V_prop', estimator='mean', errorbar='se', 
                  label='mean', ax=axes[i, 1])
    sns.pointplot(data=data, x=property + '_bin', y='V_prop', estimator='median', errorbar=None, 
                  label='median', color='green', ax=axes[i, 1])
    axes[i, 1].set_ylabel('Proportion of ct-specific genetic variance (V / [hom_g2 + V])')
    axes[i, 1].tick_params(axis='x', labelsize=fs, rotation=30)
    axes[i, 1].legend()

    # add meta regression
    line, p = plot.meta_regression(grouped, 'V_prop', 'mean')
    line2, p2 = plot.meta_regression(grouped, 'V_prop', 'median')
    axes[i, 1].text(0, 1.02, f'mean: {line}, p(slope) = {p:.3e}; median: {line2}, p(slope) = {p2:.3e}', 
                    fontsize=fs, transform=axes[i, 1].transAxes)


    # p(V)
    if 'pV' in data.columns:
        i += 1
        draw.scatter(data[property], data['pV'], ax=axes[i, 0], xlab=property, 
                    ylab='-log10 p(cell types-specific genetic variance)', 
                    linregress_label=True, coeff_determination=False, heatscatter=True, s=6)
        
        # lineplot
        data['sig_V'] = data['pV'] > -np.log10(0.05 / data.shape[0])
        grouped = data.groupby(property + '_bin')
        # print(data[[property + '_bin', 'sig_V']])

        sns.pointplot(data=data, x=property + '_bin', y='sig_V', estimator='mean', errorbar=None, 
                    ax=axes[i, 1])
        axes[i, 1].set_ylabel('Proportion of genes with p(ct-specific genetic variance) < (0.05/#genes)')
        axes[i, 1].tick_params(axis='x', labelsize=fs, rotation=30)

        # add meta regression
        line, p = plot.meta_regression(grouped, 'sig_V', 'mean')
        axes[i, 1].text(0, 1.02, f'{line}, p(slope) = {p:.3e}', 
                        fontsize=fs, transform=axes[i, 1].transAxes)

    fig.tight_layout()
    fig.savefig(png)


def main():
    # par 
    out_f, ctp_f, eds_f, dist_png, genelen_png, loeuf_png, pli_png, actenhancerno_png, eds_png, matrix_png, matrix_png2 = sys.argv[1:]

    out = np.load(out_f, allow_pickle=True).item()
    eds = pd.read_table(eds_f)
    eds = eds.loc[eds['gene_length'] < 1e6]  # TODO: drop ~10 genes with >1MB
    eds = eds.loc[eds['ActivityLinking_EnhancerNumber'] < 100]  # TODO: drop ~10 genes
    eds['gene_length (kb)'] = eds['gene_length'] / 1e3

    # h2
    h2s = util.compute_h2(out['free']['cis_hom_g2'], out['free']['cis_V'], 
                            out['free']['hom_e2'], out['free']['W'])

    # V
    Vs = np.diagonal(out['free']['cis_V'], axis1=1, axis2=2)
    gs = Vs + out['free']['cis_hom_g2'][:, np.newaxis]

    # W
    Ws = np.diagonal(out['free']['W'], axis1=1, axis2=2)
    es = Ws + out['free']['hom_e2'][:, np.newaxis]

    # 
    data = pd.DataFrame({'gene': out['gene'], 'h2': np.mean(h2s, axis=1), 
                         'V': np.mean(Vs, axis=1), 'hom_g2': out['free']['cis_hom_g2'], 'g': np.mean(gs, axis=1),
                         'W': np.mean(Ws, axis=1), 'hom_e2': out['free']['hom_e2'], 'e': np.mean(es, axis=1),
                         'var(beta)': np.var(out['free']['ct_beta'], axis=1, ddof=1),
                         })
    data['V_prop'] = data['V'] / data['g']
    data['W_prop'] = data['W'] / data['e']

    # add p values if available
    if 'p' in out.keys():
        if 'V' in out['p']['free'].keys():
            data['pV'] = (-1) * np.log10(out['p']['free']['cis_V'])
            data['pW'] = (-1) * np.log10(out['p']['free']['W'])

    # add trans if avaiable
    if 'trans_V' in out['free'].keys():
        Vts = np.diagonal(out['free']['trans_V'], axis1=1, axis2=2)
        gts = Vts + out['free']['trans_hom_g2'][:, np.newaxis]

        h2s = util.compute_h2(out['free']['cis_hom_g2'], out['free']['cis_V'], 
                              out['free']['hom_e2'] + out['free']['trans_hom_g2'], 
                              out['free']['W'] + out['free']['trans_V'])
        trans_h2s = util.compute_h2(out['free']['trans_hom_g2'], out['free']['trans_V'], 
                              out['free']['hom_e2'] + out['free']['cis_hom_g2'], 
                              out['free']['W'] + out['free']['cis_V'])

        data['h2'] = np.mean(h2s, axis=1)
        data['t-h2'] = np.mean(trans_h2s, axis=1)
        data['t-V'] = np.mean(Vts, axis=1)
        data['t-hom_g2'] = out['free']['trans_hom_g2']
        data['t-g'] = np.mean(gts, axis=1)

        data['t-V_prop'] = data['t-V'] / data['t-g']

        if 'p' in out.keys():
            if 'V' in out['p']['free'].keys():
                data['t-pV'] = (-1) * np.log10(out['p']['free']['trans_V'])
                    
    # add mean expression
    ctp = pd.read_table(ctp_f, index_col=(0, 1))
    y = ctp.mean(axis=0)
    y.name = 'Mean expression'
    data = data.merge(y.to_frame(), left_on='gene', right_index=True)

    # merge with eds
    data = data.merge(eds, left_on='gene', right_on='gene_id')

    # matrix plot
    annotations = ['gene_length (kb)', 'LOEUF', 'pLI', 'ActivityLinking_EnhancerNumber', 'EDS', 
                   'blood_connected_rank', 'combined_connected_rank']
    features = ['Mean expression', 'var(beta)', 'h2', 'g', 'hom_g2', 'V', 'V_prop', 'pV', 
                't-h2', 't-g', 't-hom_g2', 't-V', 't-V_prop', 't-pV', 
                'e', 'hom_e2', 'W', 'W_prop']
    features = [x for x in features if x in data.columns]
    plot.matrix_feature_plot(data, features, annotations, 10, matrix_png)

    # matrix plot supp
    annotations2 = ['LOEUF', 'ActivityLinking_EnhancerNumber', 'EDS', 
                    'combined_connected_rank']
    features2 = ['hom_g2', 'V', 'V_prop',
                't-hom_g2', 't-V', 't-V_prop',
                ]
    features2 = [x for x in features2 if x in data.columns]
    plot.matrix_feature_plot(data, features2, annotations2, 10, matrix_png2)


    # dist of h2 and V
    if 'trans_V' in out['free'].keys():
        fig, axes = plt.subplots(nrows=5, ncols=2, figsize=(12, 20))
    else:
        fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(12, 12))

    axes[0, 0].hist(data['hom_g2'], bins=50, density=False)
    axes[0, 0].set_xlabel('hom_g2')

    axes[0, 1].hist(data['V'], bins=50, density=False)
    axes[0, 1].text(0, 1.05, f'{(data["V"] >= 1).sum()} outliers >= 1', 
            transform=axes[0, 1].transAxes)
    axes[0, 1].set_xlabel(r'$\bar{V}$')

    axes[1, 0].hist(Vs, bins=50, density=False)
    axes[1, 0].text(0, 1.05, f'{(Vs >= 1).sum().sum()} outliers >= 1, {(Vs <= -1).sum().sum()} outliers <= -1', 
            transform=axes[1, 0].transAxes)
    axes[1, 0].set_xlabel('V')
    axes[1, 1].hist(h2s.clip(-10, 10), bins=50, density=False)
    axes[1, 1].text(0, 1.05, f'{(h2s >= 5).sum().sum()} outliers >= 5, {(h2s <= -5).sum().sum()} outliers <= -5', 
            transform=axes[1, 1].transAxes)
    axes[1, 1].set_xlabel('h2')

    axes[2, 0].hist(data['hom_e2'], bins=50, density=False)
    axes[2, 0].set_xlabel('hom_e2')
    axes[2, 1].hist(data['W'], bins=50, density=False)
    axes[2, 1].text(0, 1.05, f'{(data["W"] >= 5).sum()} outliers >= 5', 
            transform=axes[2, 1].transAxes)
    axes[2, 1].set_xlabel(r'$\bar{W}')

    if 'trans_V' in out['free'].keys():
        axes[3, 0].hist(data['t-hom_g2'], bins=50, density=False)
        axes[3, 0].set_xlabel('Trans hom_g2')
        axes[3, 1].hist(data['t-V'], bins=50, density=False)
        axes[3, 1].set_xlabel(r'Trans $\bar{V}$')

        axes[4, 0].hist(Vts.clip(-25, 25), bins=50, density=False)
        axes[4, 0].text(0, 1.05, f'{(Vts >= 10).sum().sum()} outliers >= 10, {(Vts <= -10).sum().sum()} outliers <= -10', 
                transform=axes[4, 0].transAxes)
        axes[4, 0].set_xlabel('Trans V')

        axes[4, 1].hist(trans_h2s.clip(-100, 100), bins=50, density=False)
        axes[4, 1].text(0, 1.05, f'{(trans_h2s >= 50).sum().sum()} outliers >= 50, {(trans_h2s <= -50).sum().sum()} outliers <= -50', 
                transform=axes[4, 1].transAxes)
        axes[4, 1].set_xlabel('Trans h2')

    fig.tight_layout()
    fig.savefig(dist_png)

    # plot
    scatter_bin(data, 'EDS', 10, eds_png) 

    scatter_bin(data, 'LOEUF', 10, loeuf_png)

    scatter_bin(data, 'pLI', 10, pli_png)

    scatter_bin(data, 'gene_length (kb)', 10, genelen_png)

    scatter_bin(data, 'ActivityLinking_EnhancerNumber', 10, actenhancerno_png)


if __name__ == '__main__':
    main()