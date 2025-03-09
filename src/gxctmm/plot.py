from cProfile import label
from typing import Callable, Tuple, Optional, Union, List

import re, sys
from click import option
from matplotlib import axis
import numpy as np
import pandas as pd
import statsmodels
import statsmodels.api as sm
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from scipy import stats
from statsmodels.regression.linear_model import RegressionResults
from sklearn import linear_model
from . import log


def _clip(x: Union[pd.Series, np.ndarray], cut: float=3) -> np.ndarray:
    '''
    Clip to std standard deviation
    '''
    std = np.std(x)
    return np.clip(x, -std * cut, std * cut)


def yxline(ax: mpl.axes.Axes, color: str='red', linestyle: str='--', alpha: float=0.6) -> None:
    '''
    Add y = x line to the plot
    '''
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    lim = (max(xlim[0], ylim[0]), min(xlim[1], ylim[1]))
    ax.plot(lim, lim, color=color, linestyle=linestyle, alpha=alpha)  # Plot y = x line with dashed linestyle


def scatter(x: Union[list, np.ndarray], y: Union[list, np.ndarray], 
            png_fn: Optional[str]=None, ax: Optional[mpl.axes.Axes]=None, 
            xlab: Optional[str]=None, ylab: Optional[str]=None, 
            title: Optional[str]=None, text: Optional[str]=None, 
            xyline: bool=False, xyline_color:str='red', linregress: bool=False, 
            linregress_label: bool=False, linregress_color:str='0.8', 
            coeff_determination: bool=False, 
            s: Optional[float]=None, color: Optional[str]=None, cmap: Optional[str]='cividis', 
            heatscatter: bool=False):
    '''
    Make a simple scatter plot. 
    Able to add text. 

    Parameters
    ----------
    x : x axis data
    y : y axis data
    png_fn : output png file name. optional: png_fn or ax 
    ax : axis object, optional: png_fn or ax 
    xlab : str
    ylab : str
    title : str
    text : add text to figure
    xyline: add y = x dashed line
    linregress : whether to create a linear regression line with Pearson correlation
    coeff_determination : whether to compute coefficient of determination (R^2) 
        https://scikit-learn.org/stable/modules/generated/sklearn.metrics.r2_score.html
        https://scikit-learn.org/stable/auto_examples/linear_model/plot_ols.html#sphx-glr-auto-examples-linear-model-plot-ols-py
        https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LinearRegression.html#sklearn.linear_model.LinearRegression
    s : marker size
    color : marker color
    heatscatter: add KDE 

    Returns
    -------

    Notes
    -----

    '''
    if png_fn:
        fig, ax = plt.subplots()
    elif ax:
        pass
    else:
        log.logger.error('Either png_fn or ax must be provided!\n')
    if not color:
        color = plt.rcParams['axes.prop_cycle'].by_key()['color'][0]

    x = np.array(x)
    y = np.array(y)
    if not pd.api.types.is_integer_dtype(x):
        x[(x == np.inf) | (x == -np.inf)] = np.nan
    if not pd.api.types.is_integer_dtype(y):
        y[(y == np.inf) | (y == -np.inf)] = np.nan
    removed = ( np.isnan(x) | np.isnan(y) )
    x = x[~removed]
    y = y[~removed]

    if heatscatter:
        # adopt from https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib
        # calculae the point density
        try:
            xy = np.vstack([x,y])
            z = stats.gaussian_kde(xy)(xy)
            # sort the points by density, so that the densest points are plotted last
            idx = z.argsort()
            x, y, z = x[idx], y[idx], z[idx]
            sc = ax.scatter(x, y, s=s, c=z, cmap=cmap)
        except:
            ax.scatter(x, y, s=s, color=color)
        if removed.sum() > 0:
            ax.text(0.2, 0.92, f'{removed.sum()} Inf/NA points removed', transform=ax.transAxes)
        #cbar = plt.colorbar(sc, ax=ax)
        #cbar.ax.set_ylabel('Density', rotation=270)
    else:
        ax.scatter(x, y, s=s, color=color)
        if removed.sum() > 0:
            ax.text(0.2, 0.92, f'{removed.sum()} Inf/NA points removed', transform=ax.transAxes)

    if xyline:
        # ax.plot(x, x, color=xyline_color, linestyle='--', alpha=0.6)  # Plot y = x line with dashed linestyle
        xlim = ax.get_xlim()
        ax.plot(xlim, xlim, color=xyline_color, linestyle='--', alpha=0.6)  # Plot y = x line with dashed linestyle

    if xlab:
        ax.set_xlabel(xlab)
    if ylab:
        ax.set_ylabel(ylab)
    if title:
        ax.set_title(title)
    if text:
        ax.text(0.2, 1.02, text, transform=ax.transAxes)
    if linregress:
        slope, intercept, r, p, stderr = stats.linregress(x, y)
        r, p = stats.pearsonr(x, y)
        if slope > 0:
            line = f'y={intercept:.2g}+{slope:.2g}x, r={r:.2g}, p={p:.2g}'
        else:
            line = f'y={intercept:.2g}{slope:.2g}x, r={r:.2g}, p={p:.2g}'
        if coeff_determination:
            # Create linear regression object
            regr = linear_model.LinearRegression()
            # Train the model
            regr.fit(np.array(x).reshape((-1, 1)), np.array(y).reshape((-1, 1)))
            # compute R2
            r2 = regr.score(np.array(x).reshape((-1, 1)), np.array(y).reshape((-1, 1)))
            line = line + f', R^2={r2:.2g}'
        if linregress_label:
            ax.plot(x, intercept + slope * x, linregress_color, label=line, zorder=10)
            ax.legend(fontsize='small')
        else:
            ax.plot(x, intercept + slope * x, linregress_color, zorder=10)

    if png_fn:
        fig.savefig(png_fn)


def mycolors(n: int=10, palette: str='muted', desat: float=None) -> list:
    """
    Return colors in hex format from the sns.color_palette('colorblind').

    Parameters
        n:  number of colors. Colorblind palette has 10 colors, if n > 10, cycle the palette.
        palette:    bright, muted, colorblind, pastel, deep, dark
        desat:  Proportion to desaturate each color by

    Returns:
        list of colors in hex format
    """
    colors = sns.color_palette(palette, desat=desat)
    # change the order of purple and grey
    blue, grey =colors[4], colors[7]
    colors[4] = grey
    colors[7] = blue
    colors = (colors*10)[:n]
    colors = [mpl.colors.to_hex(color) for color in colors]
    return colors


def mymarkers() -> list:
    """
    Return a list of markers for plot
    see https://matplotlib.org/stable/api/markers_api.html#module-matplotlib.markers
    """
    return ['.', 'x', 'P', 's', 'D', 'd', '1', '*']

def snsbox_get_x(x_categories: int, hue_categories: int=1, width:float =0.8) -> np.ndarray:
    """
    Get x coordinates of boxes in Seaborn boxplot.

    Parameters:
        x_categories:   number of categories on x axis
        hue_categories: number of hue categories
        width:  width of all the elements for one level of the major grouping variable. (from seaborn.boxplot with default .8)

    Returns:
        [category1:hue1, category1:hue2, category2:hue1, category2:hue2]
    """
    #print(x_categories, hue_categories, width)
    element_length = width / hue_categories
    if hue_categories == 1:
        return np.arange(x_categories)
    elif hue_categories % 2 == 0:
        shifts = np.arange(hue_categories / 2) * element_length + element_length / 2
        shifts = list(np.flip(shifts*(-1))) + list(shifts)
    else:
        shifts = (np.arange((hue_categories-1) / 2)+1) * element_length
        shifts = list(np.flip(shifts*(-1))) + [0] + list(shifts)
    shifts = [[x] for x in shifts]
    xs = np.repeat(np.atleast_2d(np.arange(x_categories)), hue_categories, axis=0)+np.array(shifts)
    xs = xs.T.flatten()
    return xs


def add_values_to_sns(ax, values):
    for xtick, value in zip(ax.get_xticks(), values):
        ax.text(xtick, value, f"{value:.3f}", ha='center', va='bottom', fontsize=10)


def sm_ols(y: Union[np.ndarray, pd.Series], X: Union[np.ndarray, pd.Series, pd.DataFrame]) -> RegressionResults:
    '''
    Returns:
        RegressionResults.params: pd.Series of estimated parameters.
        RegressionResults.pvaleus:  The two-tailed p values for the t-stats of the params.
    '''

    X = sm.add_constant(X)

    model = sm.OLS(y, X)
    results = model.fit()

    return results


def bin_gene_annotation(data: pd.DataFrame, annotation: str, num_bins: int) -> pd.DataFrame:
    data = data.copy()
    if len(data[annotation].unique()) > num_bins:
        bins, bin_labels = pd.qcut(data[annotation], q=num_bins, retbins=True)
        data[annotation + '_bin'] = bins
        grouped = data.groupby(by=annotation + '_bin', observed=True)
        medians = grouped[annotation].median()
        # format
        if annotation in ['gene_length (kb)', 'ActivityLinking_EnhancerNumber', 
                            'blood_connected_rank', 'combined_connected_rank']:
            custom_format = lambda x: format(x, '.0f')
        else:
            custom_format = lambda x: format(x, f'.{2}f') if x >= 0.01 else format(x, f'.{1}e')
        medians = medians.apply(custom_format)

        bins = bins.cat.rename_categories(medians.to_dict())
        data[annotation + '_bin'] = bins

        # print(bin_labels)
    else:
        data[annotation + '_bin'] = data[annotation]

    return data


def meta_regression(grouped_data: pd.core.groupby.DataFrameGroupBy, column: str, stat: Union[str, Callable], return_params: bool=False,
                    confounder: Optional[str]=None) -> Tuple[str, float]:
    '''
    Perform meta regression

    Parameters:
        grouped_data:   data frame grouped by a target x
        column: target y
        stat:   mean or median per group or self-defined function
        return_params:  whether to return estimated parameters
        confounder: a confounder to be included in the regression
    
    Returns:
        a tuple of 
            #. regression line
            #. p value for slope = 0
    '''

    if confounder is None:
        meta = pd.Series()
        if stat == 'mean':
            meta = grouped_data[column].mean()
        elif stat == 'median':
            meta = grouped_data[column].median()
        elif callable(stat):
            # print(grouped_data.groups.keys())
            # print(grouped_data.get_group('3'))
            # print(grouped_data.get_group('3')[column])
            meta = grouped_data[column].agg(stat)
        # print(meta) 

        X = pd.Series(meta.index.codes, index=meta.index, name=column)
        results = sm_ols(meta.to_numpy(), X)
        slope = results.params[column]
        intercept = results.params['const']

        # Print the equation of the line
        line = f"y = {slope:.4f}x + {intercept:.4f}"
        p = results.pvalues[column]
        
    else:
        meta = pd.DataFrame()
        columns = [column, confounder]
        if stat == 'mean':
            meta = grouped_data[columns].mean()
        elif stat == 'median':
            meta = grouped_data[columns].median()
        elif callable(stat):
            # print(grouped_data.groups.keys())
            # print(grouped_data.get_group('3'))
            # print(grouped_data.get_group('3')[column])
            meta = grouped_data[columns].agg(stat)
        # print(meta) 

        X = pd.DataFrame({'x': meta.index.codes, 'confounder': meta[confounder]})
        # adjust confounder
        X1 = sm.add_constant(X[['confounder']])
        model = sm.OLS(meta[column], X1).fit()
        model_resid = model.resid

        # fit residual
        results = sm_ols(model_resid, X[['x']])
        slope = results.params['x']
        intercept = results.params['const']

        # Print the equation of the line
        line = f"y = {slope:.4f}x + {intercept:.4f}"
        p = results.pvalues['x']

    if return_params:
        return line, p, slope, intercept
    else:
        return line, p


def _compute_V_prop_mean(xs: list) -> float:
    # print(xs[:5])
    data = {'a':[], 'b':[]}
    for x in xs:
        x = x.split('_')
        data['a'].append(float(x[0]))
        data['b'].append(float(x[1]))
    data = pd.DataFrame(data)

    means = data.mean()

    p = means.iloc[0] / means.sum()
    # restrict to >=0
    # if p < 0:
        # p = 0
    
    return p


def matrix_feature_plot(data: pd.DataFrame, features: list, annotations: list, num_bins: int, png: str = None,
                        h2_lower: float = -2, h2_upper: float = 2, V_upper: float = 1, 
                        V_prop_lower: float = -3, V_prop_upper: float = 3,
                        W_prop_lower: float = -3, W_prop_upper: float = 3,
                        ht2_lower: float = -2, ht2_upper: float = 2, Vt_upper: float = 1, 
                        Vt_prop_lower: float = 0, Vt_prop_upper: float = 1,
                        fs: float = 12) -> None:
    '''
    Make a matrix of plots of gene annotation (e.g. pLI) X gene feature (e.g. h2)

    Parameters:
        data:   data frame
        features:   gene feature e.g. h2, V, W
        annotations:  gene annotation e.g. gene length, pLI
        num_bins:   number of bins for each annotation
        png:    output figure file name
        h2_lower, h2_upper: cutoff for h2
        V_upper:    cutoff for V
        V_prop_lower, V_prop_upper: cutoff for V_prop
        fs: font size
    '''

    data = data.copy()

    # divide gene annotation into bins
    for annotation in annotations:
        data = bin_gene_annotation(data, annotation, num_bins)

    data.loc[(data['shared_h2'] > h2_upper) | (data['shared_h2'] < h2_lower), 'shared_h2'] = np.nan
    data.loc[data['V'] > V_upper, 'V'] = np.nan
    if 'V_prop' in data.columns:
        data.loc[(data['V_prop'] > V_prop_upper) | (data['V_prop'] < V_prop_lower), 'V_prop'] = np.nan
    if 'trans shared_h2' in data.columns:
        data.loc[(data['trans shared_h2'] > h2_upper) | (data['trans shared_h2'] < h2_lower), 'trans shared_h2'] = np.nan
    if 'trans V_prop' in data.columns:
        data.loc[(data['trans V_prop'] > V_prop_upper) | (data['trans V_prop'] < V_prop_lower), 'trans V_prop'] = np.nan
    if 'W_prop' in data.columns:
        data.loc[(data['W_prop'] > V_prop_upper) | (data['W_prop'] < V_prop_lower), 'W_prop'] = np.nan

    # plot
    ms = ['o', 's', '*']
    ylabs = {
             'Mean expression': 'Mean expression',
             'var(beta)': 'Variance of mean expression \nacross cts',
             'h2': 'Mean cis h2 across cell types', 
             'shared_h2': 'Shared cis h2 across cell types', 
             'ct_h2': 'Mean ct-specific cis h2 across cell types', 
             'h2_ratio': 'Ratio of ct-specific cis h2 to mean cis h2 across cell types', 
             'g': r'Mean cis variance (${hom_g}^2 + \bar{V}$)',
             'hom_g2': r'Shared cis variance (${hom_g}^2$)',
             'V': r'Mean ct-specific cis variance ($\bar{V}$)',
             'V_prop': 'Proportion of ct-specific cis variance\n' + r'($\bar{V} / ({hom_g}^2 + \bar{V})$)',
             'pV': 'Proportion of genes with sig cis V',
             'trans shared_h2': 'Mean trans h2 across cell types', 
             'trans g': r'Mean trans variance (${hom_g}^2 + \bar{V}$)',
             'trans hom_g2': r'Shared trans variance (${hom_g}^2$)',
             'trans V': r'Mean ct-specific trans variance ($\bar{V}$)',
             'trans V_prop': 'Proportion of ct-specific trans variance\n' + r'($\bar{V} / ({hom_g}^2 + \bar{V}$)',
             'trans pV': 'Proportion of genes with sig trans V',
             'e': r'Mean env variance (${hom_e}^2 + \bar{W}$)',
             'hom_e2': r'Shared env variance (${hom_e}^2$)',
             'W': r'Mean ct-specific env variance ($\bar{W}$)',
             'W_prop': 'Proportion of ct-specific env variance\n' + r'($\bar{W} / ({hom_e}^2 + \bar{W})$)',
             'pW': 'Proportion of genes with sig W',
             }
    cis_colors = sns.color_palette('Paired')[:2]
    trans_colors = sns.color_palette('Paired')[2:4]
    other_colors = sns.color_palette('Paired')[4:6]
    env_colors = sns.color_palette('Paired')[6:8]
    colors = {
             'Mean expression': other_colors,
             'var(beta)': other_colors,
             'shared_h2': cis_colors,
             'ct_h2': cis_colors,
             'h2_ratio': cis_colors,
             'g': cis_colors,
             'hom_g2': cis_colors,
             'V': cis_colors,
             'V_prop': cis_colors,
             'pV': cis_colors,
             'trans shared_h2': trans_colors, 
             'trans g': trans_colors,
             'trans hom_g2': trans_colors,
             'trans V': trans_colors,
             'trans V_prop': trans_colors,
             'trans pV': trans_colors,
             'e': env_colors,
             'hom_e2': env_colors,
             'W': env_colors,
             'W_prop': env_colors,
             'pW': env_colors,
             }

    
    fig, axes = plt.subplots(nrows=len(features), ncols=len(annotations), 
                             sharex='col', sharey='row', 
                             figsize=(6 * len(annotations), 4 * len(features)))

    for i, feature in enumerate(features):
        for j, annotation in enumerate(annotations):
            ax = axes[i, j]
            
            # meta regression
            grouped = data.groupby(annotation + '_bin', observed=True)

            line, p = meta_regression(grouped, feature, 'mean')
            ls = '-' if p < 0.05 else '--'

            if feature in ['V_prop']:
                pass
            else:
                sns.pointplot(data=data, x=annotation + '_bin', y=feature, estimator='mean',
                            errorbar='se', linestyles=ls, color=colors[feature][1], 
                            markers=ms[0], label='mean', ax=ax)

            if feature not in ['pV', 'trans pV']:
                line2, p2 = meta_regression(grouped, feature, 'median')
                ls2 = '-' if p2 < 0.05 else '--'                

                sns.pointplot(data=data, x=annotation + '_bin', y=feature, estimator='median',
                            errorbar=None, linestyles=ls2, color=colors[feature][0], 
                            markers=ms[1], label='median', ax=ax)
                
                if feature in ['V_prop'] or (feature == 'W_prop' and 'trans V_prop' in features):
                    # don't draw for trans analysis, since both W and t-V are too noisy and have negative prop
                    tmp_data = data.copy()
                    if feature == 'V_prop':
                        tmp_data['V_prop'] = tmp_data['V'].astype(str) + '_' + tmp_data['hom_g2'].astype(str)
                    elif feature == 'trans V_prop':
                        tmp_data['t-V_prop'] = tmp_data['t-V'].astype(str) + '_' + tmp_data['trans hom_g2'].astype(str)
                    elif feature == 'W_prop':
                        tmp_data['W_prop'] = tmp_data['W'].astype(str) + '_' + tmp_data['hom_e2'].astype(str)
                    tmp_grouped = tmp_data.groupby(annotation + '_bin', observed=True)
                    line3, p3 = meta_regression(tmp_grouped, feature, _compute_V_prop_mean)
                    ls3 = '-' if p3 < 0.05 else '--'

                    # ax.text(0.01, 0.92, f'ratio of means: p = {p3:.2e}', 
                                    # fontsize=fs, transform=ax.transAxes)
                    feature_data = tmp_grouped[feature].agg(_compute_V_prop_mean).reset_index(drop=False)
                    sns.pointplot(data=feature_data, x=annotation + '_bin', y=feature, 
                                linestyles=ls3, color='0.8', 
                                markers=ms[2], label='ratio of means', ax=ax)

                if feature in ['V_prop']:
                    ax.text(0.01, 1.02, f'p(median) = {p2:.2e}; p(ratio of means) = {p3:.2e}', 
                                    fontsize=fs, transform=ax.transAxes)
                else:
                    ax.text(0.01, 1.02, f'p(mean) = {p:.2e}; p(median) = {p2:.2e}', 
                                    fontsize=fs, transform=ax.transAxes)
            else:           
                ax.text(0.01, 1.02, f'p(mean) = {p:.2e}', 
                                fontsize=fs, transform=ax.transAxes)

            # to show mean on the top
            if feature in ['V_prop']:
                pass
            else:
                sns.pointplot(data=data, x=annotation + '_bin', y=feature, estimator='mean',
                            errorbar='se', linestyles=ls, color=colors[feature][1], 
                            markers=ms[0], ax=ax)

                
            # if i == len(features) - 1:
                # ax.tick_params(axis='x', rotation=30)
            # else:
                # ax.set_xlabel('')
            xlabel = re.sub('_', ' ', annotation + '_bin')
            ax.set_xlabel(xlabel, fontsize=18)
            ax.set_ylabel(ylabs[feature], fontsize=13)

        # axes[i, 0].set_ylabel(ylabs[feature])
    # axes[0, 0].legend()
    
    fig.tight_layout()

    if png:
        fig.savefig(png)


def matrix_feature_main_plot(data: pd.DataFrame, annotations: list, num_bins: int, 
                             png: Optional[str] = None,
                             h2_lower: float = -2, h2_upper: float = 2, 
                             V_upper: float = 1, V_prop_lower: float = -3, 
                             V_prop_upper: float = 3,
                             figsize: Optional[tuple]=None, fs: float = 12,
                             lw: float=2) -> None:

    data = data.copy()

    # divide gene annotation into bins
    for annotation in annotations:
        data = bin_gene_annotation(data, annotation, num_bins)

    data.loc[(data['shared_h2'] > h2_upper) | (data['shared_h2'] < h2_lower), 'shared_h2'] = np.nan
    data.loc[data['V'] > V_upper, 'V'] = np.nan
    if 'V_prop' in data.columns:
        data.loc[(data['V_prop'] > V_prop_upper) | (data['V_prop'] < V_prop_lower), 'V_prop'] = np.nan
    if 't-h2' in data.columns:
        data.loc[(data['trans shared_h2'] > h2_upper) | (data['trans shared_h2'] < h2_lower), 'trans shared_h2'] = np.nan
    if 't-V_prop' in data.columns:
        data.loc[(data['trans V_prop'] > V_prop_upper) | (data['trans V_prop'] < V_prop_lower), 'trans V_prop'] = np.nan
    if 'W_prop' in data.columns:
        data.loc[(data['W_prop'] > V_prop_upper) | (data['W_prop'] < V_prop_lower), 'W_prop'] = np.nan

    if figsize is None:
        figsize = (5 * len(annotations), 4 * 2)
    fig, axes = plt.subplots(nrows=2, ncols=len(annotations), 
                                sharex='col', sharey='row', figsize=figsize)
    
    for j, annotation in enumerate(annotations):
        grouped = data.groupby(annotation + '_bin', observed=True)

        line1, p1 = meta_regression(grouped, 'hom_g2', 'mean')
        ls1 = '-' if p1 < 0.05 else '--'
        line2, p2 = meta_regression(grouped, 'V', 'mean')
        ls2 = '-' if p2 < 0.05 else '--'

        ax = axes[0, j]
        sns.pointplot(data=data, x=annotation + '_bin', y='hom_g2', estimator='mean',
                    linestyles=ls1, color=sns.color_palette('bright')[0], lw=lw,
                    label=r'shared variance ($\sigma_g^2$)', ax=ax)
        sns.pointplot(data=data, x=annotation + '_bin', y='V', estimator='mean',
                    linestyles=ls2, color=sns.color_palette('bright')[1], lw=lw,
                    label=r'specific variance ($\bar{V}$)', ax=ax)
        if j != 0:
            ax.legend().set_visible(False)

        line, p = meta_regression(grouped, 'V_prop', 'median')
        ls = '-' if p < 0.05 else '--'
        ax = axes[1, j]
        sns.pointplot(data=data.loc[data['g'] > 0], x=annotation + '_bin', y='V_prop', 
                      estimator='median', linestyles=ls, color='0.5', lw=lw, ax=ax)

        xlabel = re.sub('_', ' ', annotation + '_bin')
        if annotation == 'ActivityLinking_EnhancerNumber':
            xlabel = 'Enhancer number bin'
        elif annotation == 'combined_connected_rank':
            xlabel = 'Connected rank bin'

            # set x tick labels
            xticklabels = []
            for i, label in enumerate(ax.get_xticklabels()):
                if i % 2 == 0:
                    xticklabels.append(label.get_text())
                else:
                    xticklabels.append('')
            ax.set_xticklabels(xticklabels)

        ax.set_xlabel(xlabel, fontsize=fs)

    axes[0, 0].set_ylabel('Cis variance', fontsize=fs)
    axes[1, 0].set_ylabel('Specificity: ' + r'$\bar{V}/(\sigma_g^2 + \bar{V})$', 
                          fontsize=fs)

    # add arrow
    arrow_position = (0.5, 0.6)
    axes[0, 0].annotate('', xy=(arrow_position[0] - 0.40, arrow_position[1]), xytext=arrow_position,
             arrowprops=dict(arrowstyle='->', lw=1.5),
             fontsize=10, ha='center', va='center', xycoords=axes[0, 0].transAxes)
    axes[0, 0].text(0.3, 0.62, 'more constrained', ha='center', va='bottom', 
                    fontsize=fs-1, transform=axes[0, 0].transAxes)

    if 'combined_connected_rank' in annotations:
        axes[0, 2].annotate('', xy=(arrow_position[0] - 0.40, arrow_position[1]), xytext=arrow_position,
                arrowprops=dict(arrowstyle='->', lw=1.5),
                fontsize=10, ha='center', va='center', xycoords=axes[0, 2].transAxes)
        axes[0, 2].text(0.3, 0.62, 'more connected', ha='center', va='bottom', 
                        fontsize=fs-2, transform=axes[0, 2].transAxes)

    fig.tight_layout()

    if png:
        fig.savefig(png)


def ctp_h2_plot(out: dict, ax: plt.Axes, colors: list, labels: list = ['cis', 'trans'],
               rng: np.random.Generator = np.random.default_rng(42), 
               width: float=0.2, permutation: bool=False, n_permutations=9999, 
               h2: str='ratio', model: str='free', 
               filter: Optional[np.ndarray]=None,
               ) -> Tuple[Optional[float], Optional[float]]:
    """
    Plot the cis vs trans heritability per ct.

    Parameters:
        out: output dictionary 
        ax: matplotlib axis
        colors: list of colors for shared and specific heritability
        labels: list of labels for e.g. cis vs trans
        rng: random number generator
        width: bar width
        permutation: whether to perform permutation test
        n_permutations: number of permutations to perform
        h2: 'ratio', 'mean', 'median', or 'mean_clip' across transcriptome
        model: 'free' or 'iid'
        filter: gene filter
    """     

    # extract data
    cis_hom_g2 = out['free']['hom_g2']
    trans_hom_g2 = out['free']['hom_g2_b']
    hom_e2 = out['free']['hom_e2']
    if model == 'free':
        cis_V = np.diagonal(out['free']['V'], axis1=1, axis2=2).mean(axis=1)
        trans_V = np.diagonal(out['free']['V_b'], axis1=1, axis2=2).mean(axis=1)
        W = np.diagonal(out['free']['W'], axis1=1, axis2=2).mean(axis=1)
    else:
        cis_V = out['iid']['V']
        trans_V = out['iid']['V_b']
        W = out['iid']['W']
    
    bio_var = cis_hom_g2 + cis_V + trans_hom_g2 + trans_V + hom_e2 + W

    if filter is None:
        filter = bio_var > 0
    else:
        filter = filter & (bio_var > 0)
    print(f'{filter.sum()} genes after filtering')

    cis_hom_g2 = cis_hom_g2[filter]
    trans_hom_g2 = trans_hom_g2[filter]
    hom_e2 = hom_e2[filter]
    cis_V = cis_V[filter]
    trans_V = trans_V[filter]
    W = W[filter]
    bio_var = bio_var[filter]


   # transcriptome-wide h2
    ratio_of_mean = lambda x, y: x.mean() / y.mean()
    median = lambda x, y: np.median(x / y)
    mean = lambda x, y: np.mean(x / y)
    mean_clip = lambda x, y: np.mean(np.clip((x / y), -100, 100))
    if h2 == 'ratio':
        func = ratio_of_mean
    elif h2 == 'median':
        func = median
    elif h2 =='mean':
        func = mean
    elif h2 =='mean_clip':
        func = mean_clip

    cis_bio_hom_h2 = func(cis_hom_g2, bio_var)
    cis_bio_ct_h2 = func(cis_V, bio_var)
    trans_bio_hom_h2 = func(trans_hom_g2, bio_var)
    trans_bio_ct_h2 = func(trans_V, bio_var)
    print(f"Cis shared: {cis_bio_hom_h2}; Cis specific median: {cis_bio_ct_h2};")
    print(f"Cis specific (%): {cis_bio_ct_h2 / (cis_bio_ct_h2 + cis_bio_hom_h2)}")
    print(f"Trans shared median: {trans_bio_hom_h2}; Trans specific median: {trans_bio_ct_h2};")
    print(f"Tans specific (%): {trans_bio_ct_h2 / (trans_bio_ct_h2 + trans_bio_hom_h2)}")
    cis_bio_total_h2 = cis_bio_hom_h2 + cis_bio_ct_h2
    trans_bio_total_h2 = trans_bio_hom_h2 + trans_bio_ct_h2
    print(f"Trans (%): {trans_bio_total_h2 / (cis_bio_total_h2 + trans_bio_total_h2)}")


    ## ci
    cis_bio_hom_h2_ci = stats.bootstrap((cis_hom_g2, bio_var), func, 
                                vectorized=False, paired=True, random_state=rng
                                ).confidence_interval
                                                
    cis_bio_ct_h2_ci = stats.bootstrap((cis_V, bio_var), func, 
                                    vectorized=False, paired=True,
                                    random_state=rng).confidence_interval  
    
    trans_bio_hom_h2_ci = stats.bootstrap((trans_hom_g2, bio_var), func, 
                                vectorized=False, paired=True, random_state=rng
                                ).confidence_interval
                                                
    trans_bio_ct_h2_ci = stats.bootstrap((trans_V, bio_var), func, 
                                        vectorized=False, paired=True, 
                                        random_state=rng).confidence_interval  
    

    # plot
    bottom = np.zeros(2)
    values = np.array([cis_bio_hom_h2, trans_bio_hom_h2]) 
    ci = np.array([cis_bio_hom_h2_ci, trans_bio_hom_h2_ci])
    if np.all(values > 0):
        yerr = np.abs(ci - values[:, np.newaxis])
        print(yerr)
    else:
        print(values)
        print(ci)
        ax.text(0.01, 1.01, 'Negative heritability', ha='left', va='bottom', transform=ax.transAxes)
        return None, None
    
    ax.bar(labels, values, width, yerr=yerr.T, capsize=5,
        bottom=bottom, color=colors[0], edgecolor='none', label='Shared')
    
    bottom += [cis_bio_hom_h2, trans_bio_hom_h2]
    values = np.array([cis_bio_ct_h2, trans_bio_ct_h2])
    print(bottom, values)
    ci = np.array([cis_bio_ct_h2_ci, trans_bio_ct_h2_ci])
    yerr = np.abs(ci - values[:, np.newaxis])
    ax.bar(labels, values, width, yerr=yerr.T, capsize=5,
            bottom=bottom, color=colors[1], edgecolor='none', label='Specific')
    
    # permutation test for comparing cis vs trans in median(specific h2) / (median(shared h2) + median(specific h2))
    def diff(shared_x, specific_x, shared_y, specific_y):
        x = np.median(specific_x) / (np.median(specific_x) + np.median(shared_x))
        y = np.median(specific_y) / (np.median(specific_y) + np.median(shared_y))
        return y - x

    def permutation_test(shared_x, specific_x, shared_y, specific_y, n_perm=9999):
        diff_obs = diff(shared_x, specific_x, shared_y, specific_y)
        diff_perm = np.zeros(n_perm)
        x = np.stack([shared_x, specific_x], axis=-1)
        y = np.stack([shared_y, specific_y], axis=-1)
        for i in range(n_perm):
            choice = rng.choice([True, False], size=x.shape[0])
            x_perm = np.where(choice[:, np.newaxis], x, y)
            y_perm = np.where(choice[:, np.newaxis], y, x)
            diff_perm[i] = diff(x_perm[:, 0], x_perm[:, 1], y_perm[:, 0], y_perm[:, 1])
        p = ((diff_perm >= diff_obs).sum() + 1) / (n_perm + 1)
        return diff_obs, p

    dif, p = permutation_test(cis_hom_g2 / bio_var, cis_V / bio_var, trans_hom_g2 / bio_var, trans_V / bio_var)

    return dif, p


def op_h2_plot(out: dict, P: pd.DataFrame, ax: plt.Axes, colors: list, 
               labels: list = ['cis', 'trans'],
               rng: np.random.Generator = np.random.default_rng(42), 
               width: float=0.2, permutation: bool=False, n_permutations=9999,
               h2: str='ratio') -> None:
    """
    Plot the cis vs trans heritability for OP.

    Parameters:
        out: output dictionary 
        P: pandas dataframe of cell type proportions
        ax: matplotlib axis
        colors: list of colors for shared and specific heritability
        labels: labels for e.g. cis vs trans
        rng: random number generator
        width: bar width
        permutation: whether to perform permutation test
        n_permutations: number of permutations to perform
        h2: 'ratio' or 'median' across transcriptome
    """

    # functions
    def cis_trans_diff(x, y):
        x_mean = x.mean(axis=0)
        x_ratio = x_mean[1] / x_mean.sum()
        y_mean = y.mean(axis=0)
        y_ratio = y_mean[1] / y_mean.sum()
        return x_ratio - y_ratio
        
    ratio_of_mean = lambda x, y: 100 * x.mean() / y.mean()
    median = lambda x, y: 100 * np.median(x / y)

    # cell type proportions
    pi = P.mean()
    S = P.cov()

    # vc
    beta = out['free']['ct_beta']
    vc_beta = np.diag(beta @ S @ beta.T)

    vc_cis_V = [np.trace(V @ S) + pi @ V @ pi for V in out['free']['V']]
    vc_cis_V = np.array(vc_cis_V)

    vc_trans_V = [np.trace(V @ S) + pi @ V @ pi for V in out['free']['V_b']]
    vc_trans_V = np.array(vc_trans_V)

    vc_W = [np.trace(W @ S) + pi @ W @ pi for W in out['free']['W']]

    cis_hom_g2 = out['free']['hom_g2']
    trans_hom_g2 = out['free']['hom_g2_b']
    hom_e2 = out['free']['hom_e2']

    bio_var = vc_beta + cis_hom_g2 + trans_hom_g2 + hom_e2 + vc_cis_V + vc_trans_V + vc_W

    if h2 == 'ratio':
        func = ratio_of_mean
    elif h2 == 'median':
        func = median

    # h2
    cis_bio_hom_h2 = func(cis_hom_g2, bio_var)
    cis_bio_ct_h2 = func(vc_cis_V, bio_var)
    trans_bio_hom_h2 = func(trans_hom_g2, bio_var)
    trans_bio_ct_h2 = func(vc_trans_V, bio_var)

    ## ci
    cis_bio_hom_h2_ci = stats.bootstrap((cis_hom_g2, bio_var), func, 
                                vectorized=False, paired=True, random_state=rng
                                ).confidence_interval
                                                
    cis_bio_ct_h2_ci = stats.bootstrap((vc_cis_V, bio_var), 
                                func, vectorized=False, paired=True, 
                                random_state=rng).confidence_interval  
    
    trans_bio_hom_h2_ci = stats.bootstrap((trans_hom_g2, bio_var), func, 
                                vectorized=False, paired=True, random_state=rng
                                ).confidence_interval
                                                
    trans_bio_ct_h2_ci = stats.bootstrap((vc_trans_V, bio_var), 
                                func, vectorized=False, paired=True, 
                                random_state=rng).confidence_interval  
    
    if permutation:
        ## permutation for p value cis - trans
        cis = np.stack((cis_hom_g2, vc_cis_V), axis=-1)
        trans = np.stack((trans_hom_g2, vc_trans_V), axis=-1)

        statistics = []
        for i in range(n_permutations):
            choice = rng.choice([True, False], size=cis.shape[0])
            perm_cis = np.where(choice[:, np.newaxis], cis, trans)
            perm_trans = np.where(choice[:, np.newaxis], trans, cis)
            if i == 0:
                print(cis[:5, :2])
                print(trans[:5, :2])
                print(perm_cis[:5, :2])
                print(perm_trans[:5, :2])
            statistics.append(cis_trans_diff(perm_cis, perm_trans))
        p = ((np.array(statistics) >= cis_trans_diff(cis, trans)).sum() + 1) / (n_permutations + 1)


    # plot
    bottom = np.zeros(2)
    values = np.array([cis_bio_hom_h2, trans_bio_hom_h2])
    ci = np.array([cis_bio_hom_h2_ci, trans_bio_hom_h2_ci])
    yerr = np.abs(ci - values[:, np.newaxis])
    ax.bar(labels, values, width, yerr=yerr.T, capsize=5,
        bottom=bottom, color='0.7', label='Shared')
    
    bottom += [cis_bio_hom_h2, trans_bio_hom_h2]
    values = np.array([cis_bio_ct_h2, trans_bio_ct_h2])
    ci = np.array([cis_bio_ct_h2_ci, trans_bio_ct_h2_ci])
    yerr = np.abs(ci - values[:, np.newaxis])
    ax.bar(labels, values, width, yerr=yerr.T, capsize=5,
            bottom=bottom, color=colors[0], ecolor=colors[1], label='Specific')
    # ax.text(width / 2 + 0.1, cis_bio_hom_h2 / 2, 
    #         f'{100 * (1 - cis_bio_ct_h2 / (cis_bio_hom_h2 + cis_bio_ct_h2)):.1f}%', 
    #         color='0.7', fontsize=10, ha='left', va='center')
    # ax.text(width / 2 + 0.1, cis_bio_hom_h2 + cis_bio_ct_h2 / 2, 
    #         f'{100 * (cis_bio_ct_h2 / (cis_bio_hom_h2 + cis_bio_ct_h2)):.1f}%', 
    #         color=colors[0], fontsize=10, ha='left', va='center')
    # ax.text(1 - width / 2 - 0.1, trans_bio_hom_h2 / 2, 
    #         f'{100 * (1 - trans_bio_ct_h2 / (trans_bio_hom_h2 + trans_bio_ct_h2)):.1f}%', 
    #         color='0.7', fontsize=10, ha='right', va='center')
    # ax.text(1 - width / 2 - 0.1, trans_bio_hom_h2 + trans_bio_ct_h2 / 2, 
    #         f'{100 * (trans_bio_ct_h2 / (trans_bio_hom_h2 + trans_bio_ct_h2)):.1f}%', 
    #         color=colors[0], fontsize=10, ha='right', va='center')
    
    if permutation:
        ax.text(0.95, 0.98, f'p(specificity differentiation) = {p:.2e}', 
                transform=ax.transAxes, ha='right', va='top', fontsize=10)