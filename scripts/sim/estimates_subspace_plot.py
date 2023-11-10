import os, re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from gxctmm import plot


def read_hom(res, arg):
    '''
    read Hom model estimates
    res = e.g. out['reml']['hom']
    '''
    return (pd.DataFrame({
        'hom_g2': res['hom_g2'],
        'hom_e2': res['hom_e2'],
        'arg': arg
    }))


def read_iid(res, arg):
    # read IID model estimates
    return (pd.DataFrame({
        'hom_g2': res['hom_g2'],
        'hom_e2': res['hom_e2'],
        'V': res['V'][:, 0, 0],
        'W': res['W'][:, 0, 0],
        'arg': arg
    }))


def read_free(res, arg):
    # read Free model estimates
    hom_g2, hom_e2, V, W = np.array(res['hom_g2']), np.array(res['hom_e2']), res['V'], res['W']
    C = V.shape[1]
    V = np.diagonal(V, axis1=1, axis2=2)
    W = np.diagonal(W, axis1=1, axis2=2)
    h2 = (hom_g2.reshape(-1, 1) + V) / (hom_g2.reshape(-1, 1) + V + hom_e2.reshape(-1, 1) + W)
    V, W, h2 = V.T, W.T, h2.T
    data = pd.DataFrame({'hom_g2': hom_g2, 'hom_e2': hom_e2, 'arg': arg})
    for i in range(C):
        data['V_' + str(i + 1)] = V[i]
        data['W_' + str(i + 1)] = W[i]
        data[f'h2_{i + 1}'] = h2[i]
    return (data)


def read_full(res, arg):
    # read Full model estimates
    V, W = res['V'], res['W']
    C = V.shape[1]
    V = np.diagonal(V, axis1=1, axis2=2)
    W = np.diagonal(W, axis1=1, axis2=2)
    h2 = V / (V + W)
    V, W, h2 = V.T, W.T, h2.T
    data = pd.DataFrame({})
    for i in range(C):
        data['V_' + str(i + 1)] = V[i]
        data['W_' + str(i + 1)] = W[i]
        data[f'h2_{i + 1}'] = h2[i]
    data['arg'] = arg
    return (data)


def plot_hom(data, true_hom_g2, true_hom_e2, ax, boxcolor, pointcolor):
    '''
    Not finished
    '''
    ## Hom model
    ### variance components
    sns.boxplot(x='arg', y='hom_g2', data=data, ax=ax,
                color=boxcolor)

    ### add True variances to plot
    xs = plot.snsbox_get_x(len(np.unique(data['arg'])), 1)
    ax.scatter(xs, vcs, color=pointcolor, zorder=10)

    ax.xaxis.label.set_visible(False)
    ax.set_ylabel('$\sigma_{g}^2$')
    # ax.text(-0.29, 0.95, 'Hom', fontsize=16, transform=axes[0,0].transAxes)
    # ax.text(-0.05, 1.05, '(A)', fontsize=12, transform=axes[0,0].transAxes)
    # handles, labels = axes[0,0].get_legend_handles_labels()
    # axes[0,0].legend(handles=handles, labels=labels)

    # axes[0,1].axis('off')
    # axes[0,2].axis('off')


def plot_iid(data, true_hom_g2, true_hom_e2, true_V, true_W, ax, boxcolor, pointcolor):
    '''
    Not finished
    '''
    ## IID model
    ### variance components
    sns.boxplot(x='arg', y='subject', data=summaries['iid'], ax=axes[1, 0],
                color=mycolors[0])

    ### add True variances to plot
    xs = plot.snsbox_get_x(len(np.unique(summaries['iid']['arg'])), 1)
    axes[1, 0].scatter(xs, vcs, color=pointcolor, zorder=10)

    axes[1, 0].xaxis.label.set_visible(False)
    axes[1, 0].set_ylabel('$\sigma_{hom}^2$')
    axes[1, 0].text(-0.29, 0.95, 'IID', fontsize=16, transform=axes[1, 0].transAxes)
    axes[1, 0].text(-0.05, 1.05, '(B)', fontsize=12, transform=axes[1, 0].transAxes)
    # handles, labels = axes[1,0].get_legend_handles_labels()
    # axes[1,0].legend(handles=handles, labels=labels)

    ### V
    sns.boxplot(x='arg', y='V', data=summaries['iid'], ax=axes[1, 1],
                color=mycolors[0])
    #### add true V
    trueV_ = np.array([np.diag(trueV[x]) for x in pd.unique(summaries['iid']['arg'])]).T.flatten()
    xs = plot.snsbox_get_x(len(np.unique(summaries['iid']['arg'])), 1)
    xs = list(xs) * C
    axes[1, 1].scatter(xs, trueV_, color=pointcolor, zorder=10)
    axes[1, 1].xaxis.label.set_visible(False)
    axes[1, 1].set_ylabel('V_diag (cell type-specific genetic variance)')
    axes[1, 1].text(-0.05, 1.05, '(C)', fontsize=12, transform=axes[1, 1].transAxes)

    ### W
    sns.boxplot(x='arg', y='W', data=summaries['iid'], ax=axes[1, 2],
                color=mycolors[0])
    #### add true W
    trueW_ = np.array([np.diag(trueW[x]) for x in pd.unique(summaries['iid']['arg'])]).T.flatten()
    xs = plot.snsbox_get_x(len(np.unique(summaries['iid']['arg'])), 1)
    xs = list(xs) * C
    axes[1, 2].scatter(xs, trueW_, color=pointcolor, zorder=10)
    axes[1, 2].xaxis.label.set_visible(False)
    axes[1, 2].set_ylabel('W_diag (cell type-specific noise variance)')
    axes[1, 2].text(-0.05, 1.05, '(C)', fontsize=12, transform=axes[1, 2].transAxes)



def plot_free(data, C, true_hom_g2s, true_hom_e2s, true_V, true_W, true_h2,
              boxcolor, pointcolor, colorpalette, fig_f, xlab=None):
    """
    data:   with columns of arg hom_g2, hom_e2, V_1-V_C
    boxcolor:   consistent color for g2 and  e2
    """

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(15, 10), sharex=True)
    # hom_g2
    sns.boxplot(x='arg', y='hom_g2', data=data, ax=axes[0, 0], color=boxcolor) # TODO: change to violin to compare HE vs REML

    ## add True variances
    xs = plot.snsbox_get_x(len(np.unique(data['arg'])), 1)
    axes[0, 0].scatter(xs, true_hom_g2s, color=pointcolor, zorder=10)
    # ax.xaxis.label.set_visible(False)
    axes[0, 0].set_ylabel('$\sigma_{g}^2$', fontsize=16)
    # ax.text(-0.29, 0.95, 'Free', fontsize=16, transform=axes[0,0].transAxes)
    # ax.text(-0.05, 1.05, '(D)', fontsize=12, transform=axes[2,0].transAxes)

    # hom_e2
    sns.boxplot(x='arg', y='hom_e2', data=data, ax=axes[0, 1], color=boxcolor)

    ## add True variances
    print(data)
    print(xs)
    print(true_hom_e2s)
    axes[0, 1].scatter(xs, true_hom_e2s, color=pointcolor, zorder=10)
    axes[0, 1].set_ylabel('$\sigma_{e}^2$', fontsize=16)

    # V
    V_df = pd.melt(data, id_vars=['arg'], value_vars=['V_' + str(i + 1) for i in range(C)])
    V_df['variable'] = V_df['variable'].str.replace('V', 'CT')
    V_df['value'] = plot._clip(V_df['value'], cut=3)
    sns.boxplot(x='arg', y='value', hue='variable', data=V_df, ax=axes[1, 0], palette=colorpalette)
    ## add true V
    V_diag = np.array([np.diag(true_V[x]) for x in pd.unique(V_df['arg'])]).flatten()
    xs = plot.snsbox_get_x(len(np.unique(V_df['arg'])), len(np.unique(V_df['variable'])))
    axes[1, 0].scatter(xs, V_diag, color=pointcolor, zorder=10)
    # axes[0,2].xaxis.label.set_visible(False)
    axes[1, 0].set_ylabel('V', fontsize=16)
    # axes[0,2].text(-0.05, 1.05, '(E)', fontsize=12, transform=axes[2,1].transAxes)
    handles, labels = axes[1, 0].get_legend_handles_labels()
    axes[1, 0].legend(handles=handles, labels=[r'$%s$' % (x) for x in labels])

    # W
    W_df = pd.melt(data, id_vars=['arg'], value_vars=['W_' + str(i + 1) for i in range(C)])
    W_df['value'] = plot._clip(W_df['value'], cut=3)
    sns.boxplot(x='arg', y='value', hue='variable', data=W_df, ax=axes[1, 1], palette=colorpalette)
    ## add true W
    W_diag = np.array([np.diag(true_W[x]) for x in pd.unique(W_df['arg'])]).flatten()
    xs = plot.snsbox_get_x(len(np.unique(W_df['arg'])), len(np.unique(W_df['variable'])))
    axes[1, 1].scatter(xs, W_diag, color=pointcolor, zorder=10)
    # axes[0,3].xaxis.label.set_visible(False)
    axes[1, 1].set_ylabel('W', fontsize=16)
    # axes[0,3].text(-0.05, 1.05, '(E)', fontsize=12, transform=axes[0,3].transAxes)
    axes[1, 1].legend().set_visible(False)
    # handles, labels = axes[1, 1].get_legend_handles_labels()
    # axes[1, 1].legend(handles=handles, labels=[r'$%s$' % (x) for x in labels])

    # h2
    h2_df = pd.melt(data, id_vars=['arg'], value_vars=['h2_' + str(i + 1) for i in range(C)])
    h2_df['value'] = h2_df['value'].clip(-1, 1.5)
    sns.boxplot(x='arg', y='value', hue='variable', data=h2_df, ax=axes[1, 2], palette=colorpalette)
    ## add true h2
    h2_diag = np.array([np.diag(true_h2[x]) for x in pd.unique(h2_df['arg'])]).flatten()
    xs = plot.snsbox_get_x(len(np.unique(h2_df['arg'])), len(np.unique(h2_df['variable'])))
    axes[1, 2].scatter(xs, h2_diag, color=pointcolor, zorder=10)
    axes[1, 2].set_ylabel('h2', fontsize=16)
    axes[1, 2].legend().set_visible(False)
    # handles, labels = axes[1, 2].get_legend_handles_labels()
    # axes[1, 2].legend(handles=handles, labels=[r'$h^2_{CT%s}$' % (x.split('_')[1]) for x in labels])

    if xlab:
        for ax in axes.flatten():
            ax.set_xlabel(xlab)

    fig.tight_layout()
    fig.savefig(fig_f)


def plot_full(data, C, true_V, true_W, true_h2, pointcolor, colorpalette, fig_f):
    fig, axes = plt.subplots(ncols=3, figsize=(15, 5), sharex=True)

    # V
    V_df = pd.melt(data, id_vars=['arg'], value_vars=['V_' + str(i + 1) for i in range(C)])
    sns.boxplot(x='arg', y='value', hue='variable', data=V_df, ax=axes[0], palette=colorpalette)
    ## add true V
    V_diag = np.array([np.diag(true_V[x]) for x in pd.unique(V_df['arg'])]).flatten()
    xs = plot.snsbox_get_x(len(np.unique(V_df['arg'])), len(np.unique(V_df['variable'])))
    axes[0].scatter(xs, V_diag, color=pointcolor, zorder=10)
    axes[0].set_ylabel('V', fontsize=16)
    handles, labels = axes[0].get_legend_handles_labels()
    axes[0].legend(handles=handles, labels=[r'$%s$' % (x) for x in labels])

    # W
    W_df = pd.melt(data, id_vars=['arg'], value_vars=['W_' + str(i + 1) for i in range(C)])
    sns.boxplot(x='arg', y='value', hue='variable', data=W_df, ax=axes[1], palette=colorpalette)
    ## add true W
    W_diag = np.array([np.diag(true_W[x]) for x in pd.unique(W_df['arg'])]).flatten()
    xs = plot.snsbox_get_x(len(np.unique(W_df['arg'])), len(np.unique(W_df['variable'])))
    axes[1].scatter(xs, W_diag, color=pointcolor, zorder=10)
    axes[1].set_ylabel('W', fontsize=16)
    handles, labels = axes[1].get_legend_handles_labels()
    axes[1].legend(handles=handles, labels=[r'$%s$' % (x) for x in labels])

    # h2
    h2_df = pd.melt(data, id_vars=['arg'], value_vars=['h2_' + str(i + 1) for i in range(C)])
    h2_df['value'] = h2_df['value'].clip(lower=-1, upper=2)
    sns.boxplot(x='arg', y='value', hue='variable', data=h2_df, ax=axes[2], palette=colorpalette)
    ## add true h2
    h2_diag = np.array([np.diag(true_h2[x]) for x in pd.unique(h2_df['arg'])]).flatten()
    xs = plot.snsbox_get_x(len(np.unique(h2_df['arg'])), len(np.unique(h2_df['variable'])))
    axes[2].scatter(xs, h2_diag, color=pointcolor, zorder=10)
    axes[2].set_ylabel('h2', fontsize=16)
    handles, labels = axes[2].get_legend_handles_labels()
    axes[2].legend(handles=handles, labels=[r'$h^2_{CT%s}$' % (x.split('_')[1]) for x in labels])

    fig.tight_layout()
    fig.savefig(fig_f)


def main():
    # 
    os.makedirs(os.path.dirname(snakemake.output.png), exist_ok=True)
    subspace = snakemake.params.subspace
    sim_plot_order = snakemake.params.sim_plot_order
    method = snakemake.params.get('method', False)
    print(method)
    mycolors = snakemake.params.mycolors
    pointcolor = snakemake.params.pointcolor
    colorpalette = snakemake.params.colorpalette

    # get cell type number
    out = np.load(snakemake.input.out[0], allow_pickle=True).item()
    if method:
        if 'free' in out[method].keys():
            C = out[method]['free']['V'][0].shape[0]
        elif 'full' in out[method].keys():
            C = out[method]['full']['V'][0].shape[0]
        else:
            sys.exit('Missing C!\n')
    else:
        if 'free' in out.keys():
            C = out['free']['V'][0].shape[0]
        elif 'full' in out.keys():
            C = out['full']['V'][0].shape[0]
        else:
            sys.exit('Missing C!\n')

    # collect plot order
    plot_order = np.array(
        sim_plot_order[re.sub('\d', '', snakemake.wildcards.model)][snakemake.wildcards.arg])  # e.g. free3
    args = subspace[snakemake.wildcards.arg]
    if np.any(~args.isin(plot_order)):
        sys.exit('Missing arg in plot_order!\n')
    overlap = np.intersect1d(plot_order, args)
    plot_order = plot_order[np.isin(plot_order, overlap)]

    # collect model estimates
    data = []
    for arg, out_f in zip(np.array(subspace[snakemake.wildcards.arg]), snakemake.input.out):
        out = np.load(out_f, allow_pickle=True).item()
        if snakemake.wildcards.model in ['hom', 'free']:
            if method:
                data.append(read_free(out[method]['free'], arg))
            else:
                data.append(read_free(out['free'], arg))
        elif snakemake.wildcards.model == 'full':
            if method:
                data.append(read_full(out[method][snakemake.wildcards.model], arg))
            else:
                data.append(read_full(out[snakemake.wildcards.model], arg))

    print(data)
    
    # concat and sort
    data = pd.concat(data, ignore_index=True)
    data['arg'] = pd.Categorical(data['arg'], plot_order)
    data = data.sort_values('arg').reset_index(drop=True)
    # data.head().to_csv(sys.stdout, sep='\t')

    # collect true V
    true_V = {}
    for arg, trueV_f in zip(np.array(subspace[snakemake.wildcards.arg]), snakemake.input.V):
        # read true V values
        true_V[arg] = np.loadtxt(trueV_f)

    # collect true W
    true_W = {}
    for arg, trueW_f in zip(np.array(subspace[snakemake.wildcards.arg]), snakemake.input.W):
        # read true W values
        true_W[arg] = np.loadtxt(trueW_f)

    # collect True hom_g2 hom_e2 and h2
    hom_g2s = []
    hom_e2s = []
    h2 = {}
    subspace[snakemake.wildcards.arg] = pd.Categorical(subspace[snakemake.wildcards.arg], plot_order)
    subspace = subspace.sort_values(snakemake.wildcards.arg).reset_index(drop=True)
    # print( subspace.head() )
    for index, row in subspace.iterrows():
        arg, vc = row[snakemake.wildcards.arg], row['vc']
        hom_g2, hom_e2 = float(vc.split('_')[1]), float(vc.split('_')[2])
        hom_g2s.append(hom_g2)
        hom_e2s.append(hom_e2)
        h2[arg] = (hom_g2 + true_V[arg]) / (hom_g2 + true_V[arg] + hom_e2 + true_W[arg])

    # plot
    plt.rcParams.update({'font.size': 10})
    #    #### tweak x labels
    #    xlabel = subspace_plot.get_xlabel( snakemake.wildcards.arg )
    #    for ax in axes.flatten():
    #        subspace_plot.set_xtickname( fig, ax, snakemake.wildcards.arg )
    #        if len(xlabel) < 30:
    #            ax.set_xlabel(xlabel, fontsize=16)
    #        else:
    #            ax.set_xlabel(xlabel)
    print(data.head())
    if snakemake.wildcards.model in ['hom', 'free']:
        xlab = None
        if snakemake.wildcards.arg == 'ss':
            xlab = 'sample size'
        plot_free(data, C, hom_g2s, hom_e2s, true_V, true_W, h2, mycolors[0], pointcolor,
                  colorpalette, snakemake.output.png, xlab=xlab)
    elif snakemake.wildcards.model == 'full':
        plot_full(data, C, true_V, true_W, h2, pointcolor,
                  colorpalette, snakemake.output.png)


if __name__ == '__main__':
    main()
