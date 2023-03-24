import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import plot_help

def main():
    # par

    # get cell type number
    hom = np.load(snakemake.input.out[0], allow_pickle=True).item()['ml']['hom']
    C = len(hom['beta'][0])

    ## collect subspace
    #subspace = get_subspace(wildcards.arg,
    #        overall_params.loc[overall_params['model']==wildcards.model])
    subspace = snakemake.params.subspace 

    # collect estimates from Hom IID Free Full
    summaries = {'hom':[], 'iid':[], 'free':[], 'full':[]}
    for arg_, out_f in zip(np.array(subspace[snakemake.wildcards.arg]), snakemake.input.out):
        npy = np.load(out_f, allow_pickle=True).item()
        estimates = npy['ml']
        # read Hom model estimates
        hom = estimates['hom']
        hom2, beta, celltype_main_var, cellspecific_var = hom['hom2'], hom['beta'], hom['celltype_main_var'], hom['cellspecific_var']
        beta = beta.T
        summary = pd.DataFrame({'subject':hom2, 'celltype-main':celltype_main_var,
            'celltype-genetic':np.zeros(len(hom2)), 'celltype-noise':np.zeros(len(hom2)),
            'cell-specific':cellspecific_var})
        for i in range(C):
            summary['beta_'+str(i+1)] = beta[i]
        summary['arg'] = arg_
        summaries['hom'].append(summary)

        # read IID model estimates
        iid = estimates['iid']
        hom2, beta, V, W = iid['hom2'], iid['beta'], iid['V'], iid['W']
        celltype_main_var, interxn_var, celltype_noise, cellspecific_var = iid['celltype_main_var'], iid['interxn_var'], iid['celltype_noise_var'], iid['cellspecific_var']
        V = np.array([x[0,0] for x in V]) # extract the first diagonal of V
        W = np.array([x[0,0] for x in W]) # extract the first diagonal of W
        beta = beta.T
        summary = pd.DataFrame({'subject':hom2, 'celltype-main':celltype_main_var,
            'celltype-genetic':interxn_var, 'celltype-noise':celltype_noise, 
            'cell-specific':cellspecific_var})
        for i in range(C):
            summary['beta_'+str(i+1)] = beta[i]
        summary['V'] = V
        summary['W'] = W
        summary['arg'] = arg_
        summaries['iid'].append(summary)

        # read Free model estimates
        free = estimates['free']
        hom2, beta, V, W = free['hom2'], free['beta'], free['V'], free['W']
        celltype_main_var, interxn_var, celltype_noise, cellspecific_var = free['celltype_main_var'], free['interxn_var'], free['celltype_noise_var'], free['cellspecific_var']
        V = np.array([np.diag(x) for x in V]) # extract out diagonal of V
        W = np.array([np.diag(x) for x in W]) # extract out diagonal of V
        beta, V, W = beta.T, V.T, W.T
        summary = pd.DataFrame({'subject':hom2, 'celltype-main':celltype_main_var,
            'celltype-genetic':interxn_var, 'celltype-noise':celltype_noise, 'cell-specific':cellspecific_var})
        for i in range(C):
            summary['beta_'+str(i+1)] = beta[i]
            summary['V_'+str(i+1)] = V[i]
            summary['W_'+str(i+1)] = W[i]
        summary['arg'] = arg_
        summaries['free'].append(summary)

        # read Full model estimates
        full = estimates['full']
        hom2, beta, V, W = full['hom2'], full['beta'], full['V'], full['W']
        celltype_main_var, interxn_var, celltype_noise, cellspecific_var = full['celltype_main_var'], full['interxn_var'], full['celltype_noise_var'], full['cellspecific_var']
        V_diag = np.array([np.diag(x) for x in V])
        V_tril = np.array([x[np.tril_indices(C, k=-1)].flatten() for x in V])
        W_diag = np.array([np.diag(x) for x in W])
        W_tril = np.array([x[np.tril_indices(C, k=-1)].flatten() for x in W])
        beta, V_diag, V_tril, W_diag, W_tril = beta.T, V_diag.T, V_tril.T, W_diag.T, W_tril.T
        summary = pd.DataFrame({'subject':hom2, 'celltype-main':celltype_main_var,
            'celltype-genetic':interxn_var, 'celltype-noise':celltype_noise, 'cell-specific':cellspecific_var})
        for i in range(C):
            summary['beta_'+str(i+1)] = beta[i]
            summary['V_'+str(i+1)] = V_diag[i]
            summary['W_'+str(i+1)] = W_diag[i]
        for j in range(V_tril.shape[0]):
            summary['Vlow_'+str(j+1)] = V_tril[j]
            summary['Wlow_'+str(j+1)] = W_tril[j]
        summary['arg'] = arg_
        summaries['full'].append(summary)

    # concat and sort
    for model in ['hom', 'iid', 'free', 'full']:
        summaries_ = pd.concat(summaries[model], ignore_index=True)
        plot_order_ = np.array(snakemake.params.og_plot_order[snakemake.wildcards.model][snakemake.wildcards.arg])
        plot_order_ = plot_order_[np.isin(plot_order_, summaries_['arg'])]
        summaries_['arg'] = pd.Categorical(summaries_['arg'], plot_order_)
        summaries_ = summaries_.sort_values('arg').reset_index(drop=True)
        summaries_.to_csv(sys.stdout, sep='\t')
        summaries[model] = summaries_

    # collect true beta and V and W
    truebeta = {}
    trueV = {}
    trueW = {}
    for arg_, truebeta_f, trueV_f, trueW_f in zip(
            np.array(subspace[snakemake.wildcards.arg]), snakemake.input.beta, 
            snakemake.input.V, snakemake.input.W):
        # read true beta values
        truebeta[arg_] = np.loadtxt(truebeta_f)
        # read true V values
        trueV[arg_] = np.loadtxt(trueV_f)
        # read true W values
        trueW[arg_] = np.loadtxt(trueW_f)

    # collect True variances explained
    vcs = []
    subspace_ = subspace.copy()  # may don't need to copy
    plot_order_ = np.array(snakemake.params.og_plot_order[snakemake.wildcards.model][snakemake.wildcards.arg])
    plot_order_ = plot_order_[np.isin(plot_order_, subspace_[snakemake.wildcards.arg])]
    subspace_[snakemake.wildcards.arg] = pd.Categorical(subspace_[snakemake.wildcards.arg], plot_order_)
    subspace_ = subspace_.sort_values(snakemake.wildcards.arg).reset_index(drop=True)
    print(subspace_)
    for vc in np.array(subspace_['vc']):
        print(vc)
        vcs = vcs+[float(x) for x in vc.split('_')]

    # plot
    fig, axes = plt.subplots(nrows=4, ncols=6, figsize=(24, 16), sharex=True)
    ## Hom model
    ### variance components
    var = pd.melt(summaries['hom'], id_vars=['arg'],
            value_vars=['subject', 'celltype-main', 'celltype-genetic', 'celltype-noise', 
                'cell-specific']
            )
    var.to_csv(sys.stdout, sep='\t')
    sns.boxplot(x='arg', y='value', hue='variable', data=var, ax=axes[0,0], 
            palette=snakemake.params.colorpalette)

    ### add True variances to plot
    xs = plot_help.snsbox_get_x(len(np.unique(var['arg'])), len(np.unique(var['variable'])))
    axes[0,0].scatter(xs, vcs, color=snakemake.params.pointcolor, zorder=10)

    axes[0,0].xaxis.label.set_visible(False)
    axes[0,0].set_ylabel('Proportion of variance')
    axes[0,0].text(-0.29, 0.95, 'Hom', fontsize=20, transform=axes[0,0].transAxes)
    axes[0,0].text(-0.05, 1.05, '(A)', fontsize=16, transform=axes[0,0].transAxes)
    handles, labels = axes[0,0].get_legend_handles_labels()
    axes[0,0].legend(handles=handles, labels=labels)

    ### beta
    beta_df = pd.melt(summaries['hom'], id_vars=['arg'], 
            value_vars=['beta_'+str(i+1) for i in range(C)])
    sns.boxplot(x='arg', y='value', hue='variable', data=beta_df, ax=axes[0,1], 
            palette=snakemake.params.colorpalette)
    #### add true beta
    truebeta_ = np.array([truebeta[x] for x in pd.unique(beta_df['arg'])]).flatten()
    xs = plot_help.snsbox_get_x(len(np.unique(beta_df['arg'])), 
            len(np.unique(beta_df['variable'])))
    axes[0,1].scatter(xs, truebeta_, color=snakemake.params.pointcolor, zorder=10)

    axes[0,1].xaxis.label.set_visible(False)
    axes[0,1].set_ylabel(r'$\beta$ (cell type main effect)')
    axes[0,1].text(-0.05, 1.05, '(B)', fontsize=16, transform=axes[0,1].transAxes)
    handles, labels = axes[0,1].get_legend_handles_labels()
    axes[0,1].legend(handles=handles, labels=[r'$\%s$'%(x) for x in labels])
    axes[0,2].axis('off')
    axes[0,3].axis('off')
    axes[0,4].axis('off')
    axes[0,5].axis('off')

    ## IID model
    ### variance components
    var = pd.melt(summaries['iid'], id_vars=['arg'],
            value_vars=['subject', 'celltype-main', 'celltype-genetic', 'celltype-noise', 
                'cell-specific']
            )
    var.to_csv(sys.stdout, sep='\t')
    sns.boxplot(x='arg', y='value', hue='variable', data=var, ax=axes[1,0], 
            palette=snakemake.params.colorpalette)

    ### add True variances to plot
    xs = plot_help.snsbox_get_x(len(np.unique(var['arg'])), len(np.unique(var['variable'])))
    axes[1,0].scatter(xs, vcs, color=snakemake.params.pointcolor, zorder=10)

    axes[1,0].xaxis.label.set_visible(False)
    axes[1,0].set_ylabel('Proportion of variance')
    axes[1,0].text(-0.29, 0.95, 'IID', fontsize=20, transform=axes[1,0].transAxes)
    axes[1,0].text(-0.05, 1.05, '(C)', fontsize=16, transform=axes[1,0].transAxes)
    handles, labels = axes[1,0].get_legend_handles_labels()
    axes[1,0].legend(handles=handles, labels=labels)

    ### beta
    beta_df = pd.melt(summaries['iid'], id_vars=['arg'], 
            value_vars=['beta_'+str(i+1) for i in range(C)])
    #beta_dfs.to_csv(sys.stdout, sep='\t')
    sns.boxplot(x='arg', y='value', hue='variable', data=beta_df, ax=axes[1,1], 
            palette=snakemake.params.colorpalette)
    #### add true beta
    truebeta_ = np.array([truebeta[x] for x in pd.unique(beta_df['arg'])]).flatten()
    xs = plot_help.snsbox_get_x(len(np.unique(beta_df['arg'])), 
            len(np.unique(beta_df['variable'])))
    axes[1,1].scatter(xs, truebeta_, color=snakemake.params.pointcolor, zorder=10)

    axes[1,1].xaxis.label.set_visible(False)
    axes[1,1].set_ylabel(r'$\beta$ (cell type main effect)')
    axes[1,1].text(-0.05, 1.05, '(D)', fontsize=16, transform=axes[1,1].transAxes)
    handles, labels = axes[1,1].get_legend_handles_labels()
    axes[1,1].legend(handles=handles, labels=[r'$\%s$'%(x) for x in labels])

    ### V
    sns.boxplot(x='arg', y='V', data=summaries['iid'], ax=axes[1,2], 
            color=snakemake.params.mycolors[0])
    #### add true V
    trueV_ = np.array([np.diag(trueV[x]) for x in 
        pd.unique(summaries['iid']['arg'])]).T.flatten()
    xs = plot_help.snsbox_get_x(len(np.unique(summaries['iid']['arg'])), 1)
    xs = list(xs) * C
    axes[1,2].scatter(xs, trueV_, color=snakemake.params.pointcolor, zorder=10)
    axes[1,2].xaxis.label.set_visible(False)
    axes[1,2].set_ylabel('V_diag (cell type-specific genetic variance)')
    axes[1,2].text(-0.05, 1.05, '(E)', fontsize=16, transform=axes[1,2].transAxes)

    ### W
    sns.boxplot(x='arg', y='W', data=summaries['iid'], ax=axes[1,3], 
            color=snakemake.params.mycolors[0])
    #### add true V
    trueW_ = np.array([np.diag(trueW[x]) for x in 
        pd.unique(summaries['iid']['arg'])]).T.flatten()
    xs = plot_help.snsbox_get_x(len(np.unique(summaries['iid']['arg'])), 1)
    xs = list(xs) * C
    axes[1,3].scatter(xs, trueW_, color=snakemake.params.pointcolor, zorder=10)
    axes[1,3].xaxis.label.set_visible(False)
    axes[1,3].set_ylabel('W_diag (cell type-specific noise variance)')
    axes[1,3].text(-0.05, 1.05, '(F)', fontsize=16, transform=axes[1,3].transAxes)
    axes[1,4].axis('off')
    axes[1,5].axis('off')

    ## Free model
    ### variance components
    var = pd.melt(summaries['free'], id_vars=['arg'],
            value_vars=['subject', 'celltype-main', 'celltype-genetic', 'celltype-noise', 
                'cell-specific']
            )
    var.to_csv(sys.stdout, sep='\t')
    sns.boxplot(x='arg', y='value', hue='variable', data=var, ax=axes[2,0], 
            palette=snakemake.params.colorpalette)

    #### add True variances
    xs = plot_help.snsbox_get_x(len(np.unique(var['arg'])), len(np.unique(var['variable'])))
    axes[2,0].scatter(xs, vcs, color=snakemake.params.pointcolor, zorder=10)

    axes[2,0].xaxis.label.set_visible(False)
    axes[2,0].set_ylabel('Proportion of variance')
    axes[2,0].text(-0.29, 0.95, 'Free', fontsize=20, transform=axes[2,0].transAxes)
    axes[2,0].text(-0.05, 1.05, '(G)', fontsize=16, transform=axes[2,0].transAxes)
    handles, labels = axes[2,0].get_legend_handles_labels()
    axes[2,0].legend(handles=handles, labels=labels)

    ### beta
    beta_df = pd.melt(summaries['free'], id_vars=['arg'], 
            value_vars=['beta_'+str(i+1) for i in range(C)])
    sns.boxplot(x='arg', y='value', hue='variable', data=beta_df, ax=axes[2,1], 
            palette=snakemake.params.colorpalette)
    #### add true beta
    truebeta_ = np.array([truebeta[x] for x in pd.unique(beta_df['arg'])]).flatten()
    xs = plot_help.snsbox_get_x(len(np.unique(beta_df['arg'])), 
            len(np.unique(beta_df['variable'])))
    axes[2,1].scatter(xs, truebeta_, color=snakemake.params.pointcolor, zorder=10)

    axes[2,1].xaxis.label.set_visible(False)
    axes[2,1].set_ylabel(r'$\beta$ (cell type main effect)')
    axes[2,1].text(-0.05, 1.05, '(H)', fontsize=16, transform=axes[2,1].transAxes)
    handles, labels = axes[2,1].get_legend_handles_labels()
    axes[2,1].legend(handles=handles, labels=[r'$\%s$'%(x) for x in labels])

    ### V
    V_df =  pd.melt(summaries['free'], id_vars=['arg'], 
            value_vars=['V_'+str(i+1) for i in range(C)])
    sns.boxplot(x='arg', y='value', hue='variable', data=V_df, ax=axes[2,2], 
            palette=snakemake.params.colorpalette)
    #### add true sig gam
    trueV_ = np.array([np.diag(trueV[x]) for x in pd.unique(V_df['arg'])]).flatten()
    xs = plot_help.snsbox_get_x(len(np.unique(V_df['arg'])), len(np.unique(V_df['variable'])))
    axes[2,2].scatter(xs, trueV_, color=snakemake.params.pointcolor, zorder=10)
    axes[2,2].xaxis.label.set_visible(False)
    axes[2,2].set_ylabel('V_diag (cell type-specific genetic variance')
    axes[2,2].text(-0.05, 1.05, '(I)', fontsize=16, transform=axes[2,2].transAxes)
    handles, labels = axes[2,2].get_legend_handles_labels()
    axes[2,2].legend(handles=handles, labels=[r'$%s$'%(x) for x in labels])

    ### W
    W_df =  pd.melt(summaries['free'], id_vars=['arg'], 
            value_vars=['W_'+str(i+1) for i in range(C)])
    sns.boxplot(x='arg', y='value', hue='variable', data=W_df, ax=axes[2,3], 
            palette=snakemake.params.colorpalette)
    #### add true sig gam
    trueW_ = np.array([np.diag(trueW[x]) for x in pd.unique(W_df['arg'])]).flatten()
    xs = plot_help.snsbox_get_x(len(np.unique(W_df['arg'])), len(np.unique(W_df['variable'])))
    axes[2,3].scatter(xs, trueW_, color=snakemake.params.pointcolor, zorder=10)
    axes[2,3].xaxis.label.set_visible(False)
    axes[2,3].set_ylabel('W_diag (cell type-specific noise variance')
    axes[2,3].text(-0.05, 1.05, '(J)', fontsize=16, transform=axes[2,3].transAxes)
    handles, labels = axes[2,3].get_legend_handles_labels()
    axes[2,3].legend(handles=handles, labels=[r'$%s$'%(x) for x in labels])
    axes[2,4].axis('off')
    axes[2,5].axis('off')

    ## Full model
    ### variance components
    var = pd.melt(summaries['full'], id_vars=['arg'],
            value_vars=['subject', 'celltype-main', 'celltype-genetic', 'celltype-noise', 
                'cell-specific']
            )
    sns.boxplot(x='arg', y='value', hue='variable', data=var, ax=axes[3,0], 
            palette=snakemake.params.colorpalette)

    ### add True variances to plot
    xs = plot_help.snsbox_get_x(len(np.unique(var['arg'])), len(np.unique(var['variable'])))
    axes[3,0].scatter(xs, vcs, color=snakemake.params.pointcolor, zorder=10)

    axes[3,0].set_xlabel(snakemake.wildcards.arg)
    axes[3,0].set_ylabel('Proportion of variance')
    axes[3,0].text(-0.29, 0.95, 'Full', fontsize=20, transform=axes[3,0].transAxes)
    axes[3,0].text(-0.05, 1.05, '(K)', fontsize=16, transform=axes[3,0].transAxes)
    handles, labels = axes[3,0].get_legend_handles_labels()
    axes[3,0].legend(handles=handles, labels=labels)
    ### beta
    beta_df = pd.melt(summaries['full'], id_vars=['arg'], 
            value_vars=['beta_'+str(i+1) for i in range(C)])
    sns.boxplot(x='arg', y='value', hue='variable', data=beta_df, ax=axes[3,1], 
            palette=snakemake.params.colorpalette)
    #### add true beta
    truebeta_ = np.array([truebeta[x] for x in pd.unique(beta_df['arg'])]).flatten()
    xs = plot_help.snsbox_get_x(len(np.unique(beta_df['arg'])), 
            len(np.unique(beta_df['variable'])))
    axes[3,1].scatter(xs, truebeta_, color=snakemake.params.pointcolor, zorder=10)

    axes[3,1].set_xlabel(snakemake.wildcards.arg)
    axes[3,1].set_ylabel(r'$\beta$ (main effect of cell type)')
    axes[3,1].text(-0.05, 1.05, '(L)', fontsize=16, transform=axes[3,1].transAxes)
    handles, labels = axes[3,1].get_legend_handles_labels()
    axes[3,1].legend(handles=handles, labels=[r'$\%s$'%(x) for x in labels])
    ### V
    #### diag
    V_df =  pd.melt(summaries['full'], id_vars=['arg'], 
            value_vars=['V_'+str(i+1) for i in range(C)])
    sns.boxplot(x='arg', y='value', hue='variable', data=V_df, ax=axes[3,2], 
            palette=snakemake.params.colorpalette)
    #### add true V
    trueV_ = np.array([np.diag(trueV[x]) for x in pd.unique(V_df['arg'])]).flatten()
    xs = plot_help.snsbox_get_x(len(np.unique(V_df['arg'])), len(np.unique(V_df['variable'])))
    axes[3,2].scatter(xs, trueV_, color=snakemake.params.pointcolor, zorder=10)
    axes[3,2].set_xlabel(snakemake.wildcards.arg)
    axes[3,2].set_ylabel('V_diag (cell type-specific genetic variance')
    axes[3,2].text(-0.05, 1.05, '(M)', fontsize=16, transform=axes[3,2].transAxes)
    handles, labels = axes[3,2].get_legend_handles_labels()
    axes[3,2].legend(handles=handles, labels=[r'$%s$'%(x) for x in labels])

    ### W
    #### diag
    W_df =  pd.melt(summaries['full'], id_vars=['arg'], 
            value_vars=['W_'+str(i+1) for i in range(C)])
    sns.boxplot(x='arg', y='value', hue='variable', data=W_df, ax=axes[3,3], 
            palette=snakemake.params.colorpalette)
    #### add true W
    trueW_ = np.array([np.diag(trueW[x]) for x in pd.unique(W_df['arg'])]).flatten()
    xs = plot_help.snsbox_get_x(len(np.unique(W_df['arg'])), len(np.unique(W_df['variable'])))
    axes[3,3].scatter(xs, trueW_, color=snakemake.params.pointcolor, zorder=10)
    axes[3,3].set_xlabel(snakemake.wildcards.arg)
    axes[3,3].set_ylabel('W_diag (cell type-specific noise variance')
    axes[3,3].text(-0.05, 1.05, '(O)', fontsize=16, transform=axes[3,3].transAxes)
    handles, labels = axes[3,3].get_legend_handles_labels()
    axes[3,3].legend(handles=handles, labels=[r'$%s$'%(x) for x in labels])

    #### V non-diag
    Vlow_df = pd.melt(summaries['full'], id_vars=['arg'],
            value_vars=['Vlow_'+str(i+1) for i in range(C*(C-1)//2)])
    sns.boxplot(x='arg', y='value', hue='variable', data=Vlow_df, ax=axes[3,4], 
            palette=snakemake.params.colorpalette)
    #### add true V
    trueVlow_ = np.array([trueV[x][np.tril_indices(C,k=-1)].flatten() 
        for x in pd.unique(V_df['arg'])]).flatten()
    xs = plot_help.snsbox_get_x(len(np.unique(Vlow_df['arg'])), 
            len(np.unique(Vlow_df['variable'])))
    axes[3,4].scatter(xs, trueVlow_, color=snakemake.params.pointcolor, zorder=10)
    axes[3,4].set_xlabel(snakemake.wildcards.arg)
    axes[3,4].set_ylabel('V_lowtri (cell type-specific genetic covariance')
    axes[3,4].text(-0.05, 1.05, '(N)', fontsize=16, transform=axes[3,4].transAxes)
    handles, labels = axes[3,4].get_legend_handles_labels()
    axes[3,4].legend(handles=handles, labels=labels)

    #### W non-diag
    Wlow_df = pd.melt(summaries['full'], id_vars=['arg'],
            value_vars=['Wlow_'+str(i+1) for i in range(C*(C-1)//2)])
    sns.boxplot(x='arg', y='value', hue='variable', data=Wlow_df, ax=axes[3,5], 
            palette=snakemake.params.colorpalette)
    #### add true W
    trueWlow_ = np.array([trueW[x][np.tril_indices(C,k=-1)].flatten() 
        for x in pd.unique(W_df['arg'])]).flatten()
    xs = plot_help.snsbox_get_x(len(np.unique(Wlow_df['arg'])), 
            len(np.unique(Wlow_df['variable'])))
    axes[3,5].scatter(xs, trueWlow_, color=snakemake.params.pointcolor, zorder=10)
    axes[3,5].set_xlabel(snakemake.wildcards.arg)
    axes[3,5].set_ylabel('W_lowtri (cell type-specific noise covariance')
    axes[3,5].text(-0.05, 1.05, '(P)', fontsize=16, transform=axes[3,5].transAxes)
    handles, labels = axes[3,5].get_legend_handles_labels()
    axes[3,5].legend(handles=handles, labels=labels)
    #### add dash lines
    for ax in [axes[0,0], axes[0,1], axes[1,0], axes[1,1], axes[1,2], axes[1,3], 
            axes[2,0], axes[2,1], axes[2,2], axes[2,3], axes[3,0], axes[3,1], 
            axes[3,2], axes[3,3], axes[3,4], axes[3,5]]:
        ax.axhline(0, c='0.8', ls='--', zorder=0)
    for i in range(4):
        axes[i,0].axhline(0.2, c='0.8', ls='--', zorder=0)

    #### tweak x labels
    if len(summaries['full']['arg'].values[0]) > 15:
        for ax in axes.flatten():
            ax.tick_params(axis='x', labelsize='small', labelrotation=15)
    fig.tight_layout()
    fig.savefig(snakemake.output.png)

if __name__ == '__main__':
    main()
