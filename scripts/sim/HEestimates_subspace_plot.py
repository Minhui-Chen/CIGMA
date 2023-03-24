import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import plot_help, subspace_plot

def main(): 
    # 
    os.makedirs( os.path.dirname(snakemake.output.png), exist_ok=True )
    subspace = snakemake.params.subspace 
    sim_plot_order = snakemake.params.sim_plot_order 
    mycolors = snakemake.params.mycolors 
    pointcolor = snakemake.params.pointcolor
    colorpalette = snakemake.params.colorpalette 

    # get cell type number
    free = np.load(snakemake.input.out[0], allow_pickle=True).item()['he']['free']
    C = free['V'][0].shape[0]

    # collect plot order
    plot_order_ = sim_plot_order[snakemake.wildcards.model][snakemake.wildcards.arg] 
    args = subspace[snakemake.wildcards.arg]
    plot_order_ = plot_order_ + list( args[~np.isin(args, plot_order_)] )
    plot_order_ = np.array( plot_order_ )
    plot_order_ = plot_order_[np.isin(plot_order_, subspace[snakemake.wildcards.arg])]

    # collect estimates from Free
    summaries = {'free':[]}
    for arg_, out_f in zip(np.array(subspace[snakemake.wildcards.arg]), snakemake.input.out):
        out = np.load(out_f, allow_pickle=True).item()
        estimates = out['he']
        # read Hom model estimates
        #hom = estimates['hom']
        #hom2 = hom['hom2']
        #summary = pd.DataFrame({'subject':hom2})
        #summary['arg'] = arg_
        #summaries['hom'].append(summary)

        # read IID model estimates
        #iid = estimates['iid']
        #hom2, V, W = iid['hom2'], iid['V'], iid['W']
        #V = np.array([x[0,0] for x in V]) # extract the first diagonal of V
        #W = np.array([x[0,0] for x in W]) # extract the first diagonal of W
        #summary = pd.DataFrame({'subject':hom2})
        #summary['V'] = V
        #summary['W'] = W
        #summary['arg'] = arg_
        #summaries['iid'].append(summary)

        # read Free model estimates
        free = estimates['free']
        hom_g2, hom_e2, V, W = free['hom_g2'], free['hom_e2'], free['V'], free['W']
        V = np.array([np.diag(x) for x in V]) # extract out diagonal of V
        W = np.array([np.diag(x) for x in W]) # extract out diagonal of W
        V, W = V.T, W.T
        summary = pd.DataFrame({'hom_g2':hom_g2, 'hom_e2':hom_e2})
        for i in range(C):
            summary['V_'+str(i+1)] = V[i]
            summary['W_'+str(i+1)] = W[i]
        summary['arg'] = arg_
        summaries['free'].append(summary)

    # concat and sort
    for model in summaries.keys():
        summaries_ = pd.concat(summaries[model], ignore_index=True)
        summaries_['arg'] = pd.Categorical(summaries_['arg'], plot_order_)
        summaries_ = summaries_.sort_values('arg').reset_index(drop=True)
        summaries_.to_csv(sys.stdout, sep='\t')
        summaries[model] = summaries_

    # collect true V
    trueV = {}
    for arg_, trueV_f in zip(np.array(subspace[snakemake.wildcards.arg]), snakemake.input.V):
        # read true V values
        trueV[arg_] = np.loadtxt(trueV_f)

    # collect true W
    trueW = {}
    for arg_, trueW_f in zip(np.array(subspace[snakemake.wildcards.arg]), snakemake.input.W):
        # read true W values
        trueW[arg_] = np.loadtxt(trueW_f)

    # collect True hom_g2 hom_e2
    hom_g2s = []
    hom_e2s = []
    subspace_ = subspace.copy()  # may don't need to copy
    subspace_[snakemake.wildcards.arg] = pd.Categorical(subspace_[snakemake.wildcards.arg], 
            plot_order_)
    subspace_ = subspace_.sort_values(snakemake.wildcards.arg).reset_index(drop=True)
    for vc in np.array(subspace_['vc']):
        #print( vc )
        hom_g2s.append( float(vc.split('_')[1]) )
        hom_e2s.append( float(vc.split('_')[2]) )

    # plot
    plt.rcParams.update( {'font.size' : 10} )
    fig, axes = plt.subplots(nrows=1, ncols=4, figsize=(16, 4), sharex=True, squeeze=False)
#    ## Hom model
#    ### variance components
#    sns.boxplot(x='arg', y='subject', data=summaries['hom'], ax=axes[0,0], 
#            color=mycolors[0])
#
#    ### add True variances to plot
#    xs = plot_help.snsbox_get_x(len(np.unique(summaries['hom']['arg'])), 1)
#    axes[0,0].scatter(xs, vcs, color=pointcolor, zorder=10)
#
#    axes[0,0].xaxis.label.set_visible(False)
#    axes[0,0].set_ylabel('$\sigma_{hom}^2$')
#    axes[0,0].text(-0.29, 0.95, 'Hom', fontsize=16, transform=axes[0,0].transAxes)
#    axes[0,0].text(-0.05, 1.05, '(A)', fontsize=12, transform=axes[0,0].transAxes)
#    #handles, labels = axes[0,0].get_legend_handles_labels()
#    #axes[0,0].legend(handles=handles, labels=labels)
#
#    axes[0,1].axis('off')
#    axes[0,2].axis('off')
#
#    ## IID model
#    ### variance components
#    sns.boxplot(x='arg', y='subject', data=summaries['iid'], ax=axes[1,0], 
#            color=mycolors[0])
#
#    ### add True variances to plot
#    xs = plot_help.snsbox_get_x(len(np.unique(summaries['iid']['arg'])), 1)
#    axes[1,0].scatter(xs, vcs, color=pointcolor, zorder=10)
#
#    axes[1,0].xaxis.label.set_visible(False)
#    axes[1,0].set_ylabel('$\sigma_{hom}^2$')
#    axes[1,0].text(-0.29, 0.95, 'IID', fontsize=16, transform=axes[1,0].transAxes)
#    axes[1,0].text(-0.05, 1.05, '(B)', fontsize=12, transform=axes[1,0].transAxes)
#    #handles, labels = axes[1,0].get_legend_handles_labels()
#    #axes[1,0].legend(handles=handles, labels=labels)
#
#    ### V
#    sns.boxplot(x='arg', y='V', data=summaries['iid'], ax=axes[1,1], 
#            color=mycolors[0])
#    #### add true V
#    trueV_ = np.array([np.diag(trueV[x]) for x in pd.unique(summaries['iid']['arg'])]).T.flatten()
#    xs = plot_help.snsbox_get_x(len(np.unique(summaries['iid']['arg'])), 1)
#    xs = list(xs) * C
#    axes[1,1].scatter(xs, trueV_, color=pointcolor, zorder=10)
#    axes[1,1].xaxis.label.set_visible(False)
#    axes[1,1].set_ylabel('V_diag (cell type-specific genetic variance)')
#    axes[1,1].text(-0.05, 1.05, '(C)', fontsize=12, transform=axes[1,1].transAxes)
#
#    ### W
#    sns.boxplot(x='arg', y='W', data=summaries['iid'], ax=axes[1,2], 
#            color=mycolors[0])
#    #### add true W
#    trueW_ = np.array([np.diag(trueW[x]) for x in pd.unique(summaries['iid']['arg'])]).T.flatten()
#    xs = plot_help.snsbox_get_x(len(np.unique(summaries['iid']['arg'])), 1)
#    xs = list(xs) * C
#    axes[1,2].scatter(xs, trueW_, color=pointcolor, zorder=10)
#    axes[1,2].xaxis.label.set_visible(False)
#    axes[1,2].set_ylabel('W_diag (cell type-specific noise variance)')
#    axes[1,2].text(-0.05, 1.05, '(C)', fontsize=12, transform=axes[1,2].transAxes)

    ## Free model
    ### hom_g2
    sns.boxplot(x='arg', y='hom_g2', data=summaries['free'], ax=axes[0,0], 
            color=mycolors[0])

    #### add True variances
    xs = plot_help.snsbox_get_x(len(np.unique(summaries['free']['arg'])), 1)
    axes[0,0].scatter(xs, hom_g2s, color=pointcolor, zorder=10)
    #axes[0,0].xaxis.label.set_visible(False)
    axes[0,0].set_ylabel('$\sigma_{g}^2$', fontsize=16)
    #axes[0,0].text(-0.29, 0.95, 'Free', fontsize=16, transform=axes[0,0].transAxes)
    #axes[0,0].text(-0.05, 1.05, '(D)', fontsize=12, transform=axes[2,0].transAxes)

    ### hom_e2
    sns.boxplot(x='arg', y='hom_e2', data=summaries['free'], ax=axes[0,1], 
            color=mycolors[0])

    #### add True variances
    xs = plot_help.snsbox_get_x(len(np.unique(summaries['free']['arg'])), 1)
    axes[0,1].scatter(xs, hom_e2s, color=pointcolor, zorder=10)
    #axes[0,1].xaxis.label.set_visible(False)
    axes[0,1].set_ylabel('$\sigma_{e}^2$', fontsize=16)

    ### V
    V_df =  pd.melt(summaries['free'], id_vars=['arg'], value_vars=['V_'+str(i+1) for i in range(C)])
    sns.boxplot(x='arg', y='value', hue='variable', data=V_df, ax=axes[0,2], 
            palette=colorpalette)
    #### add true V
    trueV_ = np.array([np.diag(trueV[x]) for x in pd.unique(V_df['arg'])]).flatten()
    xs = plot_help.snsbox_get_x(len(np.unique(V_df['arg'])), len(np.unique(V_df['variable'])))
    axes[0,2].scatter(xs, trueV_, color=pointcolor, zorder=10)
    #axes[0,2].xaxis.label.set_visible(False)
    axes[0,2].set_ylabel('V', fontsize=16)
    #axes[0,2].text(-0.05, 1.05, '(E)', fontsize=12, transform=axes[2,1].transAxes)
    handles, labels = axes[0,2].get_legend_handles_labels()
    axes[0,2].legend(handles=handles, labels=[r'$%s$'%(x) for x in labels])

    ### W
    W_df =  pd.melt(summaries['free'], id_vars=['arg'], value_vars=['W_'+str(i+1) for i in range(C)])
    sns.boxplot(x='arg', y='value', hue='variable', data=W_df, ax=axes[0,3], 
            palette=colorpalette)
    #### add true W
    trueW_ = np.array([np.diag(trueW[x]) for x in pd.unique(W_df['arg'])]).flatten()
    xs = plot_help.snsbox_get_x(len(np.unique(W_df['arg'])), len(np.unique(W_df['variable'])))
    axes[0,3].scatter(xs, trueW_, color=pointcolor, zorder=10)
    #axes[0,3].xaxis.label.set_visible(False)
    axes[0,3].set_ylabel('W', fontsize=16)
    #axes[0,3].text(-0.05, 1.05, '(E)', fontsize=12, transform=axes[0,3].transAxes)
    handles, labels = axes[0,3].get_legend_handles_labels()
    axes[0,3].legend(handles=handles, labels=[r'$%s$'%(x) for x in labels])

    #### x axis label
    xlabel = subspace_plot.get_xlabel( snakemake.wildcards.arg )
    for ax in axes.flatten():
        if len(xlabel) < 30:
            ax.set_xlabel(xlabel, fontsize=16)
        else:
            ax.set_xlabel(xlabel, fontsize=12)

    #### tweak x labels
    subspace_plot.set_xtickname( fig, ax, snakemake.wildcards.arg )

    fig.tight_layout()
    fig.savefig(snakemake.output.png)

if __name__ == '__main__':
    main()
