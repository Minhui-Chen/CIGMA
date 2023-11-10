import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import plot_help, subspace_plot

def main():
    # par
    os.makedirs( os.path.dirname(snakemake.output.png), exist_ok=True )
    subspace = snakemake.params.subspace
    wildcards = snakemake.wildcards 
    sim_plot_order = snakemake.params.sim_plot_order 
    mycolors = snakemake.params.mycolors 
    if snakemake.wildcards.model == 'full':
        os.system(f'touch {snakemake.output.png}')
        sys.exit()

    # select and order data
    data = pd.DataFrame({'arg':np.array(subspace[wildcards.arg]), 'out':snakemake.input.out})

    plot_order_ = sim_plot_order[wildcards.model][wildcards.arg]
    args = np.array( data['arg'] )
    plot_order_ = plot_order_ + list( args[~np.isin(args, plot_order_)] )
    plot_order_ = np.array( plot_order_ )
    plot_order_ = plot_order_[np.isin(plot_order_, data['arg'])]

    data['arg'] = pd.Categorical(data['arg'], plot_order_)
    data = data.sort_values('arg').reset_index(drop=True)
    data.to_csv(sys.stdout, sep='\t')

    he_power = {}
    for f in data['out']:
        out = np.load(f, allow_pickle=True).item()
        # wald
        he = out['wald'] # structure: he_p - model (e.g. hom_p) - statistics (e.g. hom2, V)
        for m in ['free']:
            he_m = he[m]
            if m not in he_power.keys():
                he_power[m] = {}
            for key, value in he_m.items():
                power_ = np.sum(value < 0.05, axis=0) / value.shape[0]
                if key not in he_power[m].keys():
                    he_power[m][key] = [power_]
                else:
                    he_power[m][key].append(power_)

    for m in ['free']:
        for key in he_power[m].keys():
            he_power[m][key] = np.array(he_power[m][key])
            if he_power[m][key].ndim > 1:
                he_power[m][key] = np.array(he_power[m][key]).T

    # plot
    mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=mycolors)
    plt.rcParams.update( {'font.size' : 12} )
    markers = plot_help.mymarkers()
    fig, ax = plt.subplots(ncols=1, sharex=True, sharey=True, figsize=(6, 4))
    # wald
#    ## hom
#    axes[0].plot(data['arg'], he_power['hom']['hom2'], marker=markers[0], label='$\sigma_{hom}^2$')
#    axes[0].set_title('Hom')
#    axes[0].set_xlabel(wildcards.arg)
#    axes[0].set_ylabel('Positive rate')
#    axes[0].legend()
#    ## iid
#    axes[1].plot(data['arg'], he_power['iid']['hom2'], marker=markers[0], label='$\sigma_{hom}^2$')
#    axes[1].plot(data['arg'], he_power['iid']['V'], marker=markers[0], label='$\sigma_{het}^2$')
#    axes[1].plot(data['arg'], he_power['iid']['W'], marker=markers[0], label='$\sigma_{w}^2$')
#    axes[1].set_title('IID')
#    axes[1].set_xlabel(wildcards.arg)
#    axes[1].legend()
    ## free
    #ax.plot(data['arg'], he_power['free']['sigma_g2'], marker=markers[0], label='$\sigma_{g}^2$')
    #ax.plot(data['arg'], he_power['free']['sigma_e2'], marker=markers[0], label='$\sigma_{e}^2$')
    ax.plot(data['arg'], he_power['free']['V'], marker=markers[0], label='$V$')
    #ax.plot(data['arg'], he_power['free']['W'], marker=markers[0], label='$W$')
    #ax.plot(data['arg'], he_power['free']['VW'], marker=markers[0], label='$(V, W)$')
    #for i in range(len(he_power['free']['Vi'])):
    #    axes[2].plot(data['arg'], he_power['free']['Vi'][i], marker='.', label=f'$V_{i+1}$', ls='--',
    #            color=mycolors[3+i], alpha=0.5)
    #for i in range(len(he_power['free']['Wi'])):
    #    axes[2].plot(data['arg'], he_power['free']['Wi'][i], marker='.', label=f'$W_{i+1}$', ls=':',
    #            color=mycolors[3+i], alpha=0.5)

    #axes[2].set_title('Free')
    if wildcards.model == 'free':
        ax.set_ylabel('True positive rate', fontsize=16)
    else:
        ax.set_ylabel('False positive rate', fontsize=16)
    #ax.legend()

    ax.set_ylim((0-0.02,1+0.02))
    ax.axhline(y=0.05, color='0.6', ls='--', zorder=0)

    #### tweak x labels
    subspace_plot.set_xtickname(fig, ax, wildcards.arg)
    xlabel = subspace_plot.get_xlabel(wildcards.arg)
    ax.set_xlabel( xlabel, fontsize=16 )
    fig.tight_layout()
    fig.savefig(snakemake.output.png)

if __name__ == '__main__':
    main()
