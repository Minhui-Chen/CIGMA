import sys, re, os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import plot_help, subspace_plot

def main():
    os.makedirs( os.path.dirname(snakemake.output.png), exist_ok=True )
    # select and order data
    data = pd.DataFrame({'arg':np.array(snakemake.params.subspace[snakemake.wildcards.arg]),
        'out':snakemake.input.out})

    plot_order_ = snakemake.params.sim_plot_order[snakemake.wildcards.model][snakemake.wildcards.arg]
    args = np.array( data['arg'] )
    plot_order_ = plot_order_ + list( args[~np.isin(args, plot_order_)] )
    plot_order_ = np.array( plot_order_ )
    plot_order_ = plot_order_[np.isin(plot_order_, data['arg'])]

    data['arg'] = pd.Categorical(data['arg'], plot_order_)
    data = data.sort_values('arg').reset_index(drop=True)
    data.to_csv(sys.stdout, sep='\t')

    wald_power = {}
    lrt_power = {}
    for f in data['out']:
        out = np.load(f, allow_pickle=True).item()
        # wald
#        for m in ['free']:
#            wald_m = out['reml'][m]['wald']
#            if m not in wald_power.keys():
#                wald_power[m] = {}
#            for key, value in wald_m.items():
#                power_ = np.sum(value < 0.05, axis=0) / value.shape[0]
#                if key not in wald_power[m].keys():
#                    wald_power[m][key] = [power_]
#                else:
#                    wald_power[m][key].append(power_)
        # lrt
        if 'lrt' in out['reml'].keys():
            lrt = out['reml']['lrt'] # structure: lrt - two comparing model (e.g. full_free)
            for key, value in lrt.items():
                power = np.sum(value < 0.05) / value.shape[0]
                if key not in lrt_power.keys():
                    lrt_power[key] = [power]
                else:
                    lrt_power[key].append(power)

#    for m in ['free']:
#        for key in wald_power[m].keys():
#            wald_power[m][key] = np.array(wald_power[m][key])
#            if wald_power[m][key].ndim > 1:
#                wald_power[m][key] = np.array(wald_power[m][key]).T

    # plot
    mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=snakemake.params.mycolors)
    plt.rcParams.update( {'font.size' : 12} )
    markers = plot_help.mymarkers()
    fig, ax = plt.subplots(ncols=1, sharex=True, sharey=True, figsize=(6, 4))
    # wald
#    ## hom
#    axes[0].plot(data['arg'], wald_power['hom']['hom2'], marker=markers[0], label='$\sigma_{hom}^2$')
#    axes[0].set_title('Hom')
#    axes[0].set_xlabel(snakemake.wildcards.arg)
#    axes[0].set_ylabel('Positive rate')
#    axes[0].legend()
#    ## iid
#    axes[1].plot(data['arg'], wald_power['iid']['hom2'], marker=markers[0], label='$\sigma_{hom}^2$')
#    axes[1].plot(data['arg'], wald_power['iid']['V'], marker=markers[0], label='$\sigma_{het}^2$')
#    axes[1].plot(data['arg'], wald_power['iid']['W'], marker=markers[0], label='$\sigma_{noise}^2$')
#    axes[1].set_title('IID')
#    axes[1].set_xlabel(snakemake.wildcards.arg)
#    axes[1].legend()
    ## free
    #axes[0].plot(data['arg'], wald_power['free']['sigma_g2'], marker=markers[0], label='$\sigma_{g}^2$')
    #axes[0].plot(data['arg'], wald_power['free']['sigma_e2'], marker=markers[0], label='$\sigma_{e}^2$')
    #axes[0].plot(data['arg'], wald_power['free']['V'], marker=markers[0], label='$V$')
    #axes[0].plot(data['arg'], wald_power['free']['W'], marker=markers[0], label='$W$')
    #for i in range(len(wald_power['free']['Vi'])):
    #    axes[2].plot(data['arg'], wald_power['free']['Vi'][i], marker='.', label=f'$V_{i+1}$', ls='--',
    #            color=snakemake.params.mycolors[3+i], alpha=0.5)
    #for i in range(len(wald_power['free']['Vi'])):
    #    axes[2].plot(data['arg'], wald_power['free']['Wi'][i], marker='.', label=f'$W_{i+1}$', ls=':',
    #            color=snakemake.params.mycolors[3+i], alpha=0.5)
    #axes[0].set_title('Free')
    #axes[0].set_xlabel(snakemake.wildcards.arg)
    #axes[0].legend()

    # lrt
    def format_name(x):
        x = re.sub('full','Full',re.sub('free','Free',re.sub('iid','IID',re.sub('hom','Hom',re.sub('null','Null',x)))))
        return(' vs '.join(x.split('_')[:2]))

    a, b, c, d = 0, 0, 0, 0
    for x in lrt_power.keys():
        if re.search('^hom', x):
            ax.plot(data['arg'], lrt_power[x], marker=markers[a], label=format_name(x), 
                    color=snakemake.params.mycolors[0])
            a += 1
        if re.search('^iid', x):
            ax.plot(data['arg'], lrt_power[x], marker=markers[b], label=format_name(x), 
                    color=snakemake.params.mycolors[1])
            b += 1
        if re.search('^free', x):
            print( data )
            print( lrt_power )
            print( x )
            ax.plot(data['arg'], lrt_power[x], marker=markers[c], label=format_name(x), 
                    color=snakemake.params.mycolors[2])
            c += 1
        if re.search('^full', x):
            ax.plot(data['arg'], lrt_power[x], marker=markers[d], label=format_name(x), 
                    color=snakemake.params.mycolors[3])
            d += 1
    ax.legend()
    ax.set_title('LRT')
    if snakemake.wildcards.model == 'free':
        ax.set_ylabel('True positive rate', fontsize=16)
    elif snakemake.wildcards.model == 'hom':
        ax.set_ylabel('False positive rate', fontsize=16)

    ax.set_ylim((0-0.02,1+0.02))
    ax.axhline(y=0.05, color='0.6', ls='--', zorder=0)

    #### tweak x labels
    subspace_plot.set_xtickname( fig, ax, snakemake.wildcards.arg )
    xlabel = subspace_plot.get_xlabel( snakemake.wildcards.arg )
    ax.set_xlabel( xlabel, fontsize=16 )
    fig.tight_layout()
    fig.savefig(snakemake.output.png)

if __name__ == '__main__':
    main()
