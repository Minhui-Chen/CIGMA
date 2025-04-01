import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    # read
    datas = []
    for stacked_f, trait in zip(snakemake.input.stacked, snakemake.params.traits):
        data = pd.read_table(stacked_f)
        data['trait'] = trait
        datas.append(data)
    
    data = pd.concat(datas, ignore_index=True)


    # plot
    trait_order = snakemake.params.neg_traits + snakemake.params.pos_traits
    gene_set_order = ['mean', 'var']
    fig, ax = plt.subplots(figsize=(10, 5))
    
    x = np.arange(data['trait'].nunique())
    width = 0.25

    multiplier = 0
    for gene_set in gene_set_order:
        offset = width * multiplier
        set_data = data.loc[data['Name'] == gene_set]
        if set_data.shape[0] == 0:
            continue
        set_data = set_data.set_index('trait')
        # ys = np.log10(set_data['Coefficient_P_value'][trait_order]) * (-1)
        # rects = ax.bar(x + offset, ys, width, yerr=None, label=gene_set)
        print(set_data['Coefficient'])
        print(trait_order)
        ys = set_data['Coefficient'][trait_order]
        yerrs = set_data['Coefficient_std_error'][trait_order]
        rects = ax.bar(x + offset, ys, width, yerr=yerrs, label=gene_set)
        labels = []
        for p in set_data['Coefficient_P_value'][trait_order]:
            if p > 0.05:
                labels.append('')
            elif p > 0.01:
                if gene_set == 'mean':
                    labels.append(f'{p:.3f}\n')
                else:
                    labels.append(f'{p:.3f}')
            else:
                labels.append(f'{p:.1e}')
        
        ax.bar_label(rects, labels, padding=3, fontsize=8)
        multiplier += 1
    # ax.set_ylabel(r'$-log_{10}(p)$', fontsize=22)
    ax.set_ylabel('Coefficient (tau*)', fontsize=22)
    ax.set_xticks(x + width, trait_order)
    ax.legend(loc='upper left')
    ylims = ax.get_ylim()
    ax.set_ylim(ylims[0], ylims[1] + 0.02 * (ylims[1] - ylims[0]))
    ax.axvline(len(snakemake.params.neg_traits) - 0.5 + width, linestyle=':', color='0.7')
    # ax.axhline(-np.log10(0.05), linestyle='-', color='0.3', zorder=0)
    ax.axhline(0, linestyle='-', color='0.3', zorder=0)

    fig.tight_layout()
    fig.savefig(snakemake.output.png)


if __name__ == '__main__':
     main()