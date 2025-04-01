from os import error
import sys, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import norm


def read_ldsc(data, trait, ct, ldsc_results_f):
    # get h2 and M
    log_f = re.sub('results$', 'log', ldsc_results_f)
    snp_no = 0
    h2 = 0
    for line in open(log_f):
        if 'SNPs remain' in line:
            line = line.strip().split()
            snp_no = line[-3]
            if '(' in snp_no:
                snp_no = snp_no[1:]
            snp_no = int(snp_no)
        if 'h2:' in line:
            h2 = float(line.strip().split()[-2])
    
    # print(h2, snp_no)


    df = pd.read_table(ldsc_results_f, index_col=0)
    data['trait'].append(trait)
    data['ct'].append(ct)
    data['gene set'].append('cell type-specific mean')
    data['enrichment'].append(df.loc['L2_0', 'Enrichment'])
    data['enrichment p'].append(df.loc['L2_0', 'Enrichment_p'])
    data['enrichment se'].append(df.loc['L2_0', 'Enrichment_std_error'])
    data['coefficient'].append(df.loc['L2_0', 'Coefficient'] * snp_no / h2)
    data['-log(coefficient p)'].append((-1) * np.log10(norm.sf(df.loc['L2_0', 'Coefficient_z-score'])))
    data['coefficient se'].append(df.loc['L2_0', 'Coefficient_std_error'] * snp_no / h2)

    data['trait'].append(trait)
    data['ct'].append(ct)
    data['gene set'].append('cell type-specific variance')
    data['enrichment'].append(df.loc['L2_1', 'Enrichment'])
    data['enrichment p'].append(df.loc['L2_1', 'Enrichment_p'])
    data['enrichment se'].append(df.loc['L2_1', 'Enrichment_std_error'])
    data['coefficient'].append(df.loc['L2_1', 'Coefficient'] * snp_no / h2)
    data['-log(coefficient p)'].append((-1) * np.log10(norm.sf(df.loc['L2_1', 'Coefficient_z-score'])))
    data['coefficient se'].append(df.loc['L2_1', 'Coefficient_std_error'] * snp_no / h2)

    data['trait'].append(trait)
    data['ct'].append(ct)
    data['gene set'].append('random (negative control)')
    data['enrichment'].append(df.loc['L2_2', 'Enrichment'])
    data['enrichment p'].append(df.loc['L2_2', 'Enrichment_p'])
    data['enrichment se'].append(df.loc['L2_2', 'Enrichment_std_error'])
    data['coefficient'].append(df.loc['L2_2', 'Coefficient'] * snp_no / h2)
    data['-log(coefficient p)'].append((-1) * np.log10(norm.sf(df.loc['L2_2', 'Coefficient_z-score'])))
    data['coefficient se'].append(df.loc['L2_2', 'Coefficient_std_error'] * snp_no / h2)



def main():
    # read
    P = pd.read_table(snakemake.input.P, index_col=0)
    cts = P.columns.tolist()

    data = {'ct': [], 'gene set': [], 'trait': [], 'enrichment': [], 
            'enrichment p': [], 'enrichment se': [],
            'coefficient': [], '-log(coefficient p)': [], 'coefficient se': []}
    

    for ldsc_f, trait in zip(snakemake.input.ldsc, snakemake.params.traits):
        for ct, line in zip(cts, open(ldsc_f)):
            read_ldsc(data, trait, ct, line.strip())

    # 
    data = pd.DataFrame(data)
    data.to_csv('test.txt', index=False)
    print(data.head())

    # plot
    # shift p values for mean to negative
    x = -data.loc[data['gene set'] == 'cell type-specific mean', '-log(coefficient p)']
    data.loc[data['gene set'] == 'cell type-specific mean', '-log(coefficient p)'] = x

    trait_order = snakemake.params.neg_traits + snakemake.params.pos_traits
    gene_set_order = ['random (negative control)', 'cell type-specific mean', 'cell type-specific variance']
    fig, ax = plt.subplots(figsize=(10, 8))

    sns.pointplot(x='-log(coefficient p)', y='trait', hue='ct', data=data.loc[data['gene set'] != 'random (negative control)'], 
                  order=trait_order, dodge=True, linestyles='none', errorbar=None, ax=ax,)

    sns.pointplot(x='-log(coefficient p)', y='trait', color='0.8', data=data.loc[data['gene set'] == 'random (negative control)'], 
                  order=trait_order, dodge=True, linestyles='none', errorbar=None, ax=ax,)
    
    ax.axvline(x=0, color='black', linestyle='--')
    ax.axvline(x=-2, color='0.7', linestyle='--', zorder=0)
    ax.axvline(x=2, color='0.7', linestyle='--', zorder=0)
    ax.set_xlabel('-log(coefficient p)')
    ax.set_ylabel('trait')
    ax.text(0.25, 1.02, 'Mean', ha='center', va='center', transform=ax.transAxes)
    ax.text(0.75, 1.02, 'Variance', ha='center', va='center', transform=ax.transAxes)
    xlim_max = np.max(np.abs(ax.get_xlim()))
    ax.set_xlim(-xlim_max, xlim_max)

    # Get current x-axis tick positions
    current_tick_positions, _ = plt.xticks()

    # Convert current tick positions to absolute values
    absolute_tick_positions = [abs(tick_position) for tick_position in current_tick_positions]

    # Set x-axis tick positions to absolute values
    plt.xticks(current_tick_positions, absolute_tick_positions)
    

    fig.tight_layout()
    fig.savefig(snakemake.output.png)


if __name__ == '__main__':
     main()