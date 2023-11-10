import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from gxctmm import log

def cor(V):
    std = np.sqrt(np.diag(V))
    return( V/np.outer(std,std) )


def main():
    # 
    P = pd.read_table(snakemake.input.P, index_col=0)
    cts = P.columns
    out = np.load(snakemake.input.out, allow_pickle=True).item()

    #
    Vs = out['full']['V']
    Ws = out['full']['W']

    # calculate heritability
    sig_g2s = np.array([np.diag(V).tolist() for V in Vs])
    sig_e2s = np.array([np.diag(W).tolist() for W in Ws])
    sig2s = sig_g2s + sig_e2s
    h2s = sig_g2s / sig2s 

    # calculate correlation 
    cor_Vs, cor_Ws = [], []
    for V in Vs:
        if np.any(np.diag(V) <= 0):
            continue
        else:
            cor_Vs.append( cor(V) )
    for W in Ws:
        if np.any(np.diag(W) <= 0):
            continue
        else:
            cor_Ws.append( cor(W) )
    cor_Vs = np.array( cor_Vs )
    cor_Ws = np.array( cor_Ws )
    print( cor_Vs.shape )

    # get median
    h2_median = np.median(h2s, axis=0)
    e2_median = 1 - h2_median 
    #cor_V_median = np.median(cor_Vs, axis=0)
    #print( cor_V_median )
    #cor_W_median = np.median(cor_Ws, axis=0)
    V_median = np.median(Vs, axis=0)
    W_median = np.median(Ws, axis=0)

    # plot for covariance
    fig, axes = plt.subplots(ncols=2, figsize=(20,9))

    np.fill_diagonal(V_median, h2_median)
    np.fill_diagonal(W_median, e2_median)

    # Create a boolean mask to hide the lower triangle
    mask = np.zeros_like(V_median, dtype=bool)
    mask[np.triu_indices_from(mask, k=0)] = True

    sns.heatmap(data=V_median, mask=mask, annot=False, cmap='Reds', xticklabels=cts, yticklabels=cts, ax=axes[0])
    axes[0].set_title('Genetic covariance between CTs', fontsize=16)
    sns.heatmap(data=W_median, mask=mask, annot=False, cmap='Reds', xticklabels=cts, yticklabels=cts, ax=axes[1])
    axes[1].set_title('Environmental covariance between CTs', fontsize=16)

    # 
    fig.savefig( snakemake.output.cov )

    # plot for h2
    fig, ax = plt.subplots(figsize=(10,4))

    h2s = np.clip(h2s, a_min=None, a_max=0.05)
    sns.violinplot(data=pd.DataFrame(h2s, columns=cts))
    ax.axhline(y=0, color='0.9', ls='--', zorder=0)
    ax.set_title('CT-specific heritability')

    # 
    fig.savefig( snakemake.output.h2 )

if __name__ == '__main__':
    main()

