import numpy as np
import pandas as pd
from scipy.stats import chi2
from cigma import wald


def inv_variance_meta(Ps, outs):
    '''Inverse variance meta (not used), weights are too noisy'''
    # collect number of individuals
    ns = np.array([pd.read_table(P_f, index_col=0).shape[0] for P_f in Ps])
    N = ns.sum()

    # find common genes
    out = np.load(outs[0], allow_pickle=True).item()
    C = out['free']['V'].shape[-1]
    genes = out['gene']
    for out_f in outs[1:]:
        out = np.load(out_f, allow_pickle=True).item()
        genes = genes[np.isin(genes, out['gene'])]  # intersection of gene names

    # make new output
    new_out  = {'gene': genes, 
                'free': {'r2':{}}, 
                'p': {'free': {}}
                }

    w_hom_g2 = np.zeros(len(genes))
    w_hom_e2 = np.zeros(len(genes))
    w_V = np.zeros((len(genes), C, C))
    w_W = np.zeros((len(genes), C, C))
    new_out['free']['hom_g2'] = np.zeros(len(genes))
    new_out['free']['hom_e2'] = np.zeros(len(genes))
    new_out['free']['V'] = np.zeros((len(genes), C, C))
    new_out['free']['W'] = np.zeros((len(genes), C, C))
    for i, out_f in enumerate(outs):
        out = np.load(out_f, allow_pickle=True).item()
        out_genes = out['gene']
        idx = np.isin(out_genes, genes)
        assert np.array_equal(genes, out_genes[idx])
        # hom g2
        hom_g2 = out['free']['hom_g2'][idx]
        w_hom_g2_i = 1 / out['p']['free']['var_hom_g2'][idx]
        w_hom_g2 += w_hom_g2_i
        new_out['free']['hom_g2'] += w_hom_g2_i * hom_g2

        # hom e2
        hom_e2 = out['free']['hom_e2'][idx]
        w_hom_e2_i = 1 / out['p']['free']['var_hom_e2'][idx]
        w_hom_e2 += w_hom_e2_i
        new_out['free']['hom_e2'] += w_hom_e2_i * hom_e2


        # V W
        V = np.diagonal(out['free']['V'], axis1=1, axis2=2)[idx]
        W = np.diagonal(out['free']['W'], axis1=1, axis2=2)[idx]
        w_V_i = np.array([np.linalg.inv(X) for X in out['p']['free']['var_V'][idx]])
        w_V += w_V_i
        new_out['free']['V'] += np.array([np.diag(w_V_i[k] @ V[k]) for k in range(w_V_i.shape[0])])
        w_W_i = np.array([np.linalg.inv(X) for X in out['p']['free']['var_W'][idx]])
        w_W += w_W_i
        new_out['free']['W'] += np.array([np.diag(w_W_i[k] @ W[k]) for k in range(w_W_i.shape[0])])


    # hom g2
    new_out['free']['hom_g2'] /= w_hom_g2
    new_out['p']['free']['var_hom_g2'] = 1 / w_hom_g2

    # hom e2
    new_out['free']['hom_e2'] /= w_hom_g2
    new_out['p']['free']['var_hom_e2'] = 1 / w_hom_g2

    # V W
    new_out['p']['free']['var_V'] = np.array([np.linalg.inv(X) for X in w_V])
    new_out['p']['free']['var_W'] = np.array([np.linalg.inv(X) for X in w_W])
    new_out['free']['V'] = [np.diag(new_out['p']['free']['var_V'][k] @ np.diag(new_out['free']['V'][k])) 
                            for k in range(len(genes))]
    new_out['free']['W'] = [np.diag(new_out['p']['free']['var_W'][k] @ np.diag(new_out['free']['W'][k])) 
                            for k in range(len(genes))]

    # p
    hom_g2 = new_out['free']['hom_g2']
    var_hom_g2 = new_out['p']['free']['var_hom_g2']
    p_hom_g2 = [wald.wald_test(hom_g2[i], 0, var_hom_g2[i], N - 3 * C - 2)
                for i in range(len(hom_g2))]
    new_out['p']['free']['hom_g2'] = p_hom_g2

    hom_e2 = new_out['free']['hom_e2']
    var_hom_e2 = new_out['p']['free']['var_hom_e2']
    p_hom_e2 = [wald.wald_test(hom_e2[i], 0, var_hom_e2[i], N - 3 * C - 2)
                for i in range(len(hom_e2))]
    new_out['p']['free']['hom_e2'] = p_hom_e2

    V = new_out['free']['V']
    var_V = new_out['p']['free']['var_V']
    p_V = [wald.mvwald_test(np.diag(V[i]), np.zeros(C), var_V[i], n=N, P=3 * C + 2) 
           for i in range(len(V))] # NOTE: sample size N
    new_out['p']['free']['V'] = p_V

    W = new_out['free']['W']
    var_W = new_out['p']['free']['var_W']
    p_W = [wald.mvwald_test(np.diag(W[i]), np.zeros(C), var_W[i], n=N, P=3 * C + 2) 
           for i in range(len(W))]
    new_out['p']['free']['W'] = p_W


    return new_out


def sample_size_weighted_meta(Ps, outs):
    # collect number of individuals
    ns = np.array([pd.read_table(P_f, index_col=0).shape[0] for P_f in Ps])
    ps = ns / ns.sum()

    # find common genes
    out = np.load(outs[0], allow_pickle=True).item()
    genes = out['gene']
    for out_f in outs[1:]:
        out = np.load(out_f, allow_pickle=True).item()
        genes = genes[np.isin(genes, out['gene'])]  # intersection of gene names

    # make new output
    new_out  = {'gene': genes, 
                'free': {'r2':{}}, 
                'p': {'free': {}}
                }

    # meta
    for i, out_f in enumerate(outs):
        p = ps[i]
        out = np.load(out_f, allow_pickle=True).item()
        out_genes = out['gene']
        idx = np.isin(out_genes, genes)
        # check gene order
        assert np.array_equal(genes, out_genes[idx])

        for key, val in out['free'].items():
            if key in ['nu', 'ctnu']:
                continue
            if isinstance(val, np.ndarray):
                if key not in new_out['free']:
                    new_out['free'][key] = p * val[idx]
                else:
                    new_out['free'][key] += p * val[idx]
        if 'r2' in out['free'].keys():
            for key, val in out['free']['r2'].items():
                if isinstance(val, np.ndarray):
                    if key not in new_out['free']['r2']:
                        new_out['free']['r2'][key] = p * val[idx]
                    else:
                        new_out['free']['r2'][key] += p * val[idx]
    

        # add specificity 
        V_bar = np.mean(np.diagonal(out['free']['V'], axis1=1, axis2=2), axis=1)
        hom_g2 = out['free']['hom_g2']
        specificity = V_bar / (hom_g2 + V_bar)
        if 'specificity' not in new_out['free'].keys():
            new_out['free']['specificity'] = p * specificity[idx]
        else:
            new_out['free']['specificity'] += p * specificity[idx]
    
    return new_out


def main():
    # read in the data
    P1 = pd.read_table(snakemake.input.P[0], index_col=0)
    cts = P1.columns.values

    # sanity check
    for P_f in snakemake.input.P[1:]:
        assert np.array_equal(cts, pd.read_table(P_f, index_col=0).columns.values)
        
    # find common genes
    out = np.load(snakemake.input.out[0], allow_pickle=True).item()
    genes = out['gene']
    for out_f in snakemake.input.out[1:]:
        out = np.load(out_f, allow_pickle=True).item()
        genes = genes[np.isin(genes, out['gene'])]  # intersection of gene names

    # inverse variance meta
    new_out = inv_variance_meta(snakemake.input.P, snakemake.input.out)

    # sample size meta: h2 specificity
    tmp_out = sample_size_weighted_meta(snakemake.input.P, snakemake.input.out)
    new_out['free']['shared_h2'] = tmp_out['free']['shared_h2']
    new_out['free']['specific_h2'] = tmp_out['free']['specific_h2']
    new_out['free']['specificity'] = tmp_out['free']['specificity']


    # save
    np.save(snakemake.output.out, new_out)
                

if __name__ == '__main__':
    main()