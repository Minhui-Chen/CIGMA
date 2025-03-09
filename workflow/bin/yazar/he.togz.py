import re, sys
import numpy as np
import pandas as pd


def main():
    # read
    out = np.load(snakemake.input.out, allow_pickle=True).item()

    # collect data
    data = pd.DataFrame({'gene': out['gene']})

    # cell types
    P = pd.read_table(snakemake.input.P, index_col=0)
    cts = P.columns.tolist()
    cts = [ct.replace(' ', '') for ct in cts]
    C = len(cts)

    # iid
    if 'iid' in out.keys():
        data['iid_cis_hom_g2'] = out['iid']['cis_hom_g2']
        data['iid_hom_e2'] = out['iid']['hom_e2']
        data['iid_cis_V'] = out['iid']['cis_V']
        data['iid_W'] = out['iid']['W']
        batch = out['iid']['r2']['batch']
        if len(batch.shape) > 1:
            data[[f'iid_batch_{ct}' for ct in cts]] = batch
        else:
            data['iid_batch'] = batch
        # data[[f'iid_beta_{ct}' for ct in cts]] = out['iid']['ct_beta']

        # trans
        if 'trans_hom_g2' in out['iid'].keys():
            data['iid_trans_hom_g2'] = out['iid']['trans_hom_g2']
            data['iid_trans_V'] = out['iid']['trans_V']

    # free
    if 'free' in out.keys():
        data['free_cis_hom_g2'] = out['free']['hom_g2']
        data['free_hom_e2'] = out['free']['hom_e2']
        # batch = out['free']['r2']['batch']
        # if len(batch.shape) > 1:
        #     data[[f'free_batch_{ct}' for ct in cts]] = batch
        # else:
        #     data['free_batch'] = batch
        # data['overall_nu'] = out['free']['nu'].mean(axis=1)
        # data[[f'ctnu_{ct}' for ct in cts]] = out['free']['ctnu']
        cis_V = np.diagonal(out['free']['V'], axis1=1, axis2=2)
        W = np.diagonal(out['free']['W'], axis1=1, axis2=2)
        data[[f'free_cis_V_{ct}' for ct in cts]] = cis_V
        data[[f'free_W_{ct}' for ct in cts]] = W
        # data[[f'free_h2_{ct}' for ct in cts]] = out['free']['cis_h2']
        data[[f'free_beta_{ct}' for ct in cts]] = out['free']['ct_beta']

        # trans
        if 'trans_hom_g2' in out['free'].keys():
            data['free_trans_hom_g2'] = out['free']['trans_hom_g2']
            trans_V = np.diagonal(out['free']['trans_V'], axis1=1, axis2=2)
            data[[f'free_trans_V_{ct}' for ct in cts]] = trans_V
            # data[[f'free_trans_h2_{ct}' for ct in cts]] = free['free']['trans_h2']
        
        # p
        if 'p' in out.keys():
            # se
            if 'var_hom_g2' in out['p']['free'].keys():
                data['free_se_cis_hom_g2'] = np.sqrt(out['p']['free']['var_hom_g2'])
            if 'var_hom_e2' in out['p']['free'].keys():
                data['free_se_hom_e2'] = np.sqrt(out['p']['free']['var_hom_e2'])
            if 'var_V' in out['p']['free'].keys():
                var_cis_V = np.diagonal(out['p']['free']['var_V'], axis1=1, axis2=2)
                data[[f'free_se_cis_V_{ct}' for ct in cts]] = np.sqrt(var_cis_V)
            if 'var_W' in out['p']['free'].keys():
                var_W = np.diagonal(out['p']['free']['var_W'], axis1=1, axis2=2)
                data[[f'free_se_W_{ct}' for ct in cts]] = np.sqrt(var_W)
            # p
            data[f'free_p_cis_hom_g2'] = out['p']['free']['hom_g2']
            data[f'free_p_hom_e2'] = out['p']['free']['hom_e2']
            data[f'free_p_cis_V'] = out['p']['free']['V']
            data[f'free_p_W'] = out['p']['free']['W']

    # full
    if 'full' in out.keys():
        cis_V = np.diagonal(out['full']['cis_V'], axis1=1, axis2=2)
        W = np.diagonal(out['full']['W'], axis1=1, axis2=2)
        data[[f'full_cis_V_{ct}' for ct in cts]] = cis_V
        data[[f'full_W_{ct}' for ct in cts]] = W
        # data[[f'full_beta_{ct}' for ct in cts]] = out['full']['ct_beta']
        for i in range(C-1):
            for j in range(i+1, C):
                data[f'full_V_{cts[i]}+{cts[j]}'] = out['full']['cis_V'][:, i, j]
                data[f'full_W_{cts[i]}+{cts[j]}'] = out['full']['W'][:, i, j]


    # order 
    columns = ['gene']
    columns += [x for x in data.columns if re.search('^iid', x)]
    columns += [x for x in data.columns if re.search('^free', x)]
    columns += [x for x in data.columns if re.search('^full', x)]
    columns += [x for x in data.columns if x not in columns]

    # save
    data[columns].to_csv(snakemake.output.out, index=False)


if __name__ == '__main__':
    main()