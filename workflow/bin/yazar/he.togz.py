import re
import numpy as np
import pandas as pd
from cigma import util


def main():
    # read
    out = np.load(snakemake.input.out, allow_pickle=True).item()
    # cell types
    P = pd.read_table(snakemake.input.P, index_col=0)
    cts = P.columns.tolist()
    cts = [ct.replace(' ', '') for ct in cts]
    C = len(cts)

    # collect data
    data = util.read_out(out, ['gene'])


    # free
    if 'free' in out:
        free = util.read_out(out, ['gene', 'free:hom_g2', 'free:hom_e2', 'free:v', 'free:w', 'free:V', 'free:W', 'free:ct_beta'], cts)

        # trans
        if 'hom_g2_b' in out['free'].keys():
            free_trans = util.read_out(out, ['gene', 'free:hom_g2_b', 'free:v_b', 'free:V_b'], cts)
            free = free.merge(free_trans)
        
        # p
        if 'p' in out:
            if 'hom_g2' in out['p']['free']:
                free_p = util.read_out(out, ['gene', 'p:free:hom_g2', 'p:free:hom_e2', 'p:free:V', 'p:free:W', 'p:free:vc', 'p:free:var_specificity'], cts)
                free = free.merge(free_p)
        
        data = data.merge(free)

    # full
    if 'full' in out:
        full = util.read_out(out, ['gene', 'full:v', 'full:w', 'full:V', 'full:W'], cts)
        data = data.merge(full)


    # order 
    columns = ['gene']
    columns += [x for x in data.columns if re.search('^iid', x)]
    columns += [x for x in data.columns if re.search('^free', x)]
    columns += [x for x in data.columns if re.search('^full', x)]
    columns += [x for x in data.columns if x not in columns]

    # save
    data[columns].to_csv(snakemake.output.out, index=False, float_format='%.6g')


if __name__ == '__main__':
    main()