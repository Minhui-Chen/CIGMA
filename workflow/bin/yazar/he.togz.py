import re
import sys
import numpy as np
import pandas as pd
from cigma import util


def main():
    # Get input/output files from command line or snakemake
    if len(sys.argv) > 1:
        input_out = sys.argv[1]
        input_P = sys.argv[2]
        output_out = sys.argv[3]
    else:
        try:
            input_out = snakemake.input.out
            input_P = snakemake.input.P
            output_out = snakemake.output.out
        except NameError:
            raise ValueError("Either provide command line arguments or run via snakemake")
    
    # read
    out = np.load(input_out, allow_pickle=True).item()
    # cell types
    P = pd.read_table(input_P, index_col=0)
    cts = P.columns.tolist()
    cts = [ct.replace(' ', '') for ct in cts]
    C = len(cts)

    # collect data
    data = util.read_out(out, ['gene'])


    # free
    if 'free' in out:
        free = util.read_out(out, ['gene', 'free:hom_g2', 'free:hom_e2', 'free:v', 'free:w', 'free:V', 'free:W', 'free:ct_beta'], cts)
        free['shared_h2'] = free['free:hom_g2'] / (free['free:hom_g2'] + free['free:v'] + free['free:hom_e2'] + free['free:w'])
        free['specific_h2'] = free['free:v'] / (free['free:hom_g2'] + free['free:v'] + free['free:hom_e2'] + free['free:w'])
        free['specificity'] = free['free:v'] / (free['free:v'] + free['free:hom_g2'])

        # trans
        if 'hom_g2_b' in out['free'].keys():
            free_trans = util.read_out(out, ['gene', 'free:hom_g2_b', 'free:v_b', 'free:V_b'], cts)
            free_trans['specificity_b'] = free_trans['free:v_b'] / (free_trans['free:v_b'] + free_trans['free:hom_g2_b'])
            free = free.merge(free_trans)
            bio_var = free['free:hom_g2'] + free['free:v'] + free['free:hom_g2_b'] + free['free:v_b'] + free['free:hom_e2'] + free['free:w']
            free['shared_h2'] = free['free:hom_g2'] / bio_var
            free['specific_h2'] = free['free:v'] / bio_var
            free['shared_h2_b'] = free['free:hom_g2_b'] / bio_var
            free['specific_h2_b'] = free['free:v_b'] / bio_var
        
        # p
        if 'p' in out:
            if 'hom_g2' in out['p']['free']:
                free_p = util.read_out(out, ['gene', 'p:free:hom_g2', 'p:free:hom_e2', 'p:free:V', 'p:free:W', 'p:free:vc', 
                                             'p:free:var_hom_g2', 'p:free:var_hom_e2', 'p:free:var_V', 'p:free:var_W'], cts)
                free_p['p:free:std_hom_g2'] = np.sqrt(free_p['p:free:var_hom_g2'])
                free_p = free_p.drop('p:free:var_hom_g2', axis=1)
                free_p['p:free:std_hom_e2'] = np.sqrt(free_p['p:free:var_hom_e2'])
                free_p = free_p.drop('p:free:var_hom_e2', axis=1)
                for ct in cts:
                    free_p[f'p:free:std_V_{ct}'] = np.sqrt(free_p[f'p:free:var_V_{ct}'])
                    free_p = free_p.drop(f'p:free:var_V_{ct}', axis=1)
                    free_p[f'p:free:std_W_{ct}'] = np.sqrt(free_p[f'p:free:var_W_{ct}'])
                    free_p = free_p.drop(f'p:free:var_W_{ct}', axis=1)
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
    data[columns].to_csv(output_out, index=False, float_format='%.6g')


if __name__ == '__main__':
    main()