import numpy as np
from cigma import log, fit

def main():

    outs = []
    data = np.load(snakemake.input.data, allow_pickle=True).item()
    iid_jk = snakemake.params.get('iid_jk', False)
    free_jk = snakemake.params.get('free_jk', False)

    for key in data.keys():
        log.logger.info(f'{key}')

        # read
        Y = data[key]['Y']
        K = data[key]['K']
        ctnu = data[key]['ctnu']
        # set ctnu to missing
        ctnu = np.zeros_like(ctnu)
        
        P = data[key]['P']
        if 'fixed' in data[key].keys():
            fixed = {'fixed': data[key]['fixed'][:, :-1]}
        else:
            fixed = {}
        if 'random' in data[key].keys():
            random = {'random': data[key]['random']}
        else:
            random = {}

        if snakemake.wildcards.model != 'full':
            ## IID
            iid_he, iid_he_wald = fit.iid_HE(Y, K, ctnu, P, fixed_covars=fixed, 
                                random_covars=random, jk=iid_jk)

            ## Free
            free_he, free_he_wald = fit.free_HE(Y, K, ctnu, P, fixed_covars=fixed, 
                                                random_covars=random, jk=free_jk)

            # Full
            # full_he = fit.full_HE(Y, K, ctnu, P, fixed_covars=fixed, random_covars=random)

            # save
            outs.append({'gene': key, 'iid':iid_he, 'free':free_he, 
                         's':data[key]['s'], 
                         'pi': data[key]['pi'], 'nu': data[key]['nu'].mean(), 
                         'var_y': data[key]['y'].var(), 
                         'wald': {'iid': iid_he_wald, 'free':free_he_wald}})
                        # 'full':full_he, 
        else:
            # Full
            full_he = fit.full_HE(Y, K, ctnu, P, fixed_covars=fixed, random_covars=random)

            # save
            outs.append({'gene': key, 'full':full_he, 'nu': data[key]['nu'].mean(), 
                         's':data[key]['s'], 'pi': data[key]['pi']})
        
    np.save(snakemake.output.out, outs)


if __name__ == '__main__':
    main()
