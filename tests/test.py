import numpy as np, pandas as pd
from cigma import fit
from cigma.datasets import load_sim_data

def main():

    data = load_sim_data()  # load the sample dataset
    
    Y = data['Y']  # Cell type pseudobulk expression matrix: Individuals X Cell types
    K = data['K']  # genomic relationship matrix: Individuals X Individuals
    ctnu = data['ctnu']  # cell-to-cell variance matrix \delta: Individuals X Cell types
    P = data['P']  # cell type proportion matrix: Individuals X Cell types
    sex = data['sex'] # sex (fixed covariate): Individuals X 1
    batch = data['batch'] # batch (random covariate): Individuals X batch_num

    # run CIGMA
    # out, _ = fit.free_HE(Y=Y, K=K, ctnu=ctnu, P=P, fixed_covars={'sex': sex}, random_covars={'batch': batch})  # complete in a few seconds
    out, pvalues = fit.free_HE(Y=Y, K=K, ctnu=ctnu, P=P, fixed_covars={'sex': sex}, random_covars={'batch': batch}, jk=True)  # jackknife to compute p values for cell types specific genetic effects, complete in a few seconds
    print(out)
    # output: a dictionary of estimates with keys 'hom_g2' (shared genetic variance), 'hom_e2' (shared environmental variance), 
    # 'V' (cell type-specific genetic variance matrix), 'W' (cell type-specific environmental variance matrix), 
    # 'shared_h2' (cell type-shared heritability), 'specific_h2' (cell type-specific heritability), 
    print(pvalues)
    # pvalues: a dictionary of p values for testing 'hom_g2', 'hom_e2', 'V' (null hypothesis: V[c,c] == 0 for all cell types),
    # 'vc' for each cell type (null hypothesis: V[c,c] == 0 for cell type c), and 'W' (null hypothesis: W[c,c] == 0 for all cell types)

if __name__ == '__main__':
    main()
