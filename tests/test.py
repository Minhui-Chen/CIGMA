import numpy as np, pandas as pd
from cigma import log, fit


def main():

    dataset = np.load('tests/sim.npy', allow_pickle=True).item()  # load the dataset
    data = dataset[0]  # get the first dataset
    
    Y = data['Y']  # Cell type pseudobulk expression matrix: Individuals X Cell types
    K = data['K']  # genomic relationship matrix: Individuals X Individuals
    ctnu = data['ctnu']  # cell-to-cell variance matrix \delta: Individuals X Cell types
    P = data['P']  # cell type proportion matrix: Individuals X Cell types

    # run CIGMA
    out, _ = fit.free_HE(Y=Y, K=K, ctnu=ctnu, P=P)
    print(out)

    # run CIGMA with Jackknife
    # out, p = fit.free_HE(Y=Y, K=K, ctnu=ctnu, P=P, jk=True)


if __name__ == '__main__':
    main()
