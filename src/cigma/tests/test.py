"""
Quick test script for CIGMA.

Usage:
    python -m cigma.tests.test
"""
import numpy as np
from cigma import fit
from cigma.datasets import load_sim_data


def main():
    # Load sample dataset
    data = load_sim_data()

    Y = data['Y']        # Cell type pseudobulk: individuals x cell types
    K = data['K']         # Kinship matrix: individuals x individuals
    ctnu = data['ctnu']   # Cell-to-cell variance: individuals x cell types
    P = data['P']         # Cell type proportions: individuals x cell types
    sex = data['sex']     # Fixed covariate (sex): individuals x 1
    batch = data['batch'] # Random covariate (batch): individuals x n_batches

    N, C = Y.shape
    print(f"Data loaded: {N} individuals, {C} cell types")
    print(f"  Y shape: {Y.shape}")
    print(f"  K shape: {K.shape}")
    print(f"  ctnu shape: {ctnu.shape}")
    print(f"  P shape: {P.shape}")

    # Fit CIGMA-HE with jackknife for p-values
    print("\nFitting CIGMA (Free HE with jackknife)...")
    out, pvalues = fit.free_HE(
        Y=Y, K=K, ctnu=ctnu, P=P,
        fixed_covars={'sex': sex},  # Include sex as a fixed covariate
        random_covars={'batch': batch}, # Include batch as a random covariate
        jk=True,
    )  # Finish in a few seconds

    # out, _ = fit.free_HE(
    #     Y=Y, K=K, ctnu=ctnu, P=P,
    #     fixed_covars={'sex': sex},
    #     random_covars={'batch': batch},
    # )  # If p values are not required, this is faster.

    # Print estimates
    print("\n--- Estimates ---")
    print(f"  Shared genetic variance (hom_g2):  {out['hom_g2']:.4f}")
    print(f"  Shared environ variance (hom_e2):  {out['hom_e2']:.4f}")
    print(f"  CT-specific genetic variance (V diag): {np.diag(out['V'])}")
    print(f"  CT-specific environ variance (W diag): {np.diag(out['W'])}")

    # Print p-values
    print("\n--- P-values ---")
    print("Testing if shared genetic variance (hom_g2) > 0:")
    print(f"  p-value: {pvalues['hom_g2']:.4e}")
    print("Testing if shared environ variance (hom_e2) > 0:")
    print(f"  p-value: {pvalues['hom_e2']:.4e}")
    print("Testing if CT-specific genetic variance (V) > 0:")
    print(f"  p-value: {pvalues['V']:.4e}")
    for i in range(C):
        print(f"Testing if specific genetic variance for cell type {i+1} (V[{i+1},{i+1}]) > 0:")
        print(f"  p-value: {pvalues['vc'][i]:.4e}")
    print("Testing if CT-specific environ variance (W) > 0:")
    print(f"  p-value: {pvalues['W']:.4e}")


    print("\nTest passed!")


if __name__ == '__main__':
    main()
