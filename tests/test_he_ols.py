"""
Regression tests validating the matrix-free he_ols against the dense reference
implementation (he_ols_dense) it replaces as the default.
"""
import numpy as np
import pytest

from cigma import fit
from cigma.datasets import load_sim_data


def _make_data(N, C, seed, with_Kt=False, n_random=0, p_extra_fixed=0):
    rng = np.random.default_rng(seed)

    G = rng.normal(size=(N, 25))
    G = (G - G.mean(0)) / G.std(0)
    K = G @ G.T / 25

    Kt = None
    if with_Kt:
        Gt = rng.normal(size=(N, 15))
        Gt = (Gt - Gt.mean(0)) / Gt.std(0)
        Kt = Gt @ Gt.T / 15

    Y = rng.normal(size=(N, C))
    ctnu = rng.gamma(2, 1, size=(N, C))

    X = np.zeros((N * C, C + p_extra_fixed))
    for n in range(N):
        for c in range(C):
            X[n * C + c, c] = 1.0
    for k in range(p_extra_fixed):
        col = rng.normal(size=N)
        for n in range(N):
            for c in range(C):
                X[n * C + c, C + k] = col[n]

    random_covars = {f'r{k}': rng.normal(size=(N, rng.integers(2, 4)))
                      for k in range(n_random)}

    return Y, K, X, ctnu, Kt, random_covars


@pytest.mark.parametrize('model,with_Kt', [
    ('iid', False), ('iid', True),
    ('free', False), ('free', True),
    ('full', False),
])
@pytest.mark.parametrize('n_random,random_shared', [
    (0, True), (1, True), (1, False), (2, True), (2, False),
])
def test_he_ols_matches_dense(model, with_Kt, n_random, random_shared):
    Y, K, X, ctnu, Kt, random_covars = _make_data(
        N=9, C=4, seed=hash((model, with_Kt, n_random, random_shared)) % (2**31),
        with_Kt=with_Kt, n_random=n_random)

    theta_dense = fit.he_ols_dense(Y, K, X, ctnu, model, random_covars=random_covars,
                                    Kt=Kt, random_shared=random_shared)
    theta_mf = fit.he_ols(Y, K, X, ctnu, model, random_covars=random_covars,
                           Kt=Kt, random_shared=random_shared)

    np.testing.assert_allclose(theta_mf, theta_dense, atol=1e-8, rtol=1e-6)


@pytest.mark.parametrize('model', ['iid', 'free', 'full'])
@pytest.mark.parametrize('p_extra_fixed', [0, 1, 3])
def test_he_ols_matches_dense_with_extra_fixed_covars(model, p_extra_fixed):
    Y, K, X, ctnu, Kt, random_covars = _make_data(
        N=8, C=4, seed=p_extra_fixed * 10, with_Kt=(model != 'full'),
        n_random=1, p_extra_fixed=p_extra_fixed)

    kt_arg = Kt if model != 'full' else None
    theta_dense = fit.he_ols_dense(Y, K, X, ctnu, model, random_covars=random_covars, Kt=kt_arg)
    theta_mf = fit.he_ols(Y, K, X, ctnu, model, random_covars=random_covars, Kt=kt_arg)

    np.testing.assert_allclose(theta_mf, theta_dense, atol=1e-6, rtol=1e-6)


def test_he_ols_dtype_float32():
    Y, K, X, ctnu, Kt, random_covars = _make_data(N=8, C=4, seed=999, with_Kt=True, n_random=1)

    theta_dense = fit.he_ols_dense(Y.copy(), K.copy(), X.copy(), ctnu.copy(), 'free',
                                    random_covars=random_covars, Kt=Kt, dtype='float32')
    theta_mf = fit.he_ols(Y.copy(), K.copy(), X.copy(), ctnu.copy(), 'free',
                           random_covars=random_covars, Kt=Kt, dtype='float32')

    np.testing.assert_allclose(theta_mf, theta_dense, atol=1e-3, rtol=1e-3)


def test_free_he_end_to_end_matches_dense_on_bundled_data():
    """Cross-check through the full free_HE + jackknife pipeline, not just he_ols directly."""
    data = load_sim_data()
    Y, K, ctnu, P = data['Y'], data['K'], data['ctnu'], data['P']
    sex, batch = data['sex'], data['batch']

    out_mf, p_mf = fit.free_HE(Y=Y, K=K, ctnu=ctnu, P=P, fixed_covars={'sex': sex},
                                random_covars={'batch': batch}, jk=True)

    orig_he_ols = fit.he_ols
    try:
        fit.he_ols = fit.he_ols_dense
        out_dense, p_dense = fit.free_HE(Y=Y, K=K, ctnu=ctnu, P=P, fixed_covars={'sex': sex},
                                          random_covars={'batch': batch}, jk=True)
    finally:
        fit.he_ols = orig_he_ols

    assert out_mf['hom_g2'] == pytest.approx(out_dense['hom_g2'], abs=1e-8)
    assert out_mf['hom_e2'] == pytest.approx(out_dense['hom_e2'], abs=1e-8)
    np.testing.assert_allclose(np.diag(out_mf['V']), np.diag(out_dense['V']), atol=1e-8)
    np.testing.assert_allclose(np.diag(out_mf['W']), np.diag(out_dense['W']), atol=1e-8)
    assert p_mf['hom_g2'] == pytest.approx(p_dense['hom_g2'], abs=1e-6)
