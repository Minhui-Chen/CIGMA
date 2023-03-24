import re, os, sys
import scipy
import numpy as np
import pandas as pd
import rpy2.robjects as robjects
import helper
sys.path.insert(0, 'bin')
import screml, wald 

#def hom_ML(y_f, G_f, Z_f, P_f, nu_f):
#    tmpfn = helper.generate_tmpfn()
#    os.system(f'Rscript bin/og_scml.R {y_f} {Z_f} {P_f} {nu_f} hom {tmpfn}')
#    robjects.r.load(f'{tmpfn}')
#    hom2, beta, l, hess, Vy = np.array(robjects.r('out$hom2'))[0], np.array(robjects.r('out$beta')), robjects.r('out$l')[0], np.array(robjects.r('out$hess')), np.array(robjects.r('out$sig2s'))
#    ## wald test
#    G = np.loadtxt(G_f)
#    P = np.loadtxt(P_f)
#    C = P.shape[1]
#
#    Z_designm = [G]  # hom2 is for individual, G matrix is for SNP
#    D = wald.asymptotic_dispersion_matrix(P, Z_designm, Vy)
#    ## sort analytic D matrix to change the order of beta and hom2 variance
#    neworder = [C]+list(range(C))+list(range(C+1,D.shape[0]))
#    D = D[neworder][:,neworder]
#    hom = {'hom2':hom2, 'beta':beta, 'l':l, 'D':D}
#    hom_p = {}
#    hom_p['hom2'] = wald.wald_test(hom2/G.shape[1], 0, D[0,0])  # hom2/L
#    hom_p['beta'] = [wald.wald_test(beta[i], 0, D[i+1,i+1]) for i in range(C)]
#    return(hom, hom_p)

def HE(y, Z, P, X, vs, model):
    '''
    model: hom, free, full
    '''
    N, C = P.shape
    # project out cell type main effect
    proj = np.eye(N) - X @ np.linalg.inv(X.T @ X) @ X.T
    y_p = proj @ y
    # y' y'^T - M diag(\nu) M
    Y_m = np.outer(y_p, y_p) - proj @ np.diag(vs) @ proj
    Y = Y_m.flatten('F')

    # T
    T = []
    if model != 'full':
        T.append( (proj @ Z @ proj).flatten('F') ) # sigma_g2
        T.append( proj.flatten('F') )
    if model in ['free','full']:
        for c in range(C):
            T.append( (proj @ (np.outer(P[:,c],P[:,c]) * Z) @ proj).flatten('F') ) # V diagonal
        for c in range(C):
            T.append( (proj @ np.diag(P[:,c] * P[:,c]) @ proj).flatten('F') ) # W diagonal
        if model == 'full':
            for i in range(C-1):
                for j in range(i+1,C):
                    T.append( 2 * (proj @ (np.outer(P[:,i],P[:,j]) * Z) @ proj).flatten('F') ) # V off-diagonal
            for i in range(C-1):
                for j in range(i+1,C):
                    T.append( 2 * (proj @ np.diag(P[:,i] * P[:,j]) @ proj).flatten('F') ) # W off-diagonal
    T = np.array(T).T
    
    # theta: sigma_g2, sigma_e2 in Hom; sigma_g2, sigma_e2, V diag, W diag in Free;
    # V diag, W diag, V off diag, and W off diag in Full
    theta = np.linalg.inv(T.T @ T) @ T.T @ Y

    return(theta)

def hom_HE(y_f, Z_f, P_f, nu_f):
    y = np.loadtxt(y_f) 
    Z = np.loadtxt(Z_f)
    P = np.loadtxt(P_f)
    vs = np.loadtxt(nu_f)
    N, C = P.shape # cell type number
    X = P
    n_par = 2

    sigma_g2, sigma_e2 = HE(y, Z, P, X, vs, 'hom')
    he = {'sigma_g2': sigma_g2, 'sigma_e2':sigma_e2}

    # jackknife
    #jacks = [HE(np.delete(y,i,axis=0), np.delete(np.delete(Z,i,axis=0),i,axis=1), 
    #    np.delete(P,i,axis=0), np.delete(X,i,axis=0), np.delete(vs,i,axis=0), 'hom') for i in range(N)]
    #sigma_g2_var = (len(jacks) - 1.0) * np.var([x[0] for x in jacks])
    #sigma_e2_var = (len(jacks) - 1.0) * np.var([x[1] for x in jacks])

    #p = {'sigma_g2': wald.wald_test(sigma_g2, 0, sigma_g2_var, N-n_par),
    #        'sigma_e2': wald.wald_test(sigma_e2, 0, sigma_e2_var, N-n_par)
    #        }

    return( he )

#def iid_ML(y_f, G_f, Z_f, P_f, nu_f):
#    tmpfn = helper.generate_tmpfn()
#    os.system(f'Rscript bin/og_scml.R {y_f} {Z_f} {P_f} {nu_f} iid {tmpfn}')
#    robjects.r.load(f'{tmpfn}')
#    hom2, beta, V, W, l, hess, Vy = np.array(robjects.r('out$hom2'))[0], np.array(robjects.r('out$beta')), np.array(robjects.r('out$V')), np.array(robjects.r('out$W')), robjects.r('out$l')[0], np.array(robjects.r('out$hess')), np.array(robjects.r('out$sig2s'))
#
#    ## wald test
#    G = np.loadtxt(G_f)
#    P = np.loadtxt(P_f)
#    N = P.shape[0]
#    C = P.shape[1]
#
#    Z_designm = [G, scipy.linalg.khatri_rao(G.T, P.T).T,
#        scipy.linalg.khatri_rao(np.identity(N), P.T).T]
#    D = wald.asymptotic_dispersion_matrix(P, Z_designm, Vy)
#    ## sort analytic D matrix to change the order of beta and hom2 variance
#    neworder = [C]+list(range(C))+list(range(C+1,D.shape[0]))
#    D = D[neworder][:,neworder]
#    iid = {'hom2':hom2, 'beta':beta, 'V':V, 'W':W, 'l':l, 'D':D}
#
#    iid_p = {}
#    iid_p['hom2'] = wald.wald_test(hom2/G.shape[1], 0, D[0,0])
#    iid_p['beta'] = [wald.wald_test(beta[i], 0, D[i+1,i+1]) for i in range(C)]
#    iid_p['V'] = wald.wald_test(V[0,0]/G.shape[1], 0, D[C+1,C+1])
#    iid_p['W'] = wald.wald_test(W[0,0], 0, D[C+2,C+2])
#    return(iid, iid_p)
#
#def iid_HE(y_f, Z_f, P_f, nu_f):
#    def iid_HE_(y, Z, P, vs):
#        # project out cell type main effect
#        proj = np.eye(len(y)) - P @ np.linalg.inv(P.T @ P) @ P.T
#        y_p = proj @ y
#        # momments
#        Y_m = np.outer(y_p, y_p) - proj @ np.diag(vs) @ proj
#        Y = Y_m.flatten('F')
#
#        # HE
#        x_m_1 = proj @ Z @ proj
#        x_m_2 = proj @ ((P @ P.T) * Z) @ proj
#        x_m_3 = proj @ ((P @ P.T) * np.eye(len(y))) @ proj
#        x_m = np.array([x_m_1.flatten('F'), x_m_2.flatten('F'), x_m_3.flatten('F')]).T
#        he_iid_theta = np.linalg.lstsq(x_m, Y, rcond=None)[0]
#        return(he_iid_theta)
#
#    y = np.loadtxt(y_f) 
#    Z = np.loadtxt(Z_f)
#    P = np.loadtxt(P_f)
#    vs = np.loadtxt(nu_f)
#    C = P.shape[1] # cell type number
#
#    he_iid_theta = iid_HE_(y, Z, P, vs)
#    he_iid_jacks = [iid_HE_(np.delete(y,i,axis=0), np.delete(np.delete(Z,i,axis=0),i,axis=1), 
#        np.delete(P,i,axis=0), np.delete(vs,i,axis=0)) for i in range(len(y))]
#    he_iid_var = (len(he_iid_jacks) - 1.0) * np.cov(np.array(he_iid_jacks).T, bias=True)
#    iid_he = {'hom2': he_iid_theta[0], 'V': np.eye(C) * he_iid_theta[1],
#            'W': np.eye(C) * he_iid_theta[2]}
#    iid_he_p = {'hom2': wald.wald_test(he_iid_theta[0], 0, he_iid_var[0,0]),
#            'V': wald.wald_test(he_iid_theta[1], 0, he_iid_var[1,1]),
#            'W': wald.wald_test(he_iid_theta[2], 0, he_iid_var[2,2]),
#            'VW': wald.mvwald_test(he_iid_theta[1:3], np.zeros(2), he_iid_var[1:3,1:3])}
#    return(iid_he, iid_he_p)

#def free_ML(y_f, G_f, Z_f, P_f, nu_f):
#    tmpfn = helper.generate_tmpfn()
#    os.system(f'Rscript bin/og_scml.R {y_f} {Z_f} {P_f} {nu_f} free {tmpfn}')
#    robjects.r.load(f'{tmpfn}')
#    hom2, beta, V, W, l, hess, Vy = np.array(robjects.r('out$hom2'))[0], np.array(robjects.r('out$beta')), np.array(robjects.r('out$V')), np.array(robjects.r('out$W')), robjects.r('out$l')[0], np.array(robjects.r('out$hess')), np.array(robjects.r('out$sig2s'))
#    ## wald test
#    G = np.loadtxt(G_f)
#    P = np.loadtxt(P_f)
#    C = P.shape[1]
#
#    Z_designm = [G]
#    Z_designm = Z_designm + [G * P[:,i].reshape((-1,1)) for i in range(C)]
#    Z_designm = Z_designm + [np.diag(P[:,i]) for i in range(C)]
#    D = wald.asymptotic_dispersion_matrix(P, Z_designm, Vy)
#    ## sort analytic D matrix to change the order of beta and hom2 variance
#    neworder = [C]+list(range(C))+list(range(C+1,D.shape[0]))
#    D = D[neworder][:,neworder]
#    free = {'hom2':hom2, 'beta':beta, 'V':V, 'W':W, 'l':l, 'D':D}
#    free_p = {}
#    free_p['hom2'] = wald.wald_test(hom2/G.shape[1], 0, D[0,0])
#    free_p['beta'] = [wald.wald_test(beta[i], 0, D[i+1,i+1]) for i in range(C)]
#    ### mvwald test V and W
#    free_p['V'] = wald.mvwald_test(np.diag(V)/G.shape[1], np.zeros(C), D[(C+1):(2*C+1), (C+1):(2*C+1)])
#    free_p['Vi'] = [wald.wald_test(V[i,i]/G.shape[1], 0, D[C+i+1,C+i+1]) for i in range(C)]
#    free_p['W'] = wald.mvwald_test(np.diag(W), np.zeros(C), D[(2*C+1):(3*C+1), (2*C+1):(3*C+1)])
#    free_p['Wi'] = [wald.wald_test(W[i,i], 0, D[2*C+i+1,2*C+i+1]) for i in range(C)]
#    return(free, free_p)

def free_HE(y_f, Z_f, P_f, nu_f):
    y = np.loadtxt(y_f) 
    Z = np.loadtxt(Z_f)
    P = np.loadtxt(P_f)
    vs = np.loadtxt(nu_f)
    N, C = P.shape # cell type number
    X = P
    n_par = 2 + 2 * C

    theta = HE(y, Z, P, X, vs, 'free')
    sigma_g2, sigma_e2 = theta[0], theta[1]
    V, W = np.diag(theta[2:(C+2)]), np.diag(theta[(C+2):(C*2+2)])
    jacks = [HE(np.delete(y,i,axis=0), np.delete(np.delete(Z,i,axis=0),i,axis=1),
        np.delete(P,i,axis=0), np.delete(X,i,axis=0), np.delete(vs,i,axis=0), 'free') for i in range(N)]
    covar = (len(jacks)-1.0) * np.cov(np.array(jacks).T, bias=True)
    sigma_g2_var = covar[0,0]
    sigma_e2_var = covar[1,1]
    V_var = covar[2:(C+2), 2:(C+2)]
    W_var = covar[(C+2):(C*2+2), (C+2):(C*2+2)]

    he = {'sigma_g2': sigma_g2, 'sigma_e2':sigma_e2, 'V': V, 'W': W}
    p = {'sigma_g2': wald.wald_test(sigma_g2, 0, sigma_g2_var, N-n_par),
            'sigma_e2': wald.wald_test(sigma_e2, 0, sigma_e2_var, N-n_par),
            'V': wald.mvwald_test(np.diag(V), np.zeros(C), V_var, n=N, P=n_par),
            'W': wald.mvwald_test(np.diag(W), np.zeros(C), W_var, n=N, P=n_par),
            }
    return(he, p)

#def full_ML(y_f, G_f, Z_f, P_f, nu_f):
#    tmpfn = helper.generate_tmpfn()
#    os.system(f'Rscript bin/og_scml.R {y_f} {Z_f} {P_f} {nu_f} full {tmpfn}')
#    robjects.r.load(f'{tmpfn}')
#    hom2, beta, V, W, l = np.array(robjects.r('out$hom2'))[0], np.array(robjects.r('out$beta')), np.array(robjects.r('out$V')), np.array(robjects.r('out$W')), robjects.r('out$l')[0]
#    full = {'hom2': hom2, 'beta':beta, 'V':V, 'W':W, 'l':l}
#    return(full)

def full_HE(y_f, Z_f, P_f, nu_f):
    y = np.loadtxt(y_f) 
    Z = np.loadtxt(Z_f)
    P = np.loadtxt(P_f)
    vs = np.loadtxt(nu_f)
    N, C = P.shape # cell type number
    ntril = (C-1) * C // 2
    X = P

    theta = HE(y, Z, P, X, vs, 'full')
    V, W = np.diag(theta[:C]), np.diag(theta[C:(C*2)])
    V[np.triu_indices(C,k=1)] = theta[(C*2):(C*2 + ntril)]
    V = V + V.T - np.diag(theta[:C])
    W[np.triu_indices(C,k=1)] = theta[(C*2 + ntril):(C*2 + ntril*2)]
    W = W + W.T - np.diag(theta[C:(C*2)])

    he = {'V': V, 'W': W}
    return( he )

def main():

    batch = snakemake.params.batches
    outs = [re.sub('/rep/', f'/rep{i}/', snakemake.params.out) for i in batch]
    for y_f, G_f, Z_f, P_f, nu_f, out_f in zip(
            [line.strip() for line in open(snakemake.input.y)],
            [line.strip() for line in open(snakemake.input.G)], 
            [line.strip() for line in open(snakemake.input.Z)],
            [line.strip() for line in open(snakemake.input.P)], 
            [line.strip() for line in open(snakemake.input.nu)], 
            outs):
        print(y_f, Z_f, P_f, nu_f)
        os.makedirs(os.path.dirname(out_f), exist_ok=True)
        # read
#        y = np.loadtxt(y_f) 
#        G = np.loadtxt(G_f) 
#        Z = np.loadtxt(Z_f)
#        P = np.loadtxt(P_f)
#        vs = np.loadtxt(nu_f)
#        C = P.shape[1] # cell type number
#
#        # project out cell type main effect from y
#        proj = np.eye(len(y)) - P @ np.linalg.inv(P.T @ P) @ P.T
#        y_p = proj @ y
#        # Y_ for HE estiamtes
#        Y_m = np.outer(y_p, y_p) - proj @ np.diag(vs) @ proj
#        Y_ = Y_m.flatten('F')

        # hom model
        ## ML
        #hom_ml, hom_ml_wald = hom_ML(y_f, G_f, Z_f, P_f, nu_f)

        ## HE
        #hom_he, hom_he_wald = hom_HE(y_f, Z_f, P_f, nu_f)

        # Free model
        ## ML
        #free_ml, free_ml_wald = free_ML(y_f, G_f, Z_f, P_f, nu_f)

        ## HE
        free_he, free_he_wald = free_HE(y_f, Z_f, P_f, nu_f)

        # Full model
        #full_ml = full_ML(y_f, G_f, Z_f, P_f, nu_f)

        ## HE
        full_he = full_HE(y_f, Z_f, P_f, nu_f)

        # save
        np.save(out_f,
                {'free':free_he, 'full':full_he, 
                'wald': {'free':free_he_wald} 
                }
            )
        
    with open(snakemake.output.out, 'w') as f:
        f.write('\n'.join(outs))  


if __name__ == '__main__':
    main()
