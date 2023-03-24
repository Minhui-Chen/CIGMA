library(numDeriv)
library(emdbook)
library(matrixcalc)
library(reticulate) # to use python
use_python('/apps/software/gcc-6.2.0/python/3.6.0/bin/python3')
source_python('wald.py')

# variance matrix of y
SIGMA <- function(Z, P, vs, sigma_g2=0, sigma_e2=0, V=0, W=0) {
    sigma_y2 <- sigma_g2 * Z + sigma_e2 + diag(vs)

    # add cell type-specific genetic effect
    if ( all(V == 0) ) {
        # do nothing
    } else {
        sigma_y2 <- sigma_y2 + (P %*% V %*% t(P)) * Z
    }

    # add ct-specific noise
    if ( all(W == 0) ) {
        # do nothing
    } else {
        sigma_y2 <- sigma_y2 + diag( diag(P %*% W %*% t(P)) )
    }

    return( sigma_y2 )
}
    
# loglikelihood function
LL <- function(y, Z, P, X, C, vs, sigma_g2, sigma_e2, V, W){
    sigma_y2 <- SIGMA( Z, P, vs, sigma_g2, sigma_e2, V, W )

    # inverse variance
    eval <- eigen(sigma_y2, symmetric=TRUE)$values
    if (max(eval)/(min(eval)+1e-99) > 1e8 | min(eval)<0) return(1e12)

    sigma_y2_inv <- solve(sigma_y2)

    # calculate B matrix
    m1 <- t(X) %*% sigma_y2_inv
    m2 <- m1 %*% X
    eval <- eigen(m2, symmetric=TRUE)$values
    if (max(eval)/(min(eval)+1e-99) > 1e8 | min(eval)<0) return(1e12)

    B <- sigma_y2_inv - t(m1) %*% solve(m2) %*% m1

    # calculate loglikelihood
    L <- determinant(sigma_y2, logarithm=TRUE)$modulus 
    L <- L + determinant(m2, logarithm=TRUE)$modulus
    L <- 0.5 * ( L + t(y) %*% B %*% y )
    return( L )
}

screml_hom <- function(
y, Z, P, vs
){
	C      <- ncol(P)  # cell type number
    X <- P

	sigma_g2 <- var(y) / 2
    sigma_e2 <- sigma_g2
    par  <- c( sigma_g2, sigma_e2 )

	out <- optim( par=par, fn=screml_hom_loglike, 
		y=y, Z=Z, P=P, X=X, C=C, vs=vs, method = "BFGS", hessian = TRUE)

#	hom2_ <- out$par[1]
#
#    if (hom2_ > 1 ) {
#        overwhelm_variance = TRUE
#    } else {
#        overwhelm_variance = FALSE
#    }
#
#    if (out$convergence != 0 | overwhelm_variance) {
#        for (i in 1:5) {
#            par_ <- par * rgamma(length(par), 2, scale=1/2) 
#            out_ <- optim( par=par_, fn=screml_hom_loglike, 
#                y=y, Z=Z, P=P, C=C, vs=vs, method = "BFGS", hessian = TRUE)
#
#            if (out_$convergence == 0 & (out_$value * (-1)) > (out$value * (-1)) ) {
#                out <- out_
#            }
#        }
#    }

	sigma_g2 <- out$par[1]
    #p_sigma_g2 <- wald_test(sigma_g2, 0, solve(out$hessian)[1,1], nrow(P)-2) 
    sigma_e2 <- out$par[2]
    #p_sigma_e2 <- wald_test(sigma_e2, 0, solve(out$hessian)[2,2], nrow(P)-2) 
	#sigma_y2 <- SIGMA( Z, P, vs, sigma_g2=sigma_g2, sigma_e2=sigma_e2) 
    l <- out$value * (-1)
    #stopifnot(l > -1e12)
	#print( c(sigma_g2, sigma_e2, l) )
	#print(sigma_e2)
    #print(l)

    # estimate hessian matrix
    #hess = hessian(screml_hom_loglike, x=out$par, y=y, Z=Z, P=P, C=C, vs=vs)

    #try(solve(out$hessian))
    return ( list( sigma_g2=sigma_g2, sigma_e2=sigma_e2, l=l ) )
}

screml_hom_loglike <- function(par, y, Z, P, X, C, vs){
	sigma_g2 <- par[1]
	sigma_e2 <- par[2]
    V <- matrix( rep(0,C*C), nrow=C )
    W <- matrix( rep(0,C*C), nrow=C )

    l <- LL(y, Z, P, X, C, vs, sigma_g2, sigma_e2, V, W)
    return ( l  )
}

#screml_iid <- function(
#y, Z, P, vs
#){
#	C    <- ncol(P)  # cell type number
#    pi  <- colMeans(P)
#    pd  <- scale(P, scale=False)
#    S   <- (t(pd) %*% pd ) / nrow(P)
#
#	hom2 <- var(y) / 3
#	V    <- var(y) / 3 
#    W    <- var(y) / 3
#    par  <- c( hom2, V, W )
#
#	out <- optim( par=par, fn=screml_iid_loglike, 
#		y=y, Z=Z, P=P, C=C, vs=vs, method = "BFGS", hessian = TRUE)
#
#	hom2_ <- out$par[1]
#    V_    <- out$par[2] * diag(C)
#    W_    <- out$par[3] * diag(C)
#
#    if (hom2_ > 1 | sum(diag(V_ %*% S)) + pi %*% V_ %*% pi > 1 | 
#        sum(diag(W_ %*% S)) + pi %*% W_ %*% pi > 1 ) {
#        overwhelm_variance = TRUE
#    } else {
#        overwhelm_variance = FALSE
#    }
#
#    if ( out$convergence != 0 | overwhelm_variance ) {
#        for (i in 1:5) {
#            par_ <- par * rgamma(length(par), 2, scale=1/2)
#            out_ <- optim( par=par_, fn=screml_iid_loglike, 
#                y=y, Z=Z, P=P, C=C, vs=vs, method = "BFGS", hessian = TRUE)
#
#            if (out_$convergence == 0 & (out_$value * (-1)) > (out$value * (-1)) ) {
#                out <- out_
#            }
#        }
#    }
#
#	hom2 <- out$par[1]
#    V    <- out$par[2]
#    W    <- out$par[3]
#    sig2s <- Z * hom2 + V * hadamard.prod( P %*% t(P), Z ) + W * diag(diag( P %*% t(P) )) + diag(vs)
#	V <- diag(C) * V
#    W <- diag(C) * W
#    l <- out$value * (-1)
#    stopifnot(l > -1e12)
#	#print(hom2)
#    #print(V)
#    #print(W)
#    
#    # estimate hessian
#    hess = hessian(screml_iid_loglike, x=out$par, y=y, Z=Z, P=P, C=C, vs=vs)
#
#    return ( list(hom2=hom2, V=V, W=W, l=l, hess=hess, sig2s=sig2s) )
#}
#
#screml_iid_loglike<- function(par, y, Z, P, C, vs){
#	hom2 <- par[1]
#	V <- par[2]
#    W <- par[3]
#
#	sig2s <- Z * hom2 + V * hadamard.prod( P %*% t(P), Z ) + W * diag(diag( P %*% t(P) )) + diag(vs)
#    eval   <- eigen(sig2s,symmetric=TRUE)$values
#    if (max(eval)/(min(eval)+1e-99) > 1e8 | min(eval)<0) return(1e12)
#
#	if( any( diag(sig2s) < 0 ) ) return(1e12)
#
#    sig2s_inv <- solve(sig2s)
#    A <- t(P) %*% sig2s_inv %*% P
#    M <- sig2s_inv - sig2s_inv %*% P %*% solve(A) %*% t(P) %*% sig2s_inv 
#    det_V <- determinant(sig2s, logarithm=TRUE)
#    det_A <- determinant(A, logarithm=TRUE)
#    stopifnot( (det_V$sign == 1) & (det_A$sign == 1))
#    L <- -0.5 * ( det_V$modulus  + det_A$modulus + t(y) %*% M %*% y )
#    return ( L * (-1) )
#}

screml_free <- function(
y, Z, P, vs
){
    N <- nrow(P)
	C <- ncol(P)  # cell type number
    pi <- colMeans(P)
    pd <- scale(P, scale=F)
    S <- (t(pd) %*% pd ) / nrow(P)
    X <- P

	sigma_g2 <- var(y) / 4
    sigma_e2 <- sigma_g2
	V <- rep(1,C) * as.vector( sigma_g2 )
	W <- rep(1,C) * as.vector( sigma_g2 )
    par <- c( sigma_g2, sigma_e2, V, W )

	out <- optim( par=par, fn=screml_free_loglike, 
		y=y, Z=Z, P=P, X=X, C=C, vs=vs, method = "BFGS", hessian = TRUE)

#	hom2_ <- out$par[1]
#	V_ <- diag(out$par[1+1:C])
#    W_ <- diag(out$par[(C+1)+1:C])
#
#    if (hom2_ > 1 | sum(diag(V_ %*% S)) + pi %*% V_ %*% pi > 1 | 
#        sum(diag(W_ %*% S)) + pi %*% W_ %*% pi > 1 ) {
#        overwhelm_variance = TRUE
#    } else {
#        overwhelm_variance = FALSE
#    }
#
#    if ( out$convergence != 0 | overwhelm_variance ) {
#        for (i in 1:5) {
#            par_ <- par * rgamma(length(par), 2, scale=1/2)
#            out_ <- optim( par=par_, fn=screml_free_loglike, 
#                y=y, Z=Z, P=P, C=C, vs=vs, method = "BFGS", hessian = TRUE)
#
#            if (out_$convergence == 0 & (out_$value * (-1)) > (out$value * (-1)) ) {
#                out <- out_
#            }
#        }
#    }
#
    n_par <- 10
    #dispersion <- solve(out$hessian) # singular
	sigma_g2 <- out$par[1]
    #p_sigma_g2 <- wald_test(sigma_g2, 0, dispersion[1,1], N-n_par) 
	sigma_e2 <- out$par[2]
    #p_sigma_e2 <- wald_test(sigma_e2, 0, dispersion[2,2], N-n_par) 
	V <- diag(out$par[2+1:C])
    #p_V <- mvwald_test( diag(V), rep(0,C), dispersion[2+1:C,2+1:C], n=N, P=n_par)
    W <- diag(out$par[(C+2)+1:C])
    #p_W <- mvwald_test( diag(W), rep(0,C), dispersion(out$hessian)[C+2+1:C,C+2+1:C], n=N, P=n_par)
    #sig2s <- Z * hom2 + diag(vs)
    #for (i in 1:C){
    #    sig2s <- sig2s + V[i,i] * hadamard.prod(P[,i] %*% t(P[,i]), Z) + W[i,i] * diag(diag( P[,i] %*% t(P[,i]) ))
    #}
    l <- out$value * (-1)
    #stopifnot(l > -1e12)
    #print(l)
	#print(hom2)
    #print(V)
    #print(W)

    # estimate hessian matrix
    #hess = hessian(screml_free_loglike, x=out$par, y=y, Z=Z, P=P, C=C, vs=vs)

    return ( list( sigma_g2=sigma_g2, sigma_e2=sigma_e2, V=V, W=W, l=l ) )
                  #wald=list(sigma_g2=p_sigma_g2, sigma_e2=p_sigma_e2, V=p_V, W=p_W) ) )
}

screml_free_loglike<- function(par, y, Z, P, X, C, vs){
	sigma_g2 <- par[1]
	sigma_e2 <- par[2]
	V <- diag(par[2+1:C])
    W <- diag(par[(C+2)+1:C])

    l <- LL(y, Z, P, X, C, vs, sigma_g2, sigma_e2, V, W)
    return ( l )
}

screml_full <- function(
y, Z, P, vs
){
	C <- ncol(P)  # cell type number
	ngam <- C*(C+1)/2 # number of entries in gamma matrix # should it be C + C*(C+1)/2?
    pi <- colMeans(P)
    pd <- scale(P, scale=F)
    S <- (t(pd) %*% pd ) / nrow(P)
    X <- P

	V <- diag(C)[ lower.tri(diag(C),diag=T) ] * as.vector(var(y)) / 2
	W <- diag(C)[ lower.tri(diag(C),diag=T) ] * as.vector(var(y)) / 2
    par <- c( V, W )

	out <- optim( par=par, fn=screml_full_loglike, 
		y=y, Z=Z, P=P, X=X, C=C, ngam=ngam, vs=vs, method = "BFGS", hessian = TRUE)

#    hom2_  <- out$par[1]
#    V_ <- matrix( 0, C, C )
#    V_[lower.tri(V_,diag=T)] <- out$par[1+1:ngam]
#    V_ <- V_ + t(V_)
#    W_ <-  matrix( 0, C, C )
#    W_[lower.tri(W_,diag=T)] <- out$par[(ngam+1)+1:ngam]
#    W_ <- W_ + t(W_)
#
#    if (hom2_ > 1 | sum(diag(V_ %*% S)) + pi %*% V_ %*% pi > 1 | 
#        sum(diag(W_ %*% S)) + pi %*% W_ %*% pi > 1 ) {
#        overwhelm_variance = TRUE
#    } else {
#        overwhelm_variance = FALSE
#    }
#
#    if ( out$convergence != 0 | overwhelm_variance ) {
#        for (i in 1:5) {
#            par_ <- par * rgamma(length(par), 2, scale=1/2)
#            out_ <- optim( par=par_, fn=screml_full_loglike, 
#                y=y, Z=Z, P=P, C=C, ngam=ngam, vs=vs, method = "BFGS", hessian = TRUE)
#
#            if (out_$convergence == 0 & (out_$value * (-1)) > (out$value * (-1)) ) {
#                out <- out_
#            }
#        }
#    }

    V <- matrix( 0, C, C )
    V[lower.tri(V,diag=T)] <- out$par[1:ngam]
    V <- V + t(V)
    W <-  matrix( 0, C, C )
    W[lower.tri(W,diag=T)] <- out$par[ngam+1:ngam]
    W <- W + t(W)
    l <- out$value * (-1)
    #stopifnot(l > -1e12)
	#print(beta)
    #print(V)
    #print(W)

    return ( list( V=V, W=W, l=l ))
}

screml_full_loglike <- function(par, y, Z, P, X, C, ngam, vs){
    sigma_g2 <- 0
    sigma_e2 <- 0
    V <- matrix( 0, C, C )
    V[lower.tri(V,diag=T)] <- par[1:ngam]
    V <- V + t(V)
    W <- matrix( 0, C, C )
    W[lower.tri(W,diag=T)] <- par[ngam+1:ngam]
    W <- W + t(W)

    l <- LL(y, Z, P, X, C, vs, sigma_g2, sigma_e2, V, W)
    return ( l )
}

