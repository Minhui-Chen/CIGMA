library(numDeriv)
library(emdbook)
library(matrixcalc)
library(reticulate) # to use python
use_python('/apps/software/gcc-6.2.0/python/3.6.0/bin/python3')
source_python('wald.py')
np <- import('numpy')

# variance matrix of y
SIGMA <- function(Z, ctvs, N, C, sigma_g2=0, sigma_e2=0, 
                  V=matrix(0,C,C), W=matrix(0,C,C)) 
{
    sigma_y2 <- kronecker( Z, sigma_g2 + V ) 
    sigma_y2 <- sigma_y2 + kronecker( diag(1,N), sigma_e2 + W ) 
    sigma_y2 <- sigma_y2 + diag( as.vector(t(ctvs)) )

    return( sigma_y2 )
}
    
# loglikelihood function
LL <- function(y, Z, X, N, C, ctvs, sigma_g2, sigma_e2, V, W){
    if ( N > 300 ) {
        return(1e12) # to break REML, since it take too much time
    }

    sigma_y2 <- SIGMA( Z, ctvs, N, C, sigma_g2, sigma_e2, V, W ) 

    # inverse variance
    eig <- eigen(sigma_y2, symmetric=TRUE) 
    eval <- eig$value
    evec <- eig$vectors
    if (max(eval)/(min(eval)+1e-99) > 1e8 | min(eval)<0) return(1e12)

    #sigma_y2_inv <- evec %*% diag(1/eval) %*% t(evec) # solve(sigma_y2) 

    # calculate B matrix
    m1 <- t(X) %*% evec %*% diag(1/eval) %*% t(evec)
    m2 <- m1 %*% X
    #eval <- eigen(m2, symmetric=TRUE)$values
    #if (max(eval)/(min(eval)+1e-99) > 1e8 | min(eval)<0) return(1e12)

    #B <- sigma_y2_inv - t(m1) %*% solve(m2) %*% m1

    # calculate loglikelihood
    L <-  sum( np$log(eig$value) ) # determinant(sigma_y2, logarithm=TRUE)$modulus 
    L <- L + determinant(m2, logarithm=TRUE)$modulus
    yPy <- t(y) %*% evec %*% diag(1/eval) %*% t(evec) %*% y - t(y) %*% t(m1) %*% solve(m2) %*% m1 %*% y
    L <- L + yPy
    L <- 0.5 * L
    return( L )
}

screml_hom <- function(
cty, Z, ctvs
){
	N <- nrow(cty)  
	C <- ncol(cty)  # cell type number
    y <- as.vector( t(cty) )
    X <- kronecker( rep(1,N), diag(1,C) )

	sigma_g2 <- var(y) / 2
    sigma_e2 <- sigma_g2
    par  <- c( sigma_g2, sigma_e2 )

	print( system.time( out <- optim( par=par, fn=screml_hom_loglike, 
		y=y, Z=Z, X=X, N=N, C=C, ctvs=ctvs, method = "BFGS", hessian = TRUE) ))

	sigma_g2 <- out$par[1]
    sigma_e2 <- out$par[2]
    l <- out$value * (-1)

    return ( list( sigma_g2=sigma_g2, sigma_e2=sigma_e2, l=l ) )
}

screml_hom_loglike <- function(par, y, Z, X, N, C, ctvs){
	sigma_g2 <- par[1]
	sigma_e2 <- par[2]
    V <- matrix( 0, C, C )
    W <- matrix( 0, C, C )

    l <- LL(y, Z, X, N, C, ctvs, sigma_g2, sigma_e2, V, W)
    return ( l )
}

screml_free <- function(
cty, Z, ctvs
){
    N <- nrow(cty)
	C <- ncol(cty)  # cell type number
    y <- as.vector( t(cty) )
    X <- kronecker( rep(1,N), diag(1,C) )

	sigma_g2 <- var(y) / 4
    sigma_e2 <- sigma_g2
	V <- rep(1,C) * as.vector( sigma_g2 )
	W <- rep(1,C) * as.vector( sigma_g2 )
    par <- c( sigma_g2, sigma_e2, V, W )

	print(system.time( out <- optim( par=par, fn=screml_free_loglike, 
		y=y, Z=Z, X=X, N=N, C=C, ctvs=ctvs, method = "BFGS", hessian = TRUE) ))

    #n_par <- 10
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

screml_free_loglike<- function(par, y, Z, X, N, C, ctvs){
	sigma_g2 <- par[1]
	sigma_e2 <- par[2]
	V <- diag(par[2+1:C])
    W <- diag(par[(C+2)+1:C])

    l <- LL(y, Z, X, N, C, ctvs, sigma_g2, sigma_e2, V, W) 
    return ( l )
}

screml_full <- function(
cty, Z, ctvs
){
	N <- nrow(cty)  
	C <- ncol(cty)  # cell type number
	ngam <- C*(C+1)/2 # number of entries in gamma matrix # should it be C + C*(C+1)/2?
    y <- as.vector( t(cty) )
    X <- kronecker( rep(1,N), diag(1,C) )

	V <- diag(C)[ lower.tri(diag(C),diag=T) ] * as.vector(var(y)) / 2
	W <- diag(C)[ lower.tri(diag(C),diag=T) ] * as.vector(var(y)) / 2
    par <- c( V, W )

	print( system.time( out <- optim( par=par, fn=screml_full_loglike, 
		y=y, Z=Z, X=X, N=N, C=C, ngam=ngam, ctvs=ctvs, method = "BFGS", hessian = TRUE) ))

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

screml_full_loglike <- function(par, y, Z, X, N, C, ngam, ctvs){
    sigma_g2 <- 0
    sigma_e2 <- 0
    V <- matrix( 0, C, C )
    V[lower.tri(V,diag=T)] <- par[1:ngam]
    V <- V + t(V)
    W <- matrix( 0, C, C )
    W[lower.tri(W,diag=T)] <- par[ngam+1:ngam]
    W <- W + t(W)

    l <- LL(y, Z, X, N, C, ctvs, sigma_g2, sigma_e2, V, W) 
    return ( l )
}

