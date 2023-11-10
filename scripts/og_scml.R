library(numDeriv)
library(mvtnorm)
library(matrixcalc)
args <- commandArgs(trailingOnly=TRUE)
y_f <- args[1]
Z_f <- args[2]
P_f <- args[3]
vs_f <- args[4]
model <- args[5]
out_f <- args[6]

# read data
y <- scan(y_f)
Z <- as.matrix(read.table(Z_f))
P <- as.matrix(read.table(P_f))
vs <- scan(vs_f)

screml_null <- function(
y, Z, P, vs
){

	C      <- ncol(P)  # cell type number
	beta   <- solve(t(P)%*%P) %*% ( t(P) %*% y ) # cell type effect
    par <- c( beta )

	out <- optim( par=par, fn=screml_null_loglike, 
		y=y, P=P, C=C, vs=vs, method = "BFGS", hessian = TRUE)
    
    i <- 1
    while (out$convergence != 0) {
        print(i)
        i <- i+1
        par_ <- par * rgamma(length(par), 2, scale=1/2) 
        out <- optim( par=c( par_ ), fn=screml_null_loglike, 
            y=y, P=P, C=C, vs=vs, method = "BFGS", hessian = TRUE)
    }

    beta <- out$par 
    l <- out$value * (-1)
	print(beta)
    
    # estimate hessian matrix
    hess = hessian(screml_null_loglike, x=out$par, y=y, P=P, C=C, vs=vs)

    return ( list( beta=beta, l=l, hess=hess ) )
}

screml_null_loglike<- function(par, y, P, C, vs){
	beta  <- par

	yd    <- y - P %*% beta
	sig2s <- vs

	if( any( sig2s < 0 ) ) return(1e12)

	(length(y) * log(2*pi) + sum(log( sig2s )) + sum( yd^2 / sig2s ))/2 

}

screml_hom <- function(
y, Z, P, vs
){
	C      <- ncol(P)  # cell type number
    pi  <- colMeans(P)
    pd  <- scale(P, scale=F)
    S   <- (t(pd) %*% pd ) / nrow(P)

	beta   <- solve(t(P)%*%P) %*% ( t(P) %*% y ) # cell type effect
	hom2 <- var(y - P %*% beta) 
    par  <- c( hom2, beta )

	out <- optim( par=par, fn=screml_hom_loglike, 
		y=y, Z=Z, P=P, C=C, vs=vs, method = "BFGS", hessian = TRUE)

	hom2_ <- out$par[1]
	beta_  <- out$par[1+1:C]

    if (hom2_ > 1 | beta_ %*% S %*% beta_ > 1 ) {
        overwhelm_variance = TRUE
    } else {
        overwhelm_variance = FALSE
    }

    if (out$convergence != 0 | overwhelm_variance) {
        for (i in 1:5) {
            par_ <- par * rgamma(length(par), 2, scale=1/2) 
            out_ <- optim( par=par_, fn=screml_hom_loglike, 
                y=y, Z=Z, P=P, C=C, vs=vs, method = "BFGS", hessian = TRUE)

            if (out_$convergence == 0 & (out_$value * (-1)) > (out$value * (-1)) ) {
                out <- out_
            }
        }
    }

	hom2 <- out$par[1]
	beta  <- out$par[1+1:C]
	sig2s <- Z * hom2 + diag(vs) 
    l <- out$value * (-1)
	print(hom2)
	print(beta)
    print(l)

    # estimate hessian matrix
    hess = hessian(screml_hom_loglike, x=out$par, y=y, Z=Z, P=P, C=C, vs=vs)

    #try(solve(out$hessian))
    return ( list(hom2=hom2, beta=beta, l=l, hess=hess, sig2s=sig2s) )
}

screml_hom_loglike <- function(par, y, Z, P, C, vs){
	hom2 <- par[1]
	beta <- par[1+1:C] 

	sig2s <- Z * hom2 + diag(vs) 
    eval   <- eigen(sig2s,symmetric=TRUE)$values
    if (max(eval)/(min(eval)+1e-99) > 1e8 | min(eval)<0) return(1e12)

	if( any( diag(sig2s) < 0 ) ) return(1e12)
    
    dmvnorm(y, mean=P %*% beta, sigma = sig2s, log=TRUE) * (-1)
}

screml_iid <- function(
y, Z, P, vs
){
	C    <- ncol(P)  # cell type number
    pi  <- colMeans(P)
    pd  <- scale(P, scale=F)
    S   <- (t(pd) %*% pd ) / nrow(P)

	beta <- solve(t(P)%*%P) %*% ( t(P) %*% y ) # cell type effect
	hom2 <- var(y - P %*% beta) / 3
	V    <- var(y - P %*% beta) / 3 
    W    <- var(y - P %*% beta) / 3
    par  <- c( hom2, beta, V, W )

	out <- optim( par=par, fn=screml_iid_loglike, 
		y=y, Z=Z, P=P, C=C, vs=vs, method = "BFGS", hessian = TRUE)

	hom2_ <- out$par[1]
	beta_ <- out$par[1+1:C]
    V_    <- out$par[C+2]
    W_    <- out$par[C+3]
	V_ <- diag(C) * V_
    W_ <- diag(C) * W_

    if (hom2_ > 1 | beta_ %*% S %*% beta_ > 1 | sum(diag(V_ %*% S)) + pi %*% V_ %*% pi > 1 | 
        sum(diag(W_ %*% S)) + pi %*% W_ %*% pi > 1 ) {
        overwhelm_variance = TRUE
    } else {
        overwhelm_variance = FALSE
    }

    if ( out$convergence != 0 | overwhelm_variance ) {
        for (i in 1:5) {
            par_ <- par * rgamma(length(par), 2, scale=1/2)
            out_ <- optim( par=par_, fn=screml_iid_loglike, 
                y=y, Z=Z, P=P, C=C, vs=vs, method = "BFGS", hessian = TRUE)

            if (out_$convergence == 0 & (out_$value * (-1)) > (out$value * (-1)) ) {
                out <- out_
            }
        }
    }

	hom2 <- out$par[1]
	beta <- out$par[1+1:C]
    V    <- out$par[C+2]
    W    <- out$par[C+3]
    sig2s <- Z * hom2 + V * hadamard.prod( P %*% t(P), Z ) + W * diag(diag( P %*% t(P) )) + diag(vs)
	V <- diag(C) * V
    W <- diag(C) * W
    l <- out$value * (-1)
	#print(hom2)
	#print(beta)
    #print(V)
    #print(W)
    
    # estimate hessian
    hess = hessian(screml_iid_loglike, x=out$par, y=y, Z=Z, P=P, C=C, vs=vs)

    return ( list(hom2=hom2, beta=beta, V=V, W=W, l=l, hess=hess, sig2s=sig2s) )
}

screml_iid_loglike<- function(par, y, Z, P, C, vs){
	hom2 <- par[1]
	beta    <- par[1+1:C] 
	V <- par[C+2]
    W <- par[C+3]

	sig2s <- Z * hom2 + V * hadamard.prod( P %*% t(P), Z ) + W * diag(diag( P %*% t(P) )) + diag(vs) 
    eval   <- eigen(sig2s,symmetric=TRUE)$values
    if (max(eval)/(min(eval)+1e-99) > 1e8 | min(eval)<0) return(1e12)

	if( any( diag(sig2s) < 0 ) ) return(1e12)

    dmvnorm(y, mean=P %*% beta, sigma = sig2s, log=TRUE) * (-1)
}

screml_free <- function(
y, Z, P, vs
){
	C      <- ncol(P)  # cell type number
    pi  <- colMeans(P)
    pd  <- scale(P, scale=F)
    S   <- (t(pd) %*% pd ) / nrow(P)

	beta   <- solve(t(P)%*%P) %*% ( t(P) %*% y ) # cell type effect
	hom2   <- var(y - P %*% beta) / 3
	V      <- rep(1,C) * as.vector(var(y - P %*% beta)) / 3
	W      <- rep(1,C) * as.vector(var(y - P %*% beta)) / 3
    par    <- c( hom2, beta, V, W)

	out <- optim( par=par, fn=screml_free_loglike, 
		y=y, Z=Z, P=P, C=C, vs=vs, method = "BFGS", hessian = TRUE)

	hom2_ <- out$par[1]
	beta_ <- out$par[1+1:C]
	V_ <- diag(out$par[(C+1)+1:C])
    W_ <- diag(out$par[(2*C+1)+1:C])

    if (hom2_ > 1 | beta_ %*% S %*% beta_ > 1 | sum(diag(V_ %*% S)) + pi %*% V_ %*% pi > 1 | 
        sum(diag(W_ %*% S)) + pi %*% W_ %*% pi > 1 ) {
        overwhelm_variance = TRUE
    } else {
        overwhelm_variance = FALSE
    }

    if ( out$convergence != 0 | overwhelm_variance ) {
        for (i in 1:5) {
            par_ <- par * rgamma(length(par), 2, scale=1/2)
            out_ <- optim( par=par_, fn=screml_free_loglike, 
                y=y, Z=Z, P=P, C=C, vs=vs, method = "BFGS", hessian = TRUE)

            if (out_$convergence == 0 & (out_$value * (-1)) > (out$value * (-1)) ) {
                out <- out_
            }
        }
    }

	hom2 <- out$par[1]
	beta <- out$par[1+1:C]
	V <- diag(out$par[(C+1)+1:C])
    W <- diag(out$par[(2*C+1)+1:C])
    sig2s <- Z * hom2 + diag(vs)
    for (i in 1:C){
        sig2s <- sig2s + V[i,i] * hadamard.prod(P[,i] %*% t(P[,i]), Z) + W[i,i] * diag(diag( P[,i] %*% t(P[,i]) ))
    }
    l <- out$value * (-1)
    print(l)
	print(hom2)
	print(beta)
    print(V)
    print(W)

    # estimate hessian matrix
    hess = hessian(screml_free_loglike, x=out$par, y=y, Z=Z, P=P, C=C, vs=vs)

    return ( list(hom2=hom2, beta=beta, V=V, W=W, l=l, hess=hess, sig2s=sig2s ))
}

screml_free_loglike<- function(par, y, Z, P, C, vs){
	hom2 <- par[1]
	beta <- par[1+1:C] 
	V <- diag(par[(C+1)+1:C])
    W <- diag(par[(2*C+1)+1:C])

	sig2s <- Z * hom2 + diag(vs) 
    for (i in 1:C){
        sig2s <- sig2s + V[i,i] * hadamard.prod(P[,i] %*% t(P[,i]), Z) + W[i,i] * diag(diag( P[,i] %*% t(P[,i]) ))
    }
    eval   <- eigen(sig2s,symmetric=TRUE)$values
    if (max(eval)/(min(eval)+1e-99) > 1e8 | min(eval)<0) return(1e12)

	if( any( diag(sig2s) < 0 ) ) return(1e12)

    dmvnorm(y, mean=P %*% beta, sigma = sig2s, log=TRUE) * (-1)
}

screml_full <- function(
y, Z, P, vs
){
	C      <- ncol(P)  # cell type number
	ngam   <- C*(C+1)/2 # number of entries in gamma matrix # should it be C + C*(C+1)/2?
    pi  <- colMeans(P)
    pd  <- scale(P, scale=F)
    S   <- (t(pd) %*% pd ) / nrow(P)

	beta   <- solve(t(P)%*%P) %*% ( t(P) %*% y ) # cell type effect
    hom2   <- var(y - P %*% beta) / 3
	V      <- diag(C)[ lower.tri(diag(C),diag=T) ] * as.vector(var(y - P %*% beta)) / 3
	W      <- diag(C)[ lower.tri(diag(C),diag=T) ] * as.vector(var(y - P %*% beta)) / 3
    par    <- c( hom2, beta, V, W)

	out <- optim( par=par, fn=screml_full_loglike, 
		y=y, Z=Z, P=P, C=C, ngam=ngam, vs=vs, method = "BFGS", hessian = TRUE)

    hom2_  <- out$par[1]
	beta_  <- out$par[1+1:C]
    V_ <- matrix( 0, C, C )
    V_[lower.tri(V_,diag=T)] <- out$par[(C+1)+1:ngam]
    V_ <- V_ + t(V_)
    W_ <-  matrix( 0, C, C )
    W_[lower.tri(W_,diag=T)] <- out$par[(C+ngam+1)+1:ngam]
    W_ <- W_ + t(W_)

    if (hom2_ > 1 | beta_ %*% S %*% beta_ > 1 | sum(diag(V_ %*% S)) + pi %*% V_ %*% pi > 1 | 
        sum(diag(W_ %*% S)) + pi %*% W_ %*% pi > 1 ) {
        overwhelm_variance = TRUE
    } else {
        overwhelm_variance = FALSE
    }

    if ( out$convergence != 0 | overwhelm_variance ) {
        for (i in 1:5) {
            print(i)
            par_ <- par * rgamma(length(par), 2, scale=1/2)
            out_ <- optim( par=par_, fn=screml_full_loglike, 
                y=y, Z=Z, P=P, C=C, ngam=ngam, vs=vs, method = "BFGS", hessian = TRUE)

            print(out_$convergence)

            if (out_$convergence == 0 & (out_$value * (-1)) > (out$value * (-1)) ) {
                out <- out_
            }
        }
    }

    hom2  <- out$par[1]
	beta  <- out$par[1+1:C]
    V <- matrix( 0, C, C )
    V[lower.tri(V,diag=T)] <- out$par[(C+1)+1:ngam]
    V <- V + t(V)
    W <-  matrix( 0, C, C )
    W[lower.tri(W,diag=T)] <- out$par[(C+ngam+1)+1:ngam]
    W <- W + t(W)
    l <- out$value * (-1)
	#print(beta)
    #print(V)
    #print(W)

    return ( list(hom2=hom2, beta=beta, V=V, W=W, l=l ))
}

screml_full_loglike <- function(par, y, Z, P, C, ngam, vs){
    hom2 <- par[1]
	beta <- par[1+1:C] 
    V <- matrix( 0, C, C )
    V[lower.tri(V,diag=T)] <- par[(C+1)+1:ngam]
    V <- V + t(V)
    W <- matrix( 0, C, C )
    W[lower.tri(W,diag=T)] <- par[(C+ngam+1)+1:ngam]
    W <- W + t(W)

	sig2s <- Z * hom2 + diag(vs) 
    for (p in 1:C){
        for (q in 1:C){
            sig2s <- sig2s + V[p,q]*hadamard.prod(P[,p]%*%t(P[,q]), Z) + W[p,q]*diag(diag( P[,p]%*%t(P[,q]) ))
        }
    }
    eval   <- eigen(sig2s,symmetric=TRUE)$values
    if (max(eval)/(min(eval)+1e-99) > 1e8 | min(eval)<0) return(1e12)

	if( any( diag(sig2s) < 0 ) ) return(1e12)

    dmvnorm(y, mean=P %*% beta, sigma = sig2s, log=TRUE) * (-1)
}

if (model == 'null') {
    out = screml_null(y, Z, P, vs)
    save(out, file=out_f)
} else if (model == 'hom') {
    out = screml_hom(y, Z, P, vs)
    save(out, file=out_f)
} else if (model == 'iid') {
    out = screml_iid(y, Z, P, vs)
    save(out, file=out_f)
} else if (model == 'free') {
    out = screml_free(y, Z, P, vs)
    save(out, file=out_f)
} else if (model == 'full') {
    out = screml_full(y, Z, P, vs)
    save(out, file=out_f)
}
