install.packages( setdiff(c('numDeriv','Matrix', 'optparse'), rownames(installed.packages())) )
library(optparse)
library(numDeriv)
library(Matrix)
# library(tensor)
# library(abind)

LL <- function(Y, X, K1, vs, g1, V1, e, W, K2=NULL, g2=NULL, V2=NULL,
                r2=NULL, random_MMT=NULL){
    N <- nrow( Y )
    C <- ncol( Y )
    y <- as.vector(t(Y))
    A <- matrix(1, C, C) * g1 + V1
    B <- matrix(1, C, C) * e + W

    sig2s <- kronecker(K1, A) + kronecker(diag(N), B) + diag(as.vector(t(vs)))
    if (!is.null(random_MMT)) {
        sig2s <- sig2s + r2 * random_MMT
    }
    if (!is.null(K2)) {
        A2 <- matrix(1, C, C) * g2 + V2
        sig2s <- sig2s + kronecker(K2, A2)
    }
    # print(sig2s[1:5,1:5])

    # cholesky decomposition
    R <- tryCatch({
        chol(sig2s)
    }, error = function(e) {
        # print(e)
        return(NULL)
    })
    if (is.null(R)) {
        return(1e12)
    }
    
    # calculate B matrix
    m1 <- backsolve(R, forwardsolve(t(R), X))
    m2 <- t(X) %*% m1
    eval <- eigen(m2, symmetric=TRUE, only.values=TRUE)$values
    if (max(eval) / (min(eval) + 1e-99) > 1e8 | min(eval) < 0) return(1e12)
    
    # calculate loglikelihood
    det_Vy <- 2 * sum(log(diag(R)))
    det_XVyX <- determinant(m2, logarithm = TRUE)$modulus[1]
    yBy <- t(y) %*% backsolve(R, forwardsolve(t(R), y)) - t(y) %*% m1 %*% solve(m2) %*% t(m1) %*% y
    L <- 0.5 * (det_Vy + det_XVyX + yBy)
    
    # # inv
    # e  <- eigen(sig2s,symmetric=TRUE)
    # eval <- e$values
    # evec <- e$vectors
    # if (max(eval)/(min(eval)+1e-99) > 1e8 | min(eval)<0) return(1e12)
    # sig2s_inv <- evec %*% diag(1/eval) %*% t(evec)
    # det_sig2s <- sum(log(eval))

    # F <- t(X) %*% sig2s_inv
    # B <- F %*% X
    # eval   <- eigen(B,symmetric=TRUE)$values
    # if (max(eval)/(min(eval)+1e-99) > 1e8 | min(eval)<0) return(1e12)
    # M <- sig2s_inv - t(F) %*% solve(B) %*% F
    # det_B <- determinant(B, logarithm=TRUE)$modulus
    # y <- as.vector(t(Y))
    # L <- 0.5 * ( det_sig2s + det_B + y %*% M %*% y )

    return(L)
}


LL_noChol <- function(Y, X, K1, vs, g1, V1, e, W, K2=NULL, g2=NULL, V2=NULL,
                r2=NULL, random_MMT=NULL){
    N <- nrow( Y )
    C <- ncol( Y )
    y <- as.vector(t(Y))
    A <- matrix(1, C, C) * g1 + V1
    B <- matrix(1, C, C) * e + W

    sig2s <- kronecker(K1, A) + kronecker(diag(N), B) + diag(as.vector(t(vs)))
    if (!is.null(random_MMT)) {
        sig2s <- sig2s + r2 * random_MMT
    }
    if (!is.null(K2)) {
        A2 <- matrix(1, C, C) * g2 + V2
        sig2s <- sig2s + kronecker(K2, A2)
    }
    # print(sig2s[1:5,1:5])
    
    # inv
    e  <- eigen(sig2s,symmetric=TRUE)
    eval <- e$values
    evec <- e$vectors
    if (max(eval)/(min(eval)+1e-99) > 1e8 | min(eval)<0) return(1e12)
    sig2s_inv <- evec %*% diag(1/eval) %*% t(evec)
    det_sig2s <- sum(log(eval))

    F <- t(X) %*% sig2s_inv
    B <- F %*% X
    eval   <- eigen(B,symmetric=TRUE)$values
    if (max(eval)/(min(eval)+1e-99) > 1e8 | min(eval)<0) return(1e12)
    M <- sig2s_inv - t(F) %*% solve(B) %*% F
    det_B <- determinant(B, logarithm=TRUE)$modulus
    y <- as.vector(t(Y))
    L <- 0.5 * ( det_sig2s + det_B + y %*% M %*% y )

    return(L)
}


decode_par <- function(par, C, K2=NULL, random_MMT=NULL){
	g1 <- par[1]
	V1 <- diag(par[1 + 1:C])
    e <- par[C + 2]
    W <- diag(par[C + 2 + 1:C])
    i <- 2 * C + 3
    if (!is.null(K2)) {
        g2 <- par[i]
        V2 <- diag(par[i + 1:C])
        i <- i + C + 1
    } else {
        g2 <- NULL
        V2 <- NULL
    }
    if (!is.null(random_MMT)){
        r2 <- par[i]
    } else {
        r2 <- NULL
    }

    return(list(g1=g1, V1=V1, e=e, W=W, g2=g2, V2=V2, r2=r2))
}


screml_free_loglike <- function(par, args){
    # print(par)
    Y <- args[['Y']]
    X <- args[['X']]
    K1 <- args[['K1']]
    vs <- args[['vs']]
    K2 <- args[['K2']]
    random_MMT <- args[['random_MMT']]
    C <- ncol( Y )

    par_list <- decode_par(par, C, K2, random_MMT)
    for (name in names(par_list)) {
        assign(name, par_list[[name]])
    }
    # print(r2)

    if (args[['chol']]) {
        l <- LL(Y, X, K1, vs, g1, V1, e, W, K2, g2, V2, r2, random_MMT)
    } else {
        l <- LL_noChol(Y, X, K1, vs, g1, V1, e, W, K2, g2, V2, r2, random_MMT)
    }
    # print(l)

    return( l )
}


screml_free <- function(
Y, X, K1, vs, K2=NULL, random_MMT=NULL, par=NULL, seed=NULL, maxit=NULL, 
reltol=NULL, sd=NULL, chol=TRUE
){
    N <- nrow(Y)
    C <- ncol(Y)

    y <- as.vector(t(Y))

    n_var <- 4
    n_par <- 2 * (1 + C)
    if (!is.null(random_MMT)) {
        n_var <- n_var + 1
        n_par <- n_par + 1
    }
    if (!is.null(K2)) {
        n_var <- n_var + 2
        n_par <- n_par + (1 + C)
    }

    # set initial parameters
    if (is.null(par)) {
        hom2 <- median(apply(Y, 2, var)) / n_var
        par <- rep(hom2, n_par)
    }
    if (!is.null(seed)) {
        set.seed(seed)
        if (!is.null(sd)) {
            par <- par + rnorm(length(par), 0, sd)
        } else {
            par <- par + rnorm(length(par), 0, 0.001)
        }
        # par <- par * rgamma(length(par), 2, scale=1/2)
    }
    print(par)
    flush.console()

    args <- list(Y=Y, X=X, K1=K1, vs=vs, K2=K2, random_MMT=random_MMT, chol=chol)
    if (!is.null(maxit) && !is.null(reltol)) {
        out <- optim(par=par, fn=screml_free_loglike, args=args, 
                    method='BFGS', control=list(maxit=maxit, reltol=reltol))
    } else if (!is.null(maxit) && is.null(reltol)) {
        out <- optim(par=par, fn=screml_free_loglike, args=args, 
                    method='BFGS', control=list(maxit=maxit))
    } else if (is.null(maxit) && !is.null(reltol)) {
        out <- optim(par=par, fn=screml_free_loglike, args=args, 
                    method='BFGS', control=list(reltol=reltol))
    } else {
        out <- optim(par=par, fn=screml_free_loglike, args=args, 
                    method='BFGS')
    }

    par_list <- decode_par(out$par, C, K2, random_MMT)
    for (name in names(par_list)) {
        assign(name, par_list[[name]])
    }
    l <- out$value * (-1)

    # estimate hessian matrix
    hess = hessian(screml_free_loglike, x=out$par, args=args)

    # tmp test sig2s singularity
    # A <- matrix(1, C, C) * g1 + V1
    # B <- matrix(1, C, C) * e + W
    # sig2s <- kronecker(K1, A) + kronecker(diag(N), B) + diag(as.vector(t(vs)))
    # if (!is.null(random_MMT)) {
    #     sig2s <- sig2s + r2 * random_MMT
    # }
    # print(l)
    # eval  <- eigen(sig2s,symmetric=TRUE)$values
    # print(eval[1:5])
    # print(min(eval))
    # R <- chol(sig2s)
    # print(R[1:5,1:5])

    return ( list(g1=g1, V1=V1, e=e, W=W, g2=g2, V2=V2, r2=r2, par=par,
                l=l, hess=hess, niter=out$counts[2], convergence=out$convergence, 
                message=out$message))
}



##################
# runs only when script is run by itself
if (sys.nframe() == 0){
    option_list <- list(
        make_option("--Y", type = "character", default = NULL, 
                    help = "CTP"),
        make_option("--X", type = "character", default = NULL, 
                    help = "fixed effect"),
        make_option("--K1", type = "character", default = NULL, 
                    help = "GRM1"),
        make_option("--vs", type = "character", default = NULL, 
                    help = "cell-to-cell variance ctnu"),
        make_option("--K2", type = "character", default = NULL, 
                    help = "GRM2"),
        make_option("--batch", type = "character", default = NULL, 
                    help = "batch effect design matrix"),
        make_option("--noChol", action="store_true", default = FALSE, 
                    help = "whether use Cholesky decomposition"),
        make_option("--par", type = "character", default = NULL, 
                    help = "initial parameters"),
        make_option("--seed", type = "integer", default = NULL, 
                    help = "seed for random initialization"),
        make_option("--sd", type = "numeric", default = NULL, 
                    help = "sd for normal distribution in initialization"),
        make_option("--maxit", type = "integer", default = NULL, 
                    help = "maximum number of iterations"),
        make_option("--reltol", type = "numeric", default = NULL, 
                    help = "relative tolerance for convergence"),
        make_option("--out", type = "character", default = NULL, 
                    help = "output file (.rda)")
    )
    # Parse the command line options
    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)

    Y <- as.matrix(read.table(opt$Y))
    C <- ncol(Y)
    X <- as.matrix(read.table(opt$X))
    K1 <- as.matrix(read.table(opt$K1))
    vs <- as.matrix(read.table(opt$vs))
    if (!is.null(opt$K2)) {
        K2 <- as.matrix(read.table(opt$K2))
    } else {
        K2 <- NULL
    }
    if (!is.null(opt$batch)) {
        batch <- as.matrix(read.table(opt$batch))
        random_MMT <- kronecker(batch %*% t(batch), matrix(1, C, C))
    } else {
        random_MMT <- NULL
    }
    if (!is.null(opt$par)) {
        par <- scan(opt$par)
    } else {
        par <- NULL
    }

    if (opt$noChol) {
        print("R Fitting without Cholesky decomposition")
        out <- screml_free(Y=Y, X=X, K1=K1, vs=vs, K2=K2, 
                    random_MMT=random_MMT, par=par, seed=opt$seed, sd=opt$sd,
                    maxit=opt$maxit, reltol=opt$reltol, chol=FALSE)
    } else {
        print("R Fitting with Cholesky decomposition")
        out <- screml_free(Y=Y, X=X, K1=K1, vs=vs, K2=K2, 
                    random_MMT=random_MMT, par=par, seed=opt$seed, sd=opt$sd,
                    maxit=opt$maxit, reltol=opt$reltol)
    }
    save(out, file=opt$out)
}
