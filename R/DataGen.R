##' Generate data under zero-inflated mediation models
##' 
##' @title DataGen 
##' @param distM an character value for distribution to be used for the mediator. Possible choices are 'zilonm','zinbm', or 'zipm' for zero-inflated log-normal, negative binomial, or Poisson mediators respectively.
##' @param theta vector of true parameter values
##' @param K a user supplied sequence for the number of component K in the zero-inflated mixture mediators. Default is K = 1 for zero-inflated non-mixture mediators
##' @param num_Z number of confounder variables
##' @param n number of observations
##' @param B the upper bound value B to be used in the probability mechanism of observing false zeros
##' @param x1 the first value of independent variable of interest
##' @param x2 the second value of independent variable of interest
##' @param zval the value of confounders to be conditional on in estimating effects
##' @return 
##' dat: a data frame containing variables:
##' - X: an independent variable, 
##' - Mobs: observed mediator values (with possibly false zeros)
##' - M: true mediator values, 
##' - Y: an outcome,
##' - Z: confounder variables (if any)
##' true_eff: a vector containing true effects
##' @author Zhigang Li <zhigang.li@@ufl.edu>  
##' Meilin Jiang <meilin.jiang@@ufl.edu>
##' @import stats
##' @export
##' @examples 
##' betas.tr <- c(2, 0.12, -6.6, 6.3, -3.8)
##' delta.tr <- 1.1
##' alpha0_k.tr <- c(0.4, 1.1)
##' alpha1_k.tr <- c(0.1, 0.5)
##' xi0.tr <- -1.5
##' psi_km1.tr <- c(0.6)
##' gammas.tr <- c(-1.8, 0.5)
##' eta.tr <- 1  
##' theta <- c(betas.tr, delta.tr, rbind(alpha0_k.tr,alpha1_k.tr), 
##' xi0.tr, psi_km1.tr, gammas.tr, eta.tr)
##' DataGen(distM='zilonm', theta, K=2, num_Z=0, n=200, B=20, x1=0, x2=1, zval=NULL)

DataGen <- function(distM, theta, K, num_Z = 0, n, B, x1,
    x2, zval = NULL) {
    X <- rnorm(n)
    dat_placeholder <- data.frame(NULL)
    class(dat_placeholder) <- c(distM, class(dat_placeholder))
    theta_trans <- trans(dat = dat_placeholder, theta, K,
        xval = X, num_Z, zval = NULL)

    eps <- rnorm(n, 0, theta_trans[["delta"]])

    # 1(M > 0)
    ind <- rbinom(n, 1, 1 - theta_trans[["Del_i"]])
    # ind = 1: M is from log-normal mixture dist 1 (M
    # is from the corresponding log-normal
    # distribution in the mixture)
    L_ik <- t(rmultinom(n, 1, theta_trans[["psi_k"]]))
    M <- rowSums(L_ik * apply(theta_trans[["mu_ik"]], 2,
        function(mu) rlnorm(n, mu, theta_trans[["sig"]])))

    # ind = 0: M is an excessive zero
    M[ind == 0] <- 0
    L <- apply(L_ik, 1, which.max)
    L[ind == 0] <- NA

    Y <- theta_trans[["beta0"]] + theta_trans[["beta1"]] *
        M + theta_trans[["beta2"]] * ind + theta_trans[["beta3"]] *
        X + theta_trans[["beta4"]] * X * ind + eps
    # probability of observing false zeros
    if (theta_trans[["eta"]] == Inf) {
        pf0 <- rep(0, n)
    } else {
        pf0 <- exp(-M * theta_trans[["eta"]]^2)
    }
    pf0[M > B] <- 0
    R <- rbinom(n, 1, 1 - pf0)  # R=0: false zeros
    dat <- data.frame(X, Y, M, R, Mobs = M * R, L)

    true_eff <- effects(dat = dat_placeholder, theta, x1,
        x2, K, num_Z, zval, calculate_se = F, vcovar = NULL,
        Group1 = F)$eff
    out <- list(dat = dat, true_eff = true_eff)
    return(out)
}
