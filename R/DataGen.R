##' @title DataGen
##' @description Generate data under zero-inflated mediation models and calculate the true effects
##'
##' @param distM an character value for distribution to be used for the mediator. Possible choices are 'zilonm', 'zinbm', or 'zipm' for zero-inflated log-normal, negative binomial, or Poisson mediators respectively.
##' @param theta vector of true parameter values
##' @param K a user supplied sequence for the number of component K in the zero-inflated mixture mediators. Default is K = 1 for zero-inflated non-mixture mediators
##' @param num_Z number of confounder variables
##' @param n number of observations
##' @param B the upper bound value B to be used in the probability mechanism of observing false zeros
##' @param x1 the first value of independent variable of interest
##' @param x2 the second value of independent variable of interest
##' @param zval the value of confounders to be conditional on in estimating effects
##' @param mval the fixed value of mediator to be conditional on in estimating CDE
##' @return
##' true_eff: a vector containing true effects
##' dat: a data frame containing variables:
##' - X: an independent variable,
##' - Mobs: observed mediator values (with possibly false zeros)
##' - M: true mediator values,
##' - Y: an outcome,
##' - Z: confounder variables (if any)
##' @author Zhigang Li <zhigang.li@@ufl.edu>
##' Meilin Jiang <meilin.jiang@@ufl.edu>
##' @import stats
##' @export
##' @examples
##' betas.tr <- c(2, 0.12, -6.6, 6.3, -3.8, 0)
##' delta.tr <- 1.1
##' alpha0_k.tr <- c(0.4, 1.1)
##' alpha1_k.tr <- c(0.1, 0.5)
##' alphas.tr <- rbind(alpha0_k.tr,alpha1_k.tr)
##' xi0.tr <- -1.5
##' psi_km1.tr <- c(0.6)
##' gammas.tr <- c(-1.8, 0.5)
##' eta.tr <- 1
##' theta <- c(betas.tr, delta.tr, alphas.tr,
##' xi0.tr, psi_km1.tr, gammas.tr, eta.tr)
##' out <- DataGen(distM='zilonm', theta, K=2, num_Z=0, n=200, B=20, x1=0, x2=1, zval=NULL, mval=0)
##' (true_eff <- out$true_eff)
##' dat <- out$dat


DataGen <- function(distM, theta, K, num_Z = 0, n, B, x1, x2,
    zval = NULL, mval = 0) {
    dat_placeholder <- data.frame(NULL)
    class(dat_placeholder) <- c(distM, class(dat_placeholder))

    out <- DataGen_call(dat_placeholder, theta, K, num_Z, n,
        B, x1, x2, zval, mval)
    return(out)
}

