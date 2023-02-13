##' @title Mediation Analysis for ZEro-inflated mediators
##'
##' @description A novel mediation modeling approach to address zero-inflated mediators containing true zeros and false zeros.
##'
##' @details For an independent variable \eqn{X}, a zero-inflated mediator \eqn{M} and a continuous outcome variable \eqn{Y}, the following regression equation is used to model the association between \eqn{Y} and \eqn{(X,M)}:
##' \deqn{Y_{xm1_{(m>0)}}=\beta_0+\beta_1m+\beta_2 1_{(m>0)}+\beta_3x+\beta_4x1_{(m>0)}+\beta_5xm+\epsilon}
##'
##' @param data a data frame containing variables: an independent variable X, a mediator M, an outcome Y, and confounder variables Z (if any). See example dataset: data(zinb10)
##' @param distM an optional character value for distribution to be used for the mediator. Possible choices are 'zilonm', 'zinbm', or 'zipm' for zero-inflated log-normal, negative binomial, or Poisson mediators respectively. By default, all three distributions will be fitted and the final mediation model is selected by AIC
##' @param K a user supplied sequence for the number of component K in the zero-inflated mixture mediators. Default is K = 1 for zero-inflated non-mixture mediators
##' @param selection model selection criterion to be used when more than one model is fitted. Possible choices are 'AIC' and 'BIC'. Default is 'AIC'
##' @param X name of the independent variable
##' @param M name of the mediator variable
##' @param Y name of the outcome variable
##' @param Z name(s) of confounder variables
##' @param XMint a logical vector of length 2 indicating whether to include the interaction terms between (i) X and 1(M>0) and (ii) X and M. Default is c(TRUE, FALSE)
##' @param x1 the first value of independent variable of interest
##' @param x2 the second value of independent variable of interest
##' @param zval the value of confounders to be conditional on in estimating effects
##' @param mval the fixed value of mediator to be conditional on in estimating CDE
##' @param B the upper bound value B to be used in the probability mechanism of observing false zeros
##' @param seed an optional seed number to control randomness
##' @param ncore number of cores available for parallel computing
##' @return a list containing:
##' - results_effects: a data frame for the results of estimated effects,
##' - results_parameters: a data frame for the results of model parameters,
##' - selected_model_name: a string for the distribution of M and number of components K selected in the final mediation model
##' - BIC: a numeric value for the BIC of the final mediation model,
##' - AIC: a numeric value for the AIC of the final mediation model,
##' - models: a list with all fitted models
##' - analysis2_out: a list with output from analysis2() function (used for internal check)
##' @author Zhigang Li <zhigang.li@@ufl.edu>
##' Meilin Jiang <meilin.jiang@@ufl.edu>
##' @import stats foreach doParallel MASS
##' @importFrom flexmix flexmix clusters parameters FLXMRglm
##' @importFrom pracma hessian
##' @importFrom numDeriv grad jacobian
##' @export
##' @examples
##' data(zinb10)
##' \donttest{
##' maze_out <- MAZE(data=zinb10,
##' distM=c('zilonm', 'zinbm', 'zipm'),
##' K = 1,
##' selection = 'AIC',
##' X='X', M='Mobs', Y='Y', Z=NULL,
##' XMint = c(TRUE, FALSE),
##' x1=0, x2=1, zval = NULL, mval = 0,
##' B=20, seed=1)
##' ## results of selected mediation model
##' maze_out$results_effects # indirect and direct effects
##' maze_out$results_parameters # model parameters
##' maze_out$BIC; maze_out$AIC # BIC and AIC of the selected mediation model
##' maze_out$selected_disM # distribution of the mediator in the selected mediation model
##' }

MAZE <- function(data, distM = c("zilonm", "zinbm", "zipm"),
    K = 1, selection = "AIC", X, M, Y, Z = NULL, XMint = c(TRUE,
        FALSE), x1, x2, zval = NULL, mval = 0, B = 20, seed = 1,
    ncore = 1) {
    distM_sequence <- distM
    K_sequence <- K
    limits <- 0.001
    # set.seed(seed)

    data <- as.data.frame(data)
    dat <- data.frame(X = data[, X], Y = data[, Y], Mobs = data[,
        M])
    num_Z <- length(Z)
    if (is.null(Z)) {
        Z_names <- NULL
    } else {
        dat <- data.frame(dat, data[, Z])
        Z_names <- paste0("Z", 1:num_Z)
        names(dat)[3 + 1:num_Z] <- Z_names
    }
    # complete data analysis
    dat <- dat[complete.cases(dat), ]

    out <- tryCatch({
        realanalysis(dat, distM_sequence, K_sequence, selection,
            XMint, x1, x2, zval, num_Z, Z_names, mval, limits,
            B, seed, ncore)
    }, error = function(e) {
        print(e)
        # list(results_effects = NA, results_parameters =
        # NA, BIC = Inf, AIC = Inf, e = e)
    })
    return(out)
}




