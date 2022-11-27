##' Mediation Analysis for ZEro-inflated mediators
##' 
##' A novel mediation modeling approach to address zero-inflated mediators containing true zeros and false zeros.
##' 
##' For an independent variable \eqn{X}, a zero-inflated mediator \eqn{M} and a continuous outcome variable \eqn{Y}, the following regression equation is used to model the association between $Y$ and $(X,M)$:
##' \deqn{Y_{xm1_{(m>0)}}=\beta_0+\beta_1m+\beta_2 1_{(m>0)}+\beta_3x+\beta_4x1_{(m>0)}+\beta_5xm+\epsilon}
##' 
##'  
##' @title Mediation Analysis for ZEro-inflated mediators
##' @param data a data frame containing variables an independent variable X, a mediator M, an outcome Y, and confounder variables Z (if any). See example dataset: data(zip10)
##' @param distM an optional character value for distribution to be used for the mediator. Possible choices are 'zilonm','zinbm', or 'zipm' for zero-inflated log-normal, negative binomial, or Poisson mediators respectively. By default, all three distributions will be fitted and the final mediation model is selected by AIC
##' @param X name of the independent variable
##' @param M name of the mediator variable
##' @param Y name of the outcome variable
##' @param Z name(s) of confounder variables
##' @param x1 the first value of independent variable of interest
##' @param x2 the second value of independent variable of interest
##' @param zval the value of confounders to be conditional on in estimating effects
##' @param B the upper bound value B to be used in the probability mechanism of observing false zeros
##' @param seed an optional seed number to control randomness
##' @return a list containing:  
##' - results_effects: a data frame for the results of estimated mediation effect,  
##' - results_parameters: a data frame for the results of model parameters,  
##' - BIC: a numeric value for the BIC of the final mediation model,  
##' - AIC: a numeric value for the AIC of the final mediation model,  
##' - selected_distM: a character value for the distribution of M selected in the final mediation model
##' @author Zhigang Li <zhigang.li@@ufl.edu>  
##' Meilin Jiang <meilin.jiang@@ufl.edu>
##' @import stats numDeriv 
##' @importFrom flexmix flexmix clusters parameters FLXMRglm
##' @importFrom pracma hessian
##' @importFrom MASS glm.nb
##' @export
##' @examples 
##' data(zip10)
##' \dontrun{
##' MAZE(data=zip10, distM=c('zilonm', 'zinbm', 'zipm'), 
##' X='X', M='Mobs', Y='Y', Z=NULL, x1=0, x2=1, B=20, seed=1)
##' }

MAZE <- function(data, distM = c("zilonm", "zinbm", "zipm"), X, M, Y,
    Z = NULL, x1, x2, zval = NULL, B = 20, seed = 1) {
    select_K_range = 1
    limits = 0.001
    # set.seed(seed)

    data <- as.data.frame(data)
    dat <- data.frame(X = data[, X], Y = data[, Y], Mobs = data[, M])
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

    options <- c("zilonm", "zinbm", "zipm")
    m_names <- options[options %in% distM]
    all_models <- list(NULL)
    for (i in seq_along(m_names)) {
        all_models[[i]] <- tryCatch({
            realanalysis(dat, distM = m_names[i], x1, x2, zval, num_Z,
                Z_names, select_K_range, limits, B, seed)
        }, error = function(e) {
            print(e)
            list(results_effects = NA, results_parameters = NA, BIC = Inf,
                AIC = Inf, e = e)
        })
    }
    names(all_models) <- m_names
    m_BIC <- sapply(all_models, function(t) t$BIC)
    selected_distM <- m_names[which.min(m_BIC)]
    out <- all_models[[selected_distM]]
    out$selected_distM <- selected_distM
    return(out)
}




