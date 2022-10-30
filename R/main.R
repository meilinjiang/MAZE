##' Mediation Analysis for ZEro-inflated mediators
##' 
##' A novel mediation modeling approach to address zero-inflated mediators containing true zeros and false zeros. 
##' @title Mediation Analysis for ZEro-inflated mediators
##' @param data a data frame containing variables an independent variable, a mediator, and an outcome variables. See example dataset: data(zip10)
##' @param distM an optional character value for distribution to be used for the mediator. Possible choices are 'zilonm','zinbm', or 'zipm' for zero-inflated log-normal, negative binomial, or Poisson mediators respectively. By default, all three distributions will be fitted and the final mediation model is selected by AIC
##' @param X variable name of the independent variable
##' @param M variable name of the mediator variable
##' @param Y variable name of the outcome variable
##' @param x1 the first value of independent variable of interest
##' @param x2 the second value of independent variable of interest
##' @param B the upper bound value B to be used in the probability mechanism of observing false zeros
##' @param seed an optional seed number to control randomness
##' @return a list containing:  
##' - results_parameters: a data frame for the results of model parameters,  
##' - results_effects: a data frame for the results of estimated mediation effect,  
##' - BIC: a numeric value for the BIC of the final mediation model,  
##' - AIC: a numeric value for the AIC of the final mediation model,  
##' - selected_distM: a character value for the distribution of M selected in the final mediation model
##' @author Zhigang Li <zhigang.li@@ufl.edu>  
##' Meilin Jiang <meilin.jiang@@ufl.edu>
##' @import stats flexmix numDeriv 
##' @importFrom pracma hessian
##' @importFrom MASS glm.nb
##' @export
##' @examples 
##' data(zip10)
##' MAZE(data=zip10, distM=c('zilonm', 'zinbm', 'zipm'), 
##' X='X', M='Mobs', Y='Y', x1=0, x2=1, B=20, seed=1)

MAZE <- function(data, distM = c("zilonm", "zinbm", "zipm"), X, M, Y,
    x1, x2, B = 20, seed = 1) {
    select_K_range = 1
    limits = 0.001
    # set.seed(seed)

    data <- as.data.frame(data)
    dat <- data.frame(X = data[, X], Y = data[, Y], Mobs = data[, M])
    # complete data analysis
    dat <- dat[complete.cases(dat), ]

    options <- c("zilonm", "zinbm", "zipm")
    m_names <- options[options %in% distM]
    all_models <- list(NULL)
    for (i in seq_along(m_names)) {
        all_models[[i]] <- tryCatch({
            realanalysis(dat, distM = m_names[i], x1, x2, select_K_range,
                limits, B, seed)
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




