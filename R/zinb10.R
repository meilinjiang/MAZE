#' Example dataset 'zinb10'
#'
#' An example dataset generated from the proposed model with a zero-inflated negative binomial mediator (K=1). The mediator contains 10% zero values in which half are false zeros.
#' 
#' @format A data frame with 100 rows and 3 variables:
#' \describe{
#'   \item{X}{independent variable, continuous data type}
#'   \item{Y}{outcome, contiuous data type}
#'   \item{Mobs}{observed mediator values with possibly false zeros, count data type}
#' }
"zinb10"
