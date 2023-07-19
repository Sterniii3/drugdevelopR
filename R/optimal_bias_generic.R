#' Generic function for optimizing drug development programs with bias adjustment
#' 
#' @name optimal_bias_generic
#' 
#' @param adj choose type of adjustment: \code{"multiplicative"},
#'  \code{"additive"}, \code{"both"} or \code{"all"}. When using "both",
#'  `res[1,]` contains the results using the multiplicative method and `res[2,]`
#'  contains the results using the additive method. When using "all", there are
#'  also `res[3,]` and `res[4,]`, containing the results of a multiplicative
#'  and an additive method which do not only adjust the treatment effect but
#'  also the threshold value for the decision rule.
#' @param lambdamin minimal multiplicative adjustment parameter lambda (i.e. use estimate with a retention factor)
#' @param lambdamax maximal multiplicative adjustment parameter lambda (i.e. use estimate with a retention factor)
#' @param steplambda stepsize for the adjustment parameter lambda
#' @param alphaCImin minimal additive adjustment parameter alphaCI (i.e. adjust the lower bound of the one-sided confidence interval)
#' @param alphaCImax maximal additive adjustment parameter alphaCI (i.e. adjust the lower bound of the one-sided confidence interval)
#' @param stepalphaCI stepsize for alphaCI
#' 
#' @keywords internal
optimal_bias_generic <- function(adj = "both",
                                 lambdamin = NULL, lambdamax = NULL, steplambda = NULL,
                                 alphaCImin = NULL, alphaCImax = NULL, stepalphaCI = NULL){
  # This function is only used for documentation.
  # Hence, it contains no code.
}
