#' Generic function for optimizing drug development programs with multiple endpoints
#' 
#' This function is only used for documentation generation.
#' 
#' @name optimal_multiple_generic
#' @keywords internal
#' 
#' @param rho correlation between the two endpoints
#' @param alpha one-sided significance level/family-wise error rate
#' @param n2min minimal total sample size in phase II, must be divisible by 3
#' @param n2max maximal total sample size in phase II, must be divisible by 3
#' @param stepn2 stepsize for the optimization over n2, must be divisible by 3
#' @param fixed assumed fixed treatment effect 
#' 
#' @return No return value

optimal_multiple_generic <- function(rho, alpha, n2min, n2max, stepn2, fixed){
  # This function is only used for documentation.
  # Hence, it contains no code.
}