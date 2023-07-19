#' Generic function for optimizing multi-arm programs
#' 
#' @name optimal_multiarm_generic
#' @param n2min minimal total sample size in phase II, must be divisible by 3
#' @param n2max maximal total sample size in phase II, must be divisible by 3
#' @param stepn2 stepsize for the optimization over n2, must be divisible by 3
#' @param beta type-II error rate for any pair, i.e. `1 - beta` is the (any-pair) power for calculation of the sample size for phase III
#' @param alpha one-sided significance level/family-wise error rate
#' @param strategy choose strategy: 1 (only the best promising candidate), 2 (all promising candidates) or 3 (both strategies)
#' @inheritParams optimal_generic
#' 
#' @keywords internal
optimal_multiarm_generic <- function(n2min, n2max, stepn2,
                                beta, alpha,
                                c2, c3, c02, c03,
                                K, N, S,
                                b1, b2, b3,
                                strategy, num_cl){
  # This function is only used for documentation.
  # Hence, it contains no code.
}