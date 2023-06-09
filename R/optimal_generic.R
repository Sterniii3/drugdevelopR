#' Generic function for optimizing drug development programs
#' 
#' @name optimal_generic
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: `Inf`, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: `Inf`, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: `-Inf`, e.g. no constraint
#' @param b1 expected gain for effect size category "small"
#' @param b2 expected gain for effect size category "medium"
#' @param b3 expected gain for effect size category "large"
#' @param num_cl number of clusters used for parallel computing, default: 1
#' @return NULL
#' @keywords internal
optimal_generic <- function(beta, alpha,
                            c2, c3, c02, c03,
                            K, N, S,
                            b1, b2, b3,
                            num_cl){
  # This function is only used for documentation.
  # Hence, it contains no code.
}