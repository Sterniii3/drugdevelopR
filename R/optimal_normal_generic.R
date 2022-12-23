#' @name optimal_normal_generic
#' @param w weight for
#'  \href{https://web.imbi.uni-heidelberg.de/prior/}{mixture prior distribution}
#' @param Delta1 assumed true prior treatment effect measured as the
#'  standardized difference in means, see
#'   \href{https://web.imbi.uni-heidelberg.de/prior/}{here} for details
#' @param Delta2 assumed true prior treatment effect measured as the
#'  standardized difference in means, see
#'  \href{https://web.imbi.uni-heidelberg.de/prior/}{here} for details
#' @param in1 amount of information for `Delta1` in terms of sample size, see
#'  \href{https://web.imbi.uni-heidelberg.de/prior/}{here}
#'   for details
#' @param in2 amount of information for `Delta2` in terms of sample size, see
#'  \href{https://web.imbi.uni-heidelberg.de/prior/}{here}
#'   for details
#' @param a lower boundary for the truncation of the \href{https://web.imbi.uni-heidelberg.de/prior/}{prior distribution}
#' @param b upper boundary for the truncation of the \href{https://web.imbi.uni-heidelberg.de/prior/}{prior distribution}
#' @param n2min minimal total sample size for phase II; must be an even number
#' @param n2max maximal total sample size for phase II, must be an even number
#' @param stepn2 step size for the optimization over n2; must be an even number
#' @param kappamin minimal threshold value kappa for the go/no-go decision rule
#' @param kappamax maximal threshold value  kappa for the go/no-go decision rule
#' @param stepkappa step size for the optimization over the threshold value kappa
#' @param beta type II error rate; i.e. `1 - beta` is the power for calculation of the sample size for phase III
#' @param alpha significance level
#' @param c2 variable per-patient cost for phase II in 10^5 $
#' @param c3 variable per-patient cost for phase III in 10^5 $
#' @param c02 fixed cost for phase II in 10^5 $
#' @param c03 fixed cost for phase III in 10^5 $
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param steps1 lower boundary for effect size category "small", default: 0
#' @param stepm1 lower boundary for effect size category "medium" = upper boundary for effect size category "small" default: 0.5
#' @param stepl1 lower boundary for effect size category "large" = upper boundary for effect size category "medium", default: 0.8
#' @param b1 expected gain for effect size category "small" in 10^5 $
#' @param b2 expected gain for effect size category "medium" in 10^5 $
#' @param b3 expected gain for effect size category "large" in 10^5 $
#' @param gamma to model different populations in phase II and III choose `gamma != 0`, default: 0, see
#'  \href{https://web.imbi.uni-heidelberg.de/prior/}{here}
#'   for details
#' @param fixed choose if true treatment effects are fixed or following a prior distribution, if TRUE `Delta1` is used as fixed effect
#' @param num_cl number of clusters used for parallel computing, default: 1
optimal_normal_generic <- function(w, Delta1, Delta2, in1, in2, a, b,
                                   n2min, n2max, stepn2,
                                   kappamin, kappamax, stepkappa,
                                   alpha, beta, 
                                   c2, c3, c02, c03, 
                                   K = Inf, N = Inf, S = -Inf,
                                   steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                                   b1, b2, b3,
                                   gamma = 0,  fixed = FALSE, num_cl = 1){
  # This function is only used for documentation.
  # Hence, it contains no code.
}