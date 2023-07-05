#' Generic function for optimizing programs with binary endpoints
#' 
#' @name optimal_binary_generic
#'
#' @param w weight for \href{https://web.imbi.uni-heidelberg.de/prior/}{mixture prior distribution}
#' @param p0 assumed true rate of control group, see
#'  \href{https://web.imbi.uni-heidelberg.de/prior/}{here}
#'   for details
#' @param p11 assumed true rate of treatment group, see
#'  \href{https://web.imbi.uni-heidelberg.de/prior/}{here}
#'   for details
#' @param p12 assumed true rate of treatment group, see
#'  \href{https://web.imbi.uni-heidelberg.de/prior/}{here}
#'   for details
#' @param in1 amount of information for `p11` in terms of sample size, see
#'  \href{https://web.imbi.uni-heidelberg.de/prior/}{here}
#'   for details
#' @param in2 amount of information for `p12` in terms of sample size, see
#'  \href{https://web.imbi.uni-heidelberg.de/prior/}{here}
#'   for details
#' @param n2min minimal total sample size for phase II; must be an even number
#' @param n2max maximal total sample size for phase II, must be an even number
#' @param stepn2 step size for the optimization over n2; must be an even number
#' @param rrgomin minimal threshold value for the go/no-go decision rule
#' @param rrgomax maximal threshold value for the go/no-go decision rule
#' @param steprrgo step size for the optimization over RRgo
#' @param beta type II error rate; i.e. `1 - beta` is the power for calculation of the number of events for phase III
#' @param alpha significance level
#' @param c2 variable per-patient cost for phase II in 10^5 $
#' @param c3 variable per-patient cost for phase III in 10^5 $
#' @param c02 fixed cost for phase II in 10^5 $
#' @param c03 fixed cost for phase III in 10^5 $
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param steps1 lower boundary for effect size category "small" in RR scale, default: 1
#' @param stepm1 lower boundary for effect size category "medium" in RR scale = upper boundary for effect size category "small" in RR scale, default: 0.95
#' @param stepl1 lower boundary for effect size category "large" in RR scale = upper boundary for effect size category "medium" in RR scale, default: 0.85
#' @param b1 expected gain for effect size category "small"
#' @param b2 expected gain for effect size category "medium"
#' @param b3 expected gain for effect size category "large"
#' @param gamma to model different populations in phase II and III choose `gamma != 0`, default: 0, see
#'  \href{https://web.imbi.uni-heidelberg.de/prior/}{here}
#'   for details
#' @param fixed choose if true treatment effects are fixed or random, if TRUE p11 is used as fixed effect for p1
#' @param num_cl number of clusters used for parallel computing, default: 1
#' @param progressbar logical, display a progress bar or not (default: `TRUE`)
#' @keywords internal
optimal_binary_generic <- function(w, p0, p11, p12, in1, in2,
                                   n2min, n2max, stepn2,
                                   rrgomin, rrgomax, steprrgo,
                                   alpha, beta, 
                                   c2, c3, c02, c03, 
                                   K = Inf, N = Inf, S = -Inf,
                                   steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                                   b1, b2, b3,
                                   gamma = 0, fixed = FALSE,
                                   num_cl = 1,
                                   progressbar = TRUE){
  # This function is only used for documentation.
  # Hence, it contains no code.
}