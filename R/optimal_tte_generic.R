#' @name optimal_tte_generic
#' @param w weight for mixture prior distribution, see
#'  \href{https://web.imbi.uni-heidelberg.de/prior/}{this Shiny application}
#'   for the choice of weights
#' @param hr1 first assumed true treatment effect on HR scale for \href{https://web.imbi.uni-heidelberg.de/prior/}{prior distribution}
#' @param hr2 second assumed true treatment effect on HR scale for \href{https://web.imbi.uni-heidelberg.de/prior/}{prior distribution}
#' @param id1 amount of information for hr1 in terms of number of events
#' @param id2 amount of information for hr2 in terms of number of events
#' @param d2min minimal number of events for phase II
#' @param d2max maximal number of events for phase II
#' @param stepd2 step size for the optimization over d2
#' @param hrgomin minimal threshold value for the go/no-go decision rule
#' @param hrgomax maximal threshold value for the go/no-go decision rule
#' @param stephrgo step size for the optimization over HRgo
#' @param beta type II error rate; i.e. `1 - beta` is the power for calculation of the number of events for phase III by Schoenfeld's formula (Schoenfeld 1981)
#' @param alpha significance level
#' @param xi2 assumed event rate for phase II, used for calculating the sample size of phase II via `n2 = d2/xi2`
#' @param xi3 event rate for phase III, used for calculating the sample size of phase III in analogy to `xi2`
#' @param c2 variable per-patient cost for phase II in 10^5 $.
#' @param c3 variable per-patient cost for phase III in 10^5 $.
#' @param c02 fixed cost for phase II in 10^5 $.
#' @param c03 fixed cost for phase III in 10^5 $.
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param steps1 lower boundary for effect size category "small" in HR scale, default: 1
#' @param stepm1 lower boundary for effect size category "medium" in HR scale = upper boundary for effect size category "small" in HR scale, default: 0.95
#' @param stepl1 lower boundary for effect size category "large" in HR scale = upper boundary for effect size category "medium" in HR scale, default: 0.85
#' @param b1 expected gain for effect size category "small"
#' @param b2 expected gain for effect size category "medium"
#' @param b3 expected gain for effect size category "large"
#' @param gamma to model different populations in phase II and III choose `gamma != 0`, default: 0
#' @param fixed choose if true treatment effects are fixed or random, if TRUE hr1 is used as a fixed effect and hr2 is ignored
#' @param num_cl number of clusters used for parallel computing, default: 1
optimal_tte_generic <- function(w,  hr1, hr2, id1, id2,
                                d2min, d2max, stepd2,
                                hrgomin, hrgomax, stephrgo,
                                alpha, beta, xi2, xi3,
                                c2, c3, c02, c03, 
                                K = Inf, N = Inf, S = -Inf,
                                steps1 = 1, 
                                stepm1 = 0.95, 
                                stepl1 = 0.85,
                                b1, b2, b3,
                                gamma = 0,  fixed = FALSE,
                                num_cl = 1){
  # This function is only used for documentation.
  # Hence, it contains no code.
}