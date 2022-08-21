#' @title Multiarm_04
#' @editor Johannes Cepicka
#' @editDate 2022-08-16
#' @riskAssessment
#' REQUIREMENT: ASSESSMENT


## 04. Multi-arm programs (“multiarm”)
A further extension implemented are multi-arm trials. So far, only three-arm trials with two experimental treatments and one control are implemented: Assume that two new treatments, e.g. two doses of the same drug, are tested at the same time and that the same control group is included in phase II and phase III. The aim of the phase II/III drug development program is to demonstrate efficacy for at least one of the experimental treatments. In this setting, the user should be able to supply the following parameters in addition to the parameters from the basic setting:

  *	A parameter to specify which strategy is used, i.e. if only the best candidate or all promising candidates proceed to phase III,
  *	The event rate of the control arm (in the time-to-event setting).

The best candidate is the treatment group with the highest treatment effect, the promising candidates are all candidates with a treatment effect greater than a pre-specified threshold $\kappa$

As above, possible cost or size constraints should be considered. However in contrast to the above settings, we only implemented fixed treatments effects. Treatment effects modeled on a prior distribution were yet not implemented, but the validation is easily adaptable if this feature is added in the future. Therefore, we expect from the program:

  *	04.01: Calculate the results in the multi-arm setting for normally distributed outcome variables.
  *	04.02: Calculate the results in the multi-arm setting for binary outcome variables.
  *	04.03: Calculate the results in the multi-arm setting for time-to-event outcome variables.
  *	04.04: Automatically set the internal utility function to -9999 for all sample sizes which exceed the user-defined maximum number of participants in phase II and III. This should then lead to different results for the optimal sample size that satisfy the constraints or to a result with a utility of -9999 in case the constraints cannot be satisfied.
  *	04.05: Automatically set the internal utility function to -9999 for all sample sizes whose costs exceed the user-defined maximum cost limit. This should then lead to different results for the optimal sample size that satisfy the constraints or to a result with a utility of -9999 in case the constraints cannot be satisfied.
  *	04.06: Automatically set the internal utility function to -9999 for all sample sizes whose drug development program success probability is lower than the user-defined miminum success probability. This should then lead to different results for the optimal sample size that satisfy the constraints or to a result with a utility of -9999 in case the constraints cannot be satisfied.
  *	04.07: Upon user selection, the program should calculate the results using parallel computing, i.e. using more than one core.

The program should be able to calculate the optimal sample size, the optimal threshold value and the expected utility for the possible strategies, i.e. only for the best candidate or all promising candidates (or both). Formulated as requirements this means:

  *	04.08: Upon user selection, calculate the results for the strategy where only the best candidates proceed to phase III.
  *	04.09: Upon user selection, calculate the results for the strategy where all promising candidates proceed to phase III.
As before, in addition to the main results of optimal sample size, the optimal threshold value and the expected utility, the program should be able to return the following additional data concerning the drug development program:
  *	04.10: Return the pre-specified constraint on the total costs as well as the actually calculated costs of phase II and III for the optimal sample size and threshold.
  *	04.11: Return the maximum total sample size as well as the separate sample sizes in phase II and III.
  *	04.12: Return the probability that the program proceeds to phase III.
  *	04.13: Return the probability of a successful program (for every effect size and in total).
  *	04.14: Return the selected strategy.
