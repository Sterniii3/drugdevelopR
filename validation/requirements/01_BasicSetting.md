#' @title 01. Basic setting
#' @editor Johannes Cepicka
#' @editDate 2022-08-16
#' @riskAssessment
#' 01.01: Low Risk, Very High Impact
#' 01.02: Low Risk, Very High Impact
#' 01.03: Low Risk, Very High Impact
#' 01.04: Medium Risk, Low Impact
#' 01.05: Low Risk, High Impact
#' 01.06: Medium Risk, High Impact
#' 01.07: Low Risk, Low Impact
#' 01.08: High Risk, Low Impact
#' 01.09: Low Risk, High Impact
#' 01.10: Low Risk, High Impact
#' 01.11: Low Risk, High Impact
#' 01.12: Low Risk, High Impact
#' 01.13: Low Risk, Medium Impact
#' 01.14: Low Risk, Medium Impact
#' 01.15: Low Risk, Medium Impact
#' 01.16: Low Risk, Medium Impact



## 01. Basic setting {-}

In the simplest case, the drug development program consists of a single exploratory phase II trial which is, in case of promising results, followed by one confirmatory phase III trial testing the superiority of an experimental treatment over a control treatment. The package should ask for the following user input values specific to the particular drug development program in consideration:

  *	Significance level,
  *	Assumed true treatment effects (for the intervention arm and the control arm),
  *	Event rates for phase II and III (for the time-to-event setting)
  *	Optimization set for number of participants in phase II,
  *	Optimization sets for go/no-go decision rule,
  *	Boundaries for effect size categories,
  *	Fixed and variable costs for each phase,
  *	Expected gains for each effect size,
  *	Cost constraint (maximum costs),
  *	Sample size constraint (maximum sample size),
  *	Constraint on the success probability (minimum expected probability of a successful program),
  *	A parameter to choose whether the true treatment effect is known (“fixed”) or whether the treatment effect is unknown and follows a prior distribution as specified by the following parameters:
    *	Amount of information for true treatment effects, i.e. a parameter for calculating the variance of the prior distribution: this is the sample size for the normal and the binary distribution and the number of events for the time-to-event setting
    *	Weights for prior distributions,
    *	Boundaries for truncation,
  *	Number of clusters for parallel computing.

These input values are necessary for all settings, not just for the basic setting. Input values which are only relevant in the basic setting are the following:

   *	A parameter $\gamma$ that models an offset between the treatment effects in phase II and III: the treatment effect of phase III is then calculated by adding $\gamma$ to the treatment effect of phase II
   *	A parameter to choose whether skipping phase II is an option.

Based on these input values, the program should calculate the optimal sample size and the optimal threshold value as well as the expected utility for this sample size and threshold value. As the program is required to calculate these two results in every sub-requirement within this section, we will refer to them as “the results” for the sake of simplicity. There are different possibilities how the outcome variables can be distributed. In order to fulfill the needs of different trials, the program should meet the following requirements concerning outcome variables:

  *	01.01: Calculate the results in the basic setting for normally distributed outcome variables.
  *	01.02: Calculate the results in the basic setting for binary outcome variables.
  *	01.03: Calculate the results in the basic setting for time-to-event outcome variables.

The basic setting should be adaptable to more refined use cases. One special use case may be skipping the phase II, e.g. when enough information is available to proceed directly to a confirmatory trial. In this case, no optimization for phase II sample size and threshold value should be performed, but the sample size and utility of the phase III trial should still be calculated correctly. The treatment effect for the phase III planning should then be calculated as the median of the prior distribution. This leads to the following requirement:

  *	01.04: Upon user selection, calculate the results for a setting where the phase II trial is skipped. 

Several other options for the regular setting (including both phase II and phase II should be available):
  
  *	01.05: Upon user selection, calculate the results for fixed treatment effects.
  *	01.06: Upon user selection, calculate the results for treatment effects modeled on a user-specified prior distribution depending on the distribution of the outcome variable (normal, binary, time to event).
  *	01.07: Upon user selection, calculate the results for custom population parameters (i.e. non-zero values of the offset parameter $\gamma$).
  *	01.08: Upon user selection, the program should calculate the results using parallel computing, i.e. using more than one core.
Finally, the program should be able user-specified constraints into account. In particular, the program should meet the following requirements concerning constraint violations:
  *	01.09: Automatically set the internal utility function to -9999 for all sample sizes which exceed the user-defined maximum number of participants in phase II and III. This should then lead to different results for the optimal sample size that satisfy the constraints or to a result with a utility of -9999 in case the constraints cannot be satisfied.
  *	01.10: Automatically set the internal utility function to -9999 for all sample sizes whose costs exceed the user-defined maximum cost limit. This should then lead to different results for the optimal sample size that satisfy the constraints or to a result with a utility of -9999 in case the constraints cannot be satisfied.
  *	01.11: Automatically set the internal utility function to -9999 for all sample sizes whose drug development program success probability is lower than the user-defined miminum success probability. This should then lead to different results for the optimal sample size that satisfy the constraints or to a result with a utility of -9999 in case the constraints cannot be satisfied.

In addition to the main results of optimal sample size, the optimal threshold value and the expected utility, the program should be able to return the following additional data concerning the drug development program:

  *	01.12: Return the pre-specified constraint on the total costs as well as the actually calculated costs of phase II and III for the optimal sample size and threshold.
  *	01.13: Return the total sample size as well as the separate sample sizes in phase II and III.
  *	01.14: Return the probability that the program proceeds to phase III.
  *	01.15: Return the probability of a successful program (for every effect size and in total).
  *	01.16: Return the number of events in phase II and III (in the time-to-event-setting).
