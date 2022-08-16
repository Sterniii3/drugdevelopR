#' @title BiasAdjustment_02
#' @editor Johannes Cepicka
#' @editDate 2022-08-16
#' @riskAssessment
#' REQUIREMENT: ASSESSMENT

+ Start documenting requirements here!

02. Discounting phase II results (“bias”)
As the drug development programs only continue to the next stage when preceding trials are successful, estimated treatment effects may be systematically too optimistic. The program should extend the basic setting by implementing bias adjustment . As in the basic setting, the drug development program consists of a single exploratory phase II trial which is, in case of a promising result, followed by one confirmatory phase III trial. The same time-to-event, binary, or normally distributed endpoint is used in phase II and III, respectively. In addition to the general parameters specified in the section on the basic setting, the user should be able to provide the following additional parameters:
•	A parameter to choose the bias adjustment method (additive or multiplicative adjustment),
•	Parameters to determine the size of the adjustment.
As before, the program should correctly calculate the optimal sample size, the optimal threshold value and the corresponding expected utility taking the selected adjustment method and adjustment parameter as well as all other user input parameters into account. It should be adaptable to different use cases as before. Thus, we get the following requirements:
•	02.01: Calculate the results in the bias setting for normally distributed outcome variables.
•	02.02: Calculate the results in the bias setting for binary outcome variables.
•	02.03: Calculate the results in the bias setting for time-to-event outcome variables.
•	02.04: Calculate the results for fixed treatment effects.
•	02.05: Calculate the results for treatment effects modeled on a user-specified prior distribution.
•	02.06:  Automatically set the internal utility function to -9999 for all sample sizes which exceed the user-defined maximum number of participants in phase II and III. This should then lead to different results for the optimal sample size that satisfy the constraints or to a result with a utility of -9999 in case the constraints cannot be satisfied.
•	02.07: Automatically set the internal utility function to -9999 for all sample sizes whose costs exceed the user-defined maximum cost limit. This should then lead to different results for the optimal sample size that satisfy the constraints or to a result with a utility of -9999 in case the constraints cannot be satisfied.
•	02.08: Automatically set the internal utility function to -9999 for all sample sizes whose drug development program success probability is lower than the user-defined miminum success probability. This should then lead to different results for the optimal sample size that satisfy the constraints or to a result with a utility of -9999 in case the constraints cannot be satisfied.
•	02.09: Upon user selection, the program should calculate the results using parallel computing, i.e. using more than one core.

We considered two different adjustment methods to discount (possibly) too optimistic phase II results, an additive method and a multiplicative one. Both methods adjust the estimate of the observed treatment effect of phase II. The user should be able to decide which method they want to use. If the user selects additive adjustment or multiplicative adjustment, the program should correctly adjust the treatment effect in accordance with the selected method. If the user selects the option “both”, the program should return the results for the two adjustment methods separately. If the user selects the option “all”, the program should return separate results for four different adjustment methods: the two adjustment methods named above as well as an additive and a multiplicative adjustment method that not only adjusts the treatment effect but also the threshold value for the decision rule. The trivial adjustment parameters 0 or 1 should return the results from the basic setting. Therefore, we state the following requirements:
•	02.10: Upon user selection, calculate the results using the additive adjustment method (i.e. adapting the lower bound of the confidence interval).
•	02.11: Upon user selection, calculate the results using the multiplicative adjustment method (i.e. using a retention factor).
•	02.12: Return both results of 02.09 and 02.10 if both adjustment methods are selected.
•	02.13: Upon selection of the adjustment option “all”, return both results of 02.09 and 02.10 as well as the results of an additive and a multiplicative adjustment method that not only adjust the treatment effect but also the threshold value for the decision rule.
•	02.14: Return the same results as in the basic setting if the adjustment parameters are set to zero for the additive method or set to one for the multiplicative method.

As before, in addition to the main results of optimal sample size, optimal threshold value and expected utility, the program should be able to return the following additional data concerning the drug development program:
•	02.15: Return the pre-specified constraint on the total costs as well as the actually calculated costs of phase II and III for the optimal sample size and threshold.
•	02.16: Return the maximum total sample size as well as the separate sample sizes in phase II and III.
•	02.17: Return the probability that the program proceeds to phase III.
•	02.18: Return the probability of a successful program (for every effect size and in total).
•	02.19: Return the number of events in phase II and III (in the time-to-event-setting).
•	02.20: Return the adjustment method and adjustment parameter.
