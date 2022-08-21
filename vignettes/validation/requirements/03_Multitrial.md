#' @title Multitrial_03
#' @editor Johannes Cepicka
#' @editDate 2022-08-16
#' @riskAssessment
#' REQUIREMENT: ASSESSMENT


## 03. Programs with several phase III trials (“multitrial”)

The program should also implement a framework developed for phase II/III drug development programs where several phase III trials are performed. This is of particular relevance as regulatory agencies often require statistical significance in two or more phase III trials. Different cases, defined by the number of significant trials needed for approval, should be implemented in the package. For each case, different strategies, defined by the number of phase III trials to be conducted in order to reach the goal of the case, should be implemented. For the success of the drug development program, it is necessary that the treatment effects of all phase III trials point in the same direction. For example, if we select case 3 and strategy 4, we require four phase III trials, where three need to be significant at level $\alpha$ and the treatment effect of the fourth must point in the same direction.
Hence, in addition to the parameters from the basic setting, the user should be allowed to provide the following parameters:

  *	A parameter to set the strategy, 
  *	A parameter to set the case.
The following cases and possible strategies should be implemented in the program.

| Case | Possible strategies for this case                             |
|------|---------------------------------------------------------------|
| 1    | 1, 2                                                          |
| 2    | 1 (with significance level of $\alpha^2$)*, 2, 3, 23 (=2+1)** |
| 3    | 1 (with significance level of $\alpha^3$)*, 3, 4              |


 *For cases 2 and 3, the package should also provide the strategy to only use one trial, but with adjusted significance level.

 ** For case 2, the package should also provide a 2+1 one strategy: If after conducting two trials, only one delivers a significant result and the other trial’s treatment effect points at least in the same direction, a third trial should be conducted. This strategy will be called “23” in the package.

In analogy to the other sections, we pose the following requirements:

  *	03.01: Calculate the results in the multitrial setting for normally distributed outcome variables.
  *	03.02: Calculate the results in the multitrial setting for binary outcome variables.
  *	03.03: Calculate the results in the multitrial setting for time-to-event outcome variables.
  *	03.04: Upon user selection, calculate the results in the multitrial setting for fixed treatment effects.
  *	03.05: Upon user selection, calculate the results in the multitrial setting for treatment effects modeled on a user-specified prior distribution.
  *	03.06: Automatically set the internal utility function to -9999 for all sample sizes which exceed the user-defined maximum number of participants in phase II and III. This should then lead to different results for the optimal sample size that satisfy the constraints or to a result with a utility of -9999 in case the constraints cannot be satisfied.
  *	03.07: Automatically set the internal utility function to -9999 for all sample sizes whose costs exceed the user-defined maximum cost limit. This should then lead to different results for the optimal sample size that satisfy the constraints or to a result with a utility of -9999 in case the constraints cannot be satisfied.
  *	03.08: Automatically set the internal utility function to -9999 for all sample sizes whose drug development program success probability is lower than the user-defined miminum success probability. This should then lead to different results for the optimal sample size that satisfy the constraints or to a result with a utility of -9999 in case the constraints cannot be satisfied.
  *	03.09: Upon user selection, the program should calculate the results using parallel computing, i.e. using more than one core.

For every case and possible strategy, we expect the program to calculate the optimal sample size, the optimal threshold value and the corresponding expected utility of the drug development program correctly. Methods should be implemented according to the table above. If the user sets strategy = TRUE, the program should return every strategy implemented for the specified case. (This should be the default value.) We expect an error message if a strategy is impossible for the respective case. We formulate this as follows:

  *	03.10: Calculate the results for the user-selected number of cases using the user-selected strategy.
  *	03.11: Upon user request, calculate every implemented strategy for a specific case.
  *	03.12: Return an error message if the chosen case and strategy do not match.
  *	03.13: Return the same results as in the basic setting if the input parameters match the basic setting, i.e. case = 1 and strategy = 1.
  *	03.14: Calculate the results for the 2+1 strategy (i.e. strategy 23 in the packages nomenclature).
  *	03.15: Calculate the results for cases 2 and 3 with strategy 1 (using one trial with an adjusted significance level). These results should match the results of the basic setting with manually adjusted significance level to $\alpha^2$ or $\alpha^3$, respectively.

As before, in addition to the main results of optimal sample size, the optimal threshold value and expected utility, the program should be able to return the following additional data concerning the drug development program:

  *	03.16: Return the pre-specified constraint on the total costs as well as the actually calculated costs of phase II and III for the optimal sample size and threshold.
  *	03.17: Return the maximum total sample size as well as the separate sample sizes in phase II and III.
  *	03.18: Return the probability that the program proceeds to phase III.
  *	03.19: Return the probability of a successful program (for every effect size and in total).
  *	03.20: Return the number of events in phase II and III (in the time-to-event-setting).
  *	03.21: Return the selected strategy and case.

