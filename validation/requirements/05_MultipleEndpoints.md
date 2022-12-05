#' @title 05. Multiple endpoints
#' @editor Johannes Cepicka
#' @editDate 2022-08-16
#' @riskAssessment
#' 05.01: High Risk, High Impact
#' 05.02: High Risk, High Impact
#' 05.03: High Risk, High Impact
#' 05.04: Very High Risk, High Impact
#' 05.05: Low Risk, High Impact
#' 05.06: Low Risk, High Impact
#' 05.07: Low Risk, High Impact
#' 05.08: High Risk, Low Impact
#' 05.09: Low Risk, Medium Impact
#' 05.10: Low Risk, Medium Impact
#' 05.11: Low Risk, Medium Impact
#' 05.12: Low Risk, Medium Impact
#' 05.13: High Risk, High Impact
#' 05.14: High Risk, High Impact
#' 05.15: High Risk, High Impact


## 05. Multiple endpoints {-}

The program should also provide methods for drug development programs with multiple endpoints. For now, this means that the program provides methods for two endpoints. Moreover, only normally distributed and time-to-event endpoints are implemented in the multiple endpoint setting. (Further extensions may be implemented in the future.) The definition of treatment success is different for the two endpoints:

 *	In the time-to-event setting, the drug development program is defined to be successful if it proceeds from phase II to phase III and at least one endpoint shows a statistically significant treatment effect in phase III. For example, this situation is found in oncology trials, where overall survival (OS) and progression free survival (PFS) are the two endpoints of interest.

 *	For normally distributed endpoints, the drug development program is defined to be successful if it proceeds from phase II to phase III and all endpoints show a statistically significant treatment effect in phase III. For example, this situation is found in Alzheimer’s disease trials, where a drug should show significant results in improving cognition (cognitive endpoint) as well as in improving activities of daily living (functional endpoint).

The user should be able to provide the following input values in addition to the general parameters defined in the basic setting:

  *	The correlation between the two endpoints,
  * The event rate of the control arm (in the time to event setting), and
  * The variances of the endpoints (in the normally distributed setting).

The program should correctly calculate the optimal sample size, the optimal threshold value and the corresponding expected utility for utility-based optimization of phase II/III programs with two time-to-event endpoints. We require the following:

  *	05.01: Calculate the results in the multiple endpoints setting for time-to-event outcome variables.
  *	05.02: Calculate the results in the multiple endpoints setting for normally distributed outcome variables.
  *	05.03: Upon user selection, calculate the results for fixed effects.
  *	05.04: Upon user selection, calculate the results for treatment effects modeled on a user-specified prior distribution.
  *	05.05: Automatically set the internal utility function to -9999 for all sample sizes which exceed the user-defined maximum number of participants in phase II and III. This should then lead to different results for the optimal sample size that satisfy the constraints or to a result with a utility of -9999 in case the constraints cannot be satisfied.
  *	05.06: Automatically set the internal utility function to -9999 for all sample sizes whose costs exceed the user-defined maximum cost limit. This should then lead to different results for the optimal sample size that satisfy the constraints or to a result with a utility of -9999 in case the constraints cannot be satisfied.
  *	05.07: Automatically set the internal utility function to -9999 for all sample sizes whose drug development program success probability is lower than the user-defined miminum success probability. This should then lead to different results for the optimal sample size that satisfy the constraints or to a result with a utility of -9999 in case the constraints cannot be satisfied.
  *	05.08: Upon user selection, the program should calculate the results using parallel computing, i.e. using more than one core. 


As before, in addition to the main results of optimal sample size, optimal threshold value and expected utility, the program should be able to return the following additional data concerning the drug development program:

  *	05.09: Return the pre-specified constraint on the total costs as well as the actually calculated costs of phase II and III for the optimal sample size and threshold.
  *	05.10: Return the maximum total sample size as well as the separate sample sizes in phase II and III.
  *	05.11: Return the probability that the program proceeds to phase III.
  *	05.12: Return the probability of a successful program (for every effect size and in total).

The effect size categories small, medium and large can be applied to both endpoints. In order to define an overall effect size from the two individual effect sizes, the package should implement combination rules. For normally distributed endpoints, two different combination rules should be implemented:

 *	A strict rule assigning a large overall effect in case both endpoints show an effect of large size, a small overall effect in case that at least one of the endpoints shows a small effect, and a medium overall effect otherwise.

 *	A relaxed rule assigning a large overall effect if at least one of the endpoints shows a large effect, a small effect if both endpoints show a small effect, and a medium overall effect otherwise.  

Based on this we require the following features:

  * 05.13: Upon user selection, calculate the results for above-explained strict rule regarding treatment effects. This feature is only required for normally distributed endpoints.
  *	05.14: Upon user selection, calculate the results for above-explained relaxed rule regarding treatment effects. This feature is only required for normally distributed endpoints.

On the other hand, for time-to-event endpoints, the effect size of the endpoint with larger treatment effect should be selected as overall effect size. In addition, the user should be asked to provide two triples of benefits per category. If only the less important endpoint is significant after phase III, then the smaller benefit triple should be chosen by the software. If at least the more important endpoint is significant, then the larger benefit triple should be chosen by the software. (This combination rule reflects the situation in oncological trials: If only progression-free survival is significant, then a smaller benefit can be expected compared to trials were over all survival is significant.) Based on this we require the following features:

  *	05.15: Request two benefit triples from the user. Calculate the results using above-explained combination rule regarding treatment effects and the two user-specified different benefit triples. This feature is only required for time-to-event endpoints.
