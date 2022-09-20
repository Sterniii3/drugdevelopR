#' @title 02. Bias adjustment test cases
#' @editor Johannes Cepicka
#' @editDate 2022-08-16
#' @coverage
#' 02.01: 02.03, 02.05, 02.10, 02.19
#' 02.02: 02.03, 02.05, 02.11, 02.17
#' 02.03: 02.06, 02.12, 02.16, 02.20


## 02. Bias adjustment test cases {-}

### 02.01 (shows that req. 02.03, 02.05, 02.10 and 02.19 are met): {-}
Use the function `optimal_bias()`. Supply the following input values to the function:
  
  * a significance level of 0.025,
  * a power of 0.9, i.e. $\beta$ of 0.1,
  * assumed true treatment effects of 0.69 and 0.88,
  * event rates of 0.7 for both phase II and phase III,
  * the optimization region {20, 25, …, 100} for the number of participants (events in the time-to-event setting) in phase II,
  * the optimization region {0.7, 0.72, ..., 0.9} for the threshold values,
  * boundaries of 1, 0.95 and 0.85 for the effect size categories small, medium and large,
  * expected gains of 100,000,000\$, 200,000,000\$, and 300,000,000\$ for each effect size, respectively,
  * three clusters for parallel computing,
  * fixed costs of 10,000,000\$ in phase II and of 15,000,000\$ in phase III,
  * variable costs of 75,000\$ in phase II and 100,000\$ in phase III,
  * “fixed=FALSE”, i.e. set the function to use a prior distribution,
  * weight of 0.3 for the prior distribution,
  * amount of information for prior true treatment effect given by 210 expected events in phase II and 420 events in phase III.

Furthermore set the adjustment method to "additive" and set the optimization region for the additive adjustment parameter $\alpha_CI$ to {0.3, 0.325, ..., 0.5}. Set the values for the multiplicative adjustment parameter to NULL.

Verify that the function calculates an optimal sample size of 122 in phase II and 200 in phase III (i.e. a total of 322 participants), an expected utility of 78, and an optimal threshold value of 0.78 as suggested by Stella Erdmann [2]. Furthermore, verify that one can expect 85 events in phase II and 140 events in phase III (i.e. 225 in total).

### 02.02 (shows that req. 02.03, 02.05, 02.11 and 02.17 are met): {-}
Use the function `optimal_bias`. Supply the same input values as in test case 02.01, however set the adjustment method to "mulitplicative" and set the optimization region for the multiplicative adjustment parameter $\lambda$ to {0.5, 0.55, ..., 1} and the parameters for the additive method to `NULL`. 

Verify that the function calculates an optimal sample size of 136 in phase II and 244 in phase III, i.e. a total of 380 participants), an expected utility of 99, and an optimal threshold value of 0.76 as suggested by Stella Erdmann [2].
Furthermore, verify that the probability to go to phase III is 0.38.

### 02.03 (shows that req. 02.06, 02.12, 02.16 and 02.20 are met): {-}
Use the function `optimal_bias`. Supply the same input values as in test case 02.01, however set the adjustment method to "both" and set the optimization region for the multiplicative adjustment parameter $\lambda$ to {0.5, 0.55, ..., 1} and the parameters for the additive method $\alpha_CI$ to {0.3, 0.325, ..., 0.5}. Furthermore set a constraint for the maximum sample size to be 350.

Verify that the program returns the results for both adjustment methods by returning the selected method as well as the calculated adjustment parameter. Hereby verify, that the results for the additive method are the same as in test case 02.01 as the sample size constraint is not binding and that the optimal sample size for the multiplicative method changes to 100 in phase II and 240 in phase III, (i.e a total of 340) and the expected utility changes to 98.

### 02.04 (shows that req.  are met): {-}