#' @title 02. Bias adjustment test cases
#' @editor Johannes Cepicka
#' @editDate 2022-08-16
#' @coverage
#' 02.01: 02.03, 02.05, 02.10, 02.19
#' 02.02: 02.03, 02.05, 02.11, 02.17
#' 02.03: 02.06, 02.12, 02.16, 02.20
#' 02.04: 02.13, 02.20
#' 02.05: 02.04, 02.07, 02.15
#' 02.06: 02.08, 02.18
#' 02.07: 02.14
#' 02.08: 02.09
#' 02.09: 02.01, 02.04, 02.11
#' 02.10: 02.01, 02.05, 02.13
#' 02.11: 02.02, 02.04, 02.10
#' 02.12: 02.02, 02.05, 02.12

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

Verify that the function calculates an optimal sample size of 122 in phase II and 200 in phase III (i.e. a total of 322 participants), an expected utility of 78 (in 10^5\$), and an optimal threshold value of 0.78 as suggested by Stella Erdmann [2]. Furthermore, verify that one can expect 85 events in phase II and 140 events in phase III (i.e. 225 in total).

### 02.02 (shows that req. 02.03, 02.05, 02.11 and 02.17 are met): {-}
Use the function `optimal_bias`. Supply the same input values as in test case 02.01, however set the adjustment method to "mulitplicative" and set the optimization region for the multiplicative adjustment parameter $\lambda$ to {0.5, 0.55, ..., 1} and the parameters for the additive method to `NULL`. 

Verify that the function calculates an optimal sample size of 136 in phase II and 244 in phase III, i.e. a total of 380 participants), an expected utility of 99 (in 10^5\$), and an optimal threshold value of 0.76 as suggested by Stella Erdmann [2].
Furthermore, verify that the probability to go to phase III is 0.38.

### 02.03 (shows that req. 02.06, 02.12, 02.16 and 02.20 are met): {-}
Use the function `optimal_bias`. Supply the same input values as in test case 02.01, however set the adjustment method to "both" and set the optimization region for the multiplicative adjustment parameter $\lambda$ to {0.5, 0.52, ..., 1} and the parameters for the additive method $\alpha_CI$ to {0.3, 0.325, ..., 0.5}. Furthermore set a constraint for the maximum sample size to be 350.

Verify that the program returns the results for both adjustment methods by returning the selected method as well as the calculated adjustment parameter. Hereby verify, that the results for the additive method are the same as in test case 02.01 as the sample size constraint is not binding and that the optimal sample size for the multiplicative method changes to 100 in phase II and 240 in phase III, (i.e a total of 340) and the expected utility changes to 98 (in 10^5\$).

### 02.04 (shows that req. 02.13 and 02.20 are met): {-}
Use the function `optimal_bias`. Supply the same input values as in test case 02.03 (including the optimization regions for the adjustment parameters), however set the adjustment method to "all" and let the constraint for the maximal sample size remain at a value of 350.

Verify that the program returns the results for both adjustment methods by returning the selected method as well as the calculated adjustment parameter and further returns the results of an additive and a multiplicative adjustment method that not only adjust the treatment effect but also the threshold value for the decision rule. Hereby verify that the results for the basic additive and multiplicative method the same as in test case 02.03. Moreover, verify that for the advanced method, the program returns an expected overall utility of 97.00 (in 10^5\$), an adjustment parameter of 0.78 and optimal sample sizes of 122 in phase II and 204 in phase III (i.e. an overall sample size of 326) for the advanced multiplicative method ("multipl2") and an expected overall utility of 77.40 (in 10^5\$), an adjustment parameter of 0.475 and optimal sample sizes of 136 in phase II and 206 in phase III (i.e. an overall sample size of 342) for the advanced additive method ("add2").

### 02.05 (shows that req. 02.04, 02.07 and 02.15 are met): {-}
Use the function `optimal_bias`. Supply the same input values as in test case 02.01 (including the optimization regions for the additive adjustment parameter), however set the parameter fixed to be TRUE, thus using a fixed treatment effect. Redo this, however the second time set a cost constraint of 40,000,000\$.

Verify that the expected utility changes from 853.02 (in 10^5\$) to 474.18 (in 10^5\$) and the optimal sample size changes from 144 to 44 and from 456 to 172 in phase II and III respectively due to the cost constraint. The optimal adjustment parameter remains unchanged at 0.475. Furthermore verify, that the costs in phase II and III are 208 (in 10^5\$) and 582 (in 10^5\$) without the constraint and 133 (in 10^5\$) and 263 (in 10^5\$) with the constraint, i.e. the cost constraint is met at the optimal result.

### 02.06 (shows that req. 02.08 and 02.18 are met): {-}
Use the function `optimal_bias`. Supply the same input values as in test case 02.01 (including the optimization regions for the additive adjustment parameter), however set the parameter fixed to be TRUE, thus using a fixed treatment effect. Redo this, however the second time set a constraint for the minimal success probability of 0.7.

Verify that the expected utility changes from 853.02 (in 10^5\$) to 769.80 (in 10^5\$) and the optimal sample size changes from 144 to 144 and from 456 to 710 in phase II and III respectively due to the cost constraint. The optimal adjustment parameter remains changes to 0.375. Furthermore verify, that the probability of a successful program is 0.66, with a probability of 0.07, 0.21 and 0.38 for small, medium or large treatment effects without the constraint and a probability of a successful program of 0.71, with a probability of 0.06, 0.19 and 0.46 for small, medium or large treatment effects with the constraint, i.e. the cost constraint is met at the optimal result.

### 02.07 (shows that req. 02.14 is met): {-}
Use the function `optimal_bias`. Supply the same input values as in test case 02.01, however change the setting, such that the optimization region for the additive adjustment parameter just contains the point {0.5} and the adjustment parameter for the multiplicative adjustment method just contains the point {1}. Furthermore, change the adjustment method to `"both"`. Then, use the function `optimal_tte` and supply the same input values as in test case 02.01.

Verify that both adjustment methods and the case with no bias adjustment return the same results, i.e. an expected utility of 75.8 (in 10^5\$) and optimal sample sizes of 122 participants in phase II and 210 participants in phase III (i.e. a total sample size of 332).

### 02.08 (shows that req. 02.09 is met): {-}
Use the function `optimal_bias`. Supply the same input values as in test case 02.01, however change the number of chores for parallel computing to 1. 
Verify that the computation time will increase compared to the setting in 02.01.

### 02.09 (shows that req. 02.01, 02.04 and 02.11 are met): {-}
Use the function `optimal_bias_normal()`. Supply the following input values to the function:

  * a significance level of 0.05,
  * a power of 0.9, i.e. $\beta$ of 0.1,
  * assumed true treatment effects of 0.625 and 0.325,
  * the optimization region of even numbers {20, 24, …, 400} for the number of participants in phase II,
  * the optimization region {0.02, 0.04,…, 0.4} for the threshold values,
  * boundaries of 0, 0.5 and 0.8 for the effect size categories small, medium and large,
  * expected gains of 300,000,000\$, 800,000,000\$ and 1,000,000,000\$ for each effect size, respectively,
  * three clusters for parallel computing,
  * fixed costs of 1,500,000\$ in phase II and of 2,000,000\$ in phase III,
  * variable costs of 67,500\$ in phase II and 72,000\$ in phase III,
  * “fixed=FALSE”, i.e. set the function to model the treatment effects on a prior distribution,
  * weight of 0.5 for the prior distribution,
  * 300 events in phase II and 600 events in phase III,
  * truncation values of a = 0.25 and b = 0.75.
  * three chores for parallel computing

Furthermore, use the adjustment method "multiplicative" and set the optimization region for the adjustment parameter $\lambda$ to {0.7, 0.71, ..., 0.9}.



### 02.10 (shows that req. 02.01, 02.05 and 02.13 are met): {-}
Use the function `optimal_bias_normal()`. Supply the same input values as in test case 02.09 (including the optimization regions for the multiplicative adjustment parameter), however set the parameter fixed to `"TRUE"`. Furthermore use the adjustment method "all" and provide the following optimization set for the additive adjustment parameter $\alpha_CI$: {0.25, 0.275, ..., 0.5}.

Verify that the program returns the results for both adjustment methods by returning the selected method as well as the calculated adjustment parameter and further returns the results of an additive and a multiplicative adjustment method that not only adjust the treatment effect but also the threshold value for the decision rule. 
Verify that for the basic method, the program returns an expected overall utility of 4336.42 (in 10^5\$), an adjustment parameter of 0.5 and optimal sample sizes of 102 in phase II and 544 in phase III (i.e. an overall sample size of 646) for the multiplicative method and an expected overall utility of 3871.26 (in 10^5\$), an adjustment parameter of 0.1 and optimal sample sizes of 130 in phase II and 456 in phase III (i.e. an overall sample size of 586) for the additive method.
Moreover, verify that for the advanced method, the program returns an expected overall utility of 4336.37 (in 10^5\$), an adjustment parameter of 0.5 and optimal sample sizes of 98 in phase II and 536 in phase III (i.e. an overall sample size of 634) for the advanced multiplicative method ("multipl2") and an expected overall utility of 3870.71 (in 10^5\$), an adjustment parameter of 0.1 and optimal sample sizes of 134 in phase II and 426 in phase III (i.e. an overall sample size of 560) for the advanced additive method ("add2").

### 02.11 (shows that req. 02.02, 02.04 and 02.10 are met): {-}

Use the function ` optimal_bias_binary()`. Supply the following input values to the function:

  * a significance level of 0.025,
  * a power of 0.9, i.e. $\beta$ of 0.1,
  * assumed true treatment effects of p0 = 0.6, p1 = 0.3, p2= 0.5,
  * the optimization region of all even numbers {10, 12, …, 500} for the number of participants in phase II,
  * the optimization region {0.7, 0.71, …, 0.9} for the threshold values,
  * boundaries of 1, 0.95 and 0.85 for the effect size categories small, medium and large,
  * expected gains of 100,000,000, 200,000,000, and 300,000,000 for each effect size, respectively,
  * three clusters for parallel computing,
  * fixed costs of 10,000,000\$ in phase II and of 15,000,000\$ in phase III,
  * variable costs of 75,000\$ in phase II and 100,000\$ in phase III,
  * “fixed=FALSE”, i.e. set the function to use fixed treatment effects not modeled on a prior distribution,
  * weight of 0.3 for the prior distribution,
  * 30 events in phase II and 60 events in phase III.

Furthermore, use the adjustment method "additive" and set the optimization region for the adjustment parameter $\alpha_CI$ to {0.1, 0.125, ..., 0.5}.

Verify that the function calculates an optimal sample size of 158 in phase II and 262 in phase III (i.e. a total of 420 participants), an expected utility of 708.24 (in 10^5\$), and an optimal threshold value of 0.86 as well as an optimal additive adjustment parameter of 0.4.

### 02.12 (shows that req. 02.02, 02.05 and 02.12 are met): {-}

Use the function `optimal_bias_binary()`. Supply the same input values as in test case 02.11 (including the optimization region for the additive adjustment parameter), however set the parameter fixed to `"TRUE"`. Furthermore use the adjustment method "both" and provide the following optimization set for the multiplicative adjustment parameter $\lambda$: {0.5, 0.52, ..., 1}.

Verify that the program returns the results for both adjustment methods by returning the selected method as well as the calculated adjustment parameter. Hereby verify, that that for the multiplicative method the function calculates an optimal sample size of 198 in phase II and 294 in phase III (i.e. a total of 492 participants), an expected utility of 2180.86 (in 10^5\$), and an optimal threshold value of 0.82 as well as an optimal multiplicative adjustment parameter of 0.64  and that the for the additive method the function calculates an optimal sample size of 178 in phase II and 196 in phase III (i.e. a total of 374 participants), an expected utility of 2062.42 (in 10^5\$), and an optimal threshold value of 0.78 as well as an optimal additive adjustment parameter of 0.1.
