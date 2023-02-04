#' @title 05. Multiple endpoints test cases
#' @editor Johannes Cepicka
#' @editDate 2022-08-16
#' @coverage
#' 05.01: 05.01, 05.04, 05.11, 05.15
#' 05.02: 05.01, 05.04, 05.05, 05.10
#' 05.03: 05.01, 05.03, 05.06, 05.09
#' 05.04: 05.01, 05.03, 05.07, 05.12
#' 05.05: 05.08
#' 05.06: 05.02, 05.04, 05.05 and 05.10
#' 05.07: 05.02, 05.04, 05.08 and 05.11
#' 05.08: 05.02, 05.03, 05.06 and 05.09
#' 05.09: 05.02, 05.03, 05.07 and 05.12
#' 05.10: 05.13, 05.14


##  05. Multiple endpoints test cases {-}

### 05.01 (shows that req. 05.01, 05.04, 05.11 and 05.15 are met): {-}
Use the function `optimal_multiple_tte`. Supply the following input values to the function:

  * a significance level of 0.025,
  * a power of 0.9, i.e. $\beta$ of 0.1,
  * assumed true treatment hazard ratios of 0.75 and 0.85 for endpoint 1 and endpoint 2, respectively
  * the optimization region {100, 104, …, 300} for the number of participants in phase II,
  * the optimization region {0.70, 0.72, ..., 0.86} for the threshold values,
  * boundaries of 1, 0.95 and 0.85 for the effect size categories small, medium and large,
  * expected gains of 100,000,000\$, 150,000,000\$, and 200,000,000\$ for each effect size, respectively, if only endpoint 2 shows a significant result,
  * expected gains of 100,000,000\$, 200,000,000\$, and 300,000,000\$ for each effect size, respectively, if endpoint 1 shows a significant result (independent of the significance of endpoint 2),
  * three clusters for parallel computing,
  * fixed costs of 10,000,000\$ in phase II and of 15,000,000\$ in phase III,
  * variable costs of 75,000\$ in phase II and 100,000\$ in phase III,
  * a correlation between the two endpoints of 0.6,
  * “fixed=FALSE”, i.e. set the function to model the treatment effects on a prior distribution,
  * amount of information for prior true treatment effect given by 210 and 420.
  
Verify that

### 05.02 (shows that req. 05.01, 05.04, 05.05 and 05.10 are met): {-}
Use the function `optimal_multiple_tte`. Supply the same input values as in test case 05.01, however, set a sample size constraint of 300.

Verify that

### 05.03 (shows that req. 05.01, 05.03, 05.06 and 05.09 are met): {-}
Use the function `optimal_multiple_tte`. Supply the same input values as in test case 05.01, however, set the parameter fixed to be "TRUE". Redo this, however, the second time use a maximum cost limit of 400 (in 10^5 \$).

Verify that the function returns an optimal number of participants of 228 in phase II and 443 in phase III (i.e a total of 671 participants), an optimal threshold value of 0.86 and an expected utility of 287.79 (in 10^5 \$). Furthermore, verify, that the function returns costs of 271 (in 10^5 \$) in phase II and 

### 05.04 (shows that req. 05.01, 05.03, 05.07 and 05.12 are met): {-}
Use the function `optimal_multiple_tte`. Supply the same input values as in test case 05.01, however, set the parameter fixed to be "TRUE" and set a minimum probability of a successful program of 0.7. 

Verify that

### 05.05 (shows that req. 05.08 is met): {-}
Use the function `optimal_multiple_tte`. Supply the same input values as in test case 05.01, however change the number of cores for parallel computing to 1.

Verify that the computation time will increase compared to the setting in 05.01.

### 05.06 (shows that req. 05.02, 05.04, 05.05 and 05.10 are met): {-}
Use the function `optimal_multiple_normal()`. Supply the following input values to the function:

  * a significance level of 0.05,
  * a power of 0.9, i.e. $\beta$ of 0.1,
  * assumed true treatment effects of 0.75 and 0.8 for the endpoints 1 and 2, respectively
  * the optimization region {20, 24, …, 200} for the number of participants in phase II,
  * the optimization region {0.02, 0.04,…, 0.20} for the threshold values,
  * expected gains of 100,000,000\$, 200,000,000\$ and 300,000,000\$ for each effect size, respectively,
  * three clusters for parallel computing,
  * fixed costs of 1,500,000\$ in phase II and of 2,000,000\$ in phase III,
  * variable costs of 67,500\$ in phase II and 72,000\$ in phase III,
  * `fixed=FALSE`, i.e. set the function to model the treatment effects on a prior distribution,
  * a correlation of 0.6 between the two treatment effects,
  * variances of 8 and 12 of the two treatment effects, respectively
  * sample sizes of 300 and 600 as the amount of information for the two treatment effects,
  * `relaxed=TRUE`, i.e. use the relaxed combination strategy for effect sizes.

Redo this, however, the second time set a sample size constraint of 300.
  
  Verify that

### 05.07 (shows that req. 05.02, 05.04, 05.08 and 05.11 are met): {-}
Use the function `optimal_multiple_normal`. Supply the same input values as in test case 05.06 (without sample size constraint), however, however change the number of clusters for parallel computing to 1. 

Verify that the computation time will increase compared to the setting in 05.06.

### 05.08 (shows that req. 05.02, 05.03, 05.06 and 05.09 are met): {-}
Use the function `optimal_multiple_normal`. Supply the same input values as in test case 05.06 (without sample size constraint), however, set the parameter fixed to be "TRUE". Redo this, however, the second time use a maximum cost limit of 400 (in 10^5 \$).

Verify that

### 05.09 (shows that req. 05.02, 05.03, 05.07 and 05.12 are met): {-}
Use the function `optimal_multiple_normal`. Supply the same input values as in test case 05.06 (without sample size constraint), however, set the parameter fixed to be "TRUE" and change the parameter relaxed to "FALSE".
Redo this, however, the second time set a minimum probability of a successful program of 0.7. 

Verify that

### 05.10 (shows that req. 05.13 and 05.14 are met): {-}
Use the function `optimal_multiple_normal`. Supply the same input values as in test case 05.06 (without sample size constraint), however, change the parameter relaxed to "FALSE".

Verify that
