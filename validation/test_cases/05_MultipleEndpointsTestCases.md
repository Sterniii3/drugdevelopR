#' @title 05. Multiple endpoints test cases
#' @editor Johannes Cepicka
#' @editDate 2022-08-16
#' @coverage
#' 05.01: 05.01, 05.04, 05.11, 05.15
#' 05.02: 05.01, 05.04, 05.05, 05.10
#' 05.03: 05.01, 05.03, 05.06, 05.09
#' 05.04: 05.01, 05.03, 05.07, 05.12
#' 05.05: 05.08
#' 05.06: 05.02, 05.04, 05.05, 05.10
#' 05.07: 05.02, 05.04, 05.08
#' 05.08: 05.02, 05.03, 05.06, 05.09
#' 05.09: 05.02, 05.03, 05.07, 05.12
#' 05.10: 05.11, 05.13, 05.14


##  05. Multiple endpoints test cases {-}

### 05.01 (shows that req. 05.01, 05.04, 05.11 and 05.15 are met): {-}
Use the function `optimal_multiple_tte`. Supply the following input values to the function:

  * a significance level of 0.025,
  * a power of 0.9, i.e. $\beta$ of 0.1,
  * assumed true treatment hazard ratios of 0.75 and 0.85 for endpoint 1 and endpoint 2, respectively
  * the optimization region {100, 104, …, 300} for the number of participants in phase II,
  * the optimization region {0.80, 0.82, ..., 0.9} for the threshold values,
  * boundaries of 1, 0.95 and 0.85 for the effect size categories small, medium and large,
  * expected gains of 100,000,000\$, 150,000,000\$, and 200,000,000\$ for each effect size, respectively, if only endpoint 2 shows a significant result,
  * expected gains of 100,000,000\$, 200,000,000\$, and 300,000,000\$ for each effect size, respectively, if endpoint 1 shows a significant result (independent of the significance of endpoint 2),
  * twelve clusters for parallel computing,
  * fixed costs of 10,000,000\$ in phase II and of 15,000,000\$ in phase III,
  * variable costs of 75,000\$ in phase II and 100,000\$ in phase III,
  * a correlation between the two endpoints of 0.6,
  * “fixed=FALSE”, i.e. set the function to model the treatment effects on a prior distribution,
  * amount of information for prior true treatment effect given by 210 and 420.
  
Verify that the program returns an expected utility of 457.75 (in 10^5 \$) an optimal number of participants in phase II of 296 and and optimal threshold value of 0.88. Furthermore, verify, that the program returns a probability to go to phase III of 0.75.

### 05.02 (shows that req. 05.01, 05.04, 05.05 and 05.10 are met): {-}
Use the function `optimal_multiple_tte`. Supply the same input values as in test case 05.01, however, set a sample size constraint of 500.

Verify that the program returns an expected utility of 457.75 (in 10^5 \$) an optimal number of participants of 256 in phase II and 236 in phase III (i.e. 492 in total), hence, satisfying the constraint and an optimal threshold value of 0.86. 

### 05.03 (shows that req. 05.01, 05.03, 05.06 and 05.09 are met): {-}
Use the function `optimal_multiple_tte`. Supply the same input values as in test case 05.01, however, set the parameter fixed to be "TRUE". Redo this, however, the second time use a maximum cost limit of 600 (in 10^5 \$).

Verify that the function returns an optimal number of participants of 196 in phase II and 424 in phase III (i.e. a total of 620 participants), an optimal threshold value of 0.86 and an expected utility of 161.11 (in 10^5 \$). Furthermore, verify, that the function returns costs of 247 (in 10^5 \$) in phase II and 549 (in 10^5 \$) in phase III.
With the cost constraint, the function returns an optimal number of participants of 112 in phase II and 301 in phase III (i.e. a total of 413 participants), an optimal threshold value of 0.84 and an expected utility of 137.33 (in 10^5 \$). Furthermore, verify, that the function returns costs of 184 (in 10^5 \$) in phase II and 414 (in 10^5 \$) in phase III, satisfying the constraint.

### 05.04 (shows that req. 05.01, 05.03, 05.07 and 05.12 are met): {-}
Use the function `optimal_multiple_tte`. Supply the same input values as in test case 05.01, however, set the parameter fixed to be "TRUE" and set a minimum probability of a successful program of 0.6. 

Verify that the function returns an optimal number of participants of 280 in phase II and 467 in phase III (i.e. a total of 746 participants), an optimal threshold value of 0.86 and an expected utility of 153.98 (in 10^5 \$). Furthermore, verify that the probability of a successful program is given as 0.6., satisfying the constraint and that the probability that endpoint OS is significant is 0.54.

### 05.05 (shows that req. 05.08 is met): {-}
Use the function `optimal_multiple_tte`. Supply the same input values as in test case 05.01, however change the number of cores for parallel computing to 6.

Verify that the computation time will increase compared to the setting in 05.01.

### 05.06 (shows that req. 05.02, 05.04, 05.05 and 05.10 are met): {-}
Use the function `optimal_multiple_normal()`. Supply the following input values to the function:

  * a significance level of 0.05,
  * a power of 0.9, i.e. $\beta$ of 0.1,
  * assumed true treatment effects of 0.75 and 0.8 for the endpoints 1 and 2, respectively
    * boundaries of 0, 0.5 and 0.8 for the effect size categories small, medium and large,
  * the optimization region {80, 84, …, 160} for the number of participants in phase II,
  * the optimization region {0.02, 0.04,…, 0.10} for the threshold values,
  * expected gains of 100,000,000\$, 200,000,000\$ and 300,000,000\$ for each effect size, respectively,
  * twelve clusters for parallel computing,
  * fixed costs of 10,00,000\$ in phase II and of 15,000,000\$ in phase III,
  * variable costs of 75,000\$ in phase II and 100,000\$ in phase III,
  * `fixed=FALSE`, i.e. set the function to model the treatment effects on a prior distribution,
  * a correlation of 0.5 between the two treatment effects,
  * variances of 2 and 1 of the two treatment effects, respectively
  * sample sizes of 300 and 600 as the amount of information for the two treatment effects,
  * `relaxed=TRUE`, i.e. use the relaxed combination strategy for effect sizes.

Redo this, however, the second time set a sample size constraint of 190.
  
Verify that the program returns an expected utility of 960.55 (in 10^5 \$), an optimal threshold value of 0.02 and optimal sample sizes of 108 in phase II and 85 in phase III (i.e. 193 in total). 
With the constraint, the program returns an expected utility of 959.20 (in 10^5 \$), an optimal threshold value of 0.02 and optimal sample sizes of 96 in phase II and 94 in phase III (i.e. 190 in total).

### 05.07 (shows that req. 05.02, 05.04 and 05.08 are met): {-}
Use the function `optimal_multiple_normal`. Supply the same input values as in test case 05.06 (without sample size constraint), however, however change the number of clusters for parallel computing to 6. 

Verify that the computation time will increase compared to the setting in 05.06.

### 05.08 (shows that req. 05.02, 05.03, 05.06 and 05.09 are met): {-}
Use the function `optimal_multiple_normal`. Supply the same input values as in test case 05.06 (without sample size constraint), however, set the parameter fixed to be "TRUE". Redo this, however, the second time use a maximum cost limit of 400 (in 10^5 \$).

Verify that the program returns an expected utility of 596.08 (in 10^5 \$), an optimal threshold value of 0.02, an optimal sample size in phase II of 120. Furthermore verify, that the costs are 190 (in 10^5 \$) in phase II and 217 (in 10^5 \$) in phase III, i.e. a total of 407 (in 10^5 \$).
Verify, that due to the cost constraint the program now returns an expected utility of 592.48 (in 10^5 \$), an optimal threshold value of 0.02, an optimal sample size in phase II of 104. Furthermore verify, that the costs are 178 (in 10^5 \$) in phase II and 220 (in 10^5 \$) in phase III, i.e. a total of 398 (in 10^5 \$), thus, satisfying the constraint.

### 05.09 (shows that req. 05.02, 05.03, 05.07 and 05.12 are met): {-}
Use the function `optimal_multiple_normal`. Supply the same input values as in test case 05.06 (without sample size constraint), however, set the parameter fixed to be "TRUE". Redo this and set a minimum probability of a successful program of 0.7. 

Verify that the program returns an expected utility of 596.08 (in 10^5 \$), an optimal threshold value of 0.02, an optimal sample size in phase II of 120. Furthermore, verify the probability of a successful program is 0.55, and the success probabilities for the various benefit categories are given by 0.14, 0.36 and 0.05 for small, medium and large treatment effects respectively. 
Furthermore verify, that the constraint cannot be met, within the optimization region, i.e. an expected utility of -9999 is returned, when the constraint is imposed.

### 05.10 (shows that req. 05.11, 05.13 and 05.14 are met): {-}
Use the function `optimal_multiple_normal`. Supply the same input values as in test case 05.06 (without sample size constraint), however, set parameter fixed to be "TRUE". Redo this, however the second time, set the parameter relaxed to "FALSE"

Verify that for relaxed = "TRUE" the program returns an expected utility of 596.08 (in 10^5 \$), an optimal threshold value of 0.02, an optimal sample size in phase II of 120 and a probability to go to phase III of 0.97.
For relaxed = "FALSE", the program returns an expected utility of -99.33 (in 10^5 \$), an optimal threshold value of 0.02, an optimal sample size in phase II of 96 and a probability to go to phase III of 0.96.
