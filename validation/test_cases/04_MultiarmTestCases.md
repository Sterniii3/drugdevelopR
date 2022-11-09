#' @title 04. Multiarm test cases
#' @editor Johannes Cepicka
#' @editDate 2022-08-16
#' @coverage
#' 04.01: 04.03, 04.08, 04.12
#' 04.02: 04.03, 04.09, 04.12
#' 04.03: 04.03, 04.04, 04.14
#' 04.04: 04.07
#' 04.05: 04.02, 04.08, 04.14
#' 04.06: 04.02, 04.09, 04.14
#' 04.07: 04.02, 04.06, 04.11
#' 04.08: 04.01, 04.08, 04.13
#' 04.09: 04.01, 04.09, 04.13
#' 04.10: 04.01, 04.05, 04.10


## 04. Multiarm test cases {-}

### 04.01 (shows that req. 04.03, 04.08 and 04.12 are met): {-}
Use the function `optimal_multiarm`. Supply the following input values to the function:

  * a significance level of 0.025,
  * a power of 0.9, i.e. $\beta$ of 0.1,
  * assumed true treatment hazard ratios of 0.75 and 0.85,
  * hazard ration of the control arm of 0.6
  * the optimization region {10, 11, …, 200} for the number of participants (events in the time-to-event setting) in phase II,
  * the optimization region {0.71, 0.72, ..., 0.9} for the threshold values,
  * boundaries of 1, 0.95 and 0.85 for the effect size categories small, medium and large,
  * expected gains of 100,000,000\$, 200,000,000\$, and 300,000,000\$ for each effect size, respectively,
  * three clusters for parallel computing,
  * fixed costs of 10,000,000\$ in phase II and of 15,000,000\$ in phase III,
  * variable costs of 75,000\$ in phase II and 100,000\$ in phase III,
  * strategy 1, i.e calculating the results if only the best promising treatment proceeds to phase III.
  
Verify that the program returns an optimal utility of 8.56 (in 10^5\$), an optimal sample size of 141 in phase II and 391 in phase III (i.e a total of 532) and an optimal threshold value for the go-decision of 0.82 as suggested by Stella Erdmann [2]. Furthermore, verify that the probability to go to phase III is given by 0.72.

### 04.02 (shows that req. 04.03, 04.09 and 04.12 are met): {-}
Use the function `optimal_multiarm`. Supply the same input values as in test case 04.01, however set the parameter strategy to 2, i.e. calculating the results if all promising treatments proceed to phase III.

Verify that the program returns an optimal utility of 4.36 (in 10^5\$), an optimal sample size of 79 in phase II and 426 in phase III (i.e. a total of 505) and an optimal threshold value for the go-decision of 0.78 as suggested by Stella Erdmann [2].
Furthermore, verify that the probability to go to phase III is given by 0.65.

### 04.03 (shows that req. 04.03, 04.04 and 04.14 are met): {-}
Use the function `optimal_multiarm`. Supply the same input values as in test case 04.01, however set the parameter strategy to 3, i.e. calculating the results for both strategies and set a sample size constraint of 520.

Verify that the program returns the results for both strategies, the results for strategy 2 are the same as in test case 04.02. as the constraint is not binding, however, the results strategy 1 change as follows:  The optimal utility changes to 8.34, the optimal sample size changes to 133 in phase II and 383 in phase III (i.e. a total of 516) and the optimal threshold value for the go-decision is 0.82.

### 04.04 (shows that req. 04.07 is met): {-}
Use the function `optimal_multiarm`. Supply the same input values as in test case 04.02, however change the number of chores for parallel computing to 1. 
Verify that the computation time will increase compared to the setting in 04.02.

### 04.05 (shows that req. 04.02, 04.08 and 04.14 are met): {-}
Use the function ` optimal_multiarm_binary()`. Supply the following input values to the function:

  * a significance level of 0.025,
  * a power of 0.9, i.e. $\beta$ of 0.1,
  * assumed true treatment rate of 0.5 in the control group and assumed true rates of 0.3 and 0.4 for the treatment group, 
  * the optimization region of all even numbers {10, 12, …, 400} for the number of participants in phase II,
  * the optimization region {0.7, 0.71, …, 0.9} for the threshold values,
  * boundaries of 1, 0.95 and 0.85 for the effect size categories small, medium and large,
  * expected gains of 100,000,000\$, 200,000,000\$, and 300,000,000\$ for each effect size, respectively,
  * three clusters for parallel computing,
  * fixed costs of 10,000,000\$ in phase II and of 15,000,000\$ in phase III,
  * variable costs of 75,000\$ in phase II and 100,000\$ in phase III,
  * strategy 1, i.e calculating the results if only the best promising treatment proceeds to phase III.
  
Verify that the program returns an optimal utility of 1398.06 (in 10^5\$), an optimal sample size of 370 in phase II and 289 in phase III (i.e a total of 659) and an optimal threshold value for the go-decision of 0.87. 

### 04.06 (shows that req. 04.02, 04.09 and 04.14 are met): {-}
Use the function `optimal_multiarm_binary`. Supply the same input values as in test case 04.05, however set the parameter strategy to 2, i.e. calculating the results if all promising treatments proceed to phase III.

Verify that the program returns an optimal utility of 1414.44 (in 10^5\$), an optimal sample size of 300 in phase II and 513 in phase III (i.e. a total of 505) and an optimal threshold value for the go-decision of 0.78. 

### 04.07 (shows that req. 04.02, 04.06 and 04.11 are met): {-}
Use the function `optimal_multiarm_binary`. Supply the same input values as in test case 04.05, however set the parameter strategy to 3, and set a constraint for the minimal success probability of 0.85.

Verify that the program returns the results for both strategies, the results for strategy 2 are the same as in test case 04.06. as the constraint is not binding, however, the results strategy 1 change as follows:  The optimal utility changes to -9999, indicating that the constraint can not be fulfilled, within the optimization region.

### 04.08 (shows that req. 04.01, 04.08 and 04.13 are met): {-}
Use the function `optimal_multiarm_normal()`. Supply the following input values to the function:

  * a significance level of 0.05,
  * a power of 0.9, i.e. $\beta$ of 0.1,
  * assumed true treatment effects of 0.175 and 0.225,
  * the optimization region of even numbers {10, 12, …, 500} for the number of participants in phase II,
  * the optimization region {0.02, 0.04,…, 0.3} for the threshold values,
  * boundaries of 0, 0.5 and 0.8 for the effect size categories small, medium and large,
  * expected gains of 100,000,000\$, 300,000,000\$ and 500,000,000\$ for each effect size, respectively,
  * three clusters for parallel computing,
  * fixed costs of 1,500,000\$ in phase II and of 2,000,000\$ in phase III,
  * variable costs of 67,500\$ in phase II and 72,000\$ in phase III,
  * strategy 1, i.e calculating the results if only the best promising treatment proceeds to phase III.
  
Verify that the program returns an optimal utility of 109.9 (in 10^5\$), an optimal sample size of 56 in phase II and 205 in phase III (i.e a total of 261) and an optimal threshold value for the go-decision of 0.16. Furthermore, verify, that the probability for a succesfull program is 0.32.

### 04.09 (shows that req. 04.01, 04.09 and 04.13 are met): {-}
Use the function `optimal_multiarm_normal`. Supply the same input values as in test case 04.08, however set the parameter strategy to 2, i.e. calculating the results if all promising treatments proceed to phase III.

Verify that the program returns an optimal utility of 107.09 (in 10^5\$), an optimal sample size of 30 in phase II and 247 in phase III (i.e. a total of 277) and an optimal threshold value for the go-decision of 0.20. Furthermore, verify, that the probability for a succesfull program is 0.33.

### 04.10 (shows that req. 04.01, 04.05 and 04.10 are met): {-}

Use the function `optimal_multiarm_normal`. Supply the same input values as in test case 04.08, however set the parameter strategy to 3, i.e. calculating the results for both strategies and set a cost constraint of 200.

Verify that for strategy 1 the optimal utility changes to 109.68 (in 10^5\$) and the optimal sample size changes to 46 in phase II and 190 in phase III (i.e. a total of 236). Furthermore, verify that the costs in phase II are given by 46 and the costs in phase III are given by 151 (i.e. a total of 197).
For strategy 2, verify that the optimal utility changes to 107.06 (in 10^5\$) and the optimal sample size changes to 28 in phase II and 208 in phase III (i.e. a total of 236). Furthermore, verify that the costs in phase II are given by 34 and the costs in phase III are given by 163 (i.e. a total of 197).