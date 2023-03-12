#' @title 03. Multitrial test cases
#' @editor Johannes Cepicka
#' @editDate 2022-08-16
#' @coverage
#' 03.01: 03.03, 03.05, 03.11, 03.20
#' 03.02: 03.03, 03.10, 03.15, 03.18
#' 03.03: 03.03, 03.07, 03.16
#' 03.04: 03.12
#' 03.05: 03.04, 03.14, 03.20, 03.21
#' 03.06: 03.13, 03.20
#' 03.07: 03.02, 03.05, 03.11
#' 03.08: 03.02, 03.14
#' 03.09: 03.02, 03.04, 03.21
#' 03.10: 03.06, 03.10, 03.17
#' 03.11: 03.01, 03.05, 03.11, 03.15
#' 03.12: 03.09
#' 03.13: 03.01, 03.08, 03.19
#' 03.14: 03.01, 03.04, 03.10



## 03. Multitrial test cases {-}

### 03.01 (shows that req. 03.03, 03.05, 03.11 and 03.20 are met): {-}
Use the function `optimal_multitrial`. Supply the following input values to the function:

  * a significance level of 0.025,
  * a power of 0.9, i.e. $\beta$ of 0.1,
  * assumed true treatment effects of 0.69 and 0.88,
  * event rates of 0.7 for both phase II and phase III,
  * the optimization region {100, 104, …, 300} for the number of events in phase II,
  * the optimization region {0.65, 0.71, ..., 0.8} for the threshold values,
  * expected gains of 100,000,000\$, 300,000,000\$, and 500,000,000\$ for each effect size, respectively,
  * twelve clusters for parallel computing,
  * fixed costs of 10,000,000\$ in phase II and of 15,000,000\$ in phase III,
  * variable costs of 75,000\$ in phase II and 100,000\$ in phase III,
  * `fixed=FALSE`, i.e. set the function to use a prior distribution,
  * weight of 0.3 for the prior distribution,
  * amount of information for prior true treatment effect given by 210 and 420 expected events for the two assumed treatment effects, respectively,
  * use case 2 (i.e. at least two trials have to show a significant positive treatment effect), and
  * use `strategy = TRUE`, i.e. calculating all implemented strategies for the specified case.

Verify that for strategy 1, the program returns an expected utility of -1.89, optimal sample sizes of 200 in phase II and 238 in phase III (i.e. 438 in total), 140 events in phase II and 167 events in phase III (i.e. 307 in total), and an optimal threshold value of 0.75.

For strategy 2, the program returns an expected utility of -94.31, optimal sample sizes of 172 in phase II and 172 in phase III (i.e. 344 in total), corresponding to two trials with 86 participants each, 120 events in phase II and 122 event in phase III (i.e. 242 in total), and an optimal threshold value of 0.72.

For strategy 3, the program returns an expected utility of -12.88, optimal sample sizes of 224 in phase II and 252 in phase III, (i.e. 476 in total) corresponding to three trials with 84 participants each, 156 events in phase II and 177 event in phase III (i.e. 333 in total), and an optimal threshold value of 0.72.

For strategy 23, the program returns an expected utility of 45.72, optimal sample sizes of 178 in phase II and 194 in phase III (i.e. 372 in total), 124 events in phase II and 136 events in phase III (i.e. 260 in total), and an optimal threshold value of 0.73. Furthermore, the probability, that a third trial is needed is given by 0.09.


### 03.02 (shows that req. 03.03, 03.10, 03.15 and 03.18 are met): {-}

Use the function `optimal_multitrial`. Supply the same input values as in test case 03.01, however, set the parameter case to 3 and the parameter strategy to 1. 

Verify that the program returns an optimal sample size of 160 in phase II and 130 in phase III (i.e. a total number of 290 participants), an expected utility of -148.57 and an optimal threshold value of 0.67. Furthermore, verify, that the probability to go to phase III is 0.2.

### 03.03 (shows that req. 03.03, 03.07 and 03.16 are met): {-}
Use the function `optimal_multitrial`. Supply the same input values as in test case 03.01, however set a cost constraint of 50,000,000 \$. 

Verify that for strategy 1 the program returns an expected utility of -4.21, optimal sample sizes of 172 in phase II and 212 in phase III (i.e. 384 in total), costs of 229 (in 10^5 \$) in phase II and 261 (in 10^5 \$) in phase III, and an optimal threshold value of 0.74.
For strategy 2, the results do not change compared to test case 03.01, as the cost constraint is not binding.
For strategy 3, the program returns an expected utility of -27.31, optimal sample sizes of 172 in phase II and 156 in phase III, (i.e. 328 in total) corresponding to three trials with 52 participants each, costs of 229 (in 10^5 \$) in phase II and 254 (in 10^5 \$) in phase III, and an optimal threshold value of 0.68.
For strategy 23, the results do not change compared to test case 03.01, as the cost constraint is not binding.

### 03.04 (shows that req. 03.12 is met): {-}
Use the function `optimal_multitrial`. Supply the same input values as in test case 03.01, however change the parameter case to 3 and the parameter strategy to 2. Verify that the program returns an ERROR.

### 03.05 (shows that req. 03.04, 03.14, 03.20 and 3.21 are met ): {-}

Use the function `optimal_multitrial`. Supply the same input values as in test case 03.01, however set the parameter `fixed` to be "TRUE" and change the optimization region to {200, 204, …, 400} for the number of events in phase II and to {0.8, 0.81, ..., 0.88} for the threshold values

Verify that for strategy 1, the program returns an expected utility of 1165.47, optimal sample sizes of 424 in phase II and 1044 in phase III (i.e. 1468 in total), 296 events in phase II and 731 event in phase III (i.e. 1027 in total), and an optimal threshold value of 0.86.

For strategy 2, the program returns an expected utility of 810.67, optimal sample sizes of 384 in phase II and 1040 in phase III (i.e. 1424 in total), corresponding to two trials with 520 participants each, 268 events in phase II and 728 event in phase III (i.e. 996 in total), and an optimal threshold value of 0.85.

For strategy 3, the program returns an expected utility of 1045.28, optimal sample sizes of 446 in phase II and 1398 in phase III, (i.e. 1844 in total) corresponding to three trials with 466 participants each, 312 events in phase II and 978 event in phase III (i.e. 1290 in total), and an optimal threshold value of 0.82.

For strategy 23, the program returns an expected utility of -47.41, optimal sample sizes of 286 in phase II and 454 in phase III, (i.e. 740 in total), 200 events in phase II and 318 event in phase III (i.e. 518 in total), and an optimal threshold value of 0.80.

### 03.06 (shows that req. 03.13 and 03.20 are met): {-}

Use the function `optimal_multitrial`. Supply the same input values as in test case 01.01 and set the case and the strategy to 1.

Verify, that the function returns the same results as in test case 01.01. i.e an optimal sample size of 206 in phase II and 354 in phase III (i.e. a total of 560 participants), an expected utility of 432 (in 10^5\$), and an optimal threshold value of 0.84 as suggested by Stella Erdmann [2]. Furthermore, verify that one can expect 144 events in phase II and 248 events in phase III (i.e. 392 in total).

### 03.07 (shows that req. 03.02, 03.05 and 03.11 are met): {-}
Use the function ` optimal_multitrial_binary()`. Supply the following input values to the function:

  * a significance level of 0.025,
  * a power of 0.9, i.e. $\beta$ of 0.1,
  * assumed true treatment rate of 0.6 in the control group and assumed true rates of 0.3 and 0.5 for the prior distribution of the treatment group, 
  * the optimization region of all even numbers {100, 104, …, 400} for the number of participants in phase II,
  * the optimization region {0.7, 0.71, …, 0.8} for the threshold values,
  * expected gains of 100,000,000\$, 300,000,000\$, and 500,000,000\$ for each effect size, respectively,
  * twelve clusters for parallel computing,
  * fixed costs of 10,000,000\$ in phase II and of 15,000,000\$ in phase III,
  * variable costs of 75,000\$ in phase II and 100,000\$ in phase III,
  * `fixed=FALSE`, i.e. set the function to use treatment effects modelled on a prior distribution,
  * weight of 0.3 for the prior distribution,
  * sample sizes of 30 and 60 as the amount of information for the two assumed treatment effects,
  * use case 3 (i.e. at least three trials have to show a significant positive treatment effect), and
  * use `strategy = TRUE`, hence calculating all implemented strategies for the specified case.

Verify that for strategy 1, the program returns an expected utility of 405.56 (in 10^5\$), an optimal threshold value of 0.77 and an optimal number of participants of 176 in phase II and 306 in phase III (i.e. 482 in total). Note that the value for `alpha` was automatically changed to `alpha^3` by the program.

For strategy 3, the program returns an expected utility of 1494.31 (in 10^5\$), an optimal threshold value of 0.77 and an optimal number of participants of 340 in phase II and 336 (corresponds to three trials with 112 participants) in phase III (i.e. 676 in total). 

For strategy 4, the program returns an expected utility of 1736.36 (in 10^5\$), an optimal threshold value of 0.76 and an optimal number of participants of 384 in phase II and 408 (corresponds to four trials with 102 participants) in phase III (i.e. 792 in total).

### 03.08 (shows that req. 03.02 and 03.14 are met): {-}
Use the function `optimal_multitrial_binary`. Supply the same input values as in test case 03.07, however change the case to 2 and the strategy to 23, thus, if after conducting two trials, only one delivers a significant result and the other trial’s treatment effect points at least in the same direction, a third trial is conducted.

Verify that the program returns an expected utility of 1701.93 (in 10^5\$) and an optimal number of participants of 262 in phase II and 224 in phase III.

### 03.09 (shows that req. 03.02, 03.04 and 03.21 are met): {-} 
Use the function `optimal_multitrial_binary`. Supply the same input values as in test case 03.07, however change the parameter fixed to be `"TRUE"` and the optimization region for the threshold values to {0.8, 0.81, …, 0.95} and set the case to 1 and the strategy to `"TRUE"`, hence calculating all implemented strategies for the specified case. 

Verify that for strategy 1, the program returns an expected utility of 3187.57 (in 10^5\$), an optimal threshold value of 0.91 and an optimal number of participants of 176 in phase II and 158 in phase III (i.e. 334 in total).

For strategy 2, the program returns an expected utility of 3657.70 (in 10^5\$), an optimal threshold value of 0.89 and an optimal number of participants of 240 in phase II and 284 (corresponds to two trials with 142 participants) in phase III (i.e. 524 in total).

### 03.10 (shows that req. 03.06, 03.10 and 03.17 are met): {-}
Use the function `optimal_multitrial_binary`. Supply the same input values as in test case 03.07, however change the parameter fixed to be `"TRUE"` and the optimization region for the threshold values to {0.8, 0.81, …, 0.95} and set the parameter case to 2 and the strategy to 3. Redo this, however, the second time set a sample size constraint of 600.

Verify that the expected utility changes from 2923.29 to 1795.38 (in 10^5\$), and the optimal sample size changes from 288 to 132 in phase II and from 408 to 468 in phase III (in total it changes from 656 to 600).

### 03.11 (shows that req. 03.01, 03.05, 03.11 and 03.15 are met): {-}
Use the function `optimal_multitrial_normal()`. Supply the following input values to the function:

  * a significance level of 0.05,
  * a power of 0.9, i.e. $\beta$ of 0.1,
  * assumed true treatment effects of 0.375 and 0.5,
  * the optimization region {200, 204, …, 500} for the number of participants in phase II,
  * the optimization region {0.1, 0.12,…, 0.2} for the threshold values,
  * expected gains of 300,000,000\$, 800,000,000\$ and 1,000,000,000\$ for each effect size, respectively,
  * twelve clusters for parallel computing,
  * fixed costs of 1,500,000\$ in phase II and of 2,000,000\$ in phase III,
  * variable costs of 67,500\$ in phase II and 72,000\$ in phase III,
  * `fixed=FALSE`, i.e. set the function to model the treatment effects on a prior distribution,
  * weight of 0.5 for the prior distribution,
  * sample sizes of 300 and 600 as the amount of information on which the
    two assumed treatment effects are based,
  * truncation values of a = 0.25 and b = 0.75,
  * use case 3 (i.e. at least three trials have to show a significant positive treatment effect), and
  * use `strategy = TRUE`, i.e. calculate all implemented strategies for the specified case.
  
Verify that for strategy 1, the program returns an expected utility of 1660.95 (in 10^5\$), an optimal threshold value of 0.16 and an optimal number of participants of 388 in phase II and 666 in phase III (i.e. 1054 in total). 
Moreover, the program returns a probability of a successful program of 0.81.

For strategy 3, the program returns an expected utility of 1282.14 (in 10^5\$), an optimal threshold value of 0.16 and an optimal number of participants of 332 in phase II and 702 in phase III (i.e. 1034 in total), corresponding to three trials with 234 participants each. Moreover, the program returns a probability of a successful program of 0.68.

For strategy 4, the program returns an expected utility of 1786.02 (in 10^5\$), an optimal threshold value of 0.18 and an optimal number of participants of 356 in phase II and 888 in phase III (i.e. 1244 in total), corresponding to four trials with 222 participants each. Moreover, the program returns a probability of a successful program of 0.85.

### 03.12 (shows that req. 03.09 is met): {-} 
Use the function `optimal_multitrial_normal()`. Supply the same input values as in test case 03.11, however change the number of clusters used for parallel computing from 12 to 6.
Verify that the computation time will increase compared to the setting in 03.11.

### 03.13 (shows that req. 03.01, 03.08 and 03.19 are met): {-}
Use the function `optimal_multitrial_normal()`. Supply the same input values as in test case 03.11, however change the parameter fixed to be `"TRUE"`. Redo this, however, the second time set a minimum success probability of 0.82.

Verify that for strategy 1 without the constraint, the program returns an expected utility of 1514.35 (in 10^5\$), an optimal threshold value of 0.16 and an optimal number of participants of 440 in phase II and 830 in phase III (i.e. 1270 in total). Moreover, the program returns a probability of a successful program of 0.81.

For strategy 3 without the constraint, the program returns an expected utility of 1116.60 (in 10^5\$), an optimal threshold value of 0.16 and an optimal number of participants of 364 in phase II and 882 in phase III (i.e. 1246 in total), corresponding to three trials with 294 participants each. Moreover, the program returns a probability of a successful program of 0.68.

For strategy 4 without the constraint, the program returns an expected utility of 1395.35 (in 10^5\$), an optimal threshold value of 0.18 and an optimal number of participants of 424 in phase II and 1128 in phase III (i.e. 1270 in total), corresponding to four trials with 282 participants each. Moreover, the program returns a probability of a successful program of 0.86.

For strategy 1 under the constraint, the expected utility changes to 1513.26 (in 10^5\$), the  optimal threshold value to 0.14 and an optimal number of participants of 444 in phase II and 852 in phase III (i.e. 1296 in total). Moreover, the program now returns a probability of a successful program of 0.82.

For strategy 3 under the constraint, the constraint can not be met withing the optimization region, so the program returns an expected utility of -9999.

For strategy 4 under the constraint, the program returns the same results as before as the constraint is fulfilled at the optimum for this strategy.

### 03.14 (shows that req. 03.01, 03.04 and 03.10 are met): {-}
Use the function `optimal_multitrial_normal()`. Supply the same input values as in test case 03.11, however change the parameter fixed to be `"TRUE"` and use case 2 (i.e. at least three trials have to show a significant positive treatment effect) and the strategy 3 (i.e three trials are conducted in phase III).

Verify that the program returns an expected utility of 1749.97 (in 10^5\$), an optimal threshold value of 0.16 and an optimal number of participants of 416 in phase II and 876 in phase III (i.e. 1292 in total), corresponding to three trials with 292 participants each. 
