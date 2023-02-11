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
  * the optimization region {10, 12, …, 400} for the number of events in phase II,
  * the optimization region {0.65, 0.71, ..., 0.95} for the threshold values,
  * expected gains of 100,000,000\$, 300,000,000\$, and 500,000,000\$ for each effect size, respectively,
  * three clusters for parallel computing,
  * fixed costs of 10,000,000\$ in phase II and of 15,000,000\$ in phase III,
  * variable costs of 75,000\$ in phase II and 100,000\$ in phase III,
  * `fixed=FALSE`, i.e. set the function to use a prior distribution,
  * weight of 0.3 for the prior distribution,
  * amount of information for prior true treatment effect given by 210 and 420 expected events for the two assumed treatment effects, respectively,
  * use case 2 (i.e. at least two trials have to show a significant positive treatment effect), and
  * use `strategy = TRUE`, i.ec. calculating all implemented strategies for the specified case.

Verify that for strategy 1, the program returns an expected utility of -1.55, optimal sample sizes of 180 in phase II and 238 in phase III (i.e. 418 in total), 126 events in phase II and 167 events in phase III (i.e. 293 in total), and an optimal threshold value of 0.75.

For strategy 2, the program returns an expected utility of -94.31, optimal sample sizes of 172 in phase II and 172 in phase III (i.e. 344 in total), corresponding to two trials with 86 participants each, 120 events in phase II and 122 event in phase III (i.e. 242 in total), and an optimal threshold value of 0.72.

For strategy 3, the program returns an expected utility of -11.67, optimal sample sizes of 220 in phase II and 252 in phase III, (i.e. 472 in total) corresponding to three trials with 84 participants each, 154 events in phase II and 177 event in phase III (i.e. 331 in total), and an optimal threshold value of 0.72.

For strategy 23, the program returns an expected utility of 45.84, optimal sample sizes of 180 in phase II and 194 in phase III (i.e. 374 in total), 126 events in phase II and 136 events in phase III (i.e. 262 in total), and an optimal threshold value of 0.73. Furthermore, the probability, that a third trial is needed is given by 0.09.


### 03.02 (shows that req. 03.03, 03.10, 03.15 and 03.18 are met): {-}

Use the function `optimal_multitrial`. Supply the same input values as in test case 03.01, however, set the parameter case to 3 and the parameter strategy to 1. 

Verify that the program returns an optimal sample size of 140 in phase II and 134 in phase III (i.e. a total number of 274 participants), an expected utility of -148.05 and an optimal threshold value of 0.67. Furthermore, verify, that the probability to go to phase III is 0.21.

### 03.03 (shows that req. 03.03, 03.07 and 03.16 are met): {-}
Use the function `optimal_multitrial`. Supply the same input values as in test case 03.01, however set a cost constraint of 50,000,000 \$. 

Verify that for strategy 1 the program returns an expected utility of -3.75, optimal sample sizes of 180 in phase II and 212 in phase III (i.e. 392 in total), costs of 235 (in 10^5 \$) in phase II and 261 (in 10^5 \$) in phase III, and an optimal threshold value of 0.74.
For strategy 2, the results do not change compared to test case 03.01, as the cost constraint is not binding.
For strategy 3, the program returns an expected utility of -28.05, optimal sample sizes of 206 in phase II and 150 in phase III, (i.e. 356 in total) corresponding to three trials with 50 participants each, costs of 254 (in 10^5 \$) in phase II and 243 (in 10^5 \$) in phase III, and an optimal threshold value of 0.68.
For strategy 23, the results do not change compared to test case 03.01, as the cost constraint is not binding.

### 03.04 (shows that req. 03.12 is met): {-}
Use the function `optimal_multitrial`. Supply the same input values as in test case 03.01, however change the parameter case to 3 and the parameter strategy to 2. Verify that the program returns an ERROR.

### 03.05 (shows that req. 03.04, 03.14, 03.20 and 3.21 are met ): {-}

Use the function `optimal_multitrial`. Supply the same input values as in test case 03.01, however set the parameter `fixed` to be "TRUE". 

Verify that for strategy 1, the program returns an expected utility of 1166.41, optimal sample sizes of 420 in phase II and 1044 in phase III (i.e. 1464 in total), 294 events in phase II and 731 event in phase III (i.e. 1025 in total), and an optimal threshold value of 0.86.

For strategy 2, the program returns an expected utility of 811.67, optimal sample sizes of 386 in phase II and 1040 in phase III (i.e. 1426 in total), corresponding to two trials with 520 participants each, 270 events in phase II and 728 event in phase III (i.e. 998 in total), and an optimal threshold value of 0.85.

For strategy 3, the program returns an expected utility of 1045.41, optimal sample sizes of 420 in phase II and 1386 in phase III, (i.e. 1806 in total) corresponding to three trials with 462 participants each, 294 events in phase II and 972 event in phase III (i.e. 1266 in total), and an optimal threshold value of 0.82.

For strategy 23, the program returns an expected utility of 45.84, optimal sample sizes of 180 in phase II and 194 in phase III (i.e. 374 in total), 126 events in phase II and 136 events in phase III (i.e. 262 in total), and an optimal threshold value of 0.73. Furthermore, the probability, that a third trial ("pgo3") is needed is given by 0.09 
Note, that these are the same results as in test case 03.01, as for strategy 23, there is no parameter "fixed".


### 03.06 (shows that req. 03.13 and 03.20 are met): {-}

Use the function `optimal_multitrial`. Supply the same input values as in test case 01.01 and set the case and the strategy to 1.

Verify, that the function returns the same results as in test case 01.01. i.e an optimal sample size of 206 in phase II and 354 in phase III (i.e. a total of 560 participants), an expected utility of 432 (in 10^5\$), and an optimal threshold value of 0.84 as suggested by Stella Erdmann [2]. Furthermore, verify that one can expect 144 events in phase II and 248 events in phase III (i.e. 392 in total).

### 03.07 (shows that req. 03.02, 03.05 and 03.11 are met): {-}
Use the function ` optimal_multitrial_binary()`. Supply the following input values to the function:

  * a significance level of 0.025,
  * a power of 0.9, i.e. $\beta$ of 0.1,
  * assumed true treatment rate of 0.6 in the control group and assumed true rates of 0.3 and 0.5 for the prior distribution of the treatment group, 
  * the optimization region of all even numbers {10, 12, …, 400} for the number of participants in phase II,
  * the optimization region {0.7, 0.71, …, 0.85} for the threshold values,
  * expected gains of 100,000,000\$, 300,000,000\$, and 500,000,000\$ for each effect size, respectively,
  * three clusters for parallel computing,
  * fixed costs of 10,000,000\$ in phase II and of 15,000,000\$ in phase III,
  * variable costs of 75,000\$ in phase II and 100,000\$ in phase III,
  * `fixed=FALSE`, i.e. set the function to use treatment effects modelled on a prior distribution,
  * weight of 0.3 for the prior distribution,
  * sample sizes of 30 and 60 as the amount of information for the two assumed treatment effects,
  * use case 3 (i.e. at least three trials have to show a significant positive treatment effect), and
  * use `strategy = TRUE`, hence calculating all implemented strategies for the specified case.

Verify that for strategy 1, the program returns an expected utility of 71.66 (in 10^5\$), an optimal threshold value of 0.71 and an optimal number of participants of 130 in phase II and 178 in phase III (i.e. 308 in total). Note that the value for `alpha` was automatically changed to `alpha^3` by the program.

For strategy 3, the program returns an expected utility of 585.59 (in 10^5\$), an optimal threshold value of 0.75 and an optimal number of participants of 276 in phase II and 276 (corresponds to three trials with 92 participants) in phase III (i.e. 552 in total). 

For strategy 4, the program returns an expected utility of 718.18 (in 10^5\$), an optimal threshold value of 0.73 and an optimal number of participants of 316 in phase II and 304 (corresponds to four trials with 76 participants) in phase III (i.e. 552 in total).

### 03.08 (shows that req. 03.02 and 03.14 are met): {-}
Use the function `optimal_multitrial_binary`. Supply the same input values as in test case 03.07, however change the case to 2 and the strategy to 23, thus, if after conducting two trials, only one delivers a significant result and the other trial’s treatment effect points at least in the same direction, a third trial is conducted.

Verify that the program returns an expected utility of 810.94 (in 10^5\$) and an optimal number of participants of 220 in phase II and 184 in phase III.

### 03.09 (shows that req. 03.02, 03.04 and 03.21 are met): {-} 
Use the function `optimal_multitrial_binary`. Supply the same input values as in test case 03.07, however change the parameter fixed to be `"TRUE"` and set the case to 1 and the strategy to `"TRUE"`, hence calculating all implemented strategies for the specified case. 

Verify that for strategy 1, the program returns an expected utility of 1742.40 (in 10^5\$), an optimal threshold value of 0.88 and an optimal number of participants of 162 in phase II and 160 in phase III (i.e. 322 in total).

For strategy 2, the program returns an expected utility of 1878.56 (in 10^5\$), an optimal threshold value of 0.82 and an optimal number of participants of 200 in phase II and 292 (corresponds to two trials with 146 participants) in phase III (i.e. 492 in total).

### 03.10 (shows that req. 03.06, 03.10 and 03.17 are met): {-}
Use the function `optimal_multitrial_binary`. Supply the same input values as in test case 03.07, however change the parameter fixed to be `"TRUE"` and set the parameter case to 2 and the strategy to 3. Redo this, however, the second time set a sample size constraint of 600.

Verify that the expected utility changes from 1332.94 to 1313.11 (in 10^5\$), and the optimal sample size changes from 242 to 186 in phase II and remains at 414 in phase III (in total it changes from 656 to 600).

### 03.11 (shows that req. 03.01, 03.05, 03.11 and 03.15 are met): {-}
Use the function `optimal_multitrial_normal()`. Supply the following input values to the function:

  * a significance level of 0.05,
  * a power of 0.9, i.e. $\beta$ of 0.1,
  * assumed true treatment effects of 0.375 and 0.5,
  * the optimization region {200, 204, …, 500} for the number of participants in phase II,
  * the optimization region {0.1, 0.12,…, 0.2} for the threshold values,
  * expected gains of 300,000,000\$, 800,000,000\$ and 1,000,000,000\$ for each effect size, respectively,
  * three clusters for parallel computing,
  * fixed costs of 1,500,000\$ in phase II and of 2,000,000\$ in phase III,
  * variable costs of 67,500\$ in phase II and 72,000\$ in phase III,
  * `fixed=FALSE`, i.e. set the function to model the treatment effects on a prior distribution,
  * weight of 0.5 for the prior distribution,
  * sample sizes of 300 and 600 as the amount of information on which the
    two assumed treatment effects are based,
  * truncation values of a = 0.25 and b = 0.75,
  * use case 3 (i.e. at least three trials have to show a significant positive treatment effect), and
  * use `strategy = TRUE`, i.e. calculate all implemented strategies for the specified case.
  
Verify that for strategy 1, the program returns an expected utility of 1654.08 (in 10^5\$), an optimal threshold value of 0.16 and an optimal number of participants of 364 in phase II and 672 in phase III (i.e. 1036 in total). 

For strategy 3, the program returns an expected utility of 1308.56 (in 10^5\$), an optimal threshold value of 0.16 and an optimal number of participants of 296 in phase II and 708 in phase III (i.e. 1004 in total), corresponding to three trials with 236 participants each. Moreover, the program returns a probability of a successful program of 0.68.

For strategy 4, the program returns an expected utility of 1843.1 (in 10^5\$), an optimal threshold value of 0.18 and an optimal number of participants of 342 in phase II and 888 in phase III (i.e. 1230 in total), corresponding to four trials with 222 participants each. Moreover, the program returns a probability of a successful program of 0.86.

### 03.12 (shows that req. 03.09 is met): {-} 
Use the function `optimal_multitrial_normal()`. Supply the same input values as in test case 03.11, however change the number of clusters used for parallel computing from 3 to 1.
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
