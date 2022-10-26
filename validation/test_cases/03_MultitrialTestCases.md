#' @title 03. Multitrial test cases
#' @editor Johannes Cepicka
#' @editDate 2022-08-16
#' @coverage
#' 03.01: 03.01



## 03. Multitrial test cases {-}

### 03.01 (shows that req. 03.03, 03.05, 03.11 and 03.20 are met): {-}
Use the function `optimal_multitrial`. Supply the following input values to the function:

  * a significance level of 0.025,
  * a power of 0.9, i.e. $\beta$ of 0.1,
  * assumed true treatment effects of 0.69 and 0.88,
  * event rates of 0.7 for both phase II and phase III,
  * the optimization region {10, 11, …, 400} for the number of participants (events in the time-to-event setting) in phase II,
  * the optimization region {0.71, 0.72, ..., 0.95} for the threshold values,
  * boundaries of 1, 0.95 and 0.85 for the effect size categories small, medium and large,
  * expected gains of 100,000,000\$, 300,000,000\$, and 500,000,000\$ for each effect size, respectively,
  * three clusters for parallel computing,
  * fixed costs of 10,000,000\$ in phase II and of 15,000,000\$ in phase III,
  * variable costs of 75,000\$ in phase II and 100,000\$ in phase III,
  * “fixed=FALSE”, i.e. set the function to use a prior distribution,
  * weight of 0.3 for the prior distribution,
  * amount of information for prior true treatment effect given by 210  and 420 expected events for both treatment effects.

Furthermore, use Case 2 (i.e. at least two trials have to show a significant positive treatment effect) and use "Strategy = True", hence calculating all implemented strategies for the specified case.

### 03.02 (shows that req. ): {-}

Use the function `optimal_multitrial`. Supply the same input values as in test case 03.01, however set the parameter case to 3 and the parameter strategy to 1. 

### 03.03 (shows that req. ): {-}
Use the function `optimal_multitrial`. Supply the same input values as in test case 03.01, however set a cost constraint of 50,000,000 \$. 

Verify that

### 03.04 (shows that req. ): {-}
Use the function `optimal_multitrial`. Supply the same input values as in test case 03.01, however change the parameter case to 3 and the parameter strategy to 2. Verify that the program returns an ERROR.

### 03.05 (shows that req. ): {-}

Use the function `optimal_multitrial`. Supply the same input values as in test case 03.01, however set the parameter `fixed` to be "TRUE". Furthermore, change the strategy to 23, thus, if after conducting two trials, only one delivers a significant result and the other trial’s treatment effect points at least in the same direction, a third trial is conducted.


### 03.06 (shows that req. 3.13 and 03.20 are met): {-}

Use the function `optimal_multitrial`. Supply the same input values as in test case 01.01 and set the case and the strategy to 1.

Verify, that the function returns the same results as in test case 01.01. i.e an optimal sample size of 206 in phase II and 354 in phase III (i.e. a total of 560 participants), an expected utility of 432 (in 10^5\$), and an optimal threshold value of 0.84 as suggested by Stella Erdmann [2]. Furthermore, verify that one can expect 144 events in phase II and 248 events in phase III (i.e. 392 in total).

### 03.07 (shows that req. ): {-}
Use the function ` optimal_multitrial_binary()`. Supply the following input values to the function:

  * a significance level of 0.025,
  * a power of 0.9, i.e. $\beta$ of 0.1,
  * assumed true treatment rate of 0.6 in the control group and assumed true rates of 0.3 and 0.5 for the prior distribution of the treatment group, 
  * the optimization region of all even numbers {10, 12, …, 400} for the number of participants in phase II,
  * the optimization region {0.7, 0.71, …, 0.9} for the threshold values,
  * boundaries of 1, 0.95 and 0.85 for the effect size categories small, medium and large,
  * expected gains of 100,000,000\$, 300,000,000\$, and 500,000,000\$ for each effect size, respectively,
  * three clusters for parallel computing,
  * fixed costs of 10,000,000\$ in phase II and of 15,000,000\$ in phase III,
  * variable costs of 75,000\$ in phase II and 100,000\$ in phase III,
  * “fixed=FALSE”, i.e. set the function to use fixed treatment effects not modelled on a prior distribution,
  * weight of 0.3 for the prior distribution,
  * 30 and 60 events for both treatment effects.
  
Furthermore, use Case 3 (i.e. at least three trials have to show a significant positive treatment effect) and use "Strategy = True", hence calculating all implemented strategies for the specified case.

Verify, that for strategy 1, the program returns an expected utility of 71.66 (in 10^5\$), an optimal threshold value of 0.71 and an optimal number of participants of 130 in phase II and 178 in phase III (i.e. 308 in total). Note that the value for `alpha`was automatically changed to `alpha^3`by the program. 
For strategy 3, the program returns an expected utility of 585.59 (in 10^5\$), an optimal threshold value of 0.75 and an optimal number of participants of 276 in phase II and 276 (corresponds to three trials with 92 participants) in phase III (i.e. 552 in total). 
For strategy 4, the program returns an expected utility of 718.18 (in 10^5\$), an optimal threshold value of 0.73 and an optimal number of participants of 316 in phase II and 304 (corresponds to four trials with 76 participants) in phase III (i.e. 552 in total).

### 03.08 (shows that req. ): {-}
Use the function `optimal_multitrial_binary`. Supply the same input values as in test case 03.07, however change the case to 2 and the strategy to 23, thus, if after conducting two trials, only one delivers a significant result and the other trial’s treatment effect points at least in the same direction, a third trial is conducted.

Verify that the program returns an expected utility of 810.94 (in 10^5\$) and an optimal number of participants of 220 in phase II and 184 in phase III.

### 03.09 (shows that req. ): {-} 
Use the function `optimal_multitrial_binary`. Supply the same input values as in test case 03.07, however change the parameter fixed to be `"TRUE"` and set the Case to 1 and the Strategy to `"TRUE"`, hence calculating all implemented strategies for the specified case. 

Verify, that for strategy 1, the program returns an expected utility of 1742.40 (in 10^5\$), an optimal threshold value of 0.88 and an optimal number of participants of 162 in phase II and 160 in phase III (i.e. 322 in total). For strategy 2, the program returns an expected utility of 1878.56 (in 10^5\$), an optimal threshold value of 0.82 and an optimal number of participants of 200 in phase II and 292 (corresponds to two trials with 146 participants) in phase III (i.e. 492 in total).

### 03.10 (shows that req. ): {-}
Use the function `optimal_multitrial_binary`. Supply the same input values as in test case 03.07, however change the parameter fixed to be `"TRUE"` and set the parameter Case to 2 and the parameter strategy to 3. Redo this, however, the second time set a sample size constraint of 600.

Verify, that the expected utility changes from 1332.94 to 1313.11 (in 10^5\$), and the optimal sample size changes from 242 to 186 in phase II and remains at 414 in phase III (in total it changes from 656 to 600).

### 03.11 (shows that req. ): {-}
Use the function `optimal_multitrial_normal()`. Supply the following input values to the function:

  * a significance level of 0.025,
  * a power of 0.9, i.e. $\beta$ of 0.1,
  * assumed true treatment effects of 0.375 and 0.5,
  * the optimization region of even numbers {10, 12, …, 500} for the number of participants in phase II,
  * the optimization region {0.01, 0.02,…, 0.5} for the threshold values,
  * boundaries of 0, 0.375 and 0.625 for the effect size categories small, medium and large,
  * expected gains of 62,500,000\$, 200,000,000\$ and 1,000,000,000\$ for each effect size, respectively,
  * three clusters for parallel computing,
  * fixed costs of 1,500,000\$ in phase II and of 2,000,000\$ in phase III,
  * variable costs of 67,500\$ in phase II and 72,000\$ in phase III,
  * “fixed=FALSE”, i.e. set the function to model the treatment effects on a prior distribution,
  * weight of 0.5 for the prior distribution,
  * 300 and 600 expected events for both treatment effects,
  * truncation values of a = 0 and b = 0.75.

Furthermore, use Case 3 (i.e. at least three trials have to show a significant positive treatment effect) and use "Strategy = True", hence calculating all implemented strategies for the specified case.

### 03.12 (shows that req. ): {-} 

### 03.13 (shows that req. ): {-}

### 03.14 (shows that req. ): {-}

