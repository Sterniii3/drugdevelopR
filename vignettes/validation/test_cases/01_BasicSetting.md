#' @title BasicSetting_01
#' @editor Johannes Cepicka
#' @editDate 2022-08-16
#' @coverage
#' TESTCASE: REQUIREMENT


+ Setup: DOCUMENT ANY SETUP THAT NEEDS TO BE DONE FOR TESTING

01. Test cases for the basic setting

01.01 (shows that req. 01.03, 01.06, 01.12 and 01.15 are met): Use the function `optimal_tte()`. Supply the following input values to the function:
  
  * a significance level of 0.025,
 
  * a power of 0.9, i.e. $\beta$ of 0.1,
 
  * assumed true treatment effects of 0.8 and 0.65,
  
  * event rates of 0.7 for both phase II and phase III,
 
  * the optimization region {10, 11, …, 400} for the number of participants in phase II,
 
  * the optimization region {-log(0.95), -log(0.94),…, -log(0.70)} for the threshold values,
 
  * boundaries of 1, 0.95 and 0.85 for the effect size categories small, medium and large,
 
  * expected gains of 100,000,000$, 300,000,000$, and 500,000,000$ for each effect size, respectively,
 
  * three clusters for parallel computing,
 
  * fixed costs of 10,000,000$ in phase II and of 15,000,000$ in phase III,
 
  * variable costs of 75,000$ in phase II and 100,000$ in phase III,
 
  * “fixed=FALSE”, i.e. set the function to use a prior distribution,
 
  * weight of 0.3 for the prior distribution,

  * amount of information for prior true treatment effect given by 210 expected events in phase II and 420 events in phase III.

Verify that the function calculates an optimal sample size of 206 in phase II and 354 in phase III (i.e. a total of 560 participants), an expected utility of 432, and an optimal threshold value of 0.84 as suggested by Stella Erdmann [2]. Furthermore, verify that one can expect 144 events in phase II and 248 events in phase III (i.e. 392 in total).

01.02 (shows that req. 01.03., 01.05., 01.14 and 01.16 are met): Use the function `optimal_tte()`. Supply the same input values as in test case 01.01 to the function except for the following: Set the parameter “fixed” to be TRUE, thus using fixed assumed treatment effects. Set the weight for the prior distribution to be NULL and the number of events to be NULL and NULL. 
Verify that the function calculates an optimal sample size of 240, an expected utility of 352 and an optimal threshold value of 0.88 as suggested by Stella Erdmann [2]. Furthermore, verify that the probability to go to phase III is given by 0.73 and the expected number of events in phase III and III is 168 and 546, respectively (714 in total).

01.03 (shows that req. 01.10 and 01.12 are met): Use the function `optimal_tte()`. Supply the same input values as in test case 01.01 to the function except for the following changes and additions: Set the weight for the prior distribution to be 0.6 and set a cost constraint of K=750.
Verify that the function calculates an optimal sample size of 228, an expected utility of 996 and an optimal threshold value of 0.84 as suggested by Stella Erdmann [2]. Furthermore, verify that the cost constraint is returned and that the total costs in phase II and III are 271 and 530.

01.04 (shows that req. 01.09 and 01.13 are met): Use the function `optimal_tte()`. Supply the same input values as in test case 01.01 to the function except for the following changes and additions: Set the weight for the prior distribution to be 0.6  and set a sample size constraint of N=500.
Verify that the function calculates an optimal sample size of 170 in phase II and 328 in phase III (i.e. a total sample size of 498), an expected utility of 956 and an optimal threshold value of 0.83 as suggested by Stella Erdmann [2].

01.05 (shows that req. 01.11, 01.14 and 01.15 met): Use the function `optimal_tte()`. Supply the same input values as in test case 01.01 to the function except for the following changes and additions: Set the weight for the prior distribution to be 0.6 and set a constraint on the minimal probability of a successful program of S=0.6.
Verify that the function calculates an optimal sample size of 470, an expected utility of 899 and an optimal threshold value of 0.89 as suggested by Stella Erdmann [2]. Furthermore, verify, that the probability to go phase III is given by 0.77 and the probability of a successful program is given by 0.6.

01.06 (shows that req. 01.04 is met): Use the function `optimal_tte()`. Supply the same input values as in test case 01.01 to the function except for the following changes and additions: Set the weight for the prior distribution to be 0.6 and use the option to skip phase II.
Verify that the function calculates an optimal sample size of 824, an expected utility of 1706 and an optimal threshold value of 0.76 as suggested by Stella Erdmann [2].

01.07 (shows that req. 01.07 is met): Use the function `optimal_tte()`. Supply the same input values as in test case 01.01 to the function except for the following changes and additions: Set the weight for the prior distribution to be 0.6 and use the option to model different population structures and set the parameter $\gamma$ to 0.025.
Verify that the function calculates an optimal sample size of 310, an expected utility of 1207 and an optimal threshold value of 0.86 as suggested by Stella Erdmann [2].

01.08 (shows that req. 01.02, 01.05 and 01.12 are met): Use the function ` optimal_binary()`. Supply the following input values to the function:
* a significance level of 0.025,
* a power of 0.9, i.e. $\beta$ of 0.1,
* assumed true treatment effects of p0=0.6, p1=0.5, p2= 0.3,
* the optimization region of all even numbers {10, 12, …, 500} for the number of participants in phase II,
* the optimization region {-log(0.95), -log(0.94),…, -log(0.70)} for the threshold values,
* boundaries of 1, 0.95 and 0.85 for the effect size categories small, medium and large,
* expected gains of 1000, 3000, and 5000 for each effect size, respectively,
* three clusters for parallel computing,
* fixed costs of 10,000,000$ in phase II and of 15,000,000$ in phase III,
* variable costs of 75,000$ in phase II and 100,000$ in phase III,
* “fixed=TRUE”, i.e. set the function to use fixed treatment effects not modelled on a prior distribution,
* weight of NULL for the prior distribution,
* NULL events in phase II and NULL events in phase III.
Verify that the function calculates an optimal sample size of 204, an expected utility of 299 and an optimal threshold value of 0.90 as suggested by Stella Erdmann [2]. Furthermore, verify that the programs returns the cost constraint (Inf in this case) as well as the total costs in phase II and III, 253 and 810, respectively.

01.09 (shows that req. 01.02 and 01.06 are met): Use the function `optimal_binary()`. Supply the same input values as in test case 01.08 to the function except for the following changes and additions: Set the parameter “fixed” to be FALSE. Set the weight for the prior distribution to be 0.4 and the number of events to be 30 and 60. 
Verify that the function calculates an optimal sample size of 224, an expected utility of 1542 and an optimal threshold value of 0.89 as suggested by Stella Erdmann [2].

01.10 (shows that req. 01.01 and 01.06 are met) Use the function `optimal_normal()`. Supply the following input values to the function:
* a significance level of 0.025,
* a power of 0.9, i.e. $\beta$ of 0.1,
* assumed true treatment effects of 0.625 and 0.325,
* the optimization region of even numbers {10, 12, …, 500} for the number of participants in phase II,
* the optimization region {0.01, 0.02,…, 0.5} for the threshold values,
* boundaries of 0, 0.5 and 0.8 for the effect size categories small, medium and large,
* expected gains of 62,500,000$, 200,000,000$ and 1,000,000,000$ for each effect size, respectively,
* three clusters for parallel computing,
* fixed costs of 1,500,000$ in phase II and of 2,000,000$ in phase III,
* variable costs of 67,500$ in phase II and 72,000$ in phase III,
* “fixed=FALSE”, i.e. set the function to model the treatment effects on a prior distribution,
* weight of 0.5 for the prior distribution,
* 300 events in phase II and 600 events in phase III,
* truncation values of a = 0 and b = 0.75.
Verify that the function calculates an optimal sample size of 86, an expected utility of 337 and an optimal threshold value of 0.19 as suggested by Stella Erdmann [2].

01.11 (shows that req. 01.01, 01.05 and 01.15 are met) Use the function `optimal_normal()`. Supply the same input values as in test case 01.10 to the function except for the following changes and additions: Set the parameter “fixed” to be TRUE. Set the weight for the prior distribution to be NULL, the number of events to be NULL and NULL and the truncation values to be a = NULL and b = NULL.

01.12 (shows that req. 01.08 is met) Use the `function optimal_tte`. Supply the same input values as in test case 01.02 to the function except for the following change: Set the number of cores for parallel computing to 1. Verify that the computation time will increase compared to the setting in 01.02.


+ Start documenting test case here!
