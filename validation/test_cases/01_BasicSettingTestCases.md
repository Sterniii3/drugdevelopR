#' @title 01. Basic setting test cases
#' @editor Johannes Cepicka
#' @editDate 2022-08-16
#' @coverage
#' 01.01: 01.03, 01.06, 01.12, 01.15
#' 01.02: 01.03, 01.05, 01.14, 01.16
#' 01.03: 01.10, 01.12
#' 01.04: 01.09, 01.13
#' 01.05: 01.11, 01.14, 01.15
#' 01.06: 01.04
#' 01.07: 01.07
#' 01.08: 01.02, 01.05, 01.12
#' 01.09: 01.02, 01.06
#' 01.10: 01.01, 01.06
#' 01.11: 01.01, 01.05, 01.15
#' 01.12: 01.08
#' 01.13: 01.01, 01.13, 01.14
#' 01.14: 01.02, 01.13, 01.14
#' 01.15: 01.03, 01.13, 01.14


## 01. Basic setting test cases {-}

### 01.01 (shows that req. 01.03, 01.06, 01.12 and 01.15 are met):  {-}
Use the function `optimal_tte()`. Supply the following input values to the function:
  
  * a significance level of 0.025,
  * a power of 0.9, i.e. $\beta$ of 0.1,
  * assumed true treatment effects of 0.69 and 0.88,
  * event rates of 0.7 for both phase II and phase III,
  * the optimization region {10, 11, …, 400} for the number of participants (events in the time-to-event setting) in phase II,
  * the optimization region {0.71, 0.72, ..., 0.95} for the threshold values,
  * boundaries of 1, 0.95 and 0.85 for the effect size categories small, medium and large,
  * expected gains of 100,000,000\$, 300,000,000\$, and 500,000,000\$ for each effect size, respectively,
  * twelve clusters for parallel computing,
  * fixed costs of 10,000,000\$ in phase II and of 15,000,000\$ in phase III,
  * variable costs of 75,000\$ in phase II and 100,000\$ in phase III,
  * “fixed=FALSE”, i.e. set the function to use a prior distribution,
  * weight of 0.3 for the prior distribution,
  * amount of information for prior true treatment effect given by 210  and 420 expected events for both treatment effects.

Verify that the function calculates an optimal sample size of 206 in phase II and 354 in phase III (i.e. a total of 560 participants), an expected utility of 432 (in 10^5\$), and an optimal threshold value of 0.84 as suggested by Stella Erdmann [2]. Furthermore, verify that one can expect 144 events in phase II and 248 events in phase III (i.e. 392 in total).

### 01.02 (shows that req. 01.03., 01.05., 01.14 and 01.16 are met):  {-}
Use the function `optimal_tte()`. Supply the same input values as in test case 01.01 to the function except for the following: Set the parameter “fixed” to be TRUE, thus using fixed assumed treatment effects and set the assumed true treatment effect, i.e. the hazard ratio to 0.8. As we assume fixed true treatment effects this corresponds to setting hr1 = 0.8, hr2 can be set to 0. Set the weight for the prior distribution to be NULL and the number of events to be NULL and NULL. 
Verify that the function calculates an optimal sample size of 240 in phase II, an expected utility of 352 (in 10^5\$) and an optimal threshold value of 0.88 as suggested by Stella Erdmann [2]. Furthermore, verify that the probability to go to phase III is given by 0.73 and the expected number of events in phase II and III is 168 and 546, respectively (714 in total).


### 01.03 (shows that req. 01.10 and 01.12 are met):  {-}
Use the function `optimal_tte()`. Supply the same input values as in test case 01.01 to the function except for the following changes and additions: Set the weight for the prior distribution to be 0.6 and set a cost constraint of K=750 (in 10^5\$).
Verify that the function calculates an optimal sample size of 228 in phase II, an expected utility of 996 (in 10^5\$) and an optimal threshold value of 0.84 as suggested by Stella Erdmann [2]. Furthermore, verify that the cost constraint is returned and that the total costs in phase II and III are 271 (in 10^5\$) and 478 (in 10^5\$). 

### 01.04 (shows that req. 01.09 and 01.13 are met):  {-}
Use the function `optimal_tte()`. Supply the same input values as in test case 01.01 to the function except for the following changes and additions: 
Set the weight for the prior distribution to be 0.6  and set a sample size constraint of N=500.
Verify that the function calculates an optimal sample size of 170 in phase II and 328 in phase III (i.e. a total sample size of 498), an expected utility of 956 (in 10^5\$) and an optimal threshold value of 0.83 as suggested by Stella Erdmann [2].

### 01.05 (shows that req. 01.11, 01.14 and 01.15 met):  {-}
Use the function `optimal_tte()`. Supply the same input values as in test case 01.01 to the function except for the following changes and additions: Set the weight for the prior distribution to be 0.6 and set a constraint on the minimal probability of a successful program of S=0.6.
Verify that the function calculates an optimal sample size of 470 in phase II, an expected utility of 899 (in 10^5\$) and an optimal threshold value of 0.89 as suggested by Stella Erdmann [2]. Furthermore, verify, that the probability to go phase III is given by 0.77 and the probability of a successful program is given by 0.6.

### 01.06 (shows that req. 01.04 is met):  {-}
Use the function `optimal_tte()`. Supply the same input values as in test case 01.01 to the function except for the following changes and additions: Set the weight for the prior distribution to be 0.6 and use the option to skip phase II.
Verify that the function calculates an optimal sample size of 824 in phase III (corresponding to the total sample size), an expected utility of 1706 (in 10^5\$) and an effect size of 0.76 used for planning phase III (returned as `median_prior_HR`, as the classical "optimal" threshold value `HRgo` is returned as `Inf`) as suggested by Stella Erdmann [2].

### 01.07 (shows that req. 01.07 is met):  {-}
Use the function `optimal_tte()`. Supply the same input values as in test case 01.01 to the function except for the following changes and additions: Set the weight for the prior distribution to be 0.6 and use the option to model different population structures and set the parameter $\gamma$ to 0.025.
Verify that the function calculates an optimal sample size of 310 in phase II, an expected utility of 1207 (in 10^5\$) and an optimal threshold value of 0.86 as suggested by Stella Erdmann [2].

### 01.08 (shows that req. 01.02, 01.05 and 01.12 are met):  {-}
Use the function ` optimal_binary()`. Supply the following input values to the function:

  * a significance level of 0.025,
  * a power of 0.9, i.e. $\beta$ of 0.1,
  * assumed true treatment rate of 0.6 in the control group and assumed true rates of 0.5 and NULL for the prior distribution of the treatment group, 
  * the optimization region of all even numbers {10, 12, …, 500} for the number of participants in phase II,
  * the optimization region {0.7, 0.71, …, 0.9} for the threshold values,
  * boundaries of 1, 0.95 and 0.85 for the effect size categories small, medium and large,
  * expected gains of 100,000,000\$, 300,000,000\$, and 500,000,000\$ for each effect size, respectively,
  * twelve clusters for parallel computing,
  * fixed costs of 10,000,000\$ in phase II and of 15,000,000\$ in phase III,
  * variable costs of 75,000\$ in phase II and 100,000\$ in phase III,
  * `fixed=TRUE`, i.e. set the function to use fixed treatment effects not modelled on a prior distribution,
  * weight of NULL for the prior distribution,
  * NULL and NULL events for both treatment effects.

Verify that the function calculates an optimal sample size of 204, an expected utility of 299 (in 10^5\$) and an optimal threshold value of 0.90 as suggested by Stella Erdmann [2]. Furthermore, verify that the programs returns the cost constraint (Inf in this case) as well as the total costs in phase II and III, 253 (in 10^5\$) and 769 (in 10^5\$), respectively.

### 01.09 (shows that req. 01.02 and 01.06 are met):  {-}
Use the function `optimal_binary()`. Supply the same input values as in test case 01.08 to the function except for the following changes and additions: Set the parameter “fixed” to be FALSE. Set the weight for the prior distribution to be 0.4, the assumed true rates for the prior distribution of the treatment group to be 0.3 and 0.5 and the number of events to be 30 and 60. 
Verify that the function calculates an optimal sample size of 224, an expected utility of 1542 (in 10^5\$) and an optimal threshold value of 0.89 as suggested by Stella Erdmann [2].

### 01.10 (shows that req. 01.01 and 01.06 are met): {-}
Use the function `optimal_normal()`. Supply the following input values to the function:

  * a significance level of 0.025,
  * a power of 0.9, i.e. $\beta$ of 0.1,
  * assumed true treatment effects of 0.375 and 0.5,
  * the optimization region of even numbers {10, 12, …, 500} for the number of participants in phase II,
  * the optimization region {0.01, 0.02,…, 0.5} for the threshold values,
  * boundaries of 0, 0.375 and 0.625 for the effect size categories small, medium and large,
  * expected gains of 62,500,000\$, 200,000,000\$ and 1,000,000,000\$ for each effect size, respectively,
  * twelve clusters for parallel computing,
  * fixed costs of 1,500,000\$ in phase II and of 2,000,000\$ in phase III,
  * variable costs of 67,500\$ in phase II and 72,000\$ in phase III,
  * “fixed=FALSE”, i.e. set the function to model the treatment effects on a prior distribution,
  * weight of 0.5 for the prior distribution,
  * 300 and 600 expected events for both treatment effects,
  * truncation values of a = 0 and b = 0.75.

Verify that the function calculates an optimal sample size of 86 in phase II, an expected utility of 337 (in 10^5\$) and an optimal threshold value of 0.19 as suggested by Stella Erdmann [2].

### 01.11 (shows that req. 01.01, 01.05 and 01.15 are met): {-}
Use the function `optimal_normal()`. Supply the same input values as in test case 01.10 to the function except for the following changes and additions: Set the parameter “fixed” to be TRUE and set the value for the assumed true treatment effect to 0.625 and NULL. Set the weight for the prior distribution to be NULL, the number of events to be NULL and NULL and the truncation values to be a = NULL and b = NULL.
Verify that the function calculates an optimal sample size of 78 in phase II, an expected utility of 944 (in 10^5\$)  and an optimal threshold value of 0.12 as suggested by Stella Erdmann [2]. Furthermore, verify that the probability of a successful program is given by 0.83, which is the sum of the probabilities of a small (0.51), medium (0.30) or large (0.02) treatment effect.

### 01.12 (shows that req. 01.08 is met): {-}
Use the `function optimal_tte`. Supply the same input values as in test case 01.01 to the function except for the following change: Set the number of cores for parallel computing to 6. Verify that the computation time will increase compared to the setting in 01.01.

### 01.13 (shows that req. 01.01, 01.13, 01.14 are met): {-}
Use the function `Epgo_normal()`. Supply 10 sets of the following input values:

  * threshold values $\kappa$ as parameter `kappa`,
  * sample sizes of phase II $n_2$ as parameter `n2`,
  * assumed true treatment effects $\Delta$ as parameter `Delta1`,
  * `fixed=TRUE`, i.e. set the function to use fixed treatment effects not modelled on a prior distribution, and
  * `NULL` for all other input values of the function.

Calculate the function output for all 10 parameter sets and compare the results to the results of a SAS program implementing the probability formula
$$p_{go}^{\Delta}= \mathrm{P}(\hat\Delta_2\geq \kappa|\Delta) = \Phi\left(\frac{\Delta - \kappa}{\sqrt{
		4/n_2}}\right)$$
from [@preussler2020], where $\Phi$ is the cumulative distribution function of the standard normal distribution $\mathcal N(0,1)$.

Use the function `En3_normal()`. Supply 10 sets of the following input values:

  * the same input values as above,
  * significance levels $\alpha$ as parameter `alpha`, and
  * type-II error levels $\beta$ as parameter `beta`.

Calculate the function output for all 10 parameter sets and compare the results to the results of a SAS program implementing the sample size formula
$$ \mathrm E[N_{3}^{\Delta}(\hat\Delta_2)\cdot 1_{\{ \hat\Delta_2 \geq \kappa\}}] = \int^{\infty}_{\kappa} N_{3}^{\Delta}(\hat\Delta_2) \cdot f(\hat\Delta_{2})   d\hat\Delta_{2} ,$$
where
$$N_3^{\Delta,\tau} = N_3^{\Delta,\tau}(\hat \Delta_2) = \frac{4\cdot(z_{1-\alpha}+z_{1-\beta})^2}{(\hat \Delta_2-\tau/\sigma)^2},$$
where $\tau/\sigma = 0$ (i.e. testing for superiority) and $f(\hat\Delta_{2})$ is the probability density function of $\mathcal{N}(\mu = \Delta, \sigma^2 = 4/n_2)$. This formula is the fixed case (no prior distribution) of eq. 2.8 from [@preussler2020]. (Note: R and SAS use the standard deviation for specifying normal distribution. Here, we specify using the variance.)

### 01.14 (shows that req. 01.02, 01.13, 01.14 are met): {-}
Use the function `Epgo_binary()`. Supply 10 sets of the following input values:

  * threshold risk ratios $RR_{go}$ as parameter `RRgo`,
  * sample sizes of phase II $n_2$ as parameter `n2`,
  * assumed true rate in the control group $p_0$ as parameter `p0`,
  * assumed true rate in the treatment group $p_1$ as parameter `p11`,
  * `fixed=TRUE`, i.e. set the function to use fixed treatment effects not modelled on a prior distribution, and
  * `NULL` for all other input values of the function.

Calculate the function output for all 10 parameter sets and compare the results to the results of a SAS program implementing the probability formula
$$p_{go}^{\varrho}= \Phi\left(\frac{\varrho - \kappa}{\sqrt{2/n_2 \cdot(\frac{1-p_0}{p_0}+ \frac{1-p_1}{p_1}) }}\right)$$
from [@preussler2020] where $\varrho = -\log(p_{11}/p_0)$ and $\kappa=-\log(RR_{go})$ and $\Phi$ is the cumulative distribution function of the standard normal distribution $\mathcal N(0,1)$.

Use the function `En3_binary()`. Supply 10 sets of the following input values:

  * the same input values as above,
  * significance levels $\alpha$ as parameter `alpha`, and
  * type-II error levels $\beta$ as parameter `beta`.

Calculate the function output for all 10 parameter sets and compare the results to the results of a SAS program implementing the sample size formula
$$ \mathrm E[N_{3}^{\varrho}(\hat\varrho_2)\cdot 1_{\{ \hat\varrho_2 \geq \kappa\}}] = \int^{\infty}_{\kappa} N_{3}^{\varrho}(\hat\varrho_2) \cdot f(\hat\varrho_{2})   d\hat\varrho_{2} ,$$
where
$$N_3^{\varrho}(\hat\varrho_2) = \frac{2\cdot\bigg(z_{1-\alpha}\cdot\sqrt{\frac{2 \cdot (1-p)}{p}}+z_{1-\beta}\cdot \sqrt{\frac{1-p_0}{p_0}+\frac{1-p_1}{p_1}}\bigg)^2}{\hat \varrho_2^2},$$
and $f(\hat\varrho_{2})$ is the probability density function of the normal distribution $\mathcal{N}(\mu = \varrho, \sigma^2 = 2/n_2\cdot(\frac{1-p_0}{p_0}+ \frac{1-p_1}{p_1}))$, which is the fixed case (no prior distribution) of eq. 2.8 from [@preussler2020].

### 01.15 (shows that req. 01.03, 01.13, 01.14 are met): {-}
Use the function `Epgo_tte()`. Supply 10 sets of the following input values:

  * threshold values $HR_{go}$ as parameter `HRgo`,
  * number of events in phase II $d_2$ as parameter `d2`,
  * assumed true treatment effects $HR$ as parameter `hr1`,
  * `fixed=TRUE`, i.e. set the function to use fixed treatment effects not modelled on a prior distribution, and
  * `NULL` for all other input values of the function.

Calculate the function output for all 10 parameter sets and compare the results to the results of a SAS program implementing the probability formula
$$p_{go}^{\theta}=\mathrm P(\hat\theta_2 \geq \kappa|\theta)= \Phi\left(\frac{\theta-\kappa}{\sqrt{4/d_2}}\right),$$
from [@preussler2020], where $\theta=-\log(HR)$, $\kappa=-\log(HR_{go})$ and $\Phi$ is the cumulative distribution function of the standard normal distribution $\mathcal N(0,1)$.

Use the function `Ed3_tte()`. Supply 10 sets of the following input values:

  * the same input values as above,
  * significance levels $\alpha$ as parameter `alpha`, and
  * type-II error levels $\beta$ as parameter `beta`.

Calculate the function output for all 10 parameter sets and compare the results to the results of a SAS program implementing the event number formula
$$d_3 = \int^{\infty}_{\kappa} D_{3}(\hat\theta_2) \cdot f(\hat\theta_{2}) \ d\hat\theta_{2} $$
with
$$ D_3 =\frac{4 \cdot (z_{1-\alpha}+z_{1-\beta})^{2}}{\hat\theta^{2}_{2}}$$
and $f(\hat\theta_{2})$ being probability density function of the normal distribution $\mathcal N(\mu = \theta, \sigma^2 = 4/d_2)$, which is the fixed case (no prior distribution) of eq. 2.7 from [@preussler2020].

