#' @title Expected sample size for phase III for bias adjustment programs and binary distributed outcomes
#' @description To discount for overoptimistic results in phase II when calculating the optimal sample size in phase III, 
#'it is necessary to use the functions `En3_binary_L()`, `En3_binary_L2()`, `En3_binary_R()` and `En3_binary_R2()`.
#'Each function describes a specific case:
#'- `En3_binary_L()`: calculates the optimal sample size for an additive adjustment factor (i.e. adjust the lower bound of the one-sided confidence interval), 
#'however the go-decision is not affected by the bias adjustment
#'- `En3_binary_L2()`: calculates the optimal sample size for an additive adjustment factor (i.e. adjust the lower bound of the one-sided confidence interval)
#'when the go-decision is also affected by the bias adjustment
#'- `En3_binary_R()`: calculates the optimal sample size for a multiplicative adjustment factor (i.e. use estimate with a retention factor), 
#'however the go-decision is not affected by the bias adjustment
#'- `En3_binary_R2()`: calculates the optimal sample size for a multiplicative adjustment factor (i.e. use estimate with a retention factor)
#'when the go-decision is also affected by the bias adjustment
#' @param RRgo threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param Adj adjustment parameter
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param w weight for mixture prior distribution
#' @param p0 assumed true rate of control group
#' @param p11 assumed true rate of treatment group
#' @param p12 assumed true rate of treatment group
#' @param in1 amount of information for `p11` in terms of sample size
#' @param in2 amount of information for `p12` in terms of sample size
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `p11` is used as fixed effect
#' @return The output of the functions `En3_binary_L`, `En3_binary_L2`, `En3_binary_R` and `En3_binary_R2` is the expected number of participants in phase III. 
#' @importFrom stats qnorm integrate dnorm
#' @examples res <- En3_binary_L(RRgo = 0.8, n2 = 50, Adj = 0, 
#'                              alpha = 0.025, beta = 0.1, p0 = 0.6,  w = 0.3,
#'                              p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600, 
#'                              fixed = FALSE)
#'          res <-  En3_binary_L2(RRgo = 0.8, n2 = 50, Adj = 0, 
#'                              alpha = 0.025, beta = 0.1, p0 = 0.6,  w = 0.3,
#'                              p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600, 
#'                              fixed = FALSE)
#'          res <- En3_binary_R(RRgo = 0.8, n2 = 50, Adj = 1, 
#'                              alpha = 0.025, beta = 0.1, p0 = 0.6,  w = 0.3,
#'                              p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600, 
#'                              fixed = FALSE)
#'          res <- En3_binary_R2(RRgo = 0.8, n2 = 50, Adj = 1, 
#'                              alpha = 0.025, beta = 0.1, p0 = 0.6,  w = 0.3,
#'                              p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600, 
#'                              fixed = FALSE)
#'                              
#' @name En3_bias_binary                             
#' @export 
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
En3_binary_L <- function(){}

#' @title Expected probability of a successful program for bias adjustment programs with binary distributed outcomes
#' @description To discount for overoptimistic results in phase II when calculating the optimal sample size in phase III, 
#'it is necessary to use the following functions, which each describe a specific case:
#'- `EPsProg_binary_L()`: calculates the expected probability of a successful for an additive adjustment factor (i.e. adjust the lower bound of the one-sided confidence interval), 
#'however the go-decision is not affected by the bias adjustment
#'- `EPsProg_binary_L2()`: calculates the expected probability of a successful for an additive adjustment factor (i.e. adjust the lower bound of the one-sided confidence interval)
#'when the go-decision is also affected by the bias adjustment
#'- `EPsProg_binary_R()`: calculates the expected probability of a successful for a multiplicative adjustment factor (i.e. use estimate with a retention factor), 
#'however the go-decision is not affected by the bias adjustment
#'- `EPsProg_binary_R2()`: calculates the expected probability of a successful for a multiplicative adjustment factor (i.e. use estimate with a retention factor)
#'when the go-decision is also affected by the bias adjustment
#' @param RRgo threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param Adj adjustment parameter
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param step1 lower boundary for effect size
#' @param step2 upper boundary for effect size
#' @param w weight for mixture prior distribution
#' @param p0 assumed true rate of control group
#' @param p11 assumed true rate of treatment group
#' @param p12 assumed true rate of treatment group
#' @param in1 amount of information for `p11` in terms of sample size
#' @param in2 amount of information for `p12` in terms of sample size
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `p11` is used as fixed effect
#' @return  The output of the functions `EPsProg_binary_L()`, `EPsProg_binary_L2()`, `EPsProg_binary_R()` and `EPsProg_binary_R2()` is the expected probability of a successful program.
#' @importFrom stats qnorm integrate dnorm
#' @examples res <- EPsProg_binary_L(RRgo = 0.8, n2 = 50, Adj = 0, 
#'                                 alpha = 0.025, beta = 0.1, 
#'                                 step1 = 1, step2 = 0.95, p0 = 0.6,  w = 0.3,
#'                                 p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600, 
#'                                 fixed = FALSE)
#'          res <- EPsProg_binary_L2(RRgo = 0.8, n2 = 50, Adj = 0, 
#'                                 alpha = 0.025, beta = 0.1, 
#'                                 step1 = 1, step2 = 0.95, p0 = 0.6,  w = 0.3,
#'                                 p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600, 
#'                                 fixed = FALSE)
#'          res <- EPsProg_binary_R(RRgo = 0.8, n2 = 50, Adj = 1, 
#'                                 alpha = 0.025, beta = 0.1, 
#'                                 step1 = 1, step2 = 0.95, p0 = 0.6,  w = 0.3,
#'                                 p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600, 
#'                                 fixed = FALSE)
#'          res <- EPsProg_binary_R2(RRgo = 0.8, n2 = 50, Adj = 1, 
#'                                 alpha = 0.025, beta = 0.1, 
#'                                 step1 = 1, step2 = 0.95, p0 = 0.6,  w = 0.3,
#'                                 p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600, 
#'                                 fixed = FALSE)
#' @name EPsProg_bias_binary                               
#' @export 
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
EPsProg_binary_L <- function(){}

#' @title Utility function for bias adjustment programs with binary distributed outcomes.
#' @description The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters and the expected probability of a successful program. 
#'The utility is in a further step maximized by the `optimal_bias_binary()` function.
#' @param RRgo threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param Adj adjustment parameter
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param w weight for mixture prior distribution
#' @param p0 assumed true rate of control group
#' @param p11 assumed true rate of treatment group
#' @param p12 assumed true rate of treatment group
#' @param in1 amount of information for `p11` in terms of sample size
#' @param in2 amount of information for `p12` in terms of sample size
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param steps1 lower boundary for effect size category `"small"` in RR scale, default: 1
#' @param stepm1 lower boundary for effect size category `"medium"` in RR scale = upper boundary for effect size category "small" in RR scale, default: 0.95
#' @param stepl1 lower boundary for effect size category `"large"` in RR scale = upper boundary for effect size category "medium" in RR scale, default: 0.85
#' @param b1 expected gain for effect size category `"small"`
#' @param b2 expected gain for effect size category `"medium"`
#' @param b3 expected gain for effect size category `"large"`
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `p11` is used as fixed effect
#' @return The output of the functions `utility_binary_L()`, `utility_binary_L2()`, `utility_binary_R()` and `utility_binary_R2()` is the expected utility of the program.
#' @examples res <- utility_binary_L(n2 = 50, RRgo = 0.8, Adj = 0.1, w = 0.3, 
#'                                 p0 = 0.6, p11 =  0.3, p12 = 0.5, 
#'                                 in1 = 300, in2 = 600, 
#'                                 alpha = 0.025, beta = 0.1,
#'                                 c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                 K = Inf, N = Inf, S = -Inf,
#'                                 steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
#'                                 b1 = 1000, b2 = 2000, b3 = 3000, 
#'                                 fixed = TRUE)
#'         res <- utility_binary_L2(n2 = 50, RRgo = 0.8, Adj = 0.1, w = 0.3, 
#'                                 p0 = 0.6, p11 =  0.3, p12 = 0.5, 
#'                                 in1 = 300, in2 = 600, 
#'                                 alpha = 0.025, beta = 0.1,
#'                                 c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                 K = Inf, N = Inf, S = -Inf,
#'                                 steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
#'                                 b1 = 1000, b2 = 2000, b3 = 3000, 
#'                                 fixed = TRUE)
#'         res <- utility_binary_R(n2 = 50, RRgo = 0.8, Adj = 0.9, w = 0.3, 
#'                                 p0 = 0.6, p11 =  0.3, p12 = 0.5, 
#'                                 in1 = 300, in2 = 600, 
#'                                 alpha = 0.025, beta = 0.1,
#'                                 c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                 K = Inf, N = Inf, S = -Inf,
#'                                 steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
#'                                 b1 = 1000, b2 = 2000, b3 = 3000, 
#'                                 fixed = TRUE)
#'         res <- utility_binary_R2(n2 = 50, RRgo = 0.8, Adj = 0.9, w = 0.3, 
#'                                 p0 = 0.6, p11 =  0.3, p12 = 0.5, 
#'                                 in1 = 300, in2 = 600, 
#'                                 alpha = 0.025, beta = 0.1,
#'                                 c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                 K = Inf, N = Inf, S = -Inf,
#'                                 steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
#'                                 b1 = 1000, b2 = 2000, b3 = 3000, 
#'                                 fixed = TRUE)
#' @name utility_bias_binary                                 
#' @export 
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
utility_binary_L <- function(){}

#' @title Expected probability to go to phase III for bias adjustment programs with binary distributed outcomes
#' @description In the case we do not only want do discount for overoptimistic results in phase II when calculating the sample size in phase III, 
#'but also when deciding whether to go to phase III or not the functions `Epgo_binary_L2` and `Epgo_binary_R2` are necessary.
#'The function `Epgo_binary_L2` uses an additive adjustment parameter (i.e. adjust the lower bound of the one-sided confidence interval),
#'the function `Epgo_binary_R2` uses a multiplicative adjustment parameter (i.e. use estimate with a retention factor)
#' @param RRgo threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param Adj adjustment parameter
#' @param w weight for mixture prior distribution
#' @param p0 assumed true rate of control group
#' @param p11 assumed true rate of treatment group
#' @param p12 assumed true rate of treatment group
#' @param in1 amount of information for `p11` in terms of sample size
#' @param in2 amount of information for `p12` in terms of sample size
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `p11` is used as fixed effect
#' @return The output of the functions `Epgo_normal_L2` and `Epgo_normal_R2` is the expected number of participants in phase III with conservative decision rule and sample size calculation.
#' @importFrom stats qnorm integrate dnorm
#' @examples res <- Epgo_binary_L2(RRgo = 0.8, n2 = 50, Adj = 0,  p0 = 0.6,  w = 0.3,
#'                              p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600, 
#'                              fixed = FALSE)
#'          res <- Epgo_binary_R2(RRgo = 0.8, n2 = 50, Adj = 1,  p0 = 0.6,  w = 0.3,
#'                              p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600, 
#'                              fixed = FALSE)
#' @name Epgo_bias_binary 
#' @export 
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
Epgo_binary_L2 <- function(){}

#' @title Expected sample size for phase III for bias adjustment programs and normally distributed outcomes
#' @description To discount for overoptimistic results in phase II when calculating the optimal sample size in phase III, 
#'it is necessary to use the functions `En3_normal_L()`, `En3_normal_L2()`, `En3_normal_R()` and `En3_normal_R2()`.
#'Each function describes a specific case:
#'- `En3_normal_L()`: calculates the optimal sample size for an additive adjustment factor (i.e. adjust the lower bound of the one-sided confidence interval), 
#'however the go-decision is not affected by the bias adjustment
#'- `En3_normal_L2()`: calculates the optimal sample size for an additive adjustment factor (i.e. adjust the lower bound of the one-sided confidence interval)
#'when the go-decision is also affected by the bias adjustment
#'- `En3_normal_R()`: calculates the optimal sample size for a multiplicative adjustment factor (i.e. use estimate with a retention factor), 
#'however the go-decision is not affected by the bias adjustment
#'- `En3_normal_R2()`: calculates the optimal sample size for a multiplicative adjustment factor (i.e. use estimate with a retention factor)
#'when the go-decision is also affected by the bias adjustment
#' @param kappa threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param Adj adjustment parameter
#' @param alpha significance level
#' @param beta `1-beta power` for calculation of sample size for phase III
#' @param w weight for mixture prior distribution
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param in1 amount of information for `Delta1` in terms of sample size
#' @param in2 amount of information for `Delta2` in terms of sample size
#' @param a lower boundary for the truncation
#' @param b upper boundary for the truncation
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `Delta1` is used as fixed effect
#' @return The output of the functions `En3_normal_L`, `En3_normal_L2`, `En3_normal_R` and `En3_normal_R2` is the expected number of participants in phase III.
#' @importFrom stats qnorm integrate dnorm
#' @examples res <- En3_normal_L(kappa = 0.1, n2 = 50, Adj = 0, 
#'                              alpha = 0.025, beta = 0.1, w = 0.3,
#'                              Delta1 = 0.375, Delta2 = 0.625, 
#'                              in1 = 300, in2 = 600, 
#'                              a = 0.25, b = 0.75, fixed = FALSE)
#'          res <- En3_normal_L2(kappa = 0.1, n2 = 50, Adj = 0, 
#'                              alpha = 0.025, beta = 0.1, w = 0.3,
#'                              Delta1 = 0.375, Delta2 = 0.625, 
#'                              in1 = 300, in2 = 600, 
#'                              a = 0.25, b = 0.75, fixed = TRUE)
#'          res <- En3_normal_R(kappa = 0.1, n2 = 50, Adj = 1, 
#'                              alpha = 0.025, beta = 0.1, w = 0.3,
#'                              Delta1 = 0.375, Delta2 = 0.625, 
#'                              in1 = 300, in2 = 600, 
#'                              a = 0.25, b = 0.75, fixed = FALSE)
#'          res <- En3_normal_R2(kappa = 0.1, n2 = 50, Adj = 1, 
#'                              alpha = 0.025, beta = 0.1, w = 0.3,
#'                              Delta1 = 0.375, Delta2 = 0.625, 
#'                              in1 = 300, in2 = 600, 
#'                              a = 0.25, b = 0.75, fixed = FALSE)
#' @name En3_bias_normal                           
#' @export 
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
En3_normal_L <- function(){}

#' @title Expected probability of a successful program for bias adjustment programs with normally distributed outcomes
#' @description To discount for overoptimistic results in phase II when calculating the optimal sample size in phase III, 
#'it is necessary to use the following functions, which each describe a specific case:
#'- `EPsProg_normal_L()`: calculates the expected probability of a successful for an additive adjustment factor (i.e. adjust the lower bound of the one-sided confidence interval), 
#'however the go-decision is not affected by the bias adjustment
#'- `EPsProg_normal_L2()`: calculates the expected probability of a successful for an additive adjustment factor (i.e. adjust the lower bound of the one-sided confidence interval)
#'when the go-decision is also affected by the bias adjustment
#'- `EPsProg_normal_R()`: calculates the expected probability of a successful for a multiplicative adjustment factor (i.e. use estimate with a retention factor), 
#'however the go-decision is not affected by the bias adjustment
#'- `EPsProg_normal_R2()`: calculates the expected probability of a successful for a multiplicative adjustment factor (i.e. use estimate with a retention factor)
#'when the go-decision is also affected by the bias adjustment
#' @param kappa threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param Adj adjustment parameter
#' @param alpha significance level
#' @param beta 1-beta power for calculation of sample size for phase III
#' @param step1 lower boundary for effect size
#' @param step2 upper boundary for effect size
#' @param w weight for mixture prior distribution
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param in1 amount of information for `Delta1` in terms of sample size
#' @param in2 amount of information for `Delta2` in terms of sample size
#' @param a lower boundary for the truncation
#' @param b upper boundary for the truncation
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `Delta1` is used as fixed effect
#' @return The output of the functions `EPsProg_normal_L()`, `EPsProg_normal_L2()`, `EPsProg_normal_R()` and `EPsProg_normal_R2()` is the expected probability of a successful program.
#' @importFrom stats qnorm integrate dnorm pnorm
#' @examples res <- EPsProg_normal_L(kappa = 0.1, n2 = 50, Adj = 0, 
#'                                 alpha = 0.025, beta = 0.1, w = 0.3,
#'                                 step1 = 0, step2 = 0.5,
#'                                 Delta1 = 0.375, Delta2 = 0.625, 
#'                                 in1 = 300, in2 = 600, 
#'                                 a = 0.25, b = 0.75, fixed = FALSE)
#'          res <- EPsProg_normal_L2(kappa = 0.1, n2 = 50, Adj = 0, 
#'                                 alpha = 0.025, beta = 0.1, w = 0.3,
#'                                 step1 = 0, step2 = 0.5,
#'                                 Delta1 = 0.375, Delta2 = 0.625, 
#'                                 in1 = 300, in2 = 600, 
#'                                 a = 0.25, b = 0.75, fixed = FALSE)
#'          res <- EPsProg_normal_R(kappa = 0.1, n2 = 50, Adj = 1, 
#'                                 alpha = 0.025, beta = 0.1, w = 0.3,
#'                                 step1 = 0, step2 = 0.5,
#'                                 Delta1 = 0.375, Delta2 = 0.625, 
#'                                 in1 = 300, in2 = 600, 
#'                                 a = 0.25, b = 0.75, fixed = FALSE)
#'          res <- EPsProg_normal_R2(kappa = 0.1, n2 = 50, Adj = 1, 
#'                                 alpha = 0.025, beta = 0.1, w = 0.3,
#'                                 step1 = 0, step2 = 0.5,
#'                                 Delta1 = 0.375, Delta2 = 0.625, 
#'                                 in1 = 300, in2 = 600, 
#'                                 a = 0.25, b = 0.75, fixed = FALSE)
#' @name EPsProg_bias_normal                                
#' @export 
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
EPsProg_normal_L <- function(){}

#' @title Utility function for bias adjustment programs with normally distributed outcomes.
#' @description The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters and the expected probability of a successful program. 
#'The utility is in a further step maximized by the `optimal_bias_normal()` function.
#' @param n2 total sample size for phase II; must be even number
#' @param kappa threshold value for the go/no-go decision rule
#' @param Adj adjustment parameter
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param w weight for mixture prior distribution
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param in1 amount of information for `Delta1` in terms of sample size
#' @param in2 amount of information for `Delta2` in terms of sample size
#' @param a lower boundary for the truncation
#' @param b upper boundary for the truncation
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param steps1 lower boundary for effect size category `"small"`, default: 0
#' @param stepm1 lower boundary for effect size category `"medium"` = upper boundary for effect size category "small" default: 0.5
#' @param stepl1 lower boundary for effect size category `"large"` = upper boundary for effect size category "medium", default: 0.8
#' @param b1 expected gain for effect size category `"small"`
#' @param b2 expected gain for effect size category `"medium"`
#' @param b3 expected gain for effect size category `"large"`
#' @param fixed choose if true treatment effects are fixed or random, if TRUE Delta1 is used as fixed effect
#' @return The output of the functions `utility_normal_L()`, `utility_normal_L2()`, `utility_normal_R()` and `utility_normal_R2()` is the expected utility of the program.
#' @examples res <- utility_normal_L(kappa = 0.1, n2 = 50, Adj = 0, 
#'                                 alpha = 0.025, beta = 0.1, w = 0.3,
#'                                 Delta1 = 0.375, Delta2 = 0.625, 
#'                                 in1 = 300, in2 = 600, 
#'                                 a = 0.25, b = 0.75, 
#'                                 c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                 K = Inf, N = Inf, S = -Inf, 
#'                                 steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
#'                                 b1 = 3000, b2 = 8000, b3 = 10000, 
#'                                 fixed = TRUE)
#'          res <- utility_normal_L2(kappa = 0.1, n2 = 50, Adj = 0, 
#'                                 alpha = 0.025, beta = 0.1, w = 0.3,
#'                                 Delta1 = 0.375, Delta2 = 0.625, 
#'                                 in1 = 300, in2 = 600, 
#'                                 a = 0.25, b = 0.75, 
#'                                 c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                 K = Inf, N = Inf, S = -Inf, 
#'                                 steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
#'                                 b1 = 3000, b2 = 8000, b3 = 10000, 
#'                                 fixed = TRUE)
#'          res <- utility_normal_R(kappa = 0.1, n2 = 50, Adj = 1, 
#'                                 alpha = 0.025, beta = 0.1, w = 0.3,
#'                                 Delta1 = 0.375, Delta2 = 0.625, 
#'                                 in1 = 300, in2 = 600, 
#'                                 a = 0.25, b = 0.75, 
#'                                 c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                 K = Inf, N = Inf, S = -Inf, 
#'                                 steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
#'                                 b1 = 3000, b2 = 8000, b3 = 10000, 
#'                                 fixed = TRUE)
#'          res <- utility_normal_R2(kappa = 0.1, n2 = 50, Adj = 1, 
#'                                 alpha = 0.025, beta = 0.1, w = 0.3,
#'                                 Delta1 = 0.375, Delta2 = 0.625, 
#'                                 in1 = 300, in2 = 600, 
#'                                 a = 0.25, b = 0.75, 
#'                                 c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                 K = Inf, N = Inf, S = -Inf, 
#'                                 steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
#'                                 b1 = 3000, b2 = 8000, b3 = 10000, 
#'                                 fixed = TRUE)
#' @name utility_bias_normal                               
#' @export 
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
utility_normal_L <- function(){}

#' @title Expected probability to go to phase III for bias adjustment programs with normally distributed outcomes
#' @description In the case we do not only want do discount for overoptimistic results in phase II when calculating the sample size in phase III, 
#'but also when deciding whether to go to phase III or not the functions `Epgo_normal_L2` and `Epgo_normal_R2` are necessary.
#'The function `Epgo_normal_L2` uses an additive adjustment parameter (i.e. adjust the lower bound of the one-sided confidence interval),
#'the function `Epgo_normal_R2` uses a multiplicative adjustment parameter (i.e. use estimate with a retention factor)
#' @param kappa threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param Adj adjustment parameter
#' @param w weight for mixture prior distribution
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param in1 amount of information for `Delta1` in terms of sample size
#' @param in2 amount of information for `Delta2` in terms of sample size
#' @param a lower boundary for the truncation
#' @param b upper boundary for the truncation
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `Delta1` is used as fixed effect
#' @return The output of the functions `Epgo_normal_L2` and `Epgo_normal_R2` is the expected number of participants in phase III with conservative decision rule and sample size calculation.
#' @importFrom stats qnorm integrate dnorm
#' @examples res <- Epgo_normal_L2(kappa = 0.1, n2 = 50, Adj = 0, w = 0.3,
#'                               Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#'                               a = 0.25, b = 0.75, fixed = FALSE)
#'          res <- Epgo_normal_R2(kappa = 0.1, n2 = 50, Adj = 1, w = 0.3,
#'                               Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#'                               a = 0.25, b = 0.75, fixed = FALSE)
#' @name Epgo_bias_normal                              
#' @export 
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
Epgo_normal_L2 <- function(){}

#' @title Probability to go to phase III for multiarm programs with binary distributed outcomes 
#' @description Given our parameters this function calculates the probability to go to phase III after the second phase was conducted. The considered strategies are as follows:
#'- 1. Strategy: Only best promising treatment goes to phase III
#' @param RRgo threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param p0 assumed true rate of control group
#' @param p11 assumed true rate of treatment group
#' @param p12 assumed true rate of treatment group
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)
#' @param case different cases: 1 ("nogo"), 21 (treatment 1 is promising, treatment 2 is not), 22 (treatment 2 is promising, treatment 1 is not), 31 (both treatments are promising, treatment 1 is better), 32 (both treatments are promising, treatment 2 is better)
#' @return The function pgo_binary() returns the probability to go to phase III.
#' @examples res <- pgo_binary(RRgo = 0.8 ,n2 = 50 ,p0 = 0.6, p11 =  0.3, p12 = 0.5,strategy = 3, case = 31)
#' @export 
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
pgo_binary <- function(){}

#' @title Total sample size for phase III trial with l treatments and equal allocation ratio for binary outcomes 
#' @description  Depending on the results of phase II and our strategy ,i.e. whether we proceed only with the best promising treatment (l = 1) or with all promising treatments (l = 2), this program calculates the number of participants in phase III.
#' 
#' l=1: according to Schoenfeld to guarantee power for the log rank test to detect treatment effect of phase II; 
#' 
#' l=2: according to Dunnett to guarantee y any-pair power (Horn & Vollandt)
#' @param alpha significance level
#' @param beta  1-beta power for calculation of sample size for phase III
#' @param p0 assumed true rate of control group
#' @param p11 assumed true rate of treatment group
#' @param y hat_theta_2; estimator in phase II
#' @param l number of treatments in phase III:
#'- l=1: according to Schoenfeld to guarantee power for the log rank test to detect treatment effect of phase II;  
#'- l=2: according to Dunnett to guarantee y any-pair power (Horn & Vollandt)
#' @return the function ss_binary() returns the total sample size for phase III trial with l treatments and equal allocation ratio
#' @examples res <- ss_binary(alpha = 0.05, beta = 0.1, p0 = 0.6, p11 = 0.3, y = 0.5, l = 1)
#' @export 
#' @editor Johannes Cepicka
#' @editDate 2022-05-08
ss_binary <- function(){}

#' @title Expected sample size for phase III for multiarm programs with binary distributed outcomes
#' @description Given phase II results are promising enough to get the "go"-decision to go to phase III this function now calculates the expected sample size for phase III given the cases and strategies listed below.
#'The results of this function are necessary for calculating the utility of the program, which is then in a further step maximized by the `optimal_multiarm_binary()` function
#' @param RRgo threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta  1-beta power for calculation of sample size for phase III
#' @param p0 assumed true rate of control group
#' @param p11 assumed true rate of treatment group
#' @param p12 assumed true rate of treatment group
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)
#' @param case different cases: 1 ("nogo"), 21 (treatment 1 is promising, treatment 2 is not), 22 (treatment 2 is promising, treatment 1 is not), 31 (both treatments are promising, treatment 1 is better), 32 (both treatments are promising, treatment 2 is better)
#' @return the function Ess_binary() returns the expected sample size for phase III when going to phase III
#' @examples res <- Ess_binary(RRgo = 0.8 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
#'                            p0 = 0.6, p11 =  0.3, p12 = 0.5,strategy = 3, case = 31)
#' @export 
#' @editor Johannes Cepicka
#' @editDate 2022-05-08
Ess_binary <- function(){}

#' @title Probability of a successful program for multiarm programs with binary distributed outcomes 
#' @description Given we get the "go"-decision in phase II, this functions now calculates the probability that the results of the confirmatory trial (phase III) are significant, i.e. we have a statistically relevant positive effect of the treatment.
#' @param RRgo threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be divisible by three
#' @param alpha significance level
#' @param beta  1-beta power for calculation of sample size for phase III
#' @param p0 assumed true rate of control group
#' @param p11 assumed true rate of treatment group
#' @param p12 assumed true rate of treatment group
#' @param step1 lower boundary for effect size
#' @param step2 upper boundary for effect size
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)
#' @param case different cases: 1 ("nogo"), 21 (treatment 1 is promising, treatment 2 is not), 22 (treatment 2 is promising, treatment 1 is not), 31 (both treatments are promising, treatment 1 is better), 32 (both treatments are promising, treatment 2 is better)
#' @return The function PsProg_binary() returns the probability of a successful program
#' @examples res <- PsProg_binary(RRgo = 0.8 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
#'                            p0 = 0.6, p11 =  0.3, p12 = 0.5, step1 = 1, step2 = 0.95,
#'                            strategy = 3, case = 31)
#' @export 
#' @editor Johannes Cepicka
#' @editDate 2022-05-08
PsProg_binary <- function(){}

#' @title Utility function for multiarm programs with binary distributed outcomes
#' @description The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters as on the sample size and expected probability of a successful program. 
#'The utility is in further step maximized by the `optimal_multiarm_binary()` function.
#' @param RRgo threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta  1-beta power for calculation of sample size for phase III
#' @param p0 assumed true rate of control group
#' @param p11 assumed true rate of treatment group
#' @param p12 assumed true rate of treatment group
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param steps1 lower boundary for effect size category "small" in RR scale, default: 1
#' @param stepm1 lower boundary for effect size category "medium" in RR scale = upper boundary for effect size category "small" in RR scale, default: 0.95
#' @param stepl1 lower boundary for effect size category "large" in RR scale = upper boundary for effect size category "medium" in RR scale, default: 0.85
#' @param b1 expected gain for effect size category "small"
#' @param b2 expected gain for effect size category "medium"
#' @param b3 expected gain for effect size category "large"
#' @return The output of the function `utility_multiarm_binary()` is the expected utility of the program
#' @examples res <- utility_multiarm_binary(n2 = 50, RRgo = 0.8, alpha = 0.05, beta = 0.1,
#'                            p0 = 0.6, p11 =  0.3, p12 = 0.5, strategy = 3,
#'                            c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                            K = Inf, N = Inf, S = -Inf,  
#'                            steps1 = 1, stepm1 = 0.95,   stepl1 = 0.85,
#'                            b1 = 1000, b2 = 2000, b3 = 3000)
#' @editor Johannes Cepicka
#' @editDate 2022-05-08
#' @export 
utility_multiarm_binary <- function(){}

#' @title Probability to go to phase III for multiarm programs with normally distributed outcomes 
#' @description Given our parameters this function calculates the probability to go to phase III after the second phase was conducted. The considered strategies are as follows:
#'- 1. Strategy: Only best promising treatment goes to phase III
#' @param kappa threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be an even number
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)
#' @param case different cases: 1 ("nogo"), 21 (treatment 1 is promising, treatment 2 is not), 22 (treatment 2 is promising, treatment 1 is not), 31 (both treatments are promising, treatment 1 is better), 32 (both treatments are promising, treatment 2 is better)
#' @return The function pgo_normal() returns the probability to go to phase III for multiarm programs with normally distributed outcomes 
#' @examples res <- pgo_normal(kappa = 0.1, n2 = 50, Delta1 = 0.375, Delta2 = 0.625, strategy = 3, case = 31)
#' @editor Johannes Cepicka
#' @editDate 2022-05-08
#' @export 
pgo_normal <- function(){}

#' @title Total sample size for phase III trial with l treatments and equal allocation ratio for normally distributed outcomes
#' @description  Depending on the results of phase II and our strategy ,i.e. whether we proceed only with the best promising treatment (l = 1) or with all promising treatments (l = 2), this program calculates the number of participants in phase III.
#' 
#' l=1: according to Schoenfeld to guarantee power for the log rank test to detect treatment effect of phase II; 
#' 
#' l=2: according to Dunnett to guarantee y any-pair power (Horn & Vollandt)
#' @param alpha significance level
#' @param beta  1-'beta' power for calculation of sample size for phase III
#' @param y  y_hat_theta_2; estimator in phase II
#' @param l number of treatments in phase III:
#'- l=1: according to Schoenfeld to guarantee power for the log rank test to detect treatment effect of phase II;  
#'- l=2: according to Dunnett to guarantee y any-pair power (Horn & Vollandt)
#' @return the function ss_normal() returns the total sample size for phase III trial with l treatments and equal allocation ratio
#' @examples res <- ss_normal(alpha = 0.05, beta = 0.1, y = 0.5, l = 1)
#' @editor Johannes Cepicka
#' @editDate 2022-05-08
#' @export 
ss_normal <- function(){}

#' @title Expected sample size for phase III for multiarm programs with normally distributed outcomes
#' @description Given phase II results are promising enough to get the "go"-decision to go to phase III this function now calculates the expected sample size for phase III given the cases and strategies listed below.
#'The results of this function are necessary for calculating the utility of the program, which is then in a further step maximized by the `optimal_multiarm_normal()` function
#' @param kappa threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta  1-beta power for calculation of sample size for phase III
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)
#' @param case different cases: 1 ("nogo"), 21 (treatment 1 is promising, treatment 2 is not), 22 (treatment 2 is promising, treatment 1 is not), 31 (both treatments are promising, treatment 1 is better), 32 (both treatments are promising, treatment 2 is better)
#' @return The function Ess_normal() returns the expected sample size for phase III when going to phase III when outcomes are normally distributed and we consider multiarm programs, i.e. several phase III trials with different doses or different treatments are performed
#' @examples res <- Ess_normal(kappa = 0.1 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
#'                            Delta1 = 0.375, Delta2 = 0.625, strategy = 3, case = 31)
#' @editor Johannes Cepicka
#' @editDate 2022-05-08
#' @export 
Ess_normal <- function(){}

#' @title Probability of a successful program for multiarm programs with normally distributed outcomes 
#' @description Given we get the "go"-decision in phase II, this functions now calculates the probability that the results of the confirmatory trial (phase III) are significant, i.e. we have a statistically relevant positive effect of the treatment.
#' @param kappa threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta  1-beta power for calculation of sample size for phase III
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param step1 lower boundary for effect size
#' @param step2 upper boundary for effect size
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)
#' @param case different cases: 1 ("nogo"), 21 (treatment 1 is promising, treatment 2 is not), 22 (treatment 2 is promising, treatment 1 is not), 31 (both treatments are promising, treatment 1 is better), 32 (both treatments are promising, treatment 2 is better)
#' @return The function PsProg_normal() returns the probability of a successful program.
#' @examples res <- PsProg_normal(kappa = 0.1 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
#'                            Delta1 = 0.375, Delta2 = 0.625,  step1 = 0, step2 = 0.5,
#'                            strategy = 3, case = 31)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
PsProg_normal <- function(){}

#' @title Utility function for multiarm programs with normally distributed outcomes
#' @description The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters as on the sample size and expected probability of a successful program. 
#'The utility is in further step maximized by the `optimal_multiarm_normal()` function.
#' @param kappa threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta  1-beta power for calculation of sample size for phase III
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param steps1 lower boundary for effect size category "small", default: 0
#' @param stepm1 lower boundary for effect size category "medium" = upper boundary for effect size category "small" default: 0.5
#' @param stepl1 lower boundary for effect size category "large" = upper boundary for effect size category "medium", default: 0.8
#' @param b1 expected gain for effect size category "small"
#' @param b2 expected gain for effect size category "medium"
#' @param b3 expected gain for effect size category "large"
#' @return The output of the function `utility_multiarm_normal()` is the expected utility of the program.
#' @examples res <- utility_multiarm_normal(n2 = 50, kappa = 0.8, alpha = 0.05, beta = 0.1,
#'                            Delta1 = 0.375, Delta2 = 0.625, strategy = 3,
#'                            c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                            K = Inf, N = Inf, S = -Inf,  
#'                            steps1 = 0, stepm1 = 0.5,   stepl1 = 0.8,
#'                            b1 = 1000, b2 = 2000, b3 = 3000)
#' @editor Johannes Cepicka
#' @editDate 2022-05-08
#' @export 
utility_multiarm_normal <- function(){}

#' @title Probability to go to phase III for multiarm programs with time-to-event outcomes 
#' @description Given our parameters this function calculates the probability to go to phase III after the second phase was conducted. The considered strategies are as follows:
#'- 1. Strategy: Only best promising treatment goes to phase III
#' @param HRgo threshold value for the go/no-go decision rule
#' @param n2 total sample size in phase II, must be divisible by 3
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param ec control arm event rate for phase II and III
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)
#' @param case different cases: 1 ("nogo"), 21 (treatment 1 is promising, treatment 2 is not), 22 (treatment 2 is promising, treatment 1 is not), 31 (both treatments are promising, treatment 1 is better), 32 (both treatments are promising, treatment 2 is better)
#' @return The function pgo_tte() returns the probability to go to phase III.
#' @examples res <- pgo_tte(HRgo = 0.8, n2 = 48 , ec = 0.6, hr1 = 0.7, hr2 = 0.8, strategy = 3, case = 31)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
pgo_tte <- function(){}

#' @title Total sample size for phase III trial with l treatments and equal allocation ratio for time-to-event outcomes 
#' @description  Depending on the results of phase II and our strategy ,i.e. whether we proceed only with the best promising treatment (l = 1) or with all promising treatments (l = 2), this program calculates the number of participants in phase III.
#' 
#' l=1: according to Schoenfeld to guarantee power for the log rank test to detect treatment effect of phase II; 
#' 
#' l=2: according to Dunnett to guarantee y any-pair power (Horn & Vollandt)
#' @param alpha significance level
#' @param beta  1-beta power for calculation of sample size for phase III
#' @param ec control arm event rate for phase II and III
#' @param ek event rate of arm k (either arm 1 or arm 2)
#' @param y hat_theta_2; estimator in phase II
#' @param l number of treatments in phase III:
#'- l=1: according to Schoenfeld to guarantee power for the log rank test to detect treatment effect of phase II;  
#'- l=2: according to Dunnett to guarantee y any-pair power (Horn & Vollandt)
#' @return the function ss_tte() returns the total sample size for phase III trial with l treatments and equal allocation ratio
#' @examples res <- ss_tte(alpha = 0.05, beta = 0.1, ec = 0.6, ek = 0.8, y = 0.5, l=1)
#' @editor Johannes Cepicka
#' @editDate 2022-05-08
#' @export 
ss_tte <- function(){}

#' @title Expected sample size for phase III for multiarm programs with time-to-event outcomes
#' @description Given phase II results are promising enough to get the "go"-decision to go to phase III this function now calculates the expected sample size for phase III given the cases and strategies listed below.
#'The results of this function are necessary for calculating the utility of the program, which is then in a further step maximized by the `optimal_multiarm()` function
#' @param HRgo threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta  1-beta power for calculation of sample size for phase III
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param ec control arm event rate for phase II and III
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)
#' @param case different cases: 1 ("nogo"), 21 (treatment 1 is promising, treatment 2 is not), 22 (treatment 2 is promising, treatment 1 is not), 31 (both treatments are promising, treatment 1 is better), 32 (both treatments are promising, treatment 2 is better)
#' @importFrom cubature adaptIntegrate
#' @return the function Ess_tte() returns the expected sample size for phase III when going to phase III
#' @examples res <- Ess_tte(HRgo = 0.8 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
#'                            ec = 0.6, hr1 = 0.7, hr2 = 0.8, strategy = 2, case = 21)
#' @editor Johannes Cepicka
#' @editDate 2022-05-08
#' @export 
Ess_tte <- function(){}

#' @title Probability of a successful program for multiarm programs with time-to-event outcomes 
#' @description Given we get the "go"-decision in phase II, this functions now calculates the probability that the results of the confirmatory trial (phase III) are significant, i.e. we have a statistically relevant positive effect of the treatment.
#' @param HRgo threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be divisible by three
#' @param alpha significance level
#' @param beta  1-beta power for calculation of sample size for phase III
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param ec control arm event rate for phase II and III
#' @param step1 lower boundary for effect size
#' @param step2 upper boundary for effect size
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)
#' @param case different cases: 1 ("nogo"), 21 (treatment 1 is promising, treatment 2 is not), 22 (treatment 2 is promising, treatment 1 is not), 31 (both treatments are promising, treatment 1 is better), 32 (both treatments are promising, treatment 2 is better)
#' @return The function PsProg_tte() returns the probability of a successful program
#' @examples res <- PsProg_tte(HRgo = 0.8 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
#'                            ec = 0.6, hr1 = 0.7, hr2 = 0.8, step1 = 1, step2 = 0.95,
#'                            strategy = 2, case = 21)
#' @editor Johannes Cepicka
#' @editDate 2022-05-08
#' @export 
PsProg_tte <- function(){}

#' @title Utility function for multiarm programs with time-to-event outcomes
#' @description The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters as on the sample size and expected probability of a successful program. 
#'The utility is in further step maximized by the `optimal_multiarm()` function.
#' @param HRgo threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be divisible by three
#' @param alpha significance level
#' @param beta  1-beta power for calculation of sample size for phase III
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param ec control arm event rate for phase II and III
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param steps1 lower boundary for effect size category `"small"` in HR scale, default: 1
#' @param stepm1 lower boundary for effect size category `"medium"` in HR scale = upper boundary for effect size category "small" in RR scale, default: 0.95
#' @param stepl1 lower boundary for effect size category `"large"` in HR scale = upper boundary for effect size category "medium" in RR scale, default: 0.85
#' @param b1 expected gain for effect size category `"small"`
#' @param b2 expected gain for effect size category `"medium"`
#' @param b3 expected gain for effect size category `"large"`
#' @return The output of the function `utility_multiarm()` is the expected utility of the program
#' @examples res <- utility_multiarm(n2 = 50, HRgo = 0.8, alpha = 0.05, beta = 0.1,
#'                            hr1 = 0.7, hr2 = 0.8, strategy = 3, ec = 0.6,
#'                            c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                            K = Inf, N = Inf, S = -Inf,  
#'                            steps1 = 1, stepm1 = 0.95,  stepl1 = 0.85,
#'                            b1 = 1000, b2 = 2000, b3 = 3000)
#' @editor Johannes Cepicka
#' @editDate 2022-05-08
#' @export 
utility_multiarm <- function(){}

#' @title Density for the minimum of two normally distributed random variables
#' @description The function `fmin()` will return the value of f(z), which is the value of the density function of the
#'minimum of two normally distributed random variables.
#' @details Z= min(X,Y) with X ~ N(mu1,sigma1^2), Y ~ N(mu2,sigma2^2)
#'
#'f(z)=f1(z)+f2(z)
#' @param y integral variable
#' @param mu1 mean of second endpoint 
#' @param mu2 mean of first endpoint
#' @param sigma1 standard deviation of first endpoint
#' @param sigma2 standard deviation of second endpoint
#' @param rho correlation between the two endpoints
#' @return  The function `fmin()` will return the value of f(z), which is the value of the density function of the
#'minimum of two normally distributed random variables.
#' @examples res <- fmin(y = 0.5, mu1 = 0.375, mu2 = 0.25, sigma1 = 8, sigma2 = 12, rho = 0.4 )
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
fmin <- function(){}

#' @title Density of the bivariate normal distribution
#' @param x integral variable
#' @param y integral variable
#' @param mu1 mean of second endpoint 
#' @param mu2 mean of first endpoint
#' @param sigma1 standard deviation of first endpoint
#' @param sigma2 standard deviation of second endpoint
#' @param rho correlation between the two endpoints
#' @return The Function `dbivanorm()` will return the density of a bivariate normal distribution.
#' @examples res <- dbivanorm(x = 0.5, y = 0.5, mu1 = 0.375, mu2 = 0.25, sigma1 = 8, sigma2 = 12, rho = 0.4 )
#' @name dbivanorm
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
dbivanorm <- function(){}

#' @title Probability to go to phase III for multiple endpoints with normally distributed outcomes
#' @description This function calculated the probability that we go to phase III, i.e. that results of phase II are promising enough to
#'get a successful drug development program. Successful means that both endpoints show a statistically significant positive treatment effect in phase III.
#' @param kappa threshold value for the go/no-go decision rule; vector for both endpoints 
#' @param n2 total sample size for phase II; must be even number
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param in1 amount of information for Delta1 in terms of sample size
#' @param in2 amount of information for Delta2 in terms of sample size
#' @param sigma1 standard deviation of first endpoint
#' @param sigma2 standard deviation of second endpoint
#' @param fixed choose if true treatment effects are fixed or random, if TRUE Delta1 is used as fixed effect
#' @param rho correlation between the two endpoints
#' @return The output of the function `pgo_multiple_normal()` is the probability to go to phase III.
#' @examples res <- pgo_multiple_normal(kappa = 0.1, n2 = 50,
#'                               Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#'                               sigma1 = 8, sigma2 = 4, fixed = TRUE, rho = 0.3)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
pgo_multiple_normal <- function(){}

#' @title Expected sample size for phase III for multiple endpoints with normally distributed outcomes
#' @description Given phase II results are promising enough to get the "go"-decision to go to phase III this function now calculates the expected sample size for phase III.
#'The results of this function are necessary for calculating the utility of the program, which is then in a further step maximized by the `optimal_multiple_normal()` function
#' @param kappa threshold value for the go/no-go decision rule; vector for both endpoints
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta 1-beta power for calculation of sample size for phase III
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param in1 amount of information for Delta1 in terms of sample size
#' @param in2 amount of information for Delta2 in terms of sample size
#' @param sigma1 standard deviation of first endpoint
#' @param sigma2 standard deviation of second endpoint
#' @param fixed choose if true treatment effects are fixed or random, if TRUE Delta1 is used as fixed effect
#' @param rho correlation between the two endpoints
#' @return the output of the function Ess_multiple_normal is the expected number of participants in phase III
#' @examples res <- Ess_multiple_normal(kappa = 0.1, n2 = 50, alpha = 0.025, beta = 0.1,
#'                               Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#'                               sigma1 = 8, sigma2 = 4, fixed = TRUE, rho = 0.3)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
Ess_multiple_normal <- function(){}

#' @title Probability of a successful program, when going to phase III for multiple endpoint with normally distributed outcomes
#' @description After getting the "go"-decision to go to phase III, i.e. our results of phase II are over the predefined threshold `kappa`, this function 
#'calculates the probability, that our program is successfull, i.e. that both endpoints show a statistically significant positive treatment effect in phase III.
#' @param kappa threshold value for the go/no-go decision rule; 
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param in1 amount of information for `Delta1` in terms of sample size
#' @param in2 amount of information for `Delta2` in terms of sample size
#' @param sigma1 standard deviation of first endpoint
#' @param sigma2 standard deviation of second endpoint
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `Delta1` is used as fixed effect
#' @param rho correlation between the two endpoints
#' @return The output of the function `posp_normal()` is the probability of a successful program, when going to phase III.
#' @examples res <- posp_normal(kappa = 0.1, n2 = 50, alpha = 0.025, beta = 0.1,
#'                               Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#'                               sigma1 = 8, sigma2 = 4, fixed = TRUE, rho = 0.3)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
posp_normal <- function(){}

#' @title Expected probability of a successful program for multiple endpoints and normally distributed outcomes
#' @description This function calculates the probability that our drug development program is successful.
#'Successful is defined as both endpoints showing a statistically significant positive treatment effect in phase III.
#' @param kappa threshold value for the go/no-go decision rule;
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param in1 amount of information for `Delta1` in terms of sample size
#' @param in2 amount of information for `Delta2` in terms of sample size
#' @param sigma1 standard deviation of first endpoint
#' @param sigma2 standard deviation of second endpoint
#' @param step11 lower boundary for effect size for first endpoint
#' @param step12 lower boundary for effect size for second endpoint
#' @param step21 upper boundary for effect size for first endpoint
#' @param step22 upper boundary for effect size for second endpoint
#' @param fixed choose if true treatment effects are fixed or random, if TRUE then `Delta1` is used as fixed effect
#' @param rho correlation between the two endpoints
#' @return The output of the function `EPsProg_multiple_normal()` is the expected probability of a successfull program, when going to phase III.
#' @examples res <- EPsProg_multiple_normal(kappa = 0.1, n2 = 50, alpha = 0.025, beta = 0.1,
#'                               Delta1 = 0.375, Delta2 = 0.625, sigma1 = 8, sigma2 = 4,
#'                               step11 = 0, step12 = 0, step21 = 0.5, step22 = 0.5, 
#'                               in1 = 300, in2 = 600, fixed = TRUE, rho = 0.3)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
EPsProg_multiple_normal <- function(){}

#' @title Utility function for multiple endpoints with normally distributed outcomes.
#' @description The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters and the expected probability of a successful program. 
#'The utility is in a further step maximized by the `optimal_multiple_normal()` function.
#' @param kappa threshold value for the go/no-go decision rule; vector for both endpoints
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param in1 amount of information for `Delta1` in terms of sample size
#' @param in2 amount of information for `Delta2` in terms of sample size
#' @param sigma1 standard deviation of first endpoint
#' @param sigma2 standard deviation of second endpoint
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint 
#' @param steps1 lower boundary for effect size category `"small"` in HR scale, default: 1
#' @param stepm1 lower boundary for effect size category `"medium"` in HR scale = upper boundary for effect size category `"small"` in HR scale, default: 0.95
#' @param stepl1 lower boundary for effect size category `"large"` in HR scale = upper boundary for effect size category `"medium"` in HR scale, default: 0.85
#' @param b1 expected gain for effect size category `"small"`
#' @param b2 expected gain for effect size category `"medium"`
#' @param b3 expected gain for effect size category `"large"`
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `Delta1` is used as fixed effect
#' @param rho correlation between the two endpoints
#' @param relaxed relaxed or strict decision rule
#' @return The output of the function `utility_multiple_normal()` is the expected utility of the program.
#' @examples res <- utility_multiple_normal(kappa = 0.1, n2 = 50, 
#'                               alpha = 0.025, beta = 0.1,
#'                               Delta1 = 0.375, Delta2 = 0.625, 
#'                               in1 = 300, in2 = 600, sigma1 = 8, sigma2 = 4,
#'                               c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                               K = Inf, N = Inf, S = -Inf,
#'                               steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
#'                               b1 = 1000, b2 = 2000, b3 = 3000, 
#'                               fixed = TRUE, rho = 0.3, relaxed = "TRUE")
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
utility_multiple_normal <- function(){}

#' @title Density for the maximum of two normally distributed random variables
#' @description The function `fmax()` will return the value of f(z), which is the value of the density function of the
#'maximum of two normally distributed random variables.
#' @details  Z = max(X,Y) with X ~ N(mu1,sigma1^2), Y ~ N(mu2,sigma2^2)
#'
#'f(z)=f1(-z)+f2(-z)
#' @param z integral variable
#' @param mu1 mean of second endpoint 
#' @param mu2 mean of first endpoint
#' @param sigma1 standard deviation of first endpoint
#' @param sigma2 standard deviation of second endpoint
#' @param rho correlation between the two endpoints
#' @return The function `fmax()` will return the value of f(z), which is the value of the density function of the
#'maximum of two normally distributed random variables.
#' @examples res <- fmax(z = 0.5, mu1 = 0.375, mu2 = 0.25, sigma1 = 8, sigma2 = 12, rho = 0.4 )
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
fmax <- function(){}

#' @title Probability to go to phase III for multiple endpoints in the time-to-event setting
#' @description This function calculated the probability that we go to phase III, i.e. that results of phase II are promising enough to
#'get a successful drug development program. Successful means that at least one endpoint show a statistically significant positive treatment effect in phase III.
#' @param HRgo threshold value for the go/no-go decision rule; 
#' @param n2 total sample size for phase II; must be even number
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param ec control arm event rate for phase II and III
#' @param id1 amount of information for `hr1` in terms of number of events
#' @param id2 amount of information for `hr2` in terms of number of events
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `hr1` is used as fixed effect
#' @param rho correlation between the two endpoints
#' @return The output of the function `pgo_multiple_tte()` is the probability to go to phase III.
#' @examples res <- pgo_multiple_tte(HRgo = 0.8, n2 = 50, ec = 0.6,
#'                               hr1 = 0.75, hr2 = 0.80, id1 = 300, id2 = 600, 
#'                               fixed = TRUE, rho = 0.3)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
pgo_multiple_tte <- function(){}

#' @title Expected sample size for phase III for multiple endpoints with normally distributed outcomes
#' @description Given phase II results are promising enough to get the "go"-decision to go to phase III this function now calculates the expected sample size for phase III.
#'The results of this function are necessary for calculating the utility of the program, which is then in a further step maximized by the `optimal_multiple_tte()` function
#' @param HRgo threshold value for the go/no-go decision rule; 
#' @param n2 total sample size for phase II; must be even number
#' @param beta `1-beta` power for calculation of the number of events for phase III by Schoenfeld (1981) formula
#' @param alpha one- sided significance level
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param ec control arm event rate for phase II and III
#' @param id1 amount of information for `hr1` in terms of number of events
#' @param id2 amount of information for `hr2` in terms of number of events
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `hr1` is used as fixed effect
#' @param rho correlation between the two endpoints
#' @return the output of the function `Ess_multiple_tte()` is the expected number of participants in phase III
#' @examples res <- Ess_multiple_tte(HRgo = 0.8, n2 = 50, alpha = 0.05, beta = 0.1,
#'                               ec = 0.6,hr1 = 0.75, hr2 = 0.80, 
#'                               id1 = 300, id2 = 600, 
#'                               fixed = TRUE, rho = 0.3)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
Ess_multiple_tte <- function(){}

#' @title Probabilty that effect in endpoint one larger than in endpoint two
#' @description This function calculated the probability that the treatment effect in endpoint one (or endpoint x) is larger than in endpoint two (or endpoint y), i.e. P(x>y) = P(x-y>0)
#' @details Z=X-Y is normally distributed with expectation mu_x - mu_y and variance sigma_x + sigma_y- 2 rho sdx sdy
#' @param n2 total sample size for phase II; must be even number
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param ec control arm event rate for phase II and III
#' @param id1 amount of information for `hr1` in terms of number of events
#' @param id2 amount of information for `hr2` in terms of number of events
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `hr1` is used as fixed effect
#' @param rho correlation between the two endpoints
#' @return The output of the function `pw()` is the probability that endpoint one has a better result than endpoint two
#' @examples res <- pw(n2 = 50, ec = 0.6,
#'                    hr1 = 0.75, hr2 = 0.80, id1 = 300, id2 = 600, 
#'                    fixed = FALSE, rho = 0.3)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
pw <- function(){}

#' @title Expected probability of a successful program for multiple endpoints in a time-to-event setting 
#' @description This function calculates the probability that our drug development program is successful.
#'Successful is defined as at least one endpoint showing a statistically significant positive treatment effect in phase III.
#' @param HRgo threshold value for the go/no-go decision rule; 
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param ec control arm event rate for phase II and III
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param id1 amount of information for `hr1` in terms of sample size
#' @param id2 amount of information for `hr2` in terms of sample size
#' @param step1 lower boundary for effect size
#' @param step2 upper boundary for effect size
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `hr1` is used as fixed effect
#' @param rho correlation between the two endpoints
#' @return The output of the function `EPsProg_multiple_tte()` is the expected probability of a successful program, when going to phase III.
#' @examples res <- EPsProg_multiple_tte(HRgo = 0.8, n2 = 50, alpha = 0.025, beta = 0.1,
#'                               ec = 0.6, hr1 = 0.75, hr2 = 0.80,
#'                               id1 = 300, id2 = 600, 
#'                               step1 = 1, step2 = 0.95,
#'                               fixed = TRUE, rho = 0.3)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
EPsProg_multiple_tte <- function(){}

#' @title Probability that endpoint OS significant
#' @description This function calculate the probability that the endpoint OS is statistically significant.
#'In the context of cancer research OS stands for overall survival, a positive treatment effect in this endpoints is thus sufficient for a successful program.
#' @param HRgo threshold value for the go/no-go decision rule; 
#' @param n2 total sample size for phase II; must be even number
#' @param beta 1-beta power for calculation of the number of events for phase III by Schoenfeld (1981) formula
#' @param alpha one- sided significance level
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param ec control arm event rate for phase II and III
#' @param id1 amount of information for `hr1` in terms of number of events
#' @param id2 amount of information for `hr2` in terms of number of events
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `hr1` is used as fixed effect
#' @param rho correlation between the two endpoints
#' @return The output of the function `os_tte()` is the probability that endpoint OS significant.
#' @examples res <- os_tte(HRgo = 0.8, n2 = 50, alpha = 0.05, beta = 0.1,
#'                               hr1 = 0.75, hr2 = 0.80, 
#'                               id1 = 300, id2 = 600, 
#'                               fixed = TRUE, rho = 0.3)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
os_tte <- function(){}

#' @title Utility function for multiple endpoints in a time-to-event-setting
#' @description The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters and the expected probability of a successful program. 
#'The utility is in a further step maximized by the `optimal_multiple_tte()` function. 
#'Note, that for calculating the utility of the program, two different benefit triples are necessary: 
#' - one triple for the case that the more important endpoint overall survival (OS) shows a significant positive treatment effect 
#' - one triple when only the endpoint progression-free survival (PFS) shows a significant positive treatment effect
#' @param HRgo threshold value for the go/no-go decision rule; 
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param ec control arm event rate for phase II and III
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param id1 amount of information for `hr1` in terms of sample size
#' @param id2 amount of information for `hr2` in terms of sample size
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param steps1 lower boundary for effect size category `"small"` in HR scale, default: 1
#' @param stepm1 lower boundary for effect size category `"medium"` in HR scale = upper boundary for effect size category `"small"` in HR scale, default: 0.95
#' @param stepl1 lower boundary for effect size category `"large"` in HR scale = upper boundary for effect size category `"medium"` in HR scale, default: 0.85
#' @param b11 expected gain for effect size category `"small"` if endpoint OS is significant
#' @param b21 expected gain for effect size category `"medium"`if endpoint OS is significant
#' @param b31 expected gain for effect size category `"large"` if endpoint OS is significant 
#' @param b12 expected gain for effect size category `"small"` if endpoint OS is not significant
#' @param b22 expected gain for effect size category `"medium"`if endpoint OS is not significant
#' @param b32 expected gain for effect size category `"large"` if endpoint OS is not significant
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `hr1` is used as fixed effect
#' @param rho correlation between the two endpoints
#' @return The output of the function `utility_multiple_tte()` is the expected utility of the program.
#' @examples res <- utility_multiple_tte(n2 = 50, HRgo = 0.8, alpha = 0.025, beta = 0.1,
#'                               hr1 = 0.75, hr2 = 0.80,
#'                               id1 = 300, id2 = 600, ec = 0.6,
#'                               c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,               
#'                               K = Inf, N = Inf, S = -Inf, 
#'                               steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
#'                               b11 = 1000, b21 = 2000, b31 = 3000,
#'                               b12 = 1000, b22 = 1500, b32 = 2000, 
#'                               fixed = TRUE, rho = 0.3)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23  
#' @export 
utility_multiple_tte <- function(){}

#' @title Expected probability of a successful program for multitrial programs with binary distributed outcomes
#' @description These functions calculate the expected probability of a successful program given the parameters. 
#'Each function represents a specific strategy, e.g. the function EpsProg3_binary() calculates the expected probability if three phase III trials are performed. 
#'The parameter case specifies how many of the trials have to be successful, i.e. how many trials show a significantly relevant positive treatment effect.
#' @details The following cases can be investigated by the software:
#'- Two phase III trials
#'  -  Case 1: Strategy 1/2; at least one trial significant, the treatment effect of the other one at least showing in the same direction 
#'  - Case 2: Strategy 2/2; both trials significant
#'- Three phase III trials 
#'  - Case 2: Strategy 2/3; at least two trials significant, the treatment effect of the other one at least showing in the same direction
#'  - Case 3: Strategy 3/3; all trials significant 
#'- Four phase III trials 
#'  - Case 3: Strategy 3/4; at least three trials significant, the treatment effect of the other one at least showing in the same direction
#' @param RRgo threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta 1-beta power for calculation of sample size for phase III
#' @param w weight for mixture prior distribution
#' @param p0 assumed true rate of control group
#' @param p11 assumed true rate of treatment group
#' @param p12 assumed true rate of treatment group
#' @param in1 amount of information for p11 in terms of sample size
#' @param in2 amount of information for p12 in terms of sample size
#' @param case choose case: "at least 1, 2 or 3 significant trials needed for approval"
#' @param size size category "small", "medium" or "large"
#' @param fixed choose if true treatment effects are fixed or random
#' @return The output of the function EPsProg2_binary(), EPsProg3_binary() and EPsProg4_binary() is the expected probability of a successful program when performing several phase III trials (2, 3 or 4 respectively)
#' @examples res <- EPsProg2_binary(RRgo = 0.8, n2 = 50,  alpha = 0.025, beta = 0.1, 
#'                                 p0 = 0.6,  w = 0.3, p11 =  0.3, p12 = 0.5, 
#'                                 in1 = 300, in2 = 600, case = 2, size = "small",
#'                                 fixed = FALSE)
#'          res <- EPsProg3_binary(RRgo = 0.8, n2 = 50,  alpha = 0.025, beta = 0.1, 
#'                                 p0 = 0.6,  w = 0.3, p11 =  0.3, p12 = 0.5, 
#'                                 in1 = 300, in2 = 600, case = 2, size = "small",
#'                                 fixed = FALSE)
#'          res <- EPsProg4_binary(RRgo = 0.8, n2 = 50,  alpha = 0.025, beta = 0.1, 
#'                                 p0 = 0.6,  w = 0.3, p11 =  0.3, p12 = 0.5, 
#'                                 in1 = 300, in2 = 600, case = 3, size = "small",
#'                                 fixed = FALSE)
#' @export                                 
#' @name EPsProg_multitrial_binary  
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
EPsProg2_binary <- function(){}

#' @title Utility function for multitrial programs with binary distributed outcomes
#' @description The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters and the expected probability of a successful program. 
#'The utility is in further step maximized by the `optimal_multitrial_binary()` function.
#' @param RRgo threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param w weight for mixture prior distribution
#' @param p0 assumed true rate of control group
#' @param p11 assumed true rate of treatment group
#' @param p12 assumed true rate of treatment group
#' @param in1 amount of information for `p11` in terms of sample size
#' @param in2 amount of information for `p12` in terms of sample size
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param b1 expected gain for effect size category `"small"`
#' @param b2 expected gain for effect size category `"medium"`
#' @param b3 expected gain for effect size category `"large"`
#' @param fixed choose if true treatment effects are fixed or random
#' @param case choose case: "at least 1, 2 or 3 significant trials needed for approval"
#' @return The output of the `functions utility2_binary()`, `utility3_binary()` and `utility4_binary()` is the expected utility of the program when 2, 3 or 4 phase III trials are performed.
#' @examples res <- utility2_binary(n2 = 50, RRgo = 0.8,  w = 0.3, 
#'                                 p0 = 0.6, p11 =  0.3, p12 = 0.5, 
#'                                 in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
#'                                 c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                 K = Inf, N = Inf, S = -Inf,
#'                                 b1 = 1000, b2 = 2000, b3 = 3000, 
#'                                 case = 2, fixed = TRUE)
#'          res <- utility3_binary(n2 = 50, RRgo = 0.8,  w = 0.3, 
#'                                 p0 = 0.6, p11 =  0.3, p12 = 0.5, 
#'                                 in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
#'                                 c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                 K = Inf, N = Inf, S = -Inf,
#'                                 b1 = 1000, b2 = 2000, b3 = 3000, 
#'                                 case = 2, fixed = TRUE)
#'         res <- utility4_binary(n2 = 50, RRgo = 0.8,  w = 0.3, 
#'                                 p0 = 0.6, p11 =  0.3, p12 = 0.5, 
#'                                 in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
#'                                 c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                 K = Inf, N = Inf, S = -Inf,
#'                                 b1 = 1000, b2 = 2000, b3 = 3000, 
#'                                 case = 3, fixed = TRUE)
#' @name utility_multitrial_binary 
#' @export                                
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
utility2_binary <- function(){}

#' @title Expected probability to do third phase III trial
#' @description In the setting of Case 2: Strategy 2/2( + 1); at least two trials significant (and the 
#'treatment effect of the other one at least showing in the same direction) this function calculates the probability that a third phase III trial is necessary.
#' @param RRgo threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param w weight for mixture prior distribution
#' @param p0 assumed true rate of control group
#' @param p11 assumed true rate of treatment group
#' @param p12 assumed true rate of treatment group
#' @param in1 amount of information for `p11` in terms of sample size
#' @param in2 amount of information for `p12` in terms of sample size
#' @return The output of the function `Epgo23_binary()` is the probability to a third phase III trial.
#' @examples res <- Epgo23_binary(RRgo = 0.8, n2 = 50,  p0 = 0.3, w = 0.3, alpha = 0.025, beta = 0.1,
#'                               p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600)
#' @editor Johannes Cepicka
#' @editDate 2022-05-09
#' @export 
Epgo23_binary <- function(){}

#' @title Expected probability of a successful program deciding between two or three phase III trials for a binary distributed outcome
#' @description The function `EPsProg23_binary()` calculates the expected probability of a successful program
#'with a normally distributed outcome. This function follows a special decision rule in order to determine
#'whether two or three phase III trials should be conducted. First, two phase III trials are performed. Depending
#'on their success, the decision for a third phase III trial is made:
#'- If both trials are successful, no third phase III trial will be conducted.
#'- If only one of the two trials is successful and the other trial has a treatment effect that points in the same direction,
#'a third phase III trial will be conducted with a sample size of N3 = N3(ymin), which depends on an assumed minimal clinical relevant effect (`ymin`).
#'The third trial then has to be significant at level `alpha`
#'- If only one of the two trials is successful and the treatment effect of the other points in opposite direction or 
#'if none of the two trials are successful, then no third trial is performed and the drug development development program is not successful. 
#'In the utility function, this will lead to a utility of -9999.
#' @param RRgo threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param w weight for mixture prior distribution
#' @param p0 assumed true rate of control group
#' @param p11 assumed true rate of treatment group
#' @param p12 assumed true rate of treatment group
#' @param in1 amount of information for `p11` in terms of sample size
#' @param in2 amount of information for `p12` in terms of sample size
#' @param case choose case: "at least 1, 2 or 3 significant trials needed for approval"
#' @param size size category `"small"`, `"medium"` or `"large"`
#' @param ymin assumed minimal clinical relevant effect
#' @return The output of the function `EPsProg23_binary()` is the expected probability of a successful program.
#' @examples res <- EPsProg23_binary(RRgo = 0.8, n2 = 50,  alpha = 0.025, beta = 0.1, 
#'                                 w = 0.6,  p0 = 0.3, p11 =  0.3, p12 = 0.5, 
#'                                 in1 = 300, in2 = 600, case = 2, size = "small",
#'                                 ymin = 0.5)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
EPsProg23_binary <- function(){}

#' @title Utility function for multitrial programs deciding between two or three phase III trials for a binary distributed outcome
#' @description The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters and the expected probability of a successful program. 
#'The utility is in further step maximized by the `optimal_multitrial_binary()` function.
#' @param RRgo threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param w weight for mixture prior distribution
#' @param p0 assumed true rate of control group
#' @param p11 assumed true rate of treatment group
#' @param p12 assumed true rate of treatment group
#' @param in1 amount of information for `p11` in terms of sample size
#' @param in2 amount of information for `p12` in terms of sample size
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param b1 expected gain for effect size category `"small"`
#' @param b2 expected gain for effect size category `"medium"`
#' @param b3 expected gain for effect size category `"large"`
#' @return The output of the function `utility23_binary()` is the expected utility of the program depending on whether two or three phase III trials are performed.
#' @examples #res <- utility23_binary(n2 = 50, RRgo = 0.8,  w = 0.3, 
#'          #                       p0 = 0.6, p11 =  0.3, p12 = 0.5, 
#'          #                       in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
#'          #                       c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'          #                       b1 = 1000, b2 = 2000, b3 = 3000)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
utility23_binary <- function(){}

#' @title Expected probability of a successful program for multitrial programs with normally distributed outcomes
#' @description These functions calculate the expected probability of a successful program given the parameters. 
#'Each function represents a specific strategy, e.g. the function `EpsProg3_normal()` calculates the expected probability if three phase III trials are performed. 
#'The parameter case specifies how many of the trials have to be successful, i.e. how many trials show a significantly relevant positive treatment effect.
#' @details The following cases can be investigated by the software:
#'- Two phase III trials
#'  - Case 1: Strategy 1/2; at least one trial significant, the treatment effect of the other one at least showing in the same direction 
#'  - Case 2: Strategy 2/2; both trials significant
#'- Three phase III trials 
#'  - Case 2: Strategy 2/3; at least two trials significant, the treatment effect of the other one at least showing in the same direction
#'  - Case 3: Strategy 3/3; all trials significant 
#'- Four phase III trials 
#'  - Case 3: Strategy 3/4; at least three trials significant, the treatment effect of the other one at least showing in the same direction
#' @param kappa threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param w weight for mixture prior distribution
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param in1 amount of information for `Delta1` in terms of sample size
#' @param in2 amount of information for `Delta2` in terms of sample size
#' @param a lower boundary for the truncation
#' @param b upper boundary for the truncation
#' @param case choose case: "at least 1, 2 or 3 significant trials needed for approval"
#' @param size size category `"small"`, `"medium"` or `"large"`
#' @param fixed choose if true treatment effects are fixed or random
#' @return The output of the function `EPsProg2_normal()`, `EPsProg3_normal()` and `EPsProg4_normal()` is the expected probability of a successful program when performing several phase III trials (2, 3 or 4 respectively).
#' @examples res <- EPsProg2_normal(kappa = 0.1, n2 = 50, alpha = 0.025, beta = 0.1, w = 0.3,
#'                                 Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#'                                 a = 0.25, b = 0.75, 
#'                                 case = 2, size = "small", fixed = FALSE)
#'          res <- EPsProg3_normal(kappa = 0.1, n2 = 50, alpha = 0.025, beta = 0.1, w = 0.3,
#'                                 Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#'                                 a = 0.25, b = 0.75, 
#'                                 case = 2, size = "small", fixed = TRUE)
#'          res <- EPsProg4_normal(kappa = 0.1, n2 = 50, alpha = 0.025, beta = 0.1, w = 0.3,
#'                                 Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#'                                 a = 0.25, b = 0.75, 
#'                                 case = 3, size = "small", fixed = TRUE)                      
#' @name EPsProg_multitrial_normal                                  
#' @export 
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
EPsProg2_normal <- function(){}

#' @title Utility function for multitrial programs with normally distributed outcomes
#' @description The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters and the expected probability of a successful program. 
#'The utility is in a further step maximized by the `optimal_multitrial_normal()` function.
#' @param n2 total sample size for phase II; must be even number
#' @param kappa threshold value for the go/no-go decision rule
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param w weight for mixture prior distribution
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param in1 amount of information for `Delta1` in terms of sample size
#' @param in2 amount of information for `Delta2` in terms of sample size
#' @param a lower boundary for the truncation
#' @param b upper boundary for the truncation
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param b1 expected gain for effect size category `"small"`
#' @param b2 expected gain for effect size category `"medium"`
#' @param b3 expected gain for effect size category `"large"`
#' @param case choose case: "at least 1, 2 or 3 significant trials needed for approval"
#' @param fixed choose if true treatment effects are fixed or random
#' @return The output of the functions utility2_normal(), utility3_normal() and utility4_normal() is the expected utility of the program when 2, 3 or 4 phase III trials are performed.
#' @examples res <- utility2_normal(kappa = 0.1, n2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
#'                                 Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#'                                 a = 0.25, b = 0.75, 
#'                                 c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
#'                                 K = Inf, N = Inf, S = -Inf, 
#'                                 b1 = 3000, b2 = 8000, b3 = 10000, 
#'                                 case = 2, fixed = TRUE)
#' #         res <- utility3_normal(kappa = 0.1, n2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
#' #                                Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#' #                                a = 0.25, b = 0.75, 
#' #                                c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
#' #                                K = Inf, N = Inf, S = -Inf,
#' #                                b1 = 3000, b2 = 8000, b3 = 10000, 
#' #                                case = 2, fixed = TRUE)                        
#' #         res <- utility4_normal(kappa = 0.1, n2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
#' #                                Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#' #                                a = 0.25, b = 0.75, 
#' #                                c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
#' #                                K = Inf, N = Inf, S = -Inf, 
#' #                                b1 = 3000, b2 = 8000, b3 = 10000, 
#' #                                case = 3, fixed = TRUE)
#' @name utility_multitrial_normal                           
#' @export 
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
utility2_normal <- function(){}

#' @title Expected probability to do third phase III trial
#' @description In the setting of Case 2: Strategy 2/2( + 1); at least two trials significant (and the 
#'treatment effect of the other one at least showing in the same direction) this function calculates the probability that a third phase III trial is necessary.
#' @param kappa threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param w weight for mixture prior distribution
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param in1 amount of information for `Delta1` in terms of sample size
#' @param in2 amount of information for `Delta2` in terms of sample size
#' @return The output of the function `Epgo23_normal()` is the probability to a third phase III trial.
#' @examples res <- Epgo23_normal(kappa = 0.1, n2 = 50, w = 0.3, alpha = 0.025, beta = 0.1,
#'                               Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600)
#' @export 
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
Epgo23_normal <- function(){}

#' @title Expected probability of a successful program deciding between two or three phase III trials for a normally distributed outcome
#' @description The function `EPsProg23_normal()` calculates the expected probability of a successful program
#'with a normally distributed outcome. This function follows a special decision rule in order to determine
#'whether two or three phase III trials should be conducted. First, two phase III trials are performed. Depending
#'on their success, the decision for a third phase III trial is made:
#'- If both trials are successful, no third phase III trial will be conducted.
#'- If only one of the two trials is successful and the other trial has a treatment effect that points in the same direction,
#'a third phase III trial will be conducted with a sample size of N3 = N3(ymin), which depends on an assumed minimal clinical relevant effect (`ymin`).
#'The third trial then has to be significant at level `alpha`
#'- If only one of the two trials is successful and the treatment effect of the other points in opposite direction or 
#'if none of the two trials are successful, then no third trial is performed and the drug development development program is not successful. 
#'In the utility function, this will lead to a utility of -9999.
#' @param kappa threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be an even number
#' @param alpha significance level
#' @param beta type II error rate; this means that 1 - `beta` is the power for calculating the sample size for phase III
#' @param w weight for the mixture prior distribution
#' @param Delta1 assumed true treatment effect for the standardized difference in means
#' @param Delta2 assumed true treatment effect for the standardized difference in means
#' @param in1 amount of information for `Delta1` in terms of sample size
#' @param in2 amount of information for `Delta2` in terms of sample size
#' @param a lower boundary for the truncation
#' @param b upper boundary for the truncation
#' @param case number of significant trials needed for approval; possible values are 2 and 3 for this function
#' @param size effect size category; possible values are `"small"`, `"medium"`, `"large"` and `"all"`
#' @param ymin assumed minimal clinical relevant effect
#' @return The output of the function `EPsProg23_normal()` is the expected probability of a successful program. 
#' @examples res <- EPsProg23_normal(kappa = 0.1, n2 = 50, alpha = 0.025, beta = 0.1, w = 0.3,
#'                                 Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#'                                 a = 0.25, b = 0.75, 
#'                                 case = 2, size = "small", ymin = 0.5)
#' @export 
#' @editor Johannes Cepicka, Lukas D. Sauer
#' @editDate 2022-05-06
EPsProg23_normal <- function(){}

#' @title Utility function for multitrial programs deciding between two or three phase III trials for a normally distributed outcome
#' @description The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters and the expected probability of a successful program. 
#'The utility is in a further step maximized by the `optimal_multitrial_normal()` function.
#' @param n2 total sample size for phase II; must be even number
#' @param kappa threshold value for the go/no-go decision rule
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param w weight for mixture prior distribution
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param in1 amount of information for `Delta1` in terms of sample size
#' @param in2 amount of information for `Delta2` in terms of sample size
#' @param a lower boundary for the truncation
#' @param b upper boundary for the truncation
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param b1 expected gain for effect size category `"small"`
#' @param b2 expected gain for effect size category `"medium"`
#' @param b3 expected gain for effect size category `"large"`
#' @return The output of the function utility23_normal() is the expected utility of the program depending on whether two or three phase III trials are performed.
#' @examples #res <- utility23_normal(n2 = 50, kappa = 0.2, w = 0.3,
#'      #                           Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#'      #                           a = 0.25, b = 0.75, 
#'     #                            c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
#'     #                            b1 = 3000, b2 = 8000, b3 = 10000)
#' @export 
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
utility23_normal <- function(){}

#' @title Expected probability of a successful program for multitrial programs in a time-to-event setting
#' @description These functions calculate the expected probability of a successful program given the parameters. 
#'Each function represents a specific strategy, e.g. the function EpsProg3() calculates the expected probability if three phase III trials are performed. 
#'The parameter case specifies how many of the trials have to be successful, i.e. how many trials show a significantly relevant positive treatment effect.
#' @details The following cases can be investigated by the software:
#'- Two phase III trials
#'  -  Case 1: Strategy 1/2; at least one trial significant, the treatment effect of the other one at least showing in the same direction 
#'  - Case 2: Strategy 2/2; both trials significant
#'- Three phase III trials 
#'  - Case 2: Strategy 2/3; at least two trials significant, the treatment effect of the other one at least showing in the same direction
#'  - Case 3: Strategy 3/3; all trials significant 
#'- Four phase III trials 
#'  - Case 3: Strategy 3/4; at least three trials significant, the treatment effect of the other one at least showing in the same direction
#' @param HRgo threshold value for the go/no-go decision rule
#' @param d2 total number of events in phase II
#' @param alpha significance level
#' @param beta 1-beta power for calculation of sample size for phase III
#' @param w weight for mixture prior distribution
#' @param hr1 first assumed true treatment effect on HR scale for prior distribution
#' @param hr2 second assumed true treatment effect on HR scale for prior distribution
#' @param id1 amount of information for hr1 in terms of number of events
#' @param id2 amount of information for hr2 in terms of number of events
#' @param case choose case: "at least 1, 2 or 3 significant trials needed for approval"
#' @param size size category "small", "medium" or "large"
#' @param fixed choose if true treatment effects are fixed or random
#' @return The output of the function EPsProg2(), EPsProg3() and EPsProg4() is the expected probability of a successful program when performing several phase III trials (2, 3 or 4 respectively)
#' @examples res <- EPsProg2(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, 
#'                                 w = 0.3, hr1 =  0.69, hr2 = 0.81, 
#'                                 id1 = 210, id2 = 420, case = 2, size = "small",
#'                                 fixed = FALSE)
#'          res <- EPsProg3(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, 
#'                                 w = 0.3, hr1 =  0.69, hr2 = 0.81, 
#'                                 id1 = 210, id2 = 420, case = 2, size = "small",
#'                                 fixed = TRUE)
#'          res <- EPsProg4(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, 
#'                                 w = 0.3, hr1 =  0.69, hr2 = 0.81, 
#'                                 id1 = 210, id2 = 420, case = 3, size = "small",
#'                                 fixed = TRUE)
#' @name EPsProg_multitrial                                 
#' @export 
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
EPsProg2 <- function(){}

#' @title Utility function for multitrial programs in a time-to-event setting
#' @description The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters and the expected probability of a successful program. 
#'The utility is in further step maximized by the `optimal_multitrial()` function.
#' @param HRgo threshold value for the go/no-go decision rule
#' @param d2 total number of events in phase II
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param xi2 event rate for phase II
#' @param xi3 event rate for phase III
#' @param w weight for mixture prior distribution
#' @param hr1 first assumed true treatment effect on HR scale for prior distribution
#' @param hr2 second assumed true treatment effect on HR scale for prior distribution
#' @param id1 amount of information for `hr1` in terms of number of events
#' @param id2 amount of information for `hr2` in terms of number of events
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param b1 expected gain for effect size category `"small"`
#' @param b2 expected gain for effect size category `"medium"`
#' @param b3 expected gain for effect size category `"large"`
#' @param case choose case: "at least 1, 2 or 3 significant trials needed for approval"
#' @param fixed choose if true treatment effects are fixed or random
#' @return The output of the functions `utility2()`, `utility3()` and `utility4()` is the expected utility of the program when 2, 3 or 4 phase III trials are performed.
#' @examples res <- utility2(d2 = 50, HRgo = 0.8,  w = 0.3, 
#'                                 hr1 =  0.69, hr2 = 0.81, 
#'                                 id1 = 210, id2 = 420, 
#'                                 alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,
#'                                 c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                 K = Inf, N = Inf, S = -Inf,
#'                                 b1 = 1000, b2 = 2000, b3 = 3000, 
#'                                 case = 2, fixed = TRUE)
#'          res <- utility3(d2 = 50, HRgo = 0.8,  w = 0.3, 
#'                                 hr1 =  0.69, hr2 = 0.81, 
#'                                 id1 = 210, id2 = 420, 
#'                                 alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,
#'                                 c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                 K = Inf, N = Inf, S = -Inf,
#'                                 b1 = 1000, b2 = 2000, b3 = 3000, 
#'                                 case = 2, fixed = TRUE)
#'         res <- utility4(d2 = 50, HRgo = 0.8,  w = 0.3, 
#'                                 hr1 =  0.69, hr2 = 0.81, 
#'                                 id1 = 210, id2 = 420, 
#'                                 alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,
#'                                 c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                 K = Inf, N = Inf, S = -Inf,
#'                                 b1 = 1000, b2 = 2000, b3 = 3000, 
#'                                 case = 3, fixed = TRUE)
#' @name utility_multitrial 
#' @export                                
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
utility2 <- function(){}

#' @title Expected probability to do third phase III trial
#' @description In the setting of Case 2: Strategy 2/2( + 1); at least two trials significant (and the 
#'treatment effect of the other one at least showing in the same direction) this function calculates the probability that a third phase III trial is necessary.
#' @param HRgo threshold value for the go/no-go decision rule
#' @param d2 total number of events in phase II
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param w weight for mixture prior distribution
#' @param hr1 first assumed true treatment effect on HR scale for prior distribution
#' @param hr2 second assumed true treatment effect on HR scale for prior distribution
#' @param id1 amount of information for `hr1` in terms of number of events
#' @param id2 amount of information for `hr2` in terms of number of events
#' @return The output of the function `Epgo23()` is the probability to a third phase III trial.
#' @examples res <- Epgo23(HRgo = 0.8, d2 = 50,  w = 0.3, alpha = 0.025, beta = 0.1,
#'                               hr1 =  0.69, hr2 = 0.81, id1 = 280, id2 = 420)
#' @export 
#' @editor Johannes Cepicka
#' @editDate 2022-05-09
Epgo23 <- function(){}

#' @title Expected probability of a successful program deciding between two or three phase III trials in a time-to-event setting
#' @description The function `EPsProg23()` calculates the expected probability of a successful program
#'in a time-to-event setting. This function follows a special decision rule in order to determine
#'whether two or three phase III trials should be conducted. First, two phase III trials are performed. Depending
#'on their success, the decision for a third phase III trial is made:
#'- If both trials are successful, no third phase III trial will be conducted.
#'- If only one of the two trials is successful and the other trial has a treatment effect that points in the same direction,
#'a third phase III trial will be conducted with a sample size of N3 = N3(ymin), which depends on an assumed minimal clinical relevant effect (`ymin`).
#'The third trial then has to be significant at level `alpha`
#'- If only one of the two trials is successful and the treatment effect of the other points in opposite direction or 
#'if none of the two trials are successful, then no third trial is performed and the drug development development program is not successful. 
#'In the utility function, this will lead to a utility of -9999.
#' @param HRgo threshold value for the go/no-go decision rule
#' @param d2 total number of events in phase II
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param w weight for mixture prior distribution
#' @param hr1 first assumed true treatment effect on HR scale for prior distribution
#' @param hr2 second assumed true treatment effect on HR scale for prior distribution
#' @param id1 amount of information for `hr1` in terms of number of events
#' @param id2 amount of information for `hr2` in terms of number of events
#' @param case choose case: "at least 1, 2 or 3 significant trials needed for approval"
#' @param size size category `"small"`, `"medium"` or `"large"`
#' @param ymin assumed minimal clinical relevant effect
#' @return The output of the function `EPsProg23()` is the expected probability of a successful program.
#' @examples res <- EPsProg23(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, 
#'                                  w = 0.3, hr1 =  0.69, hr2 = 0.81, 
#'                                  id1 = 280, id2 = 420, case = 2, size = "small",
#'                                  ymin = 0.5)
#' @export 
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
EPsProg23 <- function(){}

#' @title Utility function for multitrial programs deciding between two or three phase III trials in a time-to-event setting
#' @description The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters and the expected probability of a successful program. 
#'The utility is in further step maximized by the `optimal_multitrial()` function.
#' @param HRgo threshold value for the go/no-go decision rule
#' @param d2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param xi2 event rate for phase II
#' @param xi3 event rate for phase III
#' @param w weight for mixture prior distribution
#' @param hr1 first assumed true treatment effect on HR scale for prior distribution
#' @param hr2 second assumed true treatment effect on HR scale for prior distribution
#' @param id1 amount of information for `hr1` in terms of number of events
#' @param id2 amount of information for `hr2` in terms of number of events
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param b1 expected gain for effect size category `"small"`
#' @param b2 expected gain for effect size category `"medium"`
#' @param b3 expected gain for effect size category `"large"`
#' @importFrom mvtnorm pmvnorm
#' @return The output of the function `utility23()` is the expected utility of the program depending on whether two or three phase III trials are performed.
#' @examples #res <- utility23(d2 = 50, HRgo = 0.8,  w = 0.3, 
#'  #                               hr1 =  0.69, hr2 = 0.81, 
#'  #                              id1 = 280, id2 = 420, 
#'  #                               alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,
#'  #                               c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'  #                               b1 = 1000, b2 = 2000, b3 = 3000)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
utility23 <- function(){}

#' @title Optimal phase II/III drug development planning when discounting phase II results with binary endpoint
#' @description The drugdevelopR package enables planning of phase II/III drug development programs with optimal sample size allocation and go/no-go decision rules including methods for discounting of phase II results (Preussler et. al, 2020).. For binary endpoints the treatment effect is measured by the risk ratio (RR).The assumed true treatment effects can be assumed fixed or modelled by a prior distribution. The R Shiny application \href{https://web.imbi.uni-heidelberg.de/prior/}{prior} visualizes the prior distributions used in this package. Fast coputing is enabled by parallel programming.
#' @name optimal_bias_binary
#' @param w weight for mixture prior distribution
#' @param p0 assumed true rate of control group
#' @param p11 assumed true rate of treatment group
#' @param p12 assumed true rate of treatment group
#' @param in1 amount of information for p11 in terms of sample size
#' @param in2 amount of information for p12 in terms of sample size
#' @param n2min minimal total sample size for phase II; must be even number
#' @param n2max maximal total sample size for phase II, must be even number
#' @param stepn2 stepsize for the optimization over n2; must be even number
#' @param rrgomin minimal threshold value for the go/no-go decision rule
#' @param rrgomax maximal threshold value for the go/no-go decision rule
#' @param steprrgo stepsize for the optimization over RRgo
#' @param lambdamin minimal adjustment parameter lambda
#' @param lambdamax maximal adjustment parameter lambda
#' @param steplambda stepsize for the adjustment parameter lambda
#' @param alphaCImin minimal alphaCI
#' @param alphaCImax maximal alphaCI
#' @param stepalphaCI stepsize for alphaCI
#' @param adj type of adjustment, `"additive"`, `"multiplicative"` or `"both"`
#' @param beta 1-beta power for calculation of sample size for phase III
#' @param alpha significance level
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param steps1 lower boundary for effect size category "small" in RR scale, default: 1
#' @param stepm1 lower boundary for effect size category "medium" in RR scale = upper boundary for effect size category "small" in RR scale, default: 0.95
#' @param stepl1 lower boundary for effect size category "large" in RR scale = upper boundary for effect size category "medium" in RR scale, default: 0.85
#' @param b1 expected gain for effect size category "small"
#' @param b2 expected gain for effect size category "medium"
#' @param b3 expected gain for effect size category "large"
#' @param fixed choose if true treatment effects are fixed or random, if TRUE p11 is used as fixed effect for p1
#' @param num_cl number of clusters used for parallel computing, default: 1
#' @return 
#'The output of the function \code{\link{optimal_binary}} is a data.frame containing the optimization results:
#'\describe{
#'  \item{Method}{Type of adjustment: multipl. (multiplicative) or add. (additive)}
#'  \item{u}{maximal expected utility}
#'  \item{Adj}{optimal adjustment parameter (lambda or alphaCI according to Method)}
#'  \item{RRgo}{optimal threshold value for the decision rule to go to phase III}
#'  \item{n2}{total sample size for phase II}
#'  \item{n3}{total sample size for phase III; rounded to the next even natural number}
#'  \item{n}{total sample size in the program; n = n2 + n3}
#'  \item{K}{maximal costs of the program}
#'  \item{pgo}{probability to go to phase III}
#'  \item{sProg}{probability of a successful program}
#'  \item{sProg1}{probability of a successful program with "small" treatment effect in Phase III}
#'  \item{sProg2}{probability of a successful program with "medium" treatment effect in Phase III}
#'  \item{sProg3}{probability of a successful program with "large" treatment effect in Phase III }
#'  \item{K2}{expected costs for phase II}
#'  \item{K3}{expected costs for phase III}
#'}
#'and further input parameters.
#'
#'Taking cat(comment()) of the data.frame object lists the used optimization sequences, start and finish date of the optimization procedure.
#' @examples 
#'res <- optimal_bias_binary(w = 0.3,                           # define parameters for prior
#'  p0 = 0.6, p11 =  0.3, p12 = 0.5, in1 = 30, in2 = 60,   # (https://web.imbi.uni-heidelberg.de/prior/)
#'  n2min = 20, n2max = 100, stepn2 = 4,                   # define optimization set for n2
#'  rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.05,         # define optimization set for RRgo
#'  adj = "both",                                          # choose type of adjustment
#'  alpha = 0.05, beta = 0.1,                              # drug development planning parameters
#'  lambdamin = 0.2, lambdamax = 1, steplambda = 0.05,     # define optimization set for lambda
#'  alphaCImin = 0.025, alphaCImax = 0.5, stepalphaCI = 0.025, # define optimization set for alphaCI
#'  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,               # define fixed and variable costs for phase II and III,
#'  K = Inf, N = Inf, S = -Inf,                            # set maximal costs/ expected sample size for the program or minimal expected probability of a successful program
#'  steps1 = 1,                                            # define lower boundary for "small"
#'  stepm1 = 0.95,                                         # "medium"
#'  stepl1 = 0.85,                                         # and "large" treatment effect size categories as proposed by IQWiG (2016)
#'  b1 = 1000, b2 = 2000, b3 = 3000,                       # define expected benefit for a "small", "medium" and "large" treatment effect
#'  fixed = FALSE,                                         # choose if true treatment effects are fixed or random
#'  num_cl = 1)                                            # set number of cores used for parallelized computing (check maximum number possible with detectCores())
#'res
#'cat(comment(res))                                        # displays the optimization sequence, start and finish date of the optimization procedure.
#' @section drugdevelopR functions:
#'The drugdevelopR package provides the functions
#'\itemize{
#'  \item \code{\link{optimal_tte}},
#'  \item \code{\link{optimal_binary}} or
#'  \item \code{\link{optimal_normal}}
#'}
#'to plan optimal phase II/III drug development programs with
#'\itemize{
#'  \item time-to-event (treatment effect measured by hazard ratio (HR)),
#'  \item binary (treatment effect measured by risk ratio (RR)) and
#'  \item normally distributed (treatment effect measured by standardized difference in means (Delta))
#'}
#'endpoint, where the treatment effect is modelled by a \href{https://web.imbi.uni-heidelberg.de/prior/}{prior}. Optimal phase II/III drug development planning with fixed treatment effects can be done with the help of the R Shiny application \href{https://web.imbi.uni-heidelberg.de/basic/}{basic}. Extensions are 
#'\itemize{
#'  \item optimal planning of programs including methods for discounting of phase II results (function: \code{\link{optimal_bias}}, App: \href{https://web.imbi.uni-heidelberg.de/bias/}{bias}),
#'  \item optimal planning of programs with several phase III trials (function: \code{\link{optimal_multitrial}}, App: \href{https://web.imbi.uni-heidelberg.de/multitrial/}{multitrial}) and
#'  \item optimal planning of programs with multiple arms (function: \code{\link{optimal_multiarm}}, App: \href{https://web.imbi.uni-heidelberg.de/multiarm/}{multiarm}).
#'}
#' @references 
#'IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at \href{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}, assessed last 15.05.19.
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
optimal_bias_binary <- function(){}

#' @title Optimal phase II/III drug development planning with normally distributed endpoint
#' @description The drugdevelopR package enables planning of phase II/III drug development programs with optimal sample size allocation and go/no-go decision rules including methods for discounting of phase II results (Preussler et. al, 2020). For normally distributed endpoints the treatment effect is measured by the standardized difference in means (Delta). The assumed true treatment effects can be assumed fixed or modelled by a prior distribution. The R Shiny application \href{https://web.imbi.uni-heidelberg.de/prior/}{prior} visualizes the prior distributions used in this package. Fast coputing is enabled by parallel programming.
#' @name optimal_bias_normal
#' @param w weight for mixture prior distribution
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param in1 amount of information for Delta1 in terms of sample size
#' @param in2 amount of information for Delta2 in terms of sample size
#' @param a lower boundary for the truncation
#' @param b upper boundary for the truncation
#' @param n2min minimal total sample size for phase II; must be even number
#' @param n2max maximal total sample size for phase II, must be even number
#' @param stepn2 stepsize for the optimization over n2; must be even number
#' @param kappamin minimal threshold value for the go/no-go decision rule
#' @param kappamax maximal threshold value for the go/no-go decision rule
#' @param adj choose type of adjustment: "multiplicative", "additive", "both" (or "all")
#' @param stepkappa stepsize for the optimization over kappa
#' @param lambdamin minimal adjustment parameter lambda
#' @param lambdamax maximal adjustment parameter lambda
#' @param steplambda stepsize for the adjustment parameter lambda
#' @param alphaCImin minimal alphaCI
#' @param alphaCImax maximal alphaCI
#' @param stepalphaCI stepsize for alphaCI
#' @param beta 1-beta power for calculation of sample size for phase III
#' @param alpha significance level
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param steps1 lower boundary for effect size category "small", default: 0
#' @param stepm1 lower boundary for effect size category "medium" = upper boundary for effect size category "small" default: 0.5
#' @param stepl1 lower boundary for effect size category "large" = upper boundary for effect size category "medium", default: 0.8
#' @param b1 expected gain for effect size category "small"
#' @param b2 expected gain for effect size category "medium"
#' @param b3 expected gain for effect size category "large"
#' @param fixed choose if true treatment effects are fixed or random, if TRUE hr1 is used as fixed effect
#' @param num_cl number of clusters used for parallel computing, default: 1
#' @return 
#'The output of the function \code{\link{optimal_bias_normal}} is a data.frame containing the optimization results:
#'\describe{
#'  \item{Method}{Type of adjustment: multipl. (multiplicative) or add. (additive)}
#'  \item{u}{maximal expected utility}#'   
#'  \item{Adj}{optimal adjustment parameter (lambda or alphaCI according to Method)}
#'  \item{kappa}{optimal threshold value for the decision rule to go to phase III}
#'  \item{n2}{total sample size for phase II}
#'  \item{n3}{total sample size for phase III; rounded to the next even natural number}
#'  \item{n}{total sample size in the program; n = n2 + n3}
#'  \item{K}{maximal costs of the program}
#'  \item{pgo}{probability to go to phase III}
#'  \item{sProg}{probability of a successful program}
#'  \item{sProg1}{probability of a successful program with "small" treatment effect in Phase III}
#'  \item{sProg2}{probability of a successful program with "medium" treatment effect in Phase III}
#'  \item{sProg3}{probability of a successful program with "large" treatment effect in Phase III }
#'  \item{K2}{expected costs for phase II}
#'  \item{K3}{expected costs for phase III}
#'}
#'and further input parameters.
#'
#'Taking cat(comment()) of the data.frame object lists the used optimization sequences, start and finish date of the optimization procedure.
#' @examples 
#'res <- optimal_bias_normal(w=0.3,                                 # define parameters for prior
#'  Delta1 = 0.375, Delta2 = 0.625, in1=300, in2=600,          # (https://web.imbi.uni-heidelberg.de/prior/)
#'  a = 0.25, b = 0.75,
#'  n2min = 20, n2max = 100, stepn2 = 4,                       # define optimization set for n2
#'  kappamin = 0.02, kappamax = 0.2, stepkappa = 0.02,         # define optimization set for kappa
#'  adj = "both",                                              # choose type of adjustment
#'  lambdamin = 0.2, lambdamax = 1, steplambda = 0.05,         # define optimization set for lambda
#'  alphaCImin = 0.025, alphaCImax = 0.5, stepalphaCI = 0.025, # define optimization set for alphaCI
#'  alpha = 0.05, beta = 0.1,                                  # drug development planning parameters
#'  c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,                 # define fixed and variable costs for phase II and III
#'  K = Inf, N = Inf, S = -Inf,                                # set maximal costs/ expected sample size for the program or minimal expected probability of a successful program
#'  steps1 = 0,                                                # define lower boundary for "small"
#'  stepm1 = 0.5,                                              # "medium"
#'  stepl1 = 0.8,                                              # and "large" treatment effect size categories as proposed by e.g. Cohen (1988)
#'  b1 = 3000, b2 = 8000, b3 = 10000,                          # define expected benefit for a "small", "medium" and "large" treatment effect
#'  fixed = FALSE,                                             # choose if true treatment effects are fixed or random
#'  num_cl = 1)                                                # set number of cores used for parallelized computing (check maximum number possible with detectCores())
#'res
#'cat(comment(res))                                        # displays the optimization sequence, start and finish date of the optimization procedure.
#' @section drugdevelopR functions:
#'The drugdevelopR package provides the functions
#'\itemize{
#'  \item \code{\link{optimal_tte}},
#'  \item \code{\link{optimal_binary}} or
#'  \item \code{\link{optimal_normal}}
#'}
#'to plan optimal phase II/III drug development programs with
#'\itemize{
#'  \item time-to-event (treatment effect measured by hazard ratio (HR)),
#'  \item binary (treatment effect measured by risk ratio (RR)) and
#'  \item normally distributed (treatment effect measured by standardized difference in means (Delta))
#'}
#'endpoint, where the treatment effect is assumed fixed or modelled by a \href{https://web.imbi.uni-heidelberg.de/prior/}{prior}. Optimal phase II/III drug development planning with fixed treatment effects can also be done with the help of the R Shiny application \href{https://web.imbi.uni-heidelberg.de/basic/}{basic}. Extensions are 
#'\itemize{
#'  \item optimal planning of programs including methods for discounting of phase II results (function: \code{\link{optimal_bias}}, App: \href{https://web.imbi.uni-heidelberg.de/bias/}{bias}),
#'  \item optimal planning of programs with several phase III trials (function: \code{\link{optimal_multitrial}}, App: \href{https://web.imbi.uni-heidelberg.de/multitrial/}{multitrial}) and
#'  \item optimal planning of programs with multiple arms (function: \code{\link{optimal_multiarm}}, App: \href{https://web.imbi.uni-heidelberg.de/multiarm/}{multiarm}).
#'}
#' @references 
#'Cohen, J. (1988). Statistical power analysis for the behavioral sciences.
#'
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
optimal_bias_normal <- function(){}

#' @title Optimal phase II/III drug development planning with binary endpoint
#' @description The drugdevelopR package enables planning of phase II/III drug development programs with optimal sample size allocation and go/no-go decision rules. For binary endpoints the treatment effect is measured by the risk ratio (RR).The assumed true treatment effects can be assumed fixed or modelled by a prior distribution. The R Shiny application \href{https://web.imbi.uni-heidelberg.de/prior/}{prior} visualizes the prior distributions used in this package. Fast coputing is enabled by parallel programming.
#' @name optimal_multiarm_binary
#' @param p0 assumed true rate of control group
#' @param p11 assumed true rate of treatment group
#' @param p12 assumed true rate of treatment group
#' @param n2min minimal total sample size for phase II; must be even number
#' @param n2max maximal total sample size for phase II, must be even number
#' @param stepn2 stepsize for the optimization over n2; must be even number
#' @param rrgomin minimal threshold value for the go/no-go decision rule
#' @param rrgomax maximal threshold value for the go/no-go decision rule
#' @param steprrgo stepsize for the optimization over RRgo
#' @param beta 1-beta power for calculation of sample size for phase III
#' @param alpha significance level
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param steps1 lower boundary for effect size category "small" in RR scale, default: 1
#' @param stepm1 lower boundary for effect size category "medium" in RR scale = upper boundary for effect size category "small" in RR scale, default: 0.95
#' @param stepl1 lower boundary for effect size category "large" in RR scale = upper boundary for effect size category "medium" in RR scale, default: 0.85
#' @param b1 expected gain for effect size category "small"
#' @param b2 expected gain for effect size category "medium"
#' @param b3 expected gain for effect size category "large"
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)
#' @param num_cl number of clusters used for parallel computing, default: 1
#' @return 
#'The output of the function \code{\link{optimal_binary}} is a data.frame containing the optimization results:
#'\describe{
#'  \item{u}{maximal expected utility}
#'  \item{RRgo}{optimal threshold value for the decision rule to go to phase III}
#'  \item{n2}{total sample size for phase II}
#'  \item{n3}{total sample size for phase III; rounded to the next even natural number}
#'  \item{n}{total sample size in the program; n = n2 + n3}
#'  \item{K}{maximal costs of the program}
#'  \item{pgo}{probability to go to phase III}
#'  \item{sProg}{probability of a successful program}
#'  \item{sProg1}{probability of a successful program with "small" treatment effect in Phase III}
#'  \item{sProg2}{probability of a successful program with "medium" treatment effect in Phase III}
#'  \item{sProg3}{probability of a successful program with "large" treatment effect in Phase III }
#'  \item{K2}{expected costs for phase II}
#'  \item{K3}{expected costs for phase III}
#'}
#'and further input parameters.
#'
#'Taking cat(comment()) of the data.frame object lists the used optimization sequences, start and finish date of the optimization procedure.
#' @examples 
#'res <- optimal_multiarm_binary( p0 = 0.6, p11 =  0.3, p12 = 0.5, 
#'  n2min = 20, n2max = 100, stepn2 = 4,                   # define optimization set for n2
#'  rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.05,         # define optimization set for RRgo
#'  alpha = 0.05, beta = 0.1,                              # drug development planning parameters
#'  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,               # define fixed and variable costs for phase II and III,
#'  K = Inf, N = Inf, S = -Inf,                            # set maximal costs/ expected sample size for the program or minimal expected probability of a successful program
#'  steps1 = 1,                                            # define lower boundary for "small"
#'  stepm1 = 0.95,                                         # "medium"
#'  stepl1 = 0.85,                                         # and "large" treatment effect size categories as proposed by IQWiG (2016)
#'  b1 = 1000, b2 = 2000, b3 = 3000,                       # define expected benefit for a "small", "medium" and "large" treatment effect
#'  strategy = 1, num_cl = 1)                              # set number of cores used for parallelized computing (check maximum number possible with detectCores())
#'res
#'cat(comment(res))                                        # displays the optimization sequence, start and finish date of the optimization procedure.
#' @section drugdevelopR functions:
#'The drugdevelopR package provides the functions
#'\itemize{
#'  \item \code{\link{optimal_tte}},
#'  \item \code{\link{optimal_binary}} or
#'  \item \code{\link{optimal_normal}}
#'}
#'to plan optimal phase II/III drug development programs with
#'\itemize{
#'  \item time-to-event (treatment effect measured by hazard ratio (HR)),
#'  \item binary (treatment effect measured by risk ratio (RR)) and
#'  \item normally distributed (treatment effect measured by standardized difference in means (Delta))
#'}
#'endpoint, where the treatment effect is modelled by a \href{https://web.imbi.uni-heidelberg.de/prior/}{prior}. Optimal phase II/III drug development planning with fixed treatment effects can be done with the help of the R Shiny application \href{https://web.imbi.uni-heidelberg.de/basic/}{basic}. Extensions are 
#'\itemize{
#'  \item optimal planning of programs including methods for discounting of phase II results (function: \code{\link{optimal_bias}}, App: \href{https://web.imbi.uni-heidelberg.de/bias/}{bias}),
#'  \item optimal planning of programs with several phase III trials (function: \code{\link{optimal_multitrial}}, App: \href{https://web.imbi.uni-heidelberg.de/multitrial/}{multitrial}) and
#'  \item optimal planning of programs with multiple arms (function: \code{\link{optimal_multiarm}}, App: \href{https://web.imbi.uni-heidelberg.de/multiarm/}{multiarm}).
#'}
#' @references 
#'IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at \href{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}, assessed last 15.05.19.
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
optimal_multiarm_binary <- function(){}

#' @title Optimal phase II/III drug development planning with normally distributed endpoint
#' @description The drugdevelopR package enables planning of phase II/III drug development programs with optimal sample size allocation and go/no-go decision rules. For normally distributed endpoints the treatment effect is measured by the standardized difference in means (Delta). The assumed true treatment effects can be assumed fixed or modelled by a prior distribution. The R Shiny application \href{https://web.imbi.uni-heidelberg.de/prior/}{prior} visualizes the prior distributions used in this package. Fast coputing is enabled by parallel programming.
#' @name optimal_multiarm_normal
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param n2min minimal total sample size for phase II; must be even number
#' @param n2max maximal total sample size for phase II, must be even number
#' @param stepn2 stepsize for the optimization over n2; must be even number
#' @param kappamin minimal threshold value for the go/no-go decision rule
#' @param kappamax maximal threshold value for the go/no-go decision rule
#' @param stepkappa stepsize for the optimization over kappa
#' @param beta 1-beta power for calculation of sample size for phase III
#' @param alpha significance level
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param steps1 lower boundary for effect size category "small", default: 0
#' @param stepm1 lower boundary for effect size category "medium" = upper boundary for effect size category "small" default: 0.5
#' @param stepl1 lower boundary for effect size category "large" = upper boundary for effect size category "medium", default: 0.8
#' @param b1 expected gain for effect size category "small"
#' @param b2 expected gain for effect size category "medium"
#' @param b3 expected gain for effect size category "large"
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)
#' @param num_cl number of clusters used for parallel computing, default: 1
#' @return 
#'The output of the function \code{\link{optimal_normal}} is a data.frame containing the optimization results:
#'\describe{
#'  \item{u}{maximal expected utility}
#'  \item{kappa}{optimal threshold value for the decision rule to go to phase III}
#'  \item{n2}{total sample size for phase II}
#'  \item{n3}{total sample size for phase III; rounded to the next even natural number}
#'  \item{n}{total sample size in the program; n = n2 + n3}
#'  \item{K}{maximal costs of the program}
#'  \item{pgo}{probability to go to phase III}
#'  \item{sProg}{probability of a successful program}
#'  \item{sProg1}{probability of a successful program with "small" treatment effect in Phase III}
#'  \item{sProg2}{probability of a successful program with "medium" treatment effect in Phase III}
#'  \item{sProg3}{probability of a successful program with "large" treatment effect in Phase III }
#'  \item{K2}{expected costs for phase II}
#'  \item{K3}{expected costs for phase III}
#'}
#'and further input parameters.
#'
#'Taking cat(comment()) of the data.frame object lists the used optimization sequences, start and finish date of the optimization procedure.
#' @examples 
#'res <- optimal_multiarm_normal(Delta1 = 0.375, Delta2 = 0.625,     
#'  n2min = 20, n2max = 100, stepn2 = 4,                   # define optimization set for n2
#'  kappamin = 0.02, kappamax = 0.2, stepkappa = 0.02,     # define optimization set for kappa
#'  alpha = 0.05, beta = 0.1,                              # drug development planning parameters
#'  c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,             # define fixed and variable costs for phase II and III
#'  K = Inf, N = Inf, S = -Inf,                            # set maximal costs/ expected sample size for the program or minimal expected probability of a successful program
#'  steps1 = 0,                                            # define lower boundary for "small"
#'  stepm1 = 0.5,                                          # "medium"
#'  stepl1 = 0.8,                                          # and "large" treatment effect size categories as proposed by e.g. Cohen (1988)
#'  b1 = 3000, b2 = 8000, b3 = 10000,                      # define expected benefit for a "small", "medium" and "large" treatment effect
#'  strategy = 1,
#'  num_cl = 1)                                            # set number of cores used for parallelized computing (check maximum number possible with detectCores())
#'res
#'cat(comment(res))                                        # displays the optimization sequence, start and finish date of the optimization procedure.
#' @section drugdevelopR functions:
#'The drugdevelopR package provides the functions
#'\itemize{
#'  \item \code{\link{optimal_tte}},
#'  \item \code{\link{optimal_binary}} or
#'  \item \code{\link{optimal_normal}}
#'}
#'to plan optimal phase II/III drug development programs with
#'\itemize{
#'  \item time-to-event (treatment effect measured by hazard ratio (HR)),
#'  \item binary (treatment effect measured by risk ratio (RR)) and
#'  \item normally distributed (treatment effect measured by standardized difference in means (Delta))
#'}
#'endpoint, where the treatment effect is assumed fixed or modelled by a \href{https://web.imbi.uni-heidelberg.de/prior/}{prior}. Optimal phase II/III drug development planning with fixed treatment effects can also be done with the help of the R Shiny application \href{https://web.imbi.uni-heidelberg.de/basic/}{basic}. Extensions are 
#'\itemize{
#'  \item optimal planning of programs including methods for discounting of phase II results (function: \code{\link{optimal_bias}}, App: \href{https://web.imbi.uni-heidelberg.de/bias/}{bias}),
#'  \item optimal planning of programs with several phase III trials (function: \code{\link{optimal_multitrial}}, App: \href{https://web.imbi.uni-heidelberg.de/multitrial/}{multitrial}) and
#'  \item optimal planning of programs with multiple arms (function: \code{\link{optimal_multiarm}}, App: \href{https://web.imbi.uni-heidelberg.de/multiarm/}{multiarm}).
#'}
#' @references 
#'Cohen, J. (1988). Statistical power analysis for the behavioral sciences.
#'
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
optimal_multiarm_normal <- function(){}

#' @title Optimal phase II/III drug development planning for programs with multiple endpoints
#' @description The function \code{\link{optimal_multiple_normal}} of the drugdevelopR package enables planning of phase II/III drug development programs with optimal sample size allocation and go/no-go decision rules (Preussler et. al, 2019).  (planning is also possible via user friendly R Shiny App: \href{https://web.imbi.uni-heidelberg.de/multiarm/}{multiarm}). Fast coputing is enabled by parallel programming.
#' @name optimal_multiple_normal
#' @param Delta1 assumed true treatment effect on HR scale for treatment 1
#' @param Delta2 assumed true treatment effect on HR scale for treatment 2
#' @param in1 amount of information for Delta1 in terms of number of events
#' @param in2 amount of information for Delta2 in terms of number of events
#' @param sigma1 variance of endpoint 1
#' @param sigma2 variance of endpoint 1
#' @param rho correlation between the two endpoints
#' @param fixed assumed fixed treatment effect 
#' @param relaxed relaxed or strict decision rule 
#' @param n2min minimal total sample size in phase II, must be divisible by 3
#' @param n2max maximal total sample size in phase II, must be divisible by 3
#' @param stepn2 stepsize for the optimization over n2, must be divisible by 3
#' @param kappamin minimal threshold value for the go/no-go decision rule
#' @param kappamax maximal threshold value for the go/no-go decision rule
#' @param stepkappa stepsize for the optimization over HRgo
#' @param beta 1-beta (any-pair) power for calculation of the number of events for phase III
#' @param alpha one-sided significance level/ family wise error rate
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param steps1 lower boundary for effect size category "small" in HR scale, default: 1
#' @param stepm1 lower boundary for effect size category "medium" in HR scale = upper boundary for effect size category "small" in HR scale, default: 0.95
#' @param stepl1 lower boundary for effect size category "large" in HR scale = upper boundary for effect size category "medium" in HR scale, default: 0.85
#' @param b1 expected gain for effect size category "small"
#' @param b2 expected gain for effect size category "medium"
#' @param b3 expected gain for effect size category "large"
#' @param num_cl number of clusters used for parallel computing, default: 1
#' @return 
#'The output of the function \code{\link{optimal_multiple_tte}} is a data.frame containing the optimization results:
#'\describe{
#'  \item{u}{maximal expected utility}
#'  \item{HRgo}{optimal threshold value for the decision rule to go to phase III}
#'  \item{n2}{optimal total sample size in phase II}
#'  \item{n3}{total expected sample size for phase III; rounded to next natural number}
#'  \item{n}{total sample size in the program; n = n2 + n3}
#'  \item{K}{maximal costs of the program}
#'  \item{pgo}{probability to go to phase III}
#'  \item{sProg}{probability of a successful program}
#'  \item{sProg1}{probability of a successful program with "small" treatment effect in Phase III}
#'  \item{sProg2}{probability of a successful program with "medium" treatment effect in Phase III}
#'  \item{sProg3}{probability of a successful program with "large" treatment effect in Phase III }
#'  \item{K2}{expected costs for phase II}
#'  \item{K3}{expected costs for phase III}
#'  }
#'and further input parameters.
#'res
#'Taking cat(comment()) of the data.frame object lists the used optimization sequences, start and finish date of the optimization procedure.
#' @examples 
#'\dontrun{ #res <- optimal_multiple_normal(Delta1 = 0.75, Delta2 = 0.80,    # define assumed true HRs and control arm event rate
#' # in1=300, in2=600, sigma1 = 8, sigma2= 12,
#' # n2min = 30, n2max = 90, stepn2 = 6,                    # define optimization set for n2
#' # kappamin = 0.7, kappamax = 0.9, stepkappa = 0.05,         # define optimization set for HRgo
#' # alpha = 0.05, beta = 0.1,                              # drug development planning parameters
#' # c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,               # define fixed and variable costs for phase II and III
#' # K = Inf, N = Inf, S = -Inf,                            # set maximal costs/ expected sample size for the program or minimal expected probability of a successful program
#' # steps1 = 0,                                            # define lower boundary for "small"
#' # stepm1 = 0.5,                                          # "medium"
#' # stepl1 = 0.88,                                         # and "large" treatment effect size categories as proposed by IQWiG (2016)
#'#  b1 = 1000, b2 = 2000, b3 = 3000,                       # define expected benefit for a "small", "medium" and "large" treatment effect
#' # rho = 0.5, relaxed = TRUE,                             # relaxed "TRUE"
#' # fixed = TRUE,                                          #   treatment effect
#' # num_cl = 1)}                                            # set number of cores used for parallelized computing (check maximum number possible with detectCores())
#'#res
#'#cat(comment(res))                                        # displays the optimization sequence, start and finish date of the optimization procedure.
#' @section drugdevelopR functions:
#'The drugdevelopR package provides the functions
#'\itemize{
#'  \item \code{\link{optimal_tte}},
#'  \item \code{\link{optimal_binary}} and
#'  \item \code{\link{optimal_normal}}
#'}
#'to plan optimal phase II/III drug development programs with
#'\itemize{
#'  \item time-to-event (treatment effect measured by hazard ratio (HR)),
#'  \item binary (treatment effect measured by risk ratio (RR)) or
#'  \item normally distributed (treatment effect measured by standardized difference in means (Delta))
#'}
#'endpoint, where the treatment effect is modelled by a \href{https://web.imbi.uni-heidelberg.de/prior/}{prior}. Optimal phase II/III drug development planning with fixed treatment effects can be done with the help of the R Shiny application \href{https://web.imbi.uni-heidelberg.de/basic/}{basic}. Extensions are 
#'\itemize{
#'  \item optimal planning of programs including methods for discounting of phase II results (function: \code{\link{optimal_bias}}, App: \href{https://web.imbi.uni-heidelberg.de/bias/}{bias}),
#'  \item optimal planning of programs with several phase III trials (function: \code{\link{optimal_multitrial}}, App: \href{https://web.imbi.uni-heidelberg.de/multitrial/}{multitrial}) and
#'  \item optimal planning of programs with multiple arms (function: \code{\link{optimal_multiarm}}, App: \href{https://web.imbi.uni-heidelberg.de/multiarm/}{multiarm}).
#'}
#' @references 
#'Preussler, S., Kirchner, M., Goette, H., Kieser, M. (2019). Optimal Designs for Multi-Arm Phase II/III Drug Development Programs. Submitted to peer-review journal.
#'
#'IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at \href{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}, assessed last 15.05.19.
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
optimal_multiple_normal <- function(){}

#' @title Optimal phase II/III drug development planning for programs with multiple endpoints
#' @description The function \code{\link{optimal_multiple_tte}} of the drugdevelopR package enables planning of phase II/III drug development programs with optimal sample size allocation and go/no-go decision rules (Preussler et. al, 2019). 
#'Fast computing is enabled by parallel programming.
#' @name optimal_multiple_tte
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param id1 amount of information for hr1 in terms of number of events
#' @param id2 amount of information for hr2 in terms of number of events
#' @param rho correlation between the two endpoints
#' @param fixed assumed fixed treatment effect 
#' @param ec control arm event rate for phase II and III
#' @param n2min minimal number of events for phase II, must be divisible by 3
#' @param n2max maximal number of events for phase II, must be divisible by 3
#' @param stepn2 stepsize for the optimization over d2, must be divisible by 3
#' @param hrgomin minimal threshold value for the go/no-go decision rule
#' @param hrgomax maximal threshold value for the go/no-go decision rule
#' @param stephrgo stepsize for the optimization over HRgo
#' @param beta 1-beta (any-pair) power for calculation of the number of events for phase III
#' @param alpha one-sided significance level/ family wise error rate
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the maximal costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected minimal probability of a successful program, default: -Inf, e.g. no constraint
#' @param steps1 lower boundary for effect size category "small" in HR scale, default: 1
#' @param stepm1 lower boundary for effect size category "medium" in HR scale = upper boundary for effect size category "small" in HR scale, default: 0.95
#' @param stepl1 lower boundary for effect size category "large" in HR scale = upper boundary for effect size category "medium" in HR scale, default: 0.85
#' @param b11 expected gain for effect size category `"small"` if endpoint OS is significant
#' @param b21 expected gain for effect size category `"medium"`if endpoint OS is significant
#' @param b31 expected gain for effect size category `"large"` if endpoint OS is significant 
#' @param b12 expected gain for effect size category `"small"` if endpoint OS is not significant
#' @param b22 expected gain for effect size category `"medium"`if endpoint OS is not significant
#' @param b32 expected gain for effect size category `"large"` if endpoint OS is not significant
#' @param num_cl number of clusters used for parallel computing, default: 1
#' @return 
#'The output of the function \code{\link{optimal_multiple_tte}} is a data.frame containing the optimization results:
#'\describe{
#'  \item{Strategy}{Strategy, 1: "only best promising" or 2: "all promising"}
#'  \item{u}{maximal expected utility}
#'  \item{HRgo}{optimal threshold value for the decision rule to go to phase III}
#'  \item{n2}{optimal total sample size in phase II}
#'  \item{n3}{total expected sample size for phase III; rounded to next natural number}
#'  \item{n}{total sample size in the program; n = n2 + n3}
#'  \item{K}{maximal costs of the program}
#'  \item{pgo}{probability to go to phase III}
#'  \item{sProg}{probability of a successful program}
#'  \item{sProg1}{probability of a successful program with "small" treatment effect in Phase III}
#'  \item{sProg2}{probability of a successful program with "medium" treatment effect in Phase III}
#'  \item{sProg3}{probability of a successful program with "large" treatment effect in Phase III }
#'  \item{K2}{expected costs for phase II}
#'  \item{K3}{expected costs for phase III}
#'  \item{OP}{probability that one endpoint is significant}
#'  }
#'and further input parameters.
#'res
#'Taking cat(comment()) of the data.frame object lists the used optimization sequences, start and finish date of the optimization procedure.
#' @examples 
#'#res <- optimal_multiple_tte(hr1 = 0.75, hr2 = 0.80, 
#'#  ec = 0.6,                                         # define assumed true HRs and control arm event rate
#'#  id1 = 210, id2 = 420,
#'#  n2min = 30, n2max = 90, stepn2 = 6,               # define optimization set for n2
#'#  hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05,    # define optimization set for HRgo
#'#  alpha = 0.05, beta = 0.1,                         # drug development planning parameters
#'#  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,          # define fixed and variable costs for phase II and III
#'#  K = Inf, N = Inf, S = -Inf,                       # set constraints
#'#  steps1 = 1,                                       # define lower boundary for "small"
#'#  stepm1 = 0.95,                                    # "medium"
#'#  stepl1 = 0.85,                                    # and "large" treatment effect size categories as proposed by IQWiG (2016)
#'#  b11 = 1000, b21 = 2000, b31 = 3000,
#'#  b12 = 1000, b22 = 1500, b32 = 2000,               # define expected benefit for a "small", "medium" and "large" treatment effect
#'#  rho = 0.5, fixed = TRUE,                          # correlation and treatment effect
#'#  num_cl = 1)                                       # set number of cores used for parallelized computing 
#'# res
#'# cat(comment(res))                                   # displays the optimization sequence, start and finish date of the optimization 
#' @section drugdevelopR functions:
#'The drugdevelopR package provides the functions
#'\itemize{
#'  \item \code{\link{optimal_tte}},
#'  \item \code{\link{optimal_binary}} and
#'  \item \code{\link{optimal_normal}}
#'}
#'to plan optimal phase II/III drug development programs with
#'\itemize{
#'  \item time-to-event (treatment effect measured by hazard ratio (HR)),
#'  \item binary (treatment effect measured by risk ratio (RR)) or
#'  \item normally distributed (treatment effect measured by standardized difference in means (Delta))
#'}
#'endpoint, where the treatment effect is modelled by a \href{https://web.imbi.uni-heidelberg.de/prior/}{prior}. Optimal phase II/III drug development planning with fixed treatment effects can be done with the help of the R Shiny application \href{https://web.imbi.uni-heidelberg.de/basic/}{basic}. Extensions are 
#'\itemize{
#'  \item optimal planning of programs including methods for discounting of phase II results (function: \code{\link{optimal_bias}}, App: \href{https://web.imbi.uni-heidelberg.de/bias/}{bias}),
#'  \item optimal planning of programs with several phase III trials (function: \code{\link{optimal_multitrial}}, App: \href{https://web.imbi.uni-heidelberg.de/multitrial/}{multitrial}) and
#'  \item optimal planning of programs with multiple arms (function: \code{\link{optimal_multiarm}}, App: \href{https://web.imbi.uni-heidelberg.de/multiarm/}{multiarm}).
#'}
#' @references 
#'Preussler, S., Kirchner, M., Goette, H., Kieser, M. (2019). Optimal Designs for Multi-Arm Phase II/III Drug Development Programs. Submitted to peer-review journal.
#'
#'IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at \href{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}, assessed last 15.05.19.
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
optimal_multiple_tte <- function(){}

#' @title Optimal phase II/III drug development planning where several phase III trials are performed
#' @description The drugdevelopR package enables planning of phase II/III drug development programs with optimal sample size allocation and go/no-go decision rules. For binary endpoints the treatment effect is measured by the risk ratio (RR).The assumed true treatment effects can be assumed fixed or modelled by a prior distribution. The R Shiny application \href{https://web.imbi.uni-heidelberg.de/prior/}{prior} visualizes the prior distributions used in this package. Fast coputing is enabled by parallel programming.
#' @name optimal_multitrial_binary
#' @param w weight for mixture prior distribution
#' @param p0 assumed true rate of control group
#' @param p11 assumed true rate of treatment group
#' @param p12 assumed true rate of treatment group
#' @param in1 amount of information for p11 in terms of sample size
#' @param in2 amount of information for p12 in terms of sample size
#' @param n2min minimal total sample size for phase II; must be even number
#' @param n2max maximal total sample size for phase II, must be even number
#' @param stepn2 stepsize for the optimization over n2; must be even number
#' @param rrgomin minimal threshold value for the go/no-go decision rule
#' @param rrgomax maximal threshold value for the go/no-go decision rule
#' @param steprrgo stepsize for the optimization over RRgo
#' @param beta 1-beta power for calculation of sample size for phase III
#' @param alpha significance level
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param b1 expected gain for effect size category "small"
#' @param b2 expected gain for effect size category "medium"
#' @param b3 expected gain for effect size category "large"
#' @param case choose case: "at least 1, 2 or 3 significant trials needed for approval"
#' @param strategy choose strategy: "conduct 1, 2, 3 or 4 trial"; TRUE calculates all strategies of the selected Case
#' @param fixed choose if true treatment effects are fixed or random, if TRUE p11 is used as fixed effect for p1
#' @param num_cl number of clusters used for parallel computing, default: 1
#' @return 
#'The output of the function \code{\link{optimal_binary}} is a data.frame containing the optimization results:
#'\describe{
#'  \item{Case}{Case: "number of significant trials needed"}
#'  \item{Strategy}{Strategy: "number of conducted trials"}
#'  \item{u}{maximal expected utility}
#'  \item{RRgo}{optimal threshold value for the decision rule to go to phase III}
#'  \item{n2}{total sample size for phase II}
#'  \item{n3}{total sample size for phase III; rounded to the next even natural number}
#'  \item{n}{total sample size in the program; n = n2 + n3}
#'  \item{K}{maximal costs of the program}
#'  \item{pgo}{probability to go to phase III}
#'  \item{sProg}{probability of a successful program}
#'  \item{sProg1}{probability of a successful program with "small" treatment effect in Phase III}
#'  \item{sProg2}{probability of a successful program with "medium" treatment effect in Phase III}
#'  \item{sProg3}{probability of a successful program with "large" treatment effect in Phase III }
#'  \item{K2}{expected costs for phase II}
#'  \item{K3}{expected costs for phase III}
#'}
#'and further input parameters.
#'
#'Taking cat(comment()) of the data.frame object lists the used optimization sequences, start and finish date of the optimization procedure.
#' @examples 
#'res <- optimal_multitrial_binary(w = 0.3,                # define parameters for prior
#'  p0 = 0.6, p11 =  0.3, p12 = 0.5, in1 = 30, in2 = 60,   # (https://web.imbi.uni-heidelberg.de/prior/)
#'  n2min = 20, n2max = 100, stepn2 = 4,                   # define optimization set for n2
#'  rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.05,         # define optimization set for RRgo
#'  alpha = 0.05, beta = 0.1,                              # drug development planning parameters
#'  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,               # define fixed and variable costs for phase II and III,
#'  K = Inf, N = Inf, S = -Inf,                            # set constraints
#'  b1 = 1000, b2 = 2000, b3 = 3000,                       # define expected benefit for a "small", "medium" and "large" treatment effect
#'  case = 1, strategy = TRUE,                             # chose Case and Strategy                                   
#'  fixed = TRUE,                                          # choose if true treatment effects are fixed or random
#'  num_cl = 1)                                            # set number of cores used for parallelized computing 
#'res
#'cat(comment(res))                                        # displays the optimization sequence, start/finish date of procedure
#' @section drugdevelopR functions:
#'The drugdevelopR package provides the functions
#'\itemize{
#'  \item \code{\link{optimal_tte}},
#'  \item \code{\link{optimal_binary}} or
#'  \item \code{\link{optimal_normal}}
#'}
#'to plan optimal phase II/III drug development programs with
#'\itemize{
#'  \item time-to-event (treatment effect measured by hazard ratio (HR)),
#'  \item binary (treatment effect measured by risk ratio (RR)) and
#'  \item normally distributed (treatment effect measured by standardized difference in means (Delta))
#'}
#'endpoint, where the treatment effect is modelled by a \href{https://web.imbi.uni-heidelberg.de/prior/}{prior}. Optimal phase II/III drug development planning with fixed treatment effects can be done with the help of the R Shiny application \href{https://web.imbi.uni-heidelberg.de/basic/}{basic}. Extensions are 
#'\itemize{
#'  \item optimal planning of programs including methods for discounting of phase II results (function: \code{\link{optimal_bias}}, App: \href{https://web.imbi.uni-heidelberg.de/bias/}{bias}),
#'  \item optimal planning of programs with several phase III trials (function: \code{\link{optimal_multitrial}}, App: \href{https://web.imbi.uni-heidelberg.de/multitrial/}{multitrial}) and
#'  \item optimal planning of programs with multiple arms (function: \code{\link{optimal_multiarm}}, App: \href{https://web.imbi.uni-heidelberg.de/multiarm/}{multiarm}).
#'}
#' @references 
#'IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at \href{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}, assessed last 15.05.19.
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
optimal_multitrial_binary <- function(){}

#' @title Optimal phase II/III drug development planning where several phase III trials are performed
#' @description The drugdevelopR package enables planning of phase II/III drug development programs with optimal sample size allocation and go/no-go decision rules. For normally distributed endpoints the treatment effect is measured by the standardized difference in means (Delta). The assumed true treatment effects can be assumed fixed or modelled by a prior distribution. The R Shiny application \href{https://web.imbi.uni-heidelberg.de/prior/}{prior} visualizes the prior distributions used in this package. Fast coputing is enabled by parallel programming.
#' @name optimal_multitrial_normal
#' @param w weight for mixture prior distribution
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param in1 amount of information for Delta2 in terms of sample size
#' @param in2 amount of information for Delta1 in terms of sample size
#' @param a lower boundary for the truncation
#' @param b upper boundary for the truncation
#' @param n2min minimal total sample size for phase II; must be even number
#' @param n2max maximal total sample size for phase II, must be even number
#' @param stepn2 stepsize for the optimization over n2; must be even number
#' @param kappamin minimal threshold value for the go/no-go decision rule
#' @param kappamax maximal threshold value for the go/no-go decision rule
#' @param stepkappa stepsize for the optimization over kappa
#' @param beta 1-beta power for calculation of sample size for phase III
#' @param alpha significance level
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param b1 expected gain for effect size category "small"
#' @param b2 expected gain for effect size category "medium"
#' @param b3 expected gain for effect size category "large"
#' @param case choose case: "at least 1, 2 or 3 significant trials needed for approval"
#' @param strategy choose strategy: "conduct 1, 2, 3 or 4 trial"; TRUE calculates all strategies of the selected Case 
#' @param fixed choose if true treatment effects are fixed or random, if TRUE hr1 is used as fixed effect
#' @param num_cl number of clusters used for parallel computing, default: 1
#' @return 
#'The output of the function \code{\link{optimal_multitrial_normal}} is a data.frame containing the optimization results:
#'\describe{
#' \item{Case}{Case: "number of significant trials needed"}
#'  \item{Strategy}{Strategy: "number of conducted trials"}
#'  \item{u}{maximal expected utility}
#'  \item{kappa}{optimal threshold value for the decision rule to go to phase III}
#'  \item{n2}{total sample size for phase II}
#'  \item{n3}{total sample size for phase III; rounded to the next even natural number}
#'  \item{n}{total sample size in the program; n = n2 + n3}
#'  \item{K}{maximal costs of the program}
#'  \item{pgo}{probability to go to phase III}
#'  \item{sProg}{probability of a successful program}
#'  \item{sProg1}{probability of a successful program with "small" treatment effect in Phase III}
#'  \item{sProg2}{probability of a successful program with "medium" treatment effect in Phase III}
#'  \item{sProg3}{probability of a successful program with "large" treatment effect in Phase III }
#'  \item{K2}{expected costs for phase II}
#'  \item{K3}{expected costs for phase III}
#'}
#'and further input parameters.
#'
#'Taking cat(comment()) of the data.frame object lists the used optimization sequences, start and finish date of the optimization procedure.
#' @examples 
#'res <- optimal_multitrial_normal(w=0.3,                             # define parameters for prior
#'  Delta1 = 0.375, Delta2 = 0.625, in1=300, in2=600,      # (https://web.imbi.uni-heidelberg.de/prior/)
#'  a = 0.25, b = 0.75,
#'  n2min = 20, n2max = 100, stepn2 = 4,                   # define optimization set for n2
#'  kappamin = 0.02, kappamax = 0.2, stepkappa = 0.02,     # define optimization set for kappa
#'  alpha = 0.05, beta = 0.1,                              # drug development planning parameters
#'  c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,             # define fixed and variable costs for phase II and III
#'  K = Inf, N = Inf, S = -Inf,                            # set constraints
#'  b1 = 3000, b2 = 8000, b3 = 10000,                      # define expected benefit for a "small", "medium" and "large" treatment effect                                             # assume different/same population structures in phase II and III
#'  case = 1, strategy = TRUE,                             # chose Case and Strategy
#'  fixed = TRUE,                                          # choose if true treatment effects are fixed or random
#'  num_cl = 1)                                            # set number of cores used for parallelized computing
#'res
#'
#'cat(comment(res))                                        # displays the optimization sequence, start/ finish date of procedure.
#' @section drugdevelopR functions:
#'The drugdevelopR package provides the functions
#'\itemize{
#'  \item \code{\link{optimal_tte}},
#'  \item \code{\link{optimal_binary}} or
#'  \item \code{\link{optimal_normal}}
#'}
#'to plan optimal phase II/III drug development programs with
#'\itemize{
#'  \item time-to-event (treatment effect measured by hazard ratio (HR)),
#'  \item binary (treatment effect measured by risk ratio (RR)) and
#'  \item normally distributed (treatment effect measured by standardized difference in means (Delta))
#'}
#'endpoint, where the treatment effect is assumed fixed or modelled by a \href{https://web.imbi.uni-heidelberg.de/prior/}{prior}. Optimal phase II/III drug development planning with fixed treatment effects can also be done with the help of the R Shiny application \href{https://web.imbi.uni-heidelberg.de/basic/}{basic}. Extensions are 
#'\itemize{
#'  \item optimal planning of programs including methods for discounting of phase II results (function: \code{\link{optimal_bias}}, App: \href{https://web.imbi.uni-heidelberg.de/bias/}{bias}),
#'  \item optimal planning of programs with several phase III trials (function: \code{\link{optimal_multitrial}}, App: \href{https://web.imbi.uni-heidelberg.de/multitrial/}{multitrial}) and
#'  \item optimal planning of programs with multiple arms (function: \code{\link{optimal_multiarm}}, App: \href{https://web.imbi.uni-heidelberg.de/multiarm/}{multiarm}).
#'}
#' @references 
#'Cohen, J. (1988). Statistical power analysis for the behavioral sciences.
#'
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export 
optimal_multitrial_normal <- function(){}

