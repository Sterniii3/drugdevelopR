# auxiliary functions

t1 <- function(x, p0){((1-p0)/p0) + ((1-x)/x)}
t2 <- function(x, p0){sqrt(2*(1-((p0 + x)/2))/((p0 + x)/2))}
t3 <- function(x, p0){sqrt(((1-p0)/p0) + ((1-x)/x))}

# 1.1. conservative sample size calculation: use lower bound of one-sided confidence intervall
##############################################################################################

# prior distribution
# as above 

# expected probability to go to phase III
# as above

#' Expected sample size for phase III for bias adjustment programs and binary distributed outcomes
#' 
#' To discount for overoptimistic results in phase II when calculating the optimal sample size in phase III, 
#' it is necessary to use the functions `En3_binary_L()`, `En3_binary_L2()`, `En3_binary_R()` and `En3_binary_R2()`.
#' Each function describes a specific case:
#' - `En3_binary_L()`: calculates the optimal sample size for an additive adjustment factor (i.e adjust the lower bound of the one-sided confidence interval), 
#' however the go-decision is not affected by the bias adjustment
#' - `En3_binary_L2()`: calculates the optimal sample size for an additive adjustment factor (i.e adjust the lower bound of the one-sided confidence interval)
#' when the go-decision is also affected by the bias adjustment
#' - `En3_binary_R()`: calculates the optimal sample size for a multiplicative adjustment factor (i.e. use estimate with a retention factor), 
#' however the go-decision is not affected by the bias adjustment
#' - `En3_binary_R2()`: calculates the optimal sample size for a multiplicative adjustment factor (i.e. use estimate with a retention factor)
#' when the go-decision is also affected by the bias adjustment 
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
#' @return The output of the the functions `En3_binary_L`, `En3_binary_L2`, `En3_binary_R` and `En3_binary_R2` is the expected number of participants in phase III. 
#' @importFrom stats qnorm integrate dnorm
#' @examples res <- En3_binary_L(RRgo = 0.8, n2 = 50, Adj = 0, 
#'                               alpha = 0.025, beta = 0.1, p0 = 0.6,  w = 0.3,
#'                               p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600, 
#'                               fixed = FALSE)
#'           res <-  En3_binary_L2(RRgo = 0.8, n2 = 50, Adj = 0, 
#'                               alpha = 0.025, beta = 0.1, p0 = 0.6,  w = 0.3,
#'                               p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600, 
#'                               fixed = FALSE)
#'           res <- En3_binary_R(RRgo = 0.8, n2 = 50, Adj = 1, 
#'                               alpha = 0.025, beta = 0.1, p0 = 0.6,  w = 0.3,
#'                               p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600, 
#'                               fixed = FALSE)
#'           res <- En3_binary_R2(RRgo = 0.8, n2 = 50, Adj = 1, 
#'                               alpha = 0.025, beta = 0.1, p0 = 0.6,  w = 0.3,
#'                               p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600, 
#'                               fixed = FALSE)
#'                               
#' @name En3_bias_binary                             
#' @export
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
En3_binary_L <-  function(RRgo, n2, Adj, alpha, beta, p0, w, p11, p12, in1, in2, fixed){
  if(fixed){
      int   = try(integrate(function(y){
        ((2*(qnorm(1-alpha)*t2(p11, p0)+qnorm(1-beta)*t3(p11, p0))^2)/(y-qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-p11/p11))))^2) *
          dnorm(y,
                mean = -log(p11/p0),
                sd = sqrt((2/n2)*t1(p11, p0)))
      }, - log(RRgo), Inf), silent=TRUE)
    if(inherits(int ,'try-error')){
      warning(as.vector(int))
      integrated <- NA_real_
    } else {
      integrated <- int$value
    }
    return(integrated)
  }else{
    int   = try(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            ((2*(qnorm(1-alpha)*t2(x, p0)+qnorm(1-beta)*t3(x, p0))^2)/(y-qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-x/x))))^2) *
              dnorm(y,
                    mean = -log(x/p0),
                    sd = sqrt((2/n2)*t1(x, p0)))*
              prior_binary(x, w, p11, p12, in1, in2)
          }, - log(RRgo), Inf)$value
        })
      }, 0, 1), silent=TRUE)
    if(inherits(int ,'try-error')){
      warning(as.vector(int))
      integrated <- NA_real_
    } else {
      integrated <- int$value
    }
    return(integrated)
  }
}

#' Expected probability of a successful program for bias adjustment programs with binary distributed outcomes
#' 
#' To discount for overoptimistic results in phase II when calculating the optimal sample size in phase III, 
#' it is necessary to use the following functions, which each describe a specific case:
#' - `EPsProg_binary_L()`: calculates the expected probability of a successful for an additive adjustment factor (i.e adjust the lower bound of the one-sided confidence interval), 
#' however the go-decision is not affected by the bias adjustment
#' - `EPsProg_binary_L2()`: calculates the expected probability of a successful for an additive adjustment factor (i.e adjust the lower bound of the one-sided confidence interval)
#' when the go-decision is also affected by the bias adjustment
#' - `EPsProg_binary_R()`: calculates the expected probability of a successful for a multiplicative adjustment factor (i.e. use estimate with a retention factor), 
#' however the go-decision is not affected by the bias adjustment
#' - `EPsProg_binary_R2()`: calculates the expected probability of a successful for a multiplicative adjustment factor (i.e. use estimate with a retention factor)
#' when the go-decision is also affected by the bias adjustment 
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
#' @return  The output of the the functions `EPsProg_binary_L()`, `EPsProg_binary_L2()`, `EPsProg_binary_R()` and `EPsProg_binary_R2()` is the expected probability of a successful program.
#' @importFrom stats qnorm integrate dnorm
#' @examples res <- EPsProg_binary_L(RRgo = 0.8, n2 = 50, Adj = 0, 
#'                                  alpha = 0.025, beta = 0.1, 
#'                                  step1 = 1, step2 = 0.95, p0 = 0.6,  w = 0.3,
#'                                  p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600, 
#'                                  fixed = FALSE)
#'           res <- EPsProg_binary_L2(RRgo = 0.8, n2 = 50, Adj = 0, 
#'                                  alpha = 0.025, beta = 0.1, 
#'                                  step1 = 1, step2 = 0.95, p0 = 0.6,  w = 0.3,
#'                                  p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600, 
#'                                  fixed = FALSE)
#'           res <- EPsProg_binary_R(RRgo = 0.8, n2 = 50, Adj = 1, 
#'                                  alpha = 0.025, beta = 0.1, 
#'                                  step1 = 1, step2 = 0.95, p0 = 0.6,  w = 0.3,
#'                                  p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600, 
#'                                  fixed = FALSE)
#'           res <- EPsProg_binary_R2(RRgo = 0.8, n2 = 50, Adj = 1, 
#'                                  alpha = 0.025, beta = 0.1, 
#'                                  step1 = 1, step2 = 0.95, p0 = 0.6,  w = 0.3,
#'                                  p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600, 
#'                                  fixed = FALSE)
#' @name EPsProg_bias_binary                               
#' @export
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
EPsProg_binary_L <-  function(RRgo, n2, Adj, alpha, beta, step1, step2, p0, w, p11, p12, in1, in2, fixed){
  
  if(fixed){
    return(
      integrate(function(y){
        ( pnorm(qnorm(1 - alpha) -
                  log(step2)/sqrt((t1(p11, p0)*(y-qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-p11/p11))))^2)/(qnorm(1-alpha)*t2(p11, p0) +
                                                                                                          qnorm(1-beta)*t3(p11, p0))^2),
                mean = -log((p11)/p0)/sqrt((t1(p11, p0)*(y-qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-p11/p11))))^2)/(qnorm(1-alpha)*t2(p11, p0) +
                                                                                                                   qnorm(1-beta)*t3(p11, p0))^2),
                sd = 1) -
            pnorm(qnorm(1 - alpha) -
                    log(step1)/sqrt((t1(p11, p0)*y^2)/(qnorm(1-alpha)*t2(p11, p0) +
                                                         qnorm(1-beta)*t3(p11, p0))^2),
                  mean = -log((p11)/p0)/sqrt((t1(p11, p0)*(y-qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-p11/p11))))^2)/(qnorm(1-alpha)*t2(p11, p0) +
                                                                                                                     qnorm(1-beta)*t3(p11, p0))^2),
                  sd = 1) ) *
          dnorm(y,
                mean = -log(p11/p0),
                sd = sqrt((2/n2)*t1(p11, p0)))
      }, -log(RRgo), Inf)$value
    )
  }else{
    return(
      integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            ( pnorm(qnorm(1 - alpha) -
                      log(step2)/sqrt((t1(x, p0)*(y-qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-x/x))))^2)/(qnorm(1-alpha)*t2(x, p0) +
                                                                                                        qnorm(1-beta)*t3(x, p0))^2),
                    mean = -log((x)/p0)/sqrt((t1(x, p0)*(y-qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-x/x))))^2)/(qnorm(1-alpha)*t2(x, p0) +
                                                                                                               qnorm(1-beta)*t3(x, p0))^2),
                    sd = 1) -
                pnorm(qnorm(1 - alpha) -
                        log(step1)/sqrt((t1(x, p0)*(y-qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-x/x))))^2)/(qnorm(1-alpha)*t2(x, p0) +
                                                                                                          qnorm(1-beta)*t3(x, p0))^2),
                      mean = -log((x)/p0)/sqrt((t1(x, p0)*(y-qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-x/x))))^2)/(qnorm(1-alpha)*t2(x, p0) +
                                                                                                                 qnorm(1-beta)*t3(x, p0))^2),
                      sd = 1) ) *
              dnorm(y,
                    mean = -log(x/p0),
                    sd = sqrt((2/n2)*t1(x, p0))) *
              prior_binary(x, w, p11, p12, in1, in2)
          }, -log(RRgo), Inf)$value
        })
      }, 0, 1)$value
    )
  } 
  
}

#' Utility function for bias adjustment programs with binary distributed outcomes.
#' 
#' The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters and the expected probability of a successful program. 
#' The utility is in a further step maximized by the `optimal_bias_binary()` function.
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
#' @return The output of the the functions `utility_binary_L()`, `utility_binary_L2()`, `utility_binary_R()` and `utility_binary_R2()` is the expected utility of the program.
#' @examples res <- utility_binary_L(n2 = 50, RRgo = 0.8, Adj = 0.1, w = 0.3, 
#'                                  p0 = 0.6, p11 =  0.3, p12 = 0.5, 
#'                                  in1 = 300, in2 = 600, 
#'                                  alpha = 0.025, beta = 0.1,
#'                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                  K = Inf, N = Inf, S = -Inf,
#'                                  steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
#'                                  b1 = 1000, b2 = 2000, b3 = 3000, 
#'                                  fixed = TRUE)
#'          res <- utility_binary_L2(n2 = 50, RRgo = 0.8, Adj = 0.1, w = 0.3, 
#'                                  p0 = 0.6, p11 =  0.3, p12 = 0.5, 
#'                                  in1 = 300, in2 = 600, 
#'                                  alpha = 0.025, beta = 0.1,
#'                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                  K = Inf, N = Inf, S = -Inf,
#'                                  steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
#'                                  b1 = 1000, b2 = 2000, b3 = 3000, 
#'                                  fixed = TRUE)
#'          res <- utility_binary_R(n2 = 50, RRgo = 0.8, Adj = 0.9, w = 0.3, 
#'                                  p0 = 0.6, p11 =  0.3, p12 = 0.5, 
#'                                  in1 = 300, in2 = 600, 
#'                                  alpha = 0.025, beta = 0.1,
#'                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                  K = Inf, N = Inf, S = -Inf,
#'                                  steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
#'                                  b1 = 1000, b2 = 2000, b3 = 3000, 
#'                                  fixed = TRUE)
#'          res <- utility_binary_R2(n2 = 50, RRgo = 0.8, Adj = 0.9, w = 0.3, 
#'                                  p0 = 0.6, p11 =  0.3, p12 = 0.5, 
#'                                  in1 = 300, in2 = 600, 
#'                                  alpha = 0.025, beta = 0.1,
#'                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                  K = Inf, N = Inf, S = -Inf,
#'                                  steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
#'                                  b1 = 1000, b2 = 2000, b3 = 3000, 
#'                                  fixed = TRUE)
#' @name utility_bias_binary                                 
#' @export
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
utility_binary_L <-  function(n2, RRgo, Adj, w, p0, p11, p12, in1, in2,
                              alpha, beta, 
                              c2, c3, c02, c03, 
                              K, N, S,
                              steps1, stepm1, stepl1,
                              b1, b2, b3,
                              fixed){
  
  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 0
  
  
  n3  <-  En3_binary_L(RRgo = RRgo, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                          p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, fixed = fixed)
  
  if(is.na(n3)){
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))  
  }
  else{
  
  n3  <- ceiling(n3)
  
  if(round(n3/2) != n3 / 2) {n3 <- n3 + 1}
  
  if(n2+n3>N){
    
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
    
  }else{
    
    pg    <-  Epgo_binary(RRgo = RRgo, n2 = n2, p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, fixed = fixed)
    
    K2    <-  c02 + c2 * n2         # cost phase II
    K3    <-  c03 * pg + c3 * n3    # cost phase III
    
    if(K2+K3>K){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
    }else{
      # probability of a successful program; small, medium, large effect size
      prob1 <-  EPsProg_binary_L(RRgo = RRgo, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                                 step1 = steps1, step2 =  steps2,
                                 p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                                 fixed = fixed)
      prob2 <-  EPsProg_binary_L(RRgo = RRgo, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                                 step1 =  stepm1, step2 =  stepm2,
                                 p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                                 fixed = fixed)
      prob3 <-  EPsProg_binary_L(RRgo = RRgo, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                                 step1 =  stepl1, step2 = stepl2,
                                 p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                                 fixed = fixed)
      
      SP    <-  prob1 + prob2 + prob3
      
      if(SP<S){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        
      }else{
        
        G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
        
        EU    <-  - K2 - K3 + G
        
        return(c(EU, n3, SP, pg, K2, K3, prob1, prob2, prob3))
      }
    }
  }
 }
}

# 1.2. conservative decision rule and sample size calculation: 
# use lower bound of one-sided confidence intervall
##############################################################################################

# prior distribution
# as above

#' Expected probability to go to phase III for bias adjustment programs with binary distributed outcomes
#' 
#' In the case we do not only want do discount for overoptimistic results in phase II when calculating the sample size in phase III, 
#' but also when deciding whether to go to phase III or not the functions `Epgo_binary_L2` and `Epgo_binary_R2` are necessary.
#' The function `Epgo_binary_L2` uses an additive adjustment parameter (i.e adjust the lower bound of the one-sided confidence interval),
#' the function `Epgo_binary_R2` uses a multiplicative adjustment parameter (i.e. use estimate with a retention factor)
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
#' @return The output of the the functions `Epgo_normal_L2` and `Epgo_normal_R2` is the expected number of participants in phase III with conservative decision rule and sample size calculation.
#' @importFrom stats qnorm integrate dnorm
#' @examples res <- Epgo_binary_L2(RRgo = 0.8, n2 = 50, Adj = 0,  p0 = 0.6,  w = 0.3,
#'                               p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600, 
#'                               fixed = FALSE)
#'           res <- Epgo_binary_R2(RRgo = 0.8, n2 = 50, Adj = 1,  p0 = 0.6,  w = 0.3,
#'                               p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600, 
#'                               fixed = FALSE)
#' @name Epgo_bias_binary 
#' @export
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
Epgo_binary_L2 <-  function(RRgo, n2, Adj, p0, w, p11, p12, in1, in2, fixed){
  if(fixed){
    return(
      pnorm((-log(p11/p0) + log(RRgo)-qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-p11/p11))))/sqrt((2/n2)*t1(p11, p0))) 
    )
  }else{
    return(
      integrate(function(x){
        sapply(x, function(x){
          pnorm((-log(x/p0) + log(RRgo)-qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-x/x))))/sqrt((2/n2)*t1(x, p0)))  *
            prior_binary(x, w, p11, p12, in1, in2)
        })
      },  0, 1)$value
    )
  }
}

#' @rdname En3_bias_binary 
#' @export
En3_binary_L2 <-  function(RRgo, n2, Adj, alpha, beta, p0, w, p11, p12, in1, in2, fixed){
  if(fixed){
    int = try(
      integrate(function(y){
        ((2*(qnorm(1-alpha)*t2(p11, p0)+qnorm(1-beta)*t3(p11, p0))^2)/(y-qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-p11/p11))))^2) *
          dnorm(y,
                mean = -log(p11/p0),
                sd = sqrt((2/n2)*t1(p11, p0)))
      }, - log(RRgo)+qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-p11/p11))),Inf), silent=TRUE)
      if(inherits(int ,'try-error')){
        warning(as.vector(int))
        integrated <- NA_real_
      } else {
        integrated <- int$value
      }
      return(integrated)
  }else{
    int = try(
      integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            ((2*(qnorm(1-alpha)*t2(x, p0)+qnorm(1-beta)*t3(x, p0))^2)/(y-qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-x/x))))^2) *
              dnorm(y,
                    mean = -log(x/p0),
                    sd = sqrt((2/n2)*t1(x, p0)))*
              prior_binary(x, w, p11, p12, in1, in2)
          }, - log(RRgo)+qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-x/x))), Inf)$value
        })
      }, 0, 1), silent=TRUE)
    if(inherits(int ,'try-error')){
      warning(as.vector(int))
      integrated <- NA_real_
    } else {
      integrated <- int$value
    }
    return(integrated)
  }
} 

#' @rdname EPsProg_bias_binary 
#' @export
EPsProg_binary_L2 <-  function(RRgo, n2, Adj, alpha, beta, step1, step2, p0, w, p11, p12, in1, in2, fixed){
  
  if(fixed){
    return(
      integrate(function(y){
        ( pnorm(qnorm(1 - alpha) -
                  log(step2)/sqrt((t1(p11, p0)*(y-qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-p11/p11))))^2)/(qnorm(1-alpha)*t2(p11, p0) +
                                                                                                          qnorm(1-beta)*t3(p11, p0))^2),
                mean = -log((p11)/p0)/sqrt((t1(p11, p0)*(y-qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-p11/p11))))^2)/(qnorm(1-alpha)*t2(p11, p0) +
                                                                                                                 qnorm(1-beta)*t3(p11, p0))^2),
                sd = 1) -
            pnorm(qnorm(1 - alpha) -
                    log(step1)/sqrt((t1(p11, p0)*(y-qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-p11/p11))))^2)/(qnorm(1-alpha)*t2(p11, p0) +
                                                                                                            qnorm(1-beta)*t3(p11, p0))^2),
                  mean = -log((p11)/p0)/sqrt((t1(p11, p0)*(y-qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-p11/p11))))^2)/(qnorm(1-alpha)*t2(p11, p0) +
                                                                                                                     qnorm(1-beta)*t3(p11, p0))^2),
                  sd = 1) ) *
          dnorm(y,
                mean = -log(p11/p0),
                sd = sqrt((2/n2)*t1(p11, p0)))
      }, -log(RRgo) + qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-p11/p11))), Inf)$value
    )
  }else{
    return(
      integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            ( pnorm(qnorm(1 - alpha) -
                      log(step2)/sqrt((t1(x, p0)*(y-qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-x/x))))^2)/(qnorm(1-alpha)*t2(x, p0) +
                                                                                                        qnorm(1-beta)*t3(x, p0))^2),
                    mean = -log((x)/p0)/sqrt((t1(x, p0)*(y-qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-x/x))))^2)/(qnorm(1-alpha)*t2(x, p0) +
                                                                                                               qnorm(1-beta)*t3(x, p0))^2),
                    sd = 1) -
                pnorm(qnorm(1 - alpha) -
                        log(step1)/sqrt((t1(x, p0)*(y-qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-x/x))))^2)/(qnorm(1-alpha)*t2(x, p0) +
                                                                                                          qnorm(1-beta)*t3(x, p0))^2),
                      mean = -log((x)/p0)/sqrt((t1(x, p0)*(y-qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-x/x))))^2)/(qnorm(1-alpha)*t2(x, p0) +
                                                                                                                 qnorm(1-beta)*t3(x, p0))^2),
                      sd = 1) ) *
              dnorm(y,
                    mean = -log(x/p0),
                    sd = sqrt((2/n2)*t1(x, p0))) *
              prior_binary(x, w, p11, p12, in1, in2)
          }, -log(RRgo) + qnorm(1-Adj)*sqrt(2/n2*((1-p0)/p0 +(1-x/x))), Inf)$value
        })
      }, 0, 1)$value
    )
  }
  
}

#' @rdname utility_bias_binary
#' @export
utility_binary_L2 <-  function(n2, RRgo, Adj, w, p0, p11, p12, in1, in2,
                               alpha, beta, 
                               c2, c3, c02, c03, 
                               K, N, S,
                               steps1, stepm1, stepl1,
                               b1, b2, b3,
                               fixed){
  
  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 0
  
  
  n3  <-  En3_binary_L2(RRgo = RRgo, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                        p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, fixed = fixed)
  
  if(is.na(n3)){
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))  
  }
  else{
  
  n3  <- ceiling(n3)
  
  if(round(n3/2) != n3 / 2) {n3 <- n3 + 1}
  
  if(n2+n3>N){
    
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
    
  }else{
    
    pg    <-  Epgo_binary_L2(RRgo = RRgo,  Adj = Adj, n2 = n2, p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, fixed = fixed)
    
    K2    <-  c02 + c2 * n2         # cost phase II
    K3    <-  c03 * pg + c3 * n3    # cost phase III
    
    if(K2+K3>K){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
    }else{
      # probability of a successful program; small, medium, large effect size
      prob1 <-  EPsProg_binary_L2(RRgo = RRgo, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                                  step1 = steps1, step2 =  steps2,
                                  p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                                  fixed = fixed)
      prob2 <-  EPsProg_binary_L2(RRgo = RRgo, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                                  step1 =  stepm1, step2 =  stepm2,
                                  p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                                  fixed = fixed)
      prob3 <-  EPsProg_binary_L2(RRgo = RRgo, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                                  step1 =  stepl1, step2 = stepl2,
                                  p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                                  fixed = fixed)
      
      SP    <-  prob1 + prob2 + prob3
      
      if(SP<S){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        
      }else{
        
        G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
        
        EU    <-  - K2 - K3 + G
        
        return(c(EU, n3, SP, pg, K2, K3, prob1, prob2, prob3))
      }
    }
  }
 }
}

# 2.1. conservative sample size calculation: use estimate with retetion factor
##############################################################################################

# prior distribution
# as above

# expected probability to go to phase III
# as above

#' @rdname En3_bias_binary 
#' @export
En3_binary_R <-  function(RRgo, n2, Adj, alpha, beta, p0, w, p11, p12, in1, in2, fixed){
  if(fixed){
    int   = try(integrate(function(y){
        ((2*(qnorm(1-alpha)*t2(p11, p0)+qnorm(1-beta)*t3(p11, p0))^2)/(y*Adj)^2) *
          dnorm(y,
                mean = -log(p11/p0),
                sd = sqrt((2/n2)*t1(p11, p0)))
      }, - log(RRgo), Inf), silent=TRUE)
    if(inherits(int ,'try-error')){
      warning(as.vector(int))
      integrated <- NA_real_
    } else {
      integrated <- int$value
    }
    return(integrated)
  }else{
    int   = try(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            ((2*(qnorm(1-alpha)*t2(x, p0)+qnorm(1-beta)*t3(x, p0))^2)/(y*Adj)^2) *
              dnorm(y,
                    mean = -log(x/p0),
                    sd = sqrt((2/n2)*t1(x, p0)))*
              prior_binary(x, w, p11, p12, in1, in2)
          }, - log(RRgo), Inf)$value
        })
      }, 0, 1), silent=TRUE)
    if(inherits(int ,'try-error')){
      warning(as.vector(int))
      integrated <- NA_real_
    } else {
      integrated <- int$value
    }
    return(integrated)
  }
}


#' @rdname EPsProg_bias_binary 
#' @export
EPsProg_binary_R <-  function(RRgo, n2, Adj, alpha, beta, step1, step2, p0, w, p11, p12, in1, in2, fixed){
  
  if(fixed){
    return(
      integrate(function(y){
        ( pnorm(qnorm(1 - alpha) -
                  log(step2)/sqrt((t1(p11, p0)*(y*Adj)^2)/(qnorm(1-alpha)*t2(p11, p0) +
                                                             qnorm(1-beta)*t3(p11, p0))^2),
                mean = -log((p11)/p0)/sqrt((t1(p11, p0)*(y*Adj)^2)/(qnorm(1-alpha)*t2(p11, p0) +
                                                                      qnorm(1-beta)*t3(p11, p0))^2),
                sd = 1) -
            pnorm(qnorm(1 - alpha) -
                    log(step1)/sqrt((t1(p11, p0)*(y*Adj)^2)/(qnorm(1-alpha)*t2(p11, p0) +
                                                               qnorm(1-beta)*t3(p11, p0))^2),
                  mean = -log((p11)/p0)/sqrt((t1(p11, p0)*(y*Adj)^2)/(qnorm(1-alpha)*t2(p11, p0) +
                                                                        qnorm(1-beta)*t3(p11, p0))^2),
                  sd = 1) ) *
          dnorm(y,
                mean = -log(p11/p0),
                sd = sqrt((2/n2)*t1(p11, p0)))
      }, -log(RRgo), Inf)$value
    )
  }else{
    return(
      integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            ( pnorm(qnorm(1 - alpha) -
                      log(step2)/sqrt((t1(x, p0)*(y*Adj)^2)/(qnorm(1-alpha)*t2(x, p0) +
                                                               qnorm(1-beta)*t3(x, p0))^2),
                    mean = -log((x)/p0)/sqrt((t1(x, p0)*(y*Adj)^2)/(qnorm(1-alpha)*t2(x, p0) +
                                                                      qnorm(1-beta)*t3(x, p0))^2),
                    sd = 1) -
                pnorm(qnorm(1 - alpha) -
                        log(step1)/sqrt((t1(x, p0)*(y*Adj)^2)/(qnorm(1-alpha)*t2(x, p0) +
                                                                 qnorm(1-beta)*t3(x, p0))^2),
                      mean = -log((x)/p0)/sqrt((t1(x, p0)*(y*Adj)^2)/(qnorm(1-alpha)*t2(x, p0) +
                                                                        qnorm(1-beta)*t3(x, p0))^2),
                      sd = 1) ) *
              dnorm(y,
                    mean = -log(x/p0),
                    sd = sqrt((2/n2)*t1(x, p0))) *
              prior_binary(x, w, p11, p12, in1, in2)
          }, -log(RRgo), Inf)$value
        })
      }, 0, 1)$value
    )
  }
  
}

#' @rdname utility_bias_binary
#' @export
utility_binary_R <-  function(n2, RRgo, Adj, w, p0, p11, p12, in1, in2,
                              alpha, beta, 
                              c2, c3, c02, c03, 
                              K, N, S,
                              steps1, stepm1, stepl1,
                              b1, b2, b3,
                              fixed){
  
  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 0
  
  n3  <-  En3_binary_R(RRgo = RRgo, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                       p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, fixed = fixed)
  
  if(is.na(n3)){
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))  
  }
  else{
  
  n3  <- ceiling(n3)
  
  if(round(n3/2) != n3 / 2) {n3 <- n3 + 1}
  
  if(n2+n3>N){
    
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
    
  }else{
    
    pg    <-  Epgo_binary(RRgo = RRgo, n2 = n2, p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, fixed = fixed)
    
    K2    <-  c02 + c2 * n2         # cost phase II
    K3    <-  c03 * pg + c3 * n3    # cost phase III
    
    if(K2+K3>K){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
    }else{
      # probability of a successful program; small, medium, large effect size
      prob1 <-  EPsProg_binary_R(RRgo = RRgo, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                                 step1 = steps1, step2 =  steps2,
                                 p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                                 fixed = fixed)
      prob2 <-  EPsProg_binary_R(RRgo = RRgo, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                                 step1 =  stepm1, step2 =  stepm2,
                                 p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                                 fixed = fixed)
      prob3 <-  EPsProg_binary_R(RRgo = RRgo, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                                 step1 =  stepl1, step2 = stepl2,
                                 p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                                 fixed = fixed)
      
      SP    <-  prob1 + prob2 + prob3
      
      if(SP<S){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        
      }else{
        
        G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
        
        EU    <-  - K2 - K3 + G
        
        return(c(EU, n3, SP, pg, K2, K3, prob1, prob2, prob3))
      }
    }
  }
 }
}

# 2.2. conservative decision rule and sample size calculation: 
# use estimate with retetion factor
##############################################################################################

# prior distribution
# as above

#' @rdname Epgo_bias_binary 
#' @export
Epgo_binary_R2 <-  function(RRgo, n2, Adj, p0, w, p11, p12, in1, in2, fixed){
  if(fixed){
    return(
      pnorm((-log(p11/p0) + log(RRgo)/Adj)/sqrt((2/n2)*t1(p11, p0))) 
    )
  }else{
    return(
      integrate(function(x){
        sapply(x, function(x){
          pnorm((-log(x/p0) + log(RRgo)/Adj)/sqrt((2/n2)*t1(x, p0)))  *
            prior_binary(x, w, p11, p12, in1, in2)
        })
      },  0, 1)$value
    )
  }
}

#' @rdname En3_bias_binary 
#' @export
En3_binary_R2 <-  function(RRgo, n2, Adj, alpha, beta, p0, w, p11, p12, in1, in2, fixed){
  if(fixed){
    int = try(
      integrate(function(y){
        ((2*(qnorm(1-alpha)*t2(p11, p0)+qnorm(1-beta)*t3(p11, p0))^2)/(y*Adj)^2) *
          dnorm(y,
                mean = -log(p11/p0),
                sd = sqrt((2/n2)*t1(p11, p0)))
      }, - log(RRgo)/Adj, Inf), silent=TRUE)
    if(inherits(int ,'try-error')){
      warning(as.vector(int))
      integrated <- NA_real_
    } else {
      integrated <- int$value
    }
    return(integrated)
    
}else{
  int = try(
    integrate(function(x){
      sapply(x, function(x){
        integrate(function(y){
          ((2*(qnorm(1-alpha)*t2(x, p0)+qnorm(1-beta)*t3(x, p0))^2)/(y*Adj)^2) *
            dnorm(y,
                  mean = -log(x/p0),
                  sd = sqrt((2/n2)*t1(x, p0)))*
            prior_binary(x, w, p11, p12, in1, in2)
        }, - log(RRgo)/Adj, Inf)$value
      })
    }, 0, 1), silent=TRUE)
  if(inherits(int ,'try-error')){
    warning(as.vector(int))
    integrated <- NA_real_
  } else {
    integrated <- int$value
  }
  return(integrated)
}

} 

#' @rdname EPsProg_bias_binary
#' @export
EPsProg_binary_R2 <-  function(RRgo, n2, Adj, alpha, beta, step1, step2, p0, w, p11, p12, in1, in2, fixed){
  
  if(fixed){
    return(
      integrate(function(y){
        ( pnorm(qnorm(1 - alpha) -
                  log(step2)/sqrt((t1(p11, p0)*(y*Adj)^2)/(qnorm(1-alpha)*t2(p11, p0) +
                                                             qnorm(1-beta)*t3(p11, p0))^2),
                mean = -log((p11)/p0)/sqrt((t1(p11, p0)*(y*Adj)^2)/(qnorm(1-alpha)*t2(p11, p0) +
                                                                      qnorm(1-beta)*t3(p11, p0))^2),
                sd = 1) -
            pnorm(qnorm(1 - alpha) -
                    log(step1)/sqrt((t1(p11, p0)*(y*Adj)^2)/(qnorm(1-alpha)*t2(p11, p0) +
                                                               qnorm(1-beta)*t3(p11, p0))^2),
                  mean = -log((p11)/p0)/sqrt((t1(p11, p0)*(y*Adj)^2)/(qnorm(1-alpha)*t2(p11, p0) +
                                                                        qnorm(1-beta)*t3(p11, p0))^2),
                  sd = 1) ) *
          dnorm(y,
                mean = -log(p11/p0),
                sd = sqrt((2/n2)*t1(p11, p0)))
      }, -log(RRgo)/Adj, Inf)$value
    )
  }else{
    return(
      integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            ( pnorm(qnorm(1 - alpha) -
                      log(step2)/sqrt((t1(x, p0)*(y*Adj)^2)/(qnorm(1-alpha)*t2(x, p0) +
                                                               qnorm(1-beta)*t3(x, p0))^2),
                    mean = -log((x)/p0)/sqrt((t1(x, p0)*(y*Adj)^2)/(qnorm(1-alpha)*t2(x, p0) +
                                                                      qnorm(1-beta)*t3(x, p0))^2),
                    sd = 1) -
                pnorm(qnorm(1 - alpha) -
                        log(step1)/sqrt((t1(x, p0)*(y*Adj)^2)/(qnorm(1-alpha)*t2(x, p0) +
                                                                 qnorm(1-beta)*t3(x, p0))^2),
                      mean = -log((x)/p0)/sqrt((t1(x, p0)*(y*Adj)^2)/(qnorm(1-alpha)*t2(x, p0) +
                                                                        qnorm(1-beta)*t3(x, p0))^2),
                      sd = 1) ) *
              dnorm(y,
                    mean = -log(x/p0),
                    sd = sqrt((2/n2)*t1(x, p0))) *
              prior_binary(x, w, p11, p12, in1, in2)
          }, -log(RRgo)/Adj, Inf)$value
        })
      }, 0, 1)$value
    )
  }
  
}

#' @rdname utility_bias_binary
#' @export
utility_binary_R2 <-  function(n2, RRgo, Adj, w, p0, p11, p12, in1, in2,
                               alpha, beta, 
                               c2, c3, c02, c03, 
                               K, N, S,
                               steps1, stepm1, stepl1,
                               b1, b2, b3,
                               fixed){
  
 
  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 0
  
  n3  <-  En3_binary_R2(RRgo = RRgo, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                        p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, fixed = fixed)
  
  if(is.na(n3)){
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))  
  }
  else{
  
   n3  <- ceiling(n3)
  
  if(round(n3/2) != n3 / 2) {n3 <- n3 + 1}
  
  if(n2+n3>N){
    
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
    
  }else{
    
    pg    <-  Epgo_binary_R2(RRgo = RRgo, Adj = Adj,  n2 = n2, p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, fixed = fixed)
    
    K2    <-  c02 + c2 * n2         # cost phase II
    K3    <-  c03 * pg + c3 * n3    # cost phase III
    
    if(K2+K3>K){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
    }else{
      # probability of a successful program; small, medium, large effect size
      prob1 <-  EPsProg_binary_R2(RRgo = RRgo, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                                  step1 = steps1, step2 =  steps2,
                                  p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                                  fixed = fixed)
      prob2 <-  EPsProg_binary_R2(RRgo = RRgo, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                                  step1 =  stepm1, step2 =  stepm2,
                                  p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                                  fixed = fixed)
      prob3 <-  EPsProg_binary_R2(RRgo = RRgo, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                                  step1 =  stepl1, step2 = stepl2,
                                  p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                                  fixed = fixed)
      
      SP    <-  prob1 + prob2 + prob3
      
      if(SP<S){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        
      }else{
        
        G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
        
        EU    <-  - K2 - K3 + G
        
        return(c(EU, n3, SP, pg, K2, K3, prob1, prob2, prob3))
      }
    }
  }
 }
}
