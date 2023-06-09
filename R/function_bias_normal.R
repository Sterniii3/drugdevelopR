# 1.1. conservative sample size calculation: use lower bound of one-sided confidence intervall
##############################################################################################

# prior distribution
# as above

# expected probability to go to phase III
# as above

#' Expected sample size for phase III for bias adjustment programs and normally distributed outcomes
#' 
#' To discount for overoptimistic results in phase II when calculating the optimal sample size in phase III, 
#' it is necessary to use the functions `En3_normal_L()`, `En3_normal_L2()`, `En3_normal_R()` and `En3_normal_R2()`.
#' Each function describes a specific case:
#' - `En3_normal_L()`: calculates the optimal sample size for an additive adjustment factor (i.e. adjust the lower bound of the one-sided confidence interval), 
#' however the go-decision is not affected by the bias adjustment
#' - `En3_normal_L2()`: calculates the optimal sample size for an additive adjustment factor (i.e. adjust the lower bound of the one-sided confidence interval)
#' when the go-decision is also affected by the bias adjustment
#' - `En3_normal_R()`: calculates the optimal sample size for a multiplicative adjustment factor (i.e. use estimate with a retention factor), 
#' however the go-decision is not affected by the bias adjustment
#' - `En3_normal_R2()`: calculates the optimal sample size for a multiplicative adjustment factor (i.e. use estimate with a retention factor)
#' when the go-decision is also affected by the bias adjustment 
#' @param kappa threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param Adj adjustment parameter
#' @param alpha significance level
#' @param beta `1 - beta` is the power for calculation of sample size for phase III
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
#'                               alpha = 0.025, beta = 0.1, w = 0.3,
#'                               Delta1 = 0.375, Delta2 = 0.625, 
#'                               in1 = 300, in2 = 600, 
#'                               a = 0.25, b = 0.75, fixed = FALSE)
#'           res <- En3_normal_L2(kappa = 0.1, n2 = 50, Adj = 0, 
#'                               alpha = 0.025, beta = 0.1, w = 0.3,
#'                               Delta1 = 0.375, Delta2 = 0.625, 
#'                               in1 = 300, in2 = 600, 
#'                               a = 0.25, b = 0.75, fixed = TRUE)
#'           res <- En3_normal_R(kappa = 0.1, n2 = 50, Adj = 1, 
#'                               alpha = 0.025, beta = 0.1, w = 0.3,
#'                               Delta1 = 0.375, Delta2 = 0.625, 
#'                               in1 = 300, in2 = 600, 
#'                               a = 0.25, b = 0.75, fixed = FALSE)
#'           res <- En3_normal_R2(kappa = 0.1, n2 = 50, Adj = 1, 
#'                               alpha = 0.025, beta = 0.1, w = 0.3,
#'                               Delta1 = 0.375, Delta2 = 0.625, 
#'                               in1 = 300, in2 = 600, 
#'                               a = 0.25, b = 0.75, fixed = FALSE)
#' @name En3_bias_normal                           
#' @export
#' @keywords internal

En3_normal_L <-  function(kappa, n2, Adj, alpha, beta, w, Delta1, Delta2, in1, in2, a, b, fixed){
   
  if(fixed){
    int   = try(integrate(function(y){
    sapply(y,function(y){
      ( (4*(qnorm(1-alpha)+qnorm(1-beta))^2)/((y-qnorm(1-Adj)*sqrt(4/n2))^2))*
        dnorm(y,
              mean=Delta1,
              sd=sqrt(4/n2))
    })
  }, kappa,Inf),silent=TRUE)
  if(inherits(int ,'try-error')){
    warning(as.vector(int))
    integrated <- NA_real_
  } else {
    integrated <- int$value
  }
  return(integrated)
  
}else{
  int   = try(integrate(function(x){
    sapply(x,function(x){
      integrate(function(y){
        ( (4*(qnorm(1-alpha)+qnorm(1-beta))^2)/((y-qnorm(1-Adj)*sqrt(4/n2))^2))*
          dnorm(y,
                mean=x,
                sd=sqrt(4/n2))*
          prior_normal(x, w, Delta1, Delta2, in1, in2, a, b)
      }, kappa,Inf)$value
    })
  }, -Inf, Inf),silent=TRUE)
  if(inherits(int ,'try-error')){
    warning(as.vector(int))
    integrated <- NA_real_
  } else {
    integrated <- int$value
  }
  return(integrated)
}
}

#' Expected probability of a successful program for bias adjustment programs with normally distributed outcomes
#' 
#' To discount for overoptimistic results in phase II when calculating the optimal sample size in phase III, 
#' it is necessary to use the following functions, which each describe a specific case:
#' - `EPsProg_normal_L()`: calculates the expected probability of a successful for an additive adjustment factor (i.e. adjust the lower bound of the one-sided confidence interval), 
#' however the go-decision is not affected by the bias adjustment
#' - `EPsProg_normal_L2()`: calculates the expected probability of a successful for an additive adjustment factor (i.e. adjust the lower bound of the one-sided confidence interval)
#' when the go-decision is also affected by the bias adjustment
#' - `EPsProg_normal_R()`: calculates the expected probability of a successful for a multiplicative adjustment factor (i.e. use estimate with a retention factor), 
#' however the go-decision is not affected by the bias adjustment
#' - `EPsProg_normal_R2()`: calculates the expected probability of a successful for a multiplicative adjustment factor (i.e. use estimate with a retention factor)
#' when the go-decision is also affected by the bias adjustment 
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
#'                                  alpha = 0.025, beta = 0.1, w = 0.3,
#'                                  step1 = 0, step2 = 0.5,
#'                                  Delta1 = 0.375, Delta2 = 0.625, 
#'                                  in1 = 300, in2 = 600, 
#'                                  a = 0.25, b = 0.75, fixed = FALSE)
#'           res <- EPsProg_normal_L2(kappa = 0.1, n2 = 50, Adj = 0, 
#'                                  alpha = 0.025, beta = 0.1, w = 0.3,
#'                                  step1 = 0, step2 = 0.5,
#'                                  Delta1 = 0.375, Delta2 = 0.625, 
#'                                  in1 = 300, in2 = 600, 
#'                                  a = 0.25, b = 0.75, fixed = FALSE)
#'           res <- EPsProg_normal_R(kappa = 0.1, n2 = 50, Adj = 1, 
#'                                  alpha = 0.025, beta = 0.1, w = 0.3,
#'                                  step1 = 0, step2 = 0.5,
#'                                  Delta1 = 0.375, Delta2 = 0.625, 
#'                                  in1 = 300, in2 = 600, 
#'                                  a = 0.25, b = 0.75, fixed = FALSE)
#'           res <- EPsProg_normal_R2(kappa = 0.1, n2 = 50, Adj = 1, 
#'                                  alpha = 0.025, beta = 0.1, w = 0.3,
#'                                  step1 = 0, step2 = 0.5,
#'                                  Delta1 = 0.375, Delta2 = 0.625, 
#'                                  in1 = 300, in2 = 600, 
#'                                  a = 0.25, b = 0.75, fixed = FALSE)
#' @name EPsProg_bias_normal                                
#' @export
#' @keywords internal

EPsProg_normal_L <-  function(kappa, n2, Adj, alpha, beta, step1, step2, w, Delta1, Delta2, in1, in2, a, b, fixed){
  
  c = (qnorm(1 - alpha) + qnorm(1 - beta))^2
  
  if(fixed){
    return(
      integrate(function(y){
        ( pnorm(qnorm(1 - alpha) + step2/sqrt((y-qnorm(1-Adj)*sqrt(4/n2))^2/c),
                mean = (Delta1)/sqrt((y-qnorm(1-Adj)*sqrt(4/n2))^2/c),
                sd = 1) -
            pnorm(qnorm(1 - alpha) + step1/sqrt((y-qnorm(1-Adj)*sqrt(4/n2))^2/c),
                  mean = (Delta1)/sqrt((y-qnorm(1-Adj)*sqrt(4/n2))^2/c),
                  sd = 1) ) *
          dnorm(y,
                mean = Delta1,
                sd = sqrt(4/n2)) 
      }, kappa, Inf)$value
    )
  }else{
    return(
      integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            ( pnorm(qnorm(1 - alpha) + step2/sqrt((y-qnorm(1-Adj)*sqrt(4/n2))^2/c),
                    mean = (x)/sqrt((y-qnorm(1-Adj)*sqrt(4/n2))^2/c),
                    sd = 1) -
                pnorm(qnorm(1 - alpha) + step1/sqrt((y-qnorm(1-Adj)*sqrt(4/n2))^2/c),
                      mean = (x)/sqrt((y-qnorm(1-Adj)*sqrt(4/n2))^2/c),
                      sd = 1) ) *
              dnorm(y,
                    mean = x,
                    sd = sqrt(4/n2)) *
              prior_normal(x, w, Delta1, Delta2, in1, in2, a, b)
          }, kappa, Inf)$value
        })
      },  - Inf, Inf)$value
    )
  }
  
}

#' Utility function for bias adjustment programs with normally distributed outcomes.
#' 
#' The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters and the expected probability of a successful program. 
#' The utility is in a further step maximized by the `optimal_bias_normal()` function.
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
#'                                  alpha = 0.025, beta = 0.1, w = 0.3,
#'                                  Delta1 = 0.375, Delta2 = 0.625, 
#'                                  in1 = 300, in2 = 600, 
#'                                  a = 0.25, b = 0.75, 
#'                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                  K = Inf, N = Inf, S = -Inf, 
#'                                  steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
#'                                  b1 = 3000, b2 = 8000, b3 = 10000, 
#'                                  fixed = TRUE)
#'           res <- utility_normal_L2(kappa = 0.1, n2 = 50, Adj = 0, 
#'                                  alpha = 0.025, beta = 0.1, w = 0.3,
#'                                  Delta1 = 0.375, Delta2 = 0.625, 
#'                                  in1 = 300, in2 = 600, 
#'                                  a = 0.25, b = 0.75, 
#'                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                  K = Inf, N = Inf, S = -Inf, 
#'                                  steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
#'                                  b1 = 3000, b2 = 8000, b3 = 10000, 
#'                                  fixed = TRUE)
#'           res <- utility_normal_R(kappa = 0.1, n2 = 50, Adj = 1, 
#'                                  alpha = 0.025, beta = 0.1, w = 0.3,
#'                                  Delta1 = 0.375, Delta2 = 0.625, 
#'                                  in1 = 300, in2 = 600, 
#'                                  a = 0.25, b = 0.75, 
#'                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                  K = Inf, N = Inf, S = -Inf, 
#'                                  steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
#'                                  b1 = 3000, b2 = 8000, b3 = 10000, 
#'                                  fixed = TRUE)
#'           res <- utility_normal_R2(kappa = 0.1, n2 = 50, Adj = 1, 
#'                                  alpha = 0.025, beta = 0.1, w = 0.3,
#'                                  Delta1 = 0.375, Delta2 = 0.625, 
#'                                  in1 = 300, in2 = 600, 
#'                                  a = 0.25, b = 0.75, 
#'                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                  K = Inf, N = Inf, S = -Inf, 
#'                                  steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
#'                                  b1 = 3000, b2 = 8000, b3 = 10000, 
#'                                  fixed = TRUE)
#' @name utility_bias_normal                               
#' @export
#' @keywords internal

utility_normal_L <-  function(n2, kappa, Adj, w, Delta1, Delta2, in1, in2, a, b,
                            alpha, beta, 
                            c2, c3, c02, c03, 
                            K, N, S,
                            steps1, stepm1, stepl1,
                            b1, b2, b3,
                            fixed){
  
  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 1
  
  n3  <-  En3_normal_L(kappa = kappa, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                     w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                     fixed = fixed)
  
  if(is.na(n3)){
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))  
  }
  else{ 
   
  n3  <- ceiling(n3)
  
  if(round(n3/2) != n3 / 2) {n3 = n3 + 1}
  
  if(n2+n3>N){
    
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
    
  }else{
    
    pg    <-  Epgo_normal(kappa = kappa, n2 = n2, 
                          w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                          fixed = fixed)
    
    K2    <-  c02 + c2 * n2  #cost phase II
    K3    <-  c03 * pg + c3 * n3  #cost phase III
    
    if(K2+K3>K){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      #output: expected utility Eud, En3, EsP, Epgo
      
    }else{
      # probability of a successful program; small, medium, large effect size
      prob1 <-  EPsProg_normal_L(kappa = kappa, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                               step1 = steps1, step2 =  steps2,
                               w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                               fixed = fixed)
      prob2 <-  EPsProg_normal_L(kappa = kappa, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                               step1 =  stepm1, step2 =  stepm2,
                               w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                               fixed = fixed)
      prob3 <-  EPsProg_normal_L(kappa = kappa, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                               step1 =  stepl1, step2 = stepl2,
                               w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                               fixed = fixed)
      
      SP    <-  prob1 + prob2 + prob3
      
      if(SP<S){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        
      }else{
        
        G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
        
        EU    <-  - K2 - K3 + G
        
        return(c(EU, n3, SP, pg, K2, K3, prob1, prob2, prob3))
        #output: expected utility Eud, En3, EsP, Epgo, cost phase II and III
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


#' Expected probability to go to phase III for bias adjustment programs with normally distributed outcomes
#' 
#' In the case we do not only want do discount for overoptimistic results in phase II when calculating the sample size in phase III, 
#' but also when deciding whether to go to phase III or not the functions `Epgo_normal_L2` and `Epgo_normal_R2` are necessary.
#' The function `Epgo_normal_L2` uses an additive adjustment parameter (i.e. adjust the lower bound of the one-sided confidence interval),
#' the function `Epgo_normal_R2` uses a multiplicative adjustment parameter (i.e. use estimate with a retention factor)
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
#'                                Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#'                                a = 0.25, b = 0.75, fixed = FALSE)
#'           res <- Epgo_normal_R2(kappa = 0.1, n2 = 50, Adj = 1, w = 0.3,
#'                                Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#'                                a = 0.25, b = 0.75, fixed = FALSE)
#' @name Epgo_bias_normal                              
#' @export
#' @keywords internal
Epgo_normal_L2 <-  function(kappa, n2, Adj, w, Delta1, Delta2, in1, in2, a, b, fixed){
  
  if(fixed){
    return(
      pnorm((Delta1-kappa-qnorm(1-Adj)*sqrt(4/n2))/sqrt(4/n2))  
    )
  }else{
    return(
      integrate(function(x){
        sapply(x, function(x){
          pnorm((x-kappa-qnorm(1-Adj)*sqrt(4/n2))/sqrt(4/n2)) *
            prior_normal(x, w, Delta1, Delta2, in1, in2, a, b)
        })
      },  - Inf, Inf)$value 
    )
  }
}


#' @rdname En3_bias_normal 
#' @export

En3_normal_L2 <-  function(kappa, n2, Adj, alpha, beta, w, Delta1, Delta2, in1, in2, a, b, fixed){
  
  if(fixed){
    int = try(
      integrate(function(y){
        ((4*(qnorm(1-alpha)+qnorm(1-beta))^2)/(y-qnorm(1-Adj)*sqrt(4/n2))^2) *
          dnorm(y,
                mean = Delta1,
                sd = sqrt(4/n2))
      }, kappa+qnorm(1-Adj)*sqrt(4/n2), Inf), silent=TRUE)
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
            ((4*(qnorm(1-alpha)+qnorm(1-beta))^2)/(y-qnorm(1-Adj)*sqrt(4/n2))^2) *
              dnorm(y,
                    mean = x,
                    sd = sqrt(4/n2))*
              prior_normal(x, w, Delta1, Delta2, in1, in2, a, b)
          }, kappa+qnorm(1-Adj)*sqrt(4/n2), Inf)$value
        })
      },  - Inf, Inf), silent=TRUE)
          if(inherits(int ,'try-error')){
            warning(as.vector(int))
            integrated <- NA_real_
          } else {
            integrated <- int$value
          }
          return(integrated)
  }
  
}

#' @rdname EPsProg_bias_normal 
#' @export
EPsProg_normal_L2 <-  function(kappa, n2, Adj, alpha, beta, step1, step2, w, Delta1, Delta2, in1, in2, a, b, fixed){
  
  c = (qnorm(1 - alpha) + qnorm(1 - beta))^2
  
  if(fixed){
    return(
      integrate(function(y){
        ( pnorm(qnorm(1 - alpha) + step2/sqrt((y-qnorm(1-Adj)*sqrt(4/n2))^2/c),
                mean = (Delta1)/sqrt((y-qnorm(1-Adj)*sqrt(4/n2))^2/c),
                sd = 1) -
            pnorm(qnorm(1 - alpha) + step1/sqrt((y-qnorm(1-Adj)*sqrt(4/n2))^2/c),
                  mean = (Delta1)/sqrt((y-qnorm(1-Adj)*sqrt(4/n2))^2/c),
                  sd = 1) ) *
          dnorm(y,
                mean = Delta1,
                sd = sqrt(4/n2)) 
      }, kappa+qnorm(1-Adj)*sqrt(4/n2), Inf)$value
    )
  }else{
    return(
      integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            ( pnorm(qnorm(1 - alpha) + step2/sqrt((y-qnorm(1-Adj)*sqrt(4/n2))^2/c),
                    mean = (x)/sqrt((y-qnorm(1-Adj)*sqrt(4/n2))^2/c),
                    sd = 1) -
                pnorm(qnorm(1 - alpha) + step1/sqrt((y-qnorm(1-Adj)*sqrt(4/n2))^2/c),
                      mean = (x)/sqrt((y-qnorm(1-Adj)*sqrt(4/n2))^2/c),
                      sd = 1) ) *
              dnorm(y,
                    mean = x,
                    sd = sqrt(4/n2)) *
              prior_normal(x, w, Delta1, Delta2, in1, in2, a, b)
          }, kappa +qnorm(1-Adj)*sqrt(4/n2), Inf)$value
        })
      },  - Inf, Inf)$value
    )
  }
  
}

#' @rdname utility_bias_normal 
#' @keywords internal
#' @export
utility_normal_L2 <-  function(n2, kappa, Adj, w, Delta1, Delta2, in1, in2, a, b,
                            alpha, beta, 
                            c2, c3, c02, c03, 
                            K, N, S,
                            steps1, stepm1, stepl1,
                            b1, b2, b3,
                            fixed){
  
  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 1
  
  n3  <-  En3_normal_L2(kappa = kappa, Adj=Adj, n2 = n2, alpha = alpha, beta = beta,
                     w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                     fixed = fixed)
  if(is.na(n3)){
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))  
  }
  else{
  
  n3  <- ceiling(n3)
  
  if(round(n3/2) != n3 / 2) {n3 = n3 + 1}
  
  if(n2+n3>N){
    
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
    
  }else{
    
    pg    <-  Epgo_normal_L2(kappa = kappa, n2 = n2, Adj=Adj,
                          w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                          fixed = fixed)
    
    K2    <-  c02 + c2 * n2  #cost phase II
    K3    <-  c03 * pg + c3 * n3  #cost phase III
    
    if(K2+K3>K){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      #output: expected utility Eud, En3, EsP, Epgo
      
    }else{
      # probability of a successful program; small, medium, large effect size
      prob1 <-  EPsProg_normal_L2(kappa = kappa, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                               step1 = steps1, step2 =  steps2,
                               w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                               fixed = fixed)
      prob2 <-  EPsProg_normal_L2(kappa = kappa, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                               step1 =  stepm1, step2 =  stepm2,
                               w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                               fixed = fixed)
      prob3 <-  EPsProg_normal_L2(kappa = kappa, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                               step1 =  stepl1, step2 = stepl2,
                               w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                               fixed = fixed)
      
      SP    <-  prob1 + prob2 + prob3
      
      if(SP<S){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        
      }else{
        
        G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
        
        EU    <-  - K2 - K3 + G
        
        return(c(EU, n3, SP, pg, K2, K3, prob1, prob2, prob3))
        #output: expected utility Eud, En3, EsP, Epgo, cost phase II and III
      }
    }
  }
 }
}

# 2.1. conservative sample size calculation: use estimate with retention factor
##############################################################################################

# prior distribution
# as above

# expected probability to go to phase III
# as above

#' @rdname En3_bias_normal 
#' @keywords internal
#' @export

En3_normal_R <-  function(kappa, n2, Adj, alpha, beta, w, Delta1, Delta2, in1, in2, a, b, fixed){
  
  if(fixed){
    
    
    int   = try(integrate(function(y){
      sapply(y,function(y){
        ( (4*(qnorm(1-alpha)+qnorm(1-beta))^2)/((y*Adj)^2))*
          dnorm(y,
                mean=Delta1,
                sd=sqrt(4/n2)) 
      })
    }, kappa,Inf),silent=TRUE)
    if(inherits(int ,'try-error')){
      warning(as.vector(int))
      integrated <- NA_real_
    } else {
      integrated <- int$value
    }
    return(integrated)
    
  }else{
    int   = try(integrate(function(x){
      sapply(x,function(x){
        integrate(function(y){
          ( (4*(qnorm(1-alpha)+qnorm(1-beta))^2)/((y*Adj)^2))*
            dnorm(y,
                  mean=x,
                  sd=sqrt(4/n2))*
            prior_normal(x, w, Delta1, Delta2, in1, in2, a, b) 
        }, kappa,Inf)$value
      })
    }, -Inf, Inf),silent=TRUE)
    if(inherits(int ,'try-error')){
      warning(as.vector(int))
      integrated <- NA_real_
    } else {
      integrated <- int$value
    }
    return(integrated)
  }
}

#' @rdname EPsProg_bias_normal 
#' @keywords internal
#' @export
EPsProg_normal_R <-  function(kappa, n2, Adj, alpha, beta, step1, step2, w, Delta1, Delta2, in1, in2, a, b, fixed){
  
  c = (qnorm(1 - alpha) + qnorm(1 - beta))^2
  
  if(fixed){
    return(
      integrate(function(y){
        ( pnorm(qnorm(1 - alpha) + step2/sqrt((y*Adj)^2/c),
                mean = (Delta1)/sqrt((y*Adj)^2/c),
                sd = 1) -
            pnorm(qnorm(1 - alpha) + step1/sqrt((y*Adj)^2/c),
                  mean = (Delta1)/sqrt((y*Adj)^2/c),
                  sd = 1) ) *
          dnorm(y,
                mean = Delta1,
                sd = sqrt(4/n2)) 
      }, kappa, Inf)$value
    )
  }else{
    return(
      integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            ( pnorm(qnorm(1 - alpha) + step2/sqrt((y*Adj)^2/c),
                    mean = (x)/sqrt((y*Adj)^2/c),
                    sd = 1) -
                pnorm(qnorm(1 - alpha) + step1/sqrt((y*Adj)^2/c),
                      mean = (x)/sqrt((y*Adj)^2/c),
                      sd = 1) ) *
              dnorm(y,
                    mean = x,
                    sd = sqrt(4/n2)) *
              prior_normal(x, w, Delta1, Delta2, in1, in2, a, b)
          }, kappa, Inf)$value
        })
      },  - Inf, Inf)$value
    )
  }
  
}

#' @rdname utility_bias_normal 
#' @keywords internal
#' @export
utility_normal_R <-  function(n2, kappa, Adj, w, Delta1, Delta2, in1, in2, a, b,
                            alpha, beta, 
                            c2, c3, c02, c03, 
                            K, N, S,
                            steps1, stepm1, stepl1,
                            b1, b2, b3,
                            fixed){
  
  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- Inf
  
  n3  <-  En3_normal_R(kappa = kappa, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                     w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                     fixed = fixed)
  
  n3  <- ceiling(n3)
  
  if(is.na(n3)){
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))  
  }
  else{
  
  if(round(n3/2) != n3 / 2) {n3 = n3 + 1}
  
  if(n2+n3>N){
    
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
    
  }else{
    
    pg    <-  Epgo_normal(kappa = kappa, n2 = n2,
                          w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                          fixed = fixed)
    
    K2    <-  c02 + c2 * n2  #cost phase II
    K3    <-  c03 * pg + c3 * n3  #cost phase III
    
    if(K2+K3>K){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      #output: expected utility Eud, En3, EsP, Epgo
      
    }else{
      # probability of a successful program; small, medium, large effect size
      prob1 <-  EPsProg_normal_R(kappa = kappa, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                               step1 = steps1, step2 =  steps2,
                               w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                               fixed = fixed)
      prob2 <-  EPsProg_normal_R(kappa = kappa, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                               step1 =  stepm1, step2 =  stepm2,
                               w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                               fixed = fixed)
      prob3 <-  EPsProg_normal_R(kappa = kappa, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                               step1 =  stepl1, step2 = stepl2,
                               w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                               fixed = fixed)
      
      SP    <-  prob1 + prob2 + prob3
      
      if(SP<S){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        
      }else{
        
        G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
        
        EU    <-  - K2 - K3 + G
        
        return(c(EU, n3, SP, pg, K2, K3, prob1, prob2, prob3))
        #output: expected utility Eud, En3, EsP, Epgo, cost phase II and III
      }
    }
  }
 }
}

# 2.2. conservative decision rule and sample size calculation: 
# use estimate with retention factor
##############################################################################################

# prior distribution
# as above

#' @rdname Epgo_bias_normal 
#' @keywords internal
#' @export
Epgo_normal_R2 <-  function(kappa, n2, Adj, w, Delta1, Delta2, in1, in2, a, b, fixed){
  
  if(fixed){
    return(
      pnorm((Delta1 - kappa/Adj)/sqrt(4/n2))  
    )
  }else{
    return(
      integrate(function(x){
        sapply(x, function(x){
          pnorm((x-kappa/Adj)/sqrt(4/n2)) *
            prior_normal(x, w, Delta1, Delta2, in1, in2, a, b)
        })
      },  - Inf, Inf)$value
    )
  }
}

#' @rdname En3_bias_normal 
#' @keywords internal
#' @export
En3_normal_R2 <-  function(kappa, n2, Adj, alpha, beta, w, Delta1, Delta2, in1, in2, a, b, fixed){
  
  if(fixed){
     int = try(
      integrate(function(y){
        ((4*(qnorm(1-alpha)+qnorm(1-beta))^2)/(y*Adj)^2) *
          dnorm(y,
                mean = Delta1,
                sd = sqrt(4/n2))
      }, kappa/Adj, Inf),silent=TRUE)
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
            ((4*(qnorm(1-alpha)+qnorm(1-beta))^2)/(y*Adj)^2) *
              dnorm(y,
                    mean = x,
                    sd = sqrt(4/n2))*
              prior_normal(x, w, Delta1, Delta2, in1, in2, a, b)
          }, kappa/Adj, Inf)$value
        })
      },  - Inf, Inf),silent=TRUE)
    if(inherits(int ,'try-error')){
      warning(as.vector(int))
      integrated <- NA_real_
    } else {
      integrated <- int$value
    }
    return(integrated)
  }
}

#' @rdname EPsProg_bias_normal 
#' @keywords internal
#' @export
EPsProg_normal_R2 <-  function(kappa, n2, Adj, alpha, beta, step1, step2, w, Delta1, Delta2, in1, in2, a, b, fixed){
  
  c = (qnorm(1 - alpha) + qnorm(1 - beta))^2
  
  if(fixed){
    return(
      integrate(function(y){
        ( pnorm(qnorm(1 - alpha) + step2/sqrt((y*Adj)^2/c),
                mean = (Delta1)/sqrt((y*Adj)^2/c),
                sd = 1) -
            pnorm(qnorm(1 - alpha) + step1/sqrt((y*Adj)^2/c),
                  mean = (Delta1)/sqrt((y*Adj)^2/c),
                  sd = 1) ) *
          dnorm(y,
                mean = Delta1,
                sd = sqrt(4/n2)) 
      }, kappa/Adj, Inf)$value
    ) 

  
  }else{
    return(
      integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            ( pnorm(qnorm(1 - alpha) + step2/sqrt((y*Adj)^2/c),
                    mean = (x)/sqrt((y*Adj)^2/c),
                    sd = 1) -
                pnorm(qnorm(1 - alpha) + step1/sqrt((y*Adj)^2/c),
                      mean = (x)/sqrt((y*Adj)^2/c),
                      sd = 1) ) *
              dnorm(y,
                    mean = x,
                    sd = sqrt(4/n2)) *
              prior_normal(x, w, Delta1, Delta2, in1, in2, a, b)
          }, kappa/Adj, Inf)$value
        })
      },  - Inf, Inf)$value
    ) 
  }
  
}


#' @rdname utility_bias_normal 
#' @keywords internal
#' @export
utility_normal_R2 <-  function(n2, kappa, Adj,  w, Delta1, Delta2, in1, in2, a, b,
                            alpha, beta, 
                            c2, c3, c02, c03, 
                            K, N, S,
                            steps1, stepm1, stepl1,
                            b1, b2, b3,
                            fixed){
  
  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 1
  
  n3  <-  En3_normal_R2(kappa = kappa, Adj = Adj, n2 = n2, alpha = alpha, beta = beta,
                     w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                     fixed = fixed)
  
  if(is.na(n3)){
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))  
  }
  else{
     
  n3  <- ceiling(n3)
  
  if(round(n3/2) != n3 / 2) {n3 = n3 + 1}
  
  if(n2+n3>N){
    
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
    
  }else{
    
    pg    <-  Epgo_normal_R2(kappa = kappa, n2 = n2, Adj = Adj, 
                          w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                          fixed = fixed)
    
    K2    <-  c02 + c2 * n2  #cost phase II
    K3    <-  c03 * pg + c3 * n3  #cost phase III
    
    if(K2+K3>K){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      #output: expected utility Eud, En3, EsP, Epgo
      
    }else{
      # probability of a successful program; small, medium, large effect size
      prob1 <-  EPsProg_normal_R2(kappa = kappa, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                               step1 = steps1, step2 =  steps2,
                               w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                               fixed = fixed)
      prob2 <-  EPsProg_normal_R2(kappa = kappa, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                               step1 =  stepm1, step2 =  stepm2,
                               w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                               fixed = fixed)
      prob3 <-  EPsProg_normal_R2(kappa = kappa, n2 = n2, Adj = Adj, alpha = alpha, beta = beta,
                               step1 =  stepl1, step2 = stepl2,
                               w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                               fixed = fixed)
      
      SP    <-  prob1 + prob2 + prob3
      
      if(SP<S){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        
      }else{
        
        G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
        
        EU    <-  - K2 - K3 + G
        
        return(c(EU, n3, SP, pg, K2, K3, prob1, prob2, prob3))
        #output: expected utility Eud, En3, EsP, Epgo, cost phase II and III
      }
    }
  }
 }
}


