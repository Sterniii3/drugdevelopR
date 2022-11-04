
# 1.1. conservative sample size calculation: use lower bound of one-sided confidence intervall
##############################################################################################

# prior distribution
# as above

# expected probability to go to phase III
# as above

# expected number of events for phase III when going to phase III

#' Expected sample size for phase III for bias adjustment programs and time-to-event outcomes
#' 
#' To discount for overoptimistic results in phase II when calculating the optimal sample size in phase III, 
#' it is necessary to use the functions `Ed3_L()`, `Ed3_L2()`, `Ed3_R()` and `Ed3_R2()`.
#' Each function describes a specific case:
#' - `Ed3_L()`: calculates the optimal sample size for an additive adjustment factor (i.e adjust the lower bound of the one-sided confidence interval), 
#' however the go-decision is not affected by the bias adjustment
#' - `Ed3_L2()`: calculates the optimal sample size for an additive adjustment factor (i.e adjust the lower bound of the one-sided confidence interval)
#' when the go-decision is also affected by the bias adjustment
#' - `Ed3_R()`: calculates the optimal sample size for a multiplicative adjustment factor (i.e. use estimate with a retention factor), 
#' however the go-decision is not affected by the bias adjustment
#' - `Ed3_R2()`: calculates the optimal sample size for a multiplicative adjustment factor (i.e. use estimate with a retention factor)
#' when the go-decision is also affected by the bias adjustment 
#' @param HRgo threshold value for the go/no-go decision rule
#' @param d2 total events for phase II; must be even number
#' @param Adj adjustment parameter
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param w weight for mixture prior distribution
#' @param hr1 first assumed true treatment effect on HR scale for prior distribution
#' @param hr2 second assumed true treatment effect on HR scale for prior distribution
#' @param id1 amount of information for `hr1` in terms of number of events
#' @param id2 amount of information for `hr2` in terms of number of events
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `hr1` is used as fixed effect
#' @return The output of the the functions `Ed3_L`, `Ed3_L2`, `Ed3_R` and `Ed3_R2` is the expected number of participants in phase III. 
#' @examples res <-  Ed3_L(HRgo = 0.8, d2 = 50, Adj = 0.4,
#'                         alpha = 0.025, beta = 0.1, w = 0.3, 
#'                         hr1 =  0.69, hr2 = 0.81, 
#'                         id1 = 280, id2 = 420, fixed = FALSE)
#'           res <-  Ed3_L2(HRgo = 0.8, d2 = 50, Adj = 0.4,
#'                         alpha = 0.025, beta = 0.1, w = 0.3, 
#'                         hr1 =  0.69, hr2 = 0.81, 
#'                         id1 = 280, id2 = 420, fixed = FALSE)
#'           res <- Ed3_R(HRgo = 0.8, d2 = 50, Adj = 0.9,
#'                         alpha = 0.025, beta = 0.1, w = 0.3, 
#'                         hr1 =  0.69, hr2 = 0.81, 
#'                         id1 = 280, id2 = 420, fixed = FALSE)
#'           res <- Ed3_R2(HRgo = 0.8, d2 = 50, Adj = 0.9,
#'                         alpha = 0.025, beta = 0.1, w = 0.3, 
#'                         hr1 =  0.69, hr2 = 0.81, 
#'                         id1 = 280, id2 = 420, fixed = FALSE)
#'                               
#' @name Ed3_bias                            
#' @export
#' @editor Johannes Cepicka
#' @editDate 2022-11-04
Ed3_L<-function(HRgo, d2, Adj, alpha, beta, w, hr1, hr2, id1, id2, fixed){
  
  if(fixed){
    
    theta <- -log(hr1)
    
    int   = try(integrate(function(y){
      sapply(y,function(y){
        ( (4*(qnorm(1-alpha)+qnorm(1-beta))^2)/((y-qnorm(1-Adj)*sqrt(4/d2))^2))*
          dnorm(y,
                mean=theta,
                sd=sqrt(4/d2))
      })
    }, -log(HRgo),Inf),silent=TRUE)
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
            ( (4*(qnorm(1-alpha)+qnorm(1-beta))^2)/((y-qnorm(1-Adj)*sqrt(4/d2))^2))*
              dnorm(y,
                    mean=x,
                    sd=sqrt(4/d2))*
              prior_tte(x,w,hr1,hr2,id1,id2) 
        }, -log(HRgo),Inf)$value
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

#' Expected probability of a successful program for bias adjustment programs with time-to-event outcomes
#' 
#' To discount for overoptimistic results in phase II when calculating the optimal sample size in phase III, 
#' it is necessary to use the following functions, which each describe a specific case:
#' - `EPsProg_L()`: calculates the expected probability of a successful for an additive adjustment factor (i.e adjust the lower bound of the one-sided confidence interval), 
#' however the go-decision is not affected by the bias adjustment
#' - `EPsProg_L2()`: calculates the expected probability of a successful for an additive adjustment factor (i.e adjust the lower bound of the one-sided confidence interval)
#' when the go-decision is also affected by the bias adjustment
#' - `EPsProg_R()`: calculates the expected probability of a successful for a multiplicative adjustment factor (i.e. use estimate with a retention factor), 
#' however the go-decision is not affected by the bias adjustment
#' - `EPsProg_R2()`: calculates the expected probability of a successful for a multiplicative adjustment factor (i.e. use estimate with a retention factor)
#' when the go-decision is also affected by the bias adjustment 
#' @param HRgo threshold value for the go/no-go decision rule
#' @param d2 total events for phase II; must be even number
#' @param Adj adjustment parameter
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param step1 lower boundary for effect size
#' @param step2 upper boundary for effect size
#' @param w weight for mixture prior distribution
#' @param hr1 first assumed true treatment effect on HR scale for prior distribution
#' @param hr2 second assumed true treatment effect on HR scale for prior distribution
#' @param id1 amount of information for `hr1` in terms of number of events
#' @param id2 amount of information for `hr2` in terms of number of events
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `hr1` is used as fixed effect
#' @return  The output of the the functions `EPsProg_L()`, `EPsProg_L2()`, `EPsProg_R()` and `EPsProg_R2()` is the expected probability of a successful program.
#' @examples res <- EPsProg_L(HRgo = 0.8, d2 = 50, Adj = 0.4, 
#'                            alpha = 0.025, beta = 0.1, 
#'                            step1 = 1, step2 = 0.95, 
#'                            w = 0.3, hr1 = 0.69, hr2 = 0.81,
#'                            id1 = 280, id2 = 420, fixed = FALSE)
#'           res <- EPsProg_L2(HRgo = 0.8, d2 = 50, Adj = 0.4, 
#'                            alpha = 0.025, beta = 0.1, 
#'                            step1 = 1, step2 = 0.95, 
#'                            w = 0.3, hr1 = 0.69, hr2 = 0.81,
#'                            id1 = 280, id2 = 420, fixed = FALSE)
#'           res <- EPsProg_R(HRgo = 0.8, d2 = 50, Adj = 0.9, 
#'                            alpha = 0.025, beta = 0.1, 
#'                            step1 = 1, step2 = 0.95, 
#'                            w = 0.3, hr1 = 0.69, hr2 = 0.81,
#'                            id1 = 280, id2 = 420, fixed = FALSE)
#'           res <- EPsProg_R2(HRgo = 0.8, d2 = 50, Adj = 0.9, 
#'                            alpha = 0.025, beta = 0.1, 
#'                            step1 = 1, step2 = 0.95, 
#'                            w = 0.3, hr1 = 0.69, hr2 = 0.81,
#'                            id1 = 280, id2 = 420, fixed = FALSE)
#' @name EPsProg_bias                               
#' @export
#' @editor Johannes Cepicka
#' @editDate 2022-11-04
EPsProg_L<-function(HRgo, d2, Adj, alpha, beta, step1, step2, w, hr1, hr2, id1, id2, fixed){
  
  c=(qnorm(1-alpha)+qnorm(1-beta))^2
  
  if(fixed){
    
    theta <- -log(hr1)
    
    return(integrate(function(y){ 
      sapply(y,function(y){
        ( pnorm(qnorm(1-alpha)-log(step2)/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
                mean=theta/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
                sd=1)-
            pnorm(qnorm(1-alpha)-log(step1)/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
                  mean=theta/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
                  sd=1) )*
          dnorm(y,
                mean=theta,
                sd=sqrt(4/d2))
      })
    }, -log(HRgo),Inf)$value)
    
  }else{
    return(
      integrate(function(x){
        sapply(x,function(x){
          integrate(function(y){ 

              ( pnorm(qnorm(1-alpha)-log(step2)/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
                      mean=x/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
                      sd=1)-
                  pnorm(qnorm(1-alpha)-log(step1)/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
                        mean=x/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
                        sd=1) )*
                dnorm(y,
                      mean=x,
                      sd=sqrt(4/d2))*
                prior_tte(x,w,hr1,hr2,id1,id2) 

          }, -log(HRgo),Inf)$value
        })
      }, -Inf, Inf)$value
    ) 
  }
  
}

#' Utility function for bias adjustment programs with time-to-event outcomes.
#' 
#' The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters and the expected probability of a successful program. 
#' The utility is in a further step maximized by the `optimal_bias()` function.
#' @param HRgo threshold value for the go/no-go decision rule
#' @param d2 total events for phase II; must be even number
#' @param Adj adjustment parameter
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
#' @param steps1 lower boundary for effect size category `"small"` in RR scale, default: 1
#' @param stepm1 lower boundary for effect size category `"medium"` in RR scale = upper boundary for effect size category "small" in RR scale, default: 0.95
#' @param stepl1 lower boundary for effect size category `"large"` in RR scale = upper boundary for effect size category "medium" in RR scale, default: 0.85
#' @param b1 expected gain for effect size category `"small"`
#' @param b2 expected gain for effect size category `"medium"`
#' @param b3 expected gain for effect size category `"large"`
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `hr1` is used as fixed effect
#' @return The output of the the functions `utility_L()`, `utility_L2()`, `utility_R()` and `utility_R2()` is the expected utility of the program.
#' @examples res <- utility_L(d2 = 50, HRgo = 0.8, Adj = 0.4, w = 0.3, 
#'                                  hr1 =  0.69, hr2 = 0.81, 
#'                                  id1 = 280, id2 = 420, xi2 = 0.7, xi3 = 0.7,
#'                                  alpha = 0.025, beta = 0.1,
#'                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                  K = Inf, N = Inf, S = -Inf,
#'                                  steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
#'                                  b1 = 1000, b2 = 2000, b3 = 3000, 
#'                                  fixed = TRUE)
#'          res <- utility_L2(d2 = 50, HRgo = 0.8, Adj = 0.4, w = 0.3, 
#'                                  hr1 =  0.69, hr2 = 0.81, 
#'                                  id1 = 280, id2 = 420, xi2 = 0.7, xi3 = 0.7,
#'                                  alpha = 0.025, beta = 0.1,
#'                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                  K = Inf, N = Inf, S = -Inf,
#'                                  steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
#'                                  b1 = 1000, b2 = 2000, b3 = 3000, 
#'                                  fixed = TRUE)
#'          res <- utility_R(d2 = 50, HRgo = 0.8, Adj = 0.9, w = 0.3, 
#'                                  hr1 =  0.69, hr2 = 0.81, 
#'                                  id1 = 280, id2 = 420, xi2 = 0.7, xi3 = 0.7,
#'                                  alpha = 0.025, beta = 0.1,
#'                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                  K = Inf, N = Inf, S = -Inf,
#'                                  steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
#'                                  b1 = 1000, b2 = 2000, b3 = 3000, 
#'                                  fixed = TRUE)
#'          res <- utility_R2(d2 = 50, HRgo = 0.8, Adj = 0.9, w = 0.3, 
#'                                  hr1 =  0.69, hr2 = 0.81, 
#'                                  id1 = 280, id2 = 420, xi2 = 0.7, xi3 = 0.7,
#'                                  alpha = 0.025, beta = 0.1,
#'                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                  K = Inf, N = Inf, S = -Inf,
#'                                  steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
#'                                  b1 = 1000, b2 = 2000, b3 = 3000, 
#'                                  fixed = TRUE)
#' @name utility_bias                                 
#' @export
#' @editor Johannes Cepicka
#' @editDate 2022-11-04
utility_L <-  function(d2, HRgo, Adj, w, hr1, hr2, id1, id2,
                         alpha, beta, xi2, xi3,
                         c2, c3, c02, c03, 
                         K, N, S,
                         steps1, stepm1, stepl1,
                         b1, b2, b3,
                         fixed){
  
  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 0

  d3    <-  Ed3_L(HRgo = HRgo, d2 = d2, Adj = Adj,
                  alpha = alpha, beta = beta, w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                  fixed = fixed)

  if(is.na(d3)){
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))  
  }else{
    
    # sample size is rounded up to next even natural number
    n2  <- ceiling(d2 * (1/xi2))
    if(round(n2/2) != n2 / 2) {n2 <- n2 + 1}
    
    n3  <- ceiling(d3 * (1/xi3))
    if(round(n3/2) != n3 / 2) {n3 <- n3 + 1}
    
    # expected number of events is rounded to natural number for presentation
    d3  <- ceiling(d3)
    
    if(n2+n3>N){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
    }else{
      
      pg    <-  Epgo_tte(HRgo = HRgo, d2 = d2, w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                         fixed = fixed)
      
      K2    <-  c02 + c2 * n2         # cost phase II
      K3    <-  c03 * pg + c3 * n3    # cost phase III
      
      if(K2+K3>K){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        
      }else{
        # probability of a successful program; small, medium, large effect size
        prob1 <-  EPsProg_L(HRgo = HRgo, d2 = d2, Adj = Adj, alpha = alpha, beta = beta,
                            step1 = steps1, step2 =  steps2,
                            w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                            fixed = fixed)
        prob2 <-  EPsProg_L(HRgo = HRgo, d2 = d2, Adj = Adj, alpha = alpha, beta = beta,
                            step1 =  stepm1, step2 =  stepm2,
                            w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2, 
                            fixed = fixed)
        prob3 <-  EPsProg_L(HRgo = HRgo, d2 = d2, Adj = Adj, alpha = alpha, beta = beta,
                            step1 =  stepl1, step2 = stepl2,
                            w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2, 
                            fixed = fixed)
        
        SP    <-  prob1 + prob2 + prob3
        
        if(SP<S){
          
          return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
          
        }else{
          
          G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
          
          EU    <-  - K2 - K3 + G
          
          return(c(EU, d3, SP, pg, K2, K3, prob1, prob2, prob3, n2, n3))
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

#' Expected probability to go to phase III for bias adjustment programs with time-to-event outcomes
#' 
#' In the case we do not only want do discount for overoptimistic results in phase II when calculating the sample size in phase III, 
#' but also when deciding whether to go to phase III or not the functions `Epgo_L2` and `Epgo_R2` are necessary.
#' The function `Epgo_L2` uses an additive adjustment parameter (i.e adjust the lower bound of the one-sided confidence interval),
#' the function `Epgo_R2` uses a multiplicative adjustment parameter (i.e. use estimate with a retention factor)
#' @param HRgo threshold value for the go/no-go decision rule
#' @param d2 total number of events for phase II; must be even number
#' @param Adj adjustment parameter
#' @param w weight for mixture prior distribution
#' @param hr1 first assumed true treatment effect on HR scale for prior distribution
#' @param hr2 second assumed true treatment effect on HR scale for prior distribution
#' @param id1 amount of information for `hr1` in terms of number of events
#' @param id2 amount of information for `hr2` in terms of number of events
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `hr1` is used as fixed effect
#' @return The output of the the functions `Epgo_L2` and `Epgo_R2` is the expected number of participants in phase III with conservative decision rule and sample size calculation.
#' @examples res <- Epgo_binary_L2(HRgo = 0.8, d2 = 50, Adj = 0.4,  
#'                                 w = 0.3, hr1 = 0.69, hr2 = 0.81, 
#'                                 id1 = 280, id2 = 420, fixed = FALSE)
#'           res <- Epgo_binary_R2(HRgo = 0.8, d2 = 50, Adj = 0.9,  
#'                                 w = 0.3, hr1 = 0.69, hr2 = 0.81, 
#'                                 id1 = 280, id2 = 420, fixed = FALSE)
#' @name Epgo_bias 
#' @export
#' @editor Johannes Cepicka
#' @editDate 2022-11-04
Epgo_L2<-function(HRgo, d2, Adj, w, hr1, hr2, id1, id2, fixed){
  
  if(fixed){
    return(
      pnorm((-log(hr1)+log(HRgo)-qnorm(1-Adj)*sqrt(4/d2))/sqrt(4/d2))  
    )  
  }else(
    return(
      integrate(function(x){
        sapply(x,function(x){ 
          pnorm((x+log(HRgo)-qnorm(1-Adj)*sqrt(4/d2))/sqrt(4/d2))*prior_tte(x,w,hr1,hr2,id1,id2)
        })
      }, -Inf, Inf)$value   
    )
  )

}

# expected number of events for phase III when going to phase III
#' @rdname Ed3_bias 
#' @export
Ed3_L2<-function(HRgo, d2, Adj, alpha, beta, w, hr1, hr2, id1, id2, fixed){
  
  if(fixed){

    int   = try(
        integrate(function(y){
          sapply(y,function(y){
            ( (4*(qnorm(1-alpha)+qnorm(1-beta))^2)/((y-qnorm(1-Adj)*sqrt(4/d2))^2))*
              dnorm(y,
                    mean=-log(hr1),
                    sd=sqrt(4/d2)) 
          })
        }, -log(HRgo)+qnorm(1-Adj)*sqrt(4/d2),Inf),silent=TRUE)
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

              ( (4*(qnorm(1-alpha)+qnorm(1-beta))^2)/((y-qnorm(1-Adj)*sqrt(4/d2))^2))*
                dnorm(y,
                      mean=x,
                      sd=sqrt(4/d2))*
                prior_tte(x,w,hr1,hr2,id1,id2) 

          }, -log(HRgo)+qnorm(1-Adj)*sqrt(4/d2),Inf)$value
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

# expected probability of a successful program
#' @rdname EPsProg_bias 
#' @export
EPsProg_L2<-function(HRgo, d2, Adj, alpha, beta, step1, step2, w, hr1, hr2, id1, id2, fixed){
  
  c=(qnorm(1-alpha)+qnorm(1-beta))^2
  
  if(fixed){
    return(

          integrate(function(y){ 
            
            ( pnorm(qnorm(1-alpha)-log(step2)/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
                    mean=-log(hr1)/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
                    sd=1)-
                pnorm(qnorm(1-alpha)-log(step1)/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
                      mean=-log(hr1)/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
                      sd=1) )*
              dnorm(y,
                    mean=-log(hr1),
                    sd=sqrt(4/d2)) 
            
          }, -log(HRgo)+qnorm(1-Adj)*sqrt(4/d2),Inf)$value
  
    )
  }else{
    return(
      integrate(function(x){
        sapply(x,function(x){
          integrate(function(y){ 

              ( pnorm(qnorm(1-alpha)-log(step2)/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
                      mean=x/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
                      sd=1)-
                  pnorm(qnorm(1-alpha)-log(step1)/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
                        mean=x/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
                        sd=1) )*
                dnorm(y,
                      mean=x,
                      sd=sqrt(4/d2))*
                prior_tte(x,w,hr1,hr2,id1,id2) 

          }, -log(HRgo)+qnorm(1-Adj)*sqrt(4/d2),Inf)$value
        })
      }, -Inf, Inf)$value
    )
  }
  
}


# Utility function
#' @rdname utility_bias
#' @export
utility_L2 <-  function(d2, HRgo, Adj, w, hr1, hr2, id1, id2,
                       alpha, beta, xi2, xi3,
                       c2, c3, c02, c03, 
                       K, N, S,
                       steps1, stepm1, stepl1,
                       b1, b2, b3,
                       fixed){
  
  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 0
  
  d3    <-  Ed3_L2(HRgo = HRgo, d2 = d2, Adj = Adj,
                  alpha = alpha, beta = beta, w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                  fixed = fixed)
  if(is.na(d3)){
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))  
  }else{
    # sample size is rounded up to next even natural number
    n2  <- ceiling(d2 * (1/xi2))
    if(round(n2/2) != n2 / 2) {n2 <- n2 + 1}
    
    n3  <- ceiling(d3 * (1/xi3))
    if(round(n3/2) != n3 / 2) {n3 <- n3 + 1}
    
    # expected number of events is rounded to natural number for presentation
    d3  <- ceiling(d3)
    
    if(n2+n3>N){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
    }else{
      
      pg    <-  Epgo_L2(HRgo = HRgo, d2 = d2, Adj = Adj, w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                         fixed = fixed)
      
      K2    <-  c02 + c2 * n2         # cost phase II
      K3    <-  c03 * pg + c3 * n3    # cost phase III
      
      if(K2+K3>K){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        
      }else{
        # probability of a successful program; small, medium, large effect size
        prob1 <-  EPsProg_L2(HRgo = HRgo, d2 = d2, Adj = Adj, alpha = alpha, beta = beta,
                            step1 = steps1, step2 =  steps2,
                            w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                            fixed = fixed)
        prob2 <-  EPsProg_L2(HRgo = HRgo, d2 = d2, Adj = Adj, alpha = alpha, beta = beta,
                            step1 =  stepm1, step2 =  stepm2,
                            w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2, 
                            fixed = fixed)
        prob3 <-  EPsProg_L2(HRgo = HRgo, d2 = d2, Adj = Adj, alpha = alpha, beta = beta,
                            step1 =  stepl1, step2 = stepl2,
                            w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2, 
                            fixed = fixed)
        
        SP    <-  prob1 + prob2 + prob3
        
        if(SP<S){
          
          return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
          
        }else{
          
          G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
          
          EU    <-  - K2 - K3 + G
          
          return(c(EU, d3, SP, pg, K2, K3, prob1, prob2, prob3, n2, n3))
          #output: expected utility Eud, En3, EsP, Epgo, cost phase II and III
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

# expected number of events for phase III when going to phase III
#' @rdname Ed3_bias 
#' @export
Ed3_R<-function(HRgo, d2, Adj, alpha, beta, w, hr1, hr2, id1, id2, fixed){
  
  if(fixed){
    
    theta <- -log(hr1)
    
    int   = try(integrate(function(y){
      sapply(y,function(y){
        ( (4*(qnorm(1-alpha)+qnorm(1-beta))^2)/((y*Adj)^2))*
          dnorm(y,
                mean=theta,
                sd=sqrt(4/d2)) 
      })
    }, -log(HRgo),Inf),silent=TRUE)
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
                    sd=sqrt(4/d2))*
              prior_tte(x,w,hr1,hr2,id1,id2) 
        }, -log(HRgo),Inf)$value
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

# expected probability of a successful program
#' @rdname EPsProg_bias 
#' @export
EPsProg_R<-function(HRgo, d2, Adj, alpha, beta, step1, step2, w, hr1, hr2, id1, id2, fixed){
  
  c=(qnorm(1-alpha)+qnorm(1-beta))^2
  
  if(fixed){
    
    theta <- -log(hr1)
    
    return(    integrate(function(y){ 
      sapply(y,function(y){
        ( pnorm(qnorm(1-alpha)-log(step2)/(sqrt((y*Adj)^2/c)),
                mean=theta/(sqrt((y*Adj)^2/c)),
                sd=1)-
            pnorm(qnorm(1-alpha)-log(step1)/(sqrt((y*Adj)^2/c)),
                  mean=theta/(sqrt((y*Adj)^2/c)),
                  sd=1) )*
          dnorm(y,
                mean=theta,
                sd=sqrt(4/d2))
      })
    }, -log(HRgo),Inf)$value)
    
  }else{
    return(
      integrate(function(x){
        sapply(x,function(x){
          integrate(function(y){ 

              ( pnorm(qnorm(1-alpha)-log(step2)/(sqrt((y*Adj)^2/c)),
                      mean=x/(sqrt((y*Adj)^2/c)),
                      sd=1)-
                  pnorm(qnorm(1-alpha)-log(step1)/(sqrt((y*Adj)^2/c)),
                        mean=x/(sqrt((y*Adj)^2/c)),
                        sd=1) )*
                dnorm(y,
                      mean=x,
                      sd=sqrt(4/d2))*
                prior_tte(x,w,hr1,hr2,id1,id2) 

          }, -log(HRgo),Inf)$value
        })
      }, -Inf, Inf)$value
    ) 
  }
  
}

# Utility function
#' @rdname utility_bias
#' @export
utility_R <-  function(d2, HRgo, Adj, w, hr1, hr2, id1, id2,
                       alpha, beta, xi2, xi3,
                       c2, c3, c02, c03, 
                       K, N, S,
                       steps1, stepm1, stepl1,
                       b1, b2, b3,
                       fixed){
  
  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 0
  
  d3    <-  Ed3_R(HRgo = HRgo, d2 = d2, Adj = Adj,
                  alpha = alpha, beta = beta, w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                  fixed = fixed)
  
  if(is.na(d3)){
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))  
  }else{
  
    # sample size is rounded up to next even natural number
    n2  <- ceiling(d2 * (1/xi2))
    if(round(n2/2) != n2 / 2) {n2 <- n2 + 1}
    
    n3  <- ceiling(d3 * (1/xi3))
    if(round(n3/2) != n3 / 2) {n3 <- n3 + 1}
    
    # expected number of events is rounded to natural number for presentation
    d3  <- ceiling(d3)
    
    if(n2+n3>N){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
    }else{
      
      pg    <-  Epgo_tte(HRgo = HRgo, d2 = d2, w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                         fixed = fixed)
      
      K2    <-  c02 + c2 * n2         # cost phase II
      K3    <-  c03 * pg + c3 * n3    # cost phase III
      
      if(K2+K3>K){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        
      }else{
        # probability of a successful program; small, medium, large effect size
        prob1 <-  EPsProg_R(HRgo = HRgo, d2 = d2, Adj = Adj, alpha = alpha, beta = beta,
                            step1 = steps1, step2 =  steps2,
                            w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                            fixed = fixed)
        prob2 <-  EPsProg_R(HRgo = HRgo, d2 = d2, Adj = Adj, alpha = alpha, beta = beta,
                            step1 =  stepm1, step2 =  stepm2,
                            w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2, 
                            fixed = fixed)
        prob3 <-  EPsProg_R(HRgo = HRgo, d2 = d2, Adj = Adj, alpha = alpha, beta = beta,
                            step1 =  stepl1, step2 = stepl2,
                            w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2, 
                            fixed = fixed)
        
        SP    <-  prob1 + prob2 + prob3
        
        if(SP<S){
          
          return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
          
        }else{
          
          G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
          
          EU    <-  - K2 - K3 + G
          
          return(c(EU, d3, SP, pg, K2, K3, prob1, prob2, prob3, n2, n3))
          #output: expected utility Eud, En3, EsP, Epgo, cost phase II and III
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

# expected probability to go to phase III
#' @rdname Epgo_bias
#' @export
Epgo_R2<-function(HRgo, d2, Adj, w, hr1, hr2, id1, id2, fixed){

  if(fixed){
    return(
      pnorm((-log(hr1)+(log(HRgo)/Adj))/sqrt(4/d2))
    )
  }else{
    return(
      integrate(function(x){
        sapply(x,function(x){ 
          pnorm((x+(log(HRgo)/Adj))/sqrt(4/d2))*prior_tte(x,w,hr1,hr2,id1,id2)
        })
      }, -Inf, Inf)$value 
    )
  }
}

# expected number of events for phase III when going to phase III
#' @rdname Ed3_bias 
#' @export
Ed3_R2<-function(HRgo, d2, Adj, alpha, beta, w, hr1, hr2, id1, id2, fixed){
  
  if(fixed){
    
    theta <- -log(hr1)
    
    int   = try(integrate(function(y){
          
          ( (4*(qnorm(1-alpha)+qnorm(1-beta))^2)/((y*Adj)^2))*
            dnorm(y,
                  mean=theta,
                  sd=sqrt(4/d2))
          
        }, -log(HRgo)/Adj,Inf),silent=TRUE)
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
                    sd=sqrt(4/d2))*
              prior_tte(x,w,hr1,hr2,id1,id2) 

        }, -log(HRgo)/Adj,Inf)$value
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

# expected probability of a successful program
#' @rdname EPsProg_bias 
#' @export
EPsProg_R2<-function(HRgo, d2, Adj, alpha, beta, step1, step2, w, hr1, hr2, id1, id2, fixed){
  
  c=(qnorm(1-alpha)+qnorm(1-beta))^2
  
  if(fixed){
      return(
        integrate(function(y){ 
          ( pnorm(qnorm(1-alpha)-log(step2)/(sqrt((y*Adj)^2/c)),
                  mean=-log(hr1)/(sqrt((y*Adj)^2/c)),
                  sd=1)-
              pnorm(qnorm(1-alpha)-log(step1)/(sqrt((y*Adj)^2/c)),
                    mean=-log(hr1)/(sqrt((y*Adj)^2/c)),
                    sd=1) )*
            dnorm(y,
                  mean=-log(hr1),
                  sd=sqrt(4/d2))
        }, -log(HRgo)/Adj,Inf)$value

    )
    
  }else{
    return(
      integrate(function(x){
        sapply(x,function(x){
          integrate(function(y){ 
              ( pnorm(qnorm(1-alpha)-log(step2)/(sqrt((y*Adj)^2/c)),
                      mean=x/(sqrt((y*Adj)^2/c)),
                      sd=1)-
                  pnorm(qnorm(1-alpha)-log(step1)/(sqrt((y*Adj)^2/c)),
                        mean=x/(sqrt((y*Adj)^2/c)),
                        sd=1) )*
                dnorm(y,
                      mean=x,
                      sd=sqrt(4/d2))*
                prior_tte(x,w,hr1,hr2,id1,id2) 
          }, -log(HRgo)/Adj,Inf)$value
        })
      }, -Inf, Inf)$value
    ) 
  }
  
}

# Utility function
#' @rdname utility_bias
#' @export
utility_R2 <-  function(d2, HRgo, Adj, w, hr1, hr2, id1, id2,
                       alpha, beta, xi2, xi3,
                       c2, c3, c02, c03, 
                       K, N, S,
                       steps1, stepm1, stepl1,
                       b1, b2, b3,
                       fixed){
  
  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 0
  
  d3    <-  Ed3_R2(HRgo = HRgo, d2 = d2, Adj = Adj,
                  alpha = alpha, beta = beta, w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                  fixed = fixed)
  
  if(is.na(d3)){
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))  
  }else{
  
    # sample size is rounded up to next even natural number
    n2  <- ceiling(d2 * (1/xi2))
    if(round(n2/2) != n2 / 2) {n2 <- n2 + 1}
    
    n3  <- ceiling(d3 * (1/xi3))
    if(round(n3/2) != n3 / 2) {n3 <- n3 + 1}
    
    # expected number of events is rounded to natural number for presentation
    d3  <- ceiling(d3)
    
    if(n2+n3>N){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
    }else{
      
      pg    <-  Epgo_R2(HRgo = HRgo, d2 = d2, Adj = Adj, w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                         fixed = fixed)
      
      K2    <-  c02 + c2 * n2         # cost phase II
      K3    <-  c03 * pg + c3 * n3    # cost phase III
      
      if(K2+K3>K){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        
      }else{
        # probability of a successful program; small, medium, large effect size
        prob1 <-  EPsProg_R2(HRgo = HRgo, d2 = d2, Adj = Adj, alpha = alpha, beta = beta,
                            step1 = steps1, step2 =  steps2,
                            w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                            fixed = fixed)
        prob2 <-  EPsProg_R2(HRgo = HRgo, d2 = d2, Adj = Adj, alpha = alpha, beta = beta,
                            step1 =  stepm1, step2 =  stepm2,
                            w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2, 
                            fixed = fixed)
        prob3 <-  EPsProg_R2(HRgo = HRgo, d2 = d2, Adj = Adj, alpha = alpha, beta = beta,
                            step1 =  stepl1, step2 = stepl2,
                            w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2, 
                            fixed = fixed)
        
        SP    <-  prob1 + prob2 + prob3
        
        if(SP<S){
          
          return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
          
        }else{
          
          G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
          
          EU    <-  - K2 - K3 + G
          
          return(c(EU, d3, SP, pg, K2, K3, prob1, prob2, prob3, n2, n3))
          #output: expected utility Eud, En3, EsP, Epgo, cost phase II and III
        } 
        
      } 
      
    }
  }
}
