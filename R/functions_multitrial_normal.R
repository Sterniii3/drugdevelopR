########################
# Two phase III trials #
########################
# Case 1: Strategy 1/2; at least one trial significant, the treatment effect of the other one at least showing in the same direction 
# Case 2: Strategy 2/2; both trials significant 

#' Expected probability of a successful program for multitrial programs with normally distributed outcomes
#' 
#' These functions calculate the expected probability of a successful program given the parameters. 
#' Each function represents a specific strategy, e.g. the function `EpsProg3_normal()` calculates the expected probability if three phase III trials are performed. 
#' The parameter case specifies how many of the trials have to be successful, i.e. how many trials show a significantly relevant positive treatment effect.
#' 
#' The following cases can be investigated by the software:
#' - Two phase III trials
#'   - Case 1: Strategy 1/2; at least one trial significant, the treatment effect of the other one at least showing in the same direction 
#'   - Case 2: Strategy 2/2; both trials significant
#' - Three phase III trials 
#'   - Case 2: Strategy 2/3; at least two trials significant, the treatment effect of the other one at least showing in the same direction
#'   - Case 3: Strategy 3/3; all trials significant 
#' - Four phase III trials 
#'   - Case 3: Strategy 3/4; at least three trials significant, the treatment effect of the other one at least showing in the same direction
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
#' @return The output of the the function `EPsProg2_normal()`, `EPsProg3_normal()` and `EPsProg4_normal()` is the expected probability of a successful program when performing several phase III trials (2, 3 or 4 respectively).
#' @examples res <- EPsProg2_normal(kappa = 0.1, n2 = 50, alpha = 0.025, beta = 0.1, w = 0.3,
#'                                  Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#'                                  a = 0.25, b = 0.75, 
#'                                  case = 2, size = "small", fixed = FALSE)
#'           res <- EPsProg3_normal(kappa = 0.1, n2 = 50, alpha = 0.025, beta = 0.1, w = 0.3,
#'                                  Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#'                                  a = 0.25, b = 0.75, 
#'                                  case = 2, size = "small", fixed = TRUE)
#'           res <- EPsProg4_normal(kappa = 0.1, n2 = 50, alpha = 0.025, beta = 0.1, w = 0.3,
#'                                  Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#'                                  a = 0.25, b = 0.75, 
#'                                  case = 3, size = "small", fixed = TRUE)                      
#' @name EPsProg_multitrial_normal                                  
#' @export
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
EPsProg2_normal <-  function(kappa, n2, alpha, beta, w, Delta1, Delta2, in1, in2, a, b, case, size, fixed){
  
  SIGMA <-  diag(2)
  c     <-  (qnorm(1 - alpha) + qnorm(1 - beta))^2
  
  if(fixed){
    
    if(case == 1){
      if(size == "small"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( pmvnorm(lower = c(0, 
                                0), 
                      upper = c(qnorm(1 - alpha) + 0.5/sqrt(y^2/c),
                                qnorm(1 - alpha) + 0.5/sqrt(y^2/c)), 
                      mean = c((Delta1)/sqrt(y^2/c), 
                               (Delta1)/sqrt(y^2/c)), 
                      sigma = SIGMA) - 
                pmvnorm(lower = c(0, 
                                  0), 
                        upper = c(qnorm(1 - alpha), 
                                  qnorm(1 - alpha)), 
                        mean = c((Delta1)/sqrt(y^2/c), 
                                 (Delta1)/sqrt(y^2/c)), 
                        sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = Delta1, 
                    sd = sqrt(4/n2))  
          })
        },  kappa, Inf)$value)  
      }
      if(size == "large"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( pmvnorm(lower = c(0, 
                                0), 
                      upper = c(Inf, 
                                Inf), 
                      mean = c((Delta1)/sqrt(y^2/c), 
                               (Delta1)/sqrt(y^2/c)), 
                      sigma = SIGMA) - 
                pmvnorm(lower = c(0, 
                                  0), 
                        upper = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                  qnorm(1 - alpha) + 0.8/sqrt(y^2/c)), 
                        mean = c((Delta1)/sqrt(y^2/c), 
                                 (Delta1)/sqrt(y^2/c)), 
                        sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = Delta1, 
                    sd = sqrt(4/n2)) 
          })
        },  kappa, Inf)$value)  
      }
      if(size == "all"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( pmvnorm(lower = c(0, 
                                0), 
                      upper = c(Inf, 
                                Inf), 
                      mean = c((Delta1)/sqrt(y^2/c), 
                               (Delta1)/sqrt(y^2/c)), 
                      sigma = SIGMA) - 
                pmvnorm(lower = c(0, 
                                  0), 
                        upper = c(qnorm(1 - alpha), 
                                  qnorm(1 - alpha)), 
                        mean = c((Delta1)/sqrt(y^2/c), 
                                 (Delta1)/sqrt(y^2/c)), 
                        sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = Delta1, 
                    sd = sqrt(4/n2))
          })
        },  kappa, Inf)$value)     
      }
    }
    if(case == 2){
      if(size == "small"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                qnorm(1 - alpha)), 
                      upper = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                qnorm(1 - alpha) + 0.8/sqrt(y^2/c)), 
                      mean = c((Delta1)/sqrt(y^2/c), 
                               (Delta1)/sqrt(y^2/c)), 
                      sigma = SIGMA) - 
                pmvnorm(lower = c(qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                  qnorm(1 - alpha) + 0.5/sqrt(y^2/c)), 
                        upper = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                  qnorm(1 - alpha) + 0.8/sqrt(y^2/c)), 
                        mean = c((Delta1)/sqrt(y^2/c), 
                                 (Delta1)/sqrt(y^2/c)), 
                        sigma = SIGMA)) * 
              dnorm(y, 
                    mean = Delta1, 
                    sd = sqrt(4/n2)) 
          })
        },  kappa, Inf)$value) 
      }
      if(size == "large"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( pmvnorm(lower = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                qnorm(1 - alpha) + 0.8/sqrt(y^2/c)), 
                      upper = c(Inf, 
                                Inf), 
                      mean = c((Delta1)/sqrt(y^2/c), 
                               (Delta1)/sqrt(y^2/c)), 
                      sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = Delta1, 
                    sd = sqrt(4/n2)) 
          })
        },  kappa, Inf)$value)    
      }
      if(size == "all"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                qnorm(1 - alpha)), 
                      upper = c(Inf, 
                                Inf), 
                      mean = c((Delta1)/sqrt(y^2/c), 
                               (Delta1)/sqrt(y^2/c)), 
                      sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = Delta1, 
                    sd = sqrt(4/n2)) 
          })
        }, kappa, Inf)$value)    
      }
    }
  }else{
    
    if(case == 1){
      if(size == "small"){
        return(  integrate(function(x){
          sapply(x, function(x){
            integrate(function(y){
              sapply(y, function(y){
                ( pmvnorm(lower = c(0, 
                                    0), 
                          upper = c(qnorm(1 - alpha) + 0.5/sqrt(y^2/c),
                                    qnorm(1 - alpha) + 0.5/sqrt(y^2/c)), 
                          mean = c(x/sqrt(y^2/c), 
                                   x/sqrt(y^2/c)), 
                          sigma = SIGMA) - 
                    pmvnorm(lower = c(0, 
                                      0), 
                            upper = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha)), 
                            mean = c(x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c)), 
                            sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(4/n2)) * 
                  prior_normal(x, w, Delta1, Delta2, in1, in2, a, b)  
              })
            },  kappa, Inf)$value  
          })
        },  - Inf, Inf)$value)
      }
      if(size == "large"){
        return(  integrate(function(x){
          sapply(x, function(x){
            integrate(function(y){
              sapply(y, function(y){
                ( pmvnorm(lower = c(0, 
                                    0), 
                          upper = c(Inf, 
                                    Inf), 
                          mean = c(x/sqrt(y^2/c), 
                                   x/sqrt(y^2/c)), 
                          sigma = SIGMA) - 
                    pmvnorm(lower = c(0, 
                                      0), 
                            upper = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                      qnorm(1 - alpha) + 0.8/sqrt(y^2/c)), 
                            mean = c(x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c)), 
                            sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(4/n2)) * 
                  prior_normal(x, w, Delta1, Delta2, in1, in2, a, b)  
              })
            },  kappa, Inf)$value  
          })
        },  - Inf, Inf)$value)  
      }
      if(size == "all"){
        return(  integrate(function(x){
          sapply(x, function(x){
            integrate(function(y){
              sapply(y, function(y){
                ( pmvnorm(lower = c(0, 
                                    0), 
                          upper = c(Inf, 
                                    Inf), 
                          mean = c(x/sqrt(y^2/c), 
                                   x/sqrt(y^2/c)), 
                          sigma = SIGMA) - 
                    pmvnorm(lower = c(0, 
                                      0), 
                            upper = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha)), 
                            mean = c(x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c)), 
                            sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(4/n2)) * 
                  prior_normal(x, w, Delta1, Delta2, in1, in2, a, b)  
              })
            },  kappa, Inf)$value  
          })
        },  - Inf, Inf)$value)     
      }
    }
    if(case == 2){
      if(size == "small"){
        return(integrate(function(x){
          sapply(x, function(x){
            integrate(function(y){
              sapply(y, function(y){
                ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                    qnorm(1 - alpha)), 
                          upper = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                    qnorm(1 - alpha) + 0.8/sqrt(y^2/c)), 
                          mean = c(x/sqrt(y^2/c), 
                                   x/sqrt(y^2/c)), 
                          sigma = SIGMA) - 
                    pmvnorm(lower = c(qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                      qnorm(1 - alpha) + 0.5/sqrt(y^2/c)), 
                            upper = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                      qnorm(1 - alpha) + 0.8/sqrt(y^2/c)), 
                            mean = c(x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c)), 
                            sigma = SIGMA)) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(4/n2)) * 
                  prior_normal(x, w, Delta1, Delta2, in1, in2, a, b) 
              })
            },  kappa, Inf)$value
          })
        },  - Inf, Inf)$value) 
      }
      if(size == "large"){
        return(integrate(function(x){
          sapply(x, function(x){
            integrate(function(y){
              sapply(y, function(y){
                ( pmvnorm(lower = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                    qnorm(1 - alpha) + 0.8/sqrt(y^2/c)), 
                          upper = c(Inf, 
                                    Inf), 
                          mean = c(x/sqrt(y^2/c), 
                                   x/sqrt(y^2/c)), 
                          sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(4/n2)) * 
                  prior_normal(x, w, Delta1, Delta2, in1, in2, a, b) 
              })
            },  kappa, Inf)$value
          })
        },  - Inf, Inf)$value)    
      }
      if(size == "all"){
        return(integrate(function(x){
          sapply(x, function(x){
            integrate(function(y){
              sapply(y, function(y){
                ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                    qnorm(1 - alpha)), 
                          upper = c(Inf, 
                                    Inf), 
                          mean = c(x/sqrt(y^2/c), 
                                   x/sqrt(y^2/c)), 
                          sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(4/n2)) * 
                  prior_normal(x, w, Delta1, Delta2, in1, in2, a, b) 
              })
            },  kappa, Inf)$value
          })
        },  - Inf, Inf)$value)    
      }
    }
    
  }
  
  
  
}

#' Utility function for multitrial programs with normally distributed outcomes
#' 
#' The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters and the expected probability of a successful program. 
#' The utility is in a further step maximized by the `optimal_multitrial_normal()` function.
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
#' @return The output of the the functions utility2_normal(), utility3_normal() and utility4_normal() is the expected utility of the program when 2, 3 or 4 phase III trials are performed.
#' @examples res <- utility2_normal(kappa = 0.1, n2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
#'                                  Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#'                                  a = 0.25, b = 0.75, 
#'                                  c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
#'                                  K = Inf, N = Inf, S = -Inf, 
#'                                  b1 = 3000, b2 = 8000, b3 = 10000, 
#'                                  case = 2, fixed = TRUE)
#'  #         res <- utility3_normal(kappa = 0.1, n2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
#'  #                                Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#'  #                                a = 0.25, b = 0.75, 
#'  #                                c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
#'  #                                K = Inf, N = Inf, S = -Inf,
#'  #                                b1 = 3000, b2 = 8000, b3 = 10000, 
#'  #                                case = 2, fixed = TRUE)                        
#'  #         res <- utility4_normal(kappa = 0.1, n2 = 50,  alpha = 0.025, beta = 0.1, w = 0.3,
#'  #                                Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#'  #                                a = 0.25, b = 0.75, 
#'  #                                c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
#'  #                                K = Inf, N = Inf, S = -Inf, 
#'  #                                b1 = 3000, b2 = 8000, b3 = 10000, 
#'  #                                case = 3, fixed = TRUE)
#' @name utility_multitrial_normal                           
#' @export
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
utility2_normal <-  function(n2, kappa, w, Delta1, Delta2, in1, in2, a, b,
                      alpha, beta, 
                      c2, c3, c02, c03, 
                      K, N, S,
                      b1, b2, b3,
                      case, fixed){ 
  
  
  n3  <-  En3_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                     w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                     fixed = fixed)
  
  n3  <- ceiling(n3)
  
  if(round(n3/2) != n3 / 2) {n3 = n3 + 1}
  
  if(n2+ 2*n3>N){
    
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
    
  }else{
    
    pg    <-  Epgo_normal(kappa = kappa, n2 = n2,
                          w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                          fixed = fixed)
    
    K2    <-  c02 + c2 * n2  #cost phase II
    K3    <-  c03 * pg + c3 * n3  #cost phase III
    
    if(K2+2*K3>K){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
    }else{
      # probability of a successful program; small, medium, large effect size
      prob1 <-  EPsProg2_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                         w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                         case = case, size = "small", fixed = fixed)
      prob3 <-  EPsProg2_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                         w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                         case = case, size = "large", fixed = fixed)
      prob2 <-  EPsProg2_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                         w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                         case = case, size = "all", fixed = fixed) - prob1 - prob3
      
      SP    <-  prob1 + prob2 + prob3
      
        
        if(SP<S){
          
          return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
          
        }else{
          
          G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
          
          EU    <-  - K2 - 2*K3 + G
          
          return(c(EU, 2*n3, SP, pg, K2, 2*K3, prob1, prob2, prob3))
          #output: expected utility Eud, En3, EsP, Epgo, cost phase II and III
        }
      }
    
  }
  
}



##########################
# Three phase III trials #
##########################
# Case 2: Strategy 2/3; at least two trials significant, the treatment effect 
# of the other one at least showing in the same direction
# Case 3: Strategy 3/3; all trials significant

#' @rdname EPsProg_multitrial_normal 
#' @export
EPsProg3_normal <-  function(kappa, n2, alpha, beta, w, Delta1, Delta2, in1, in2, a, b, case, size, fixed){
  
  SIGMA <-  diag(3)
  c     <-  (qnorm(1 - alpha) + qnorm(1 - beta))^2
  
  if(fixed){
    
    
    if(case == 2){
      if(size == "small"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( 3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    0), 
                          upper = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                    qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                    qnorm(1 - alpha) + 0.5/sqrt(y^2/c)), 
                          mean = c(Delta1/sqrt(y^2/c), 
                                   Delta1/sqrt(y^2/c), 
                                   Delta1/sqrt(y^2/c)), 
                          sigma = SIGMA) - 
                2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha), 
                                      qnorm(1 - alpha)), 
                            upper = c(qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                      qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                      qnorm(1 - alpha) + 0.5/sqrt(y^2/c)), 
                            mean = c(Delta1/sqrt(y^2/c), 
                                     Delta1/sqrt(y^2/c), 
                                     Delta1/sqrt(y^2/c)), 
                            sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = Delta1, 
                    sd = sqrt(4/n2)) 
          })
        },  kappa, Inf)$value)
      }
      if(size == "large"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( 3 * pmvnorm(lower = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                    qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                    0), 
                          upper = c(Inf, 
                                    Inf, 
                                    Inf), 
                          mean = c(Delta1/sqrt(y^2/c), 
                                   Delta1/sqrt(y^2/c), 
                                   Delta1/sqrt(y^2/c)), 
                          sigma = SIGMA) - 
                2 * pmvnorm(lower = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                      qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                      qnorm(1 - alpha) + 0.8/sqrt(y^2/c)), 
                            upper = c(Inf, 
                                      Inf, 
                                      Inf), 
                            mean = c(Delta1/sqrt(y^2/c), 
                                     Delta1/sqrt(y^2/c), 
                                     Delta1/sqrt(y^2/c)), 
                            sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = Delta1, 
                    sd = sqrt(4/n2)) 
          })
        },  kappa, Inf)$value)
      }
      if(size == "all"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( 3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    0), 
                          upper = c(Inf, 
                                    Inf, 
                                    Inf), 
                          mean = c(Delta1/sqrt(y^2/c), 
                                   Delta1/sqrt(y^2/c), 
                                   Delta1/sqrt(y^2/c)), 
                          sigma = SIGMA) - 
                2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha), 
                                      qnorm(1 - alpha)), 
                            upper = c(Inf, 
                                      Inf, 
                                      Inf), 
                            mean = c(Delta1/sqrt(y^2/c), 
                                     Delta1/sqrt(y^2/c), 
                                     Delta1/sqrt(y^2/c)), 
                            sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = Delta1, 
                    sd = sqrt(4/n2))  
          })
        },  kappa, Inf)$value)
      }
    }
    if(case == 3){
      if(size == "small"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( 3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    qnorm(1 - alpha)), 
                          upper = c(qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                    qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                    qnorm(1 - alpha) + 0.5/sqrt(y^2/c)), 
                          mean = c(Delta1/sqrt(y^2/c), 
                                   Delta1/sqrt(y^2/c), 
                                   Delta1/sqrt(y^2/c)), 
                          sigma = SIGMA) - 
                2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha), 
                                      qnorm(1 - alpha)), 
                            upper = c(qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                      qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                      qnorm(1 - alpha) + 0.5/sqrt(y^2/c)), 
                            mean = c(Delta1/sqrt(y^2/c), 
                                     Delta1/sqrt(y^2/c), 
                                     Delta1/sqrt(y^2/c)), 
                            sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = Delta1, 
                    sd = sqrt(4/n2))
          })
        },  kappa, Inf)$value)
      }
      if(size == "large"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( 3 * pmvnorm(lower = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                    qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                    qnorm(1 - alpha)), 
                          upper = c(Inf, 
                                    Inf, 
                                    Inf), 
                          mean = c(Delta1/sqrt(y^2/c), 
                                   Delta1/sqrt(y^2/c), 
                                   Delta1/sqrt(y^2/c)), 
                          sigma = SIGMA) - 
                2 * pmvnorm(lower = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                      qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                      qnorm(1 - alpha) + 0.8/sqrt(y^2/c)), 
                            upper = c(Inf, 
                                      Inf, 
                                      Inf), 
                            mean = c(Delta1/sqrt(y^2/c), 
                                     Delta1/sqrt(y^2/c), 
                                     Delta1/sqrt(y^2/c)), 
                            sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = Delta1, 
                    sd = sqrt(4/n2)) 
          })
        },  kappa, Inf)$value)
      }
      if(size == "all"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                qnorm(1 - alpha), 
                                qnorm(1 - alpha)), 
                      upper = c(Inf, 
                                Inf, 
                                Inf), 
                      mean = c(Delta1/sqrt(y^2/c), 
                               Delta1/sqrt(y^2/c), 
                               Delta1/sqrt(y^2/c)), 
                      sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = Delta1, 
                    sd = sqrt(4/n2)) 
          })
        },  kappa, Inf)$value)
      }
    }
  }else{
    
    if(case == 2){
      if(size == "small"){
        return(integrate(function(x){
          sapply(x, function(x){
            integrate(function(y){
              sapply(y, function(y){
                ( 3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                        qnorm(1 - alpha), 
                                        0), 
                              upper = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                        qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                        qnorm(1 - alpha) + 0.5/sqrt(y^2/c)), 
                              mean = c(x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c)), 
                              sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                          qnorm(1 - alpha), 
                                          qnorm(1 - alpha)), 
                                upper = c(qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                          qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                          qnorm(1 - alpha) + 0.5/sqrt(y^2/c)), 
                                mean = c(x/sqrt(y^2/c), 
                                         x/sqrt(y^2/c), 
                                         x/sqrt(y^2/c)), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(4/n2)) * 
                  prior_normal(x, w, Delta1, Delta2, in1, in2, a, b)
              })
            },  kappa, Inf)$value   
          })
        },  0, Inf)$value)
      }
      if(size == "large"){
        return(integrate(function(x){
          sapply(x, function(x){
            integrate(function(y){
              sapply(y, function(y){
                ( 3 * pmvnorm(lower = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                        qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                        0), 
                              upper = c(Inf, 
                                        Inf, 
                                        Inf), 
                              mean = c(x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c)), 
                              sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                          qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                          qnorm(1 - alpha) + 0.8/sqrt(y^2/c)), 
                                upper = c(Inf, 
                                          Inf, 
                                          Inf), 
                                mean = c(x/sqrt(y^2/c), 
                                         x/sqrt(y^2/c), 
                                         x/sqrt(y^2/c)), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(4/n2)) * 
                  prior_normal(x, w, Delta1, Delta2, in1, in2, a, b) 
              })
            },  kappa, Inf)$value
          })
        },  0, Inf)$value)
      }
      if(size == "all"){
        return(integrate(function(x){
          sapply(x, function(x){
            integrate(function(y){
              sapply(y, function(y){
                ( 3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                        qnorm(1 - alpha), 
                                        0), 
                              upper = c(Inf, 
                                        Inf, 
                                        Inf), 
                              mean = c(x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c)), 
                              sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                          qnorm(1 - alpha), 
                                          qnorm(1 - alpha)), 
                                upper = c(Inf, 
                                          Inf, 
                                          Inf), 
                                mean = c(x/sqrt(y^2/c), 
                                         x/sqrt(y^2/c), 
                                         x/sqrt(y^2/c)), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(4/n2)) * 
                  prior_normal(x, w, Delta1, Delta2, in1, in2, a, b) 
              })
            },  kappa, Inf)$value
          })
        },  0, Inf)$value)
      }
    }
    if(case == 3){
      if(size == "small"){
        return(integrate(function(x){
          sapply(x, function(x){
            integrate(function(y){
              sapply(y, function(y){
                ( 3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                        qnorm(1 - alpha), 
                                        qnorm(1 - alpha)), 
                              upper = c(qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                        qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                        qnorm(1 - alpha) + 0.8/sqrt(y^2/c)), 
                              mean = c(x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c)), 
                              sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                          qnorm(1 - alpha), 
                                          qnorm(1 - alpha)), 
                                upper = c(qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                          qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                          qnorm(1 - alpha) + 0.5/sqrt(y^2/c)), 
                                mean = c(x/sqrt(y^2/c), 
                                         x/sqrt(y^2/c), 
                                         x/sqrt(y^2/c)), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(4/n2)) * 
                  prior_normal(x, w, Delta1, Delta2, in1, in2, a, b) 
              })
            },  kappa, Inf)$value
          })
        },  0, Inf)$value)
      }
      if(size == "large"){
        return(integrate(function(x){
          sapply(x, function(x){
            integrate(function(y){
              sapply(y, function(y){
                ( 3 * pmvnorm(lower = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                        qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                        qnorm(1 - alpha)), 
                              upper = c(Inf, 
                                        Inf, 
                                        Inf), 
                              mean = c(x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c)), 
                              sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                          qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                          qnorm(1 - alpha) + 0.8/sqrt(y^2/c)), 
                                upper = c(Inf, 
                                          Inf, 
                                          Inf), 
                                mean = c(x/sqrt(y^2/c), 
                                         x/sqrt(y^2/c), 
                                         x/sqrt(y^2/c)), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(4/n2)) * 
                  prior_normal(x, w, Delta1, Delta2, in1, in2, a, b) 
              })
            },  kappa, Inf)$value
          })
        },  0, Inf)$value)
      }
      if(size == "all"){
        return(integrate(function(x){
          sapply(x, function(x){
            integrate(function(y){
              sapply(y, function(y){
                ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    qnorm(1 - alpha)), 
                          upper = c(Inf, 
                                    Inf, 
                                    Inf), 
                          mean = c(x/sqrt(y^2/c), 
                                   x/sqrt(y^2/c), 
                                   x/sqrt(y^2/c)), 
                          sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(4/n2)) * 
                  prior_normal(x, w, Delta1, Delta2, in1, in2, a, b) 
              })
            },  kappa, Inf)$value
          })
        },  0, Inf)$value)
      }
    }  
    
  }
  
}

#' @rdname utility_multitrial_normal 
#' @export
utility3_normal <-  function(n2, kappa, w, Delta1, Delta2, in1, in2, a, b,
                             alpha, beta, 
                             c2, c3, c02, c03, 
                             K, N, S,
                             b1, b2, b3,
                             case, fixed){ 
  
  n3  <-  En3_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                     w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                     fixed = fixed)
  
  n3  <- ceiling(n3)
  
  if(round(n3/2) != n3 / 2) {n3 = n3 + 1}
  
  if(n2+ 3*n3>N){
    
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
    
  }else{
    
    pg    <-  Epgo_normal(kappa = kappa, n2 = n2,
                          w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                          fixed = fixed)
    
    K2    <-  c02 + c2 * n2  #cost phase II
    K3    <-  c03 * pg + c3 * n3  #cost phase III
    
    if(K2+3*K3>K){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
    }else{
      # probability of a successful program; small, medium, large effect size
      prob1 <-  EPsProg3_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                                w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                                case = case, size = "small", fixed = fixed)
      prob3 <-  EPsProg3_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                                w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                                case = case, size = "large", fixed = fixed)
      prob2 <-  EPsProg3_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                                w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                                case = case, size = "all", fixed = fixed) - prob1 - prob3
      
      SP    <-  prob1 + prob2 + prob3
      
        
        if(SP<S){
          
          return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
          
        }else{
          
          G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
          
          EU    <-  - K2 - 3*K3 + G
          
          return(c(EU, 3*n3, SP, pg, K2, 3*K3, prob1, prob2, prob3))
          #output: expected utility Eud, En3, EsP, Epgo, cost phase II and III
        }
      }
    
  }
  
}





#########################
# Four phase III trials #
#########################
# Case 3: Strategy 3/4; at least three trials significant, the treatment effect 
# of the other one at least showing in the same direction

#' @rdname EPsProg_multitrial_normal 
#' @export
EPsProg4_normal <-  function(kappa, n2, alpha, beta, w, Delta1, Delta2, in1, in2, a, b, case, size,fixed){
  
  SIGMA <-  diag(4)
  c     <-  (qnorm(1 - alpha) + qnorm(1 - beta))^2
  
  if(fixed){
    
    if(size == "small"){
      return(integrate(function(y){
        sapply(y, function(y){
          ( 4 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                  qnorm(1 - alpha), 
                                  qnorm(1 - alpha), 
                                  0), 
                        upper = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                  qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                  qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                  qnorm(1 - alpha) + 0.5/sqrt(y^2/c)), 
                        mean = c(Delta1/sqrt(y^2/c), 
                                 Delta1/sqrt(y^2/c), 
                                 Delta1/sqrt(y^2/c), 
                                 Delta1/sqrt(y^2/c)), 
                        sigma = SIGMA)  - 
              3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    qnorm(1 - alpha)), 
                          upper = c(qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                    qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                    qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                    qnorm(1 - alpha) + 0.5/sqrt(y^2/c)), 
                          mean = c(Delta1/sqrt(y^2/c), 
                                   Delta1/sqrt(y^2/c), 
                                   Delta1/sqrt(y^2/c), 
                                   Delta1/sqrt(y^2/c)), 
                          sigma = SIGMA) ) * 
            dnorm(y, 
                  mean = Delta1, 
                  sd = sqrt(4/n2)) 
        })
      },  kappa, Inf)$value)
    }
    if(size == "large"){
      return(integrate(function(y){
        sapply(y, function(y){
          ( 4 * pmvnorm(lower = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                  qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                  qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                  0), 
                        upper = c(Inf, 
                                  Inf, 
                                  Inf, 
                                  Inf), 
                        mean = c(Delta1/sqrt(y^2/c), 
                                 Delta1/sqrt(y^2/c), 
                                 Delta1/sqrt(y^2/c), 
                                 Delta1/sqrt(y^2/c)), 
                        sigma = SIGMA)  - 
              3 * pmvnorm(lower = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                    qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                    qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                    qnorm(1 - alpha) + 0.8/sqrt(y^2/c)), 
                          upper = c(Inf, 
                                    Inf, 
                                    Inf, 
                                    Inf), 
                          mean = c(Delta1/sqrt(y^2/c), 
                                   Delta1/sqrt(y^2/c), 
                                   Delta1/sqrt(y^2/c), 
                                   Delta1/sqrt(y^2/c)), 
                          sigma = SIGMA) ) * 
            dnorm(y, 
                  mean = Delta1, 
                  sd = sqrt(4/n2)) 
        })
      },  kappa, Inf)$value)
    }
    if(size == "all"){
      return(integrate(function(y){
        sapply(y, function(y){
          ( 4 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                  qnorm(1 - alpha), 
                                  qnorm(1 - alpha), 
                                  0), 
                        upper = c(Inf, 
                                  Inf, 
                                  Inf, 
                                  Inf), 
                        mean = c(Delta1/sqrt(y^2/c), 
                                 Delta1/sqrt(y^2/c), 
                                 Delta1/sqrt(y^2/c), 
                                 Delta1/sqrt(y^2/c)), 
                        sigma = SIGMA) - 
              3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    qnorm(1 - alpha)), 
                          upper = c(Inf, 
                                    Inf, 
                                    Inf, 
                                    Inf), 
                          mean = c(Delta1/sqrt(y^2/c), 
                                   Delta1/sqrt(y^2/c), 
                                   Delta1/sqrt(y^2/c), 
                                   Delta1/sqrt(y^2/c)), 
                          sigma = SIGMA) ) * 
            dnorm(y, 
                  mean = Delta1, 
                  sd = sqrt(4/n2)) 
        })
      },  kappa, Inf)$value)
    } 
  }else{
    
    if(size == "small"){
      return(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            sapply(y, function(y){
              ( 4 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha), 
                                      qnorm(1 - alpha), 
                                      0), 
                            upper = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                      qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                      qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                      qnorm(1 - alpha) + 0.5/sqrt(y^2/c)), 
                            mean = c(x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c)), 
                            sigma = SIGMA)  - 
                  3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                        qnorm(1 - alpha), 
                                        qnorm(1 - alpha), 
                                        qnorm(1 - alpha)), 
                              upper = c(qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                        qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                        qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                        qnorm(1 - alpha) + 0.5/sqrt(y^2/c)), 
                              mean = c(x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c)), 
                              sigma = SIGMA) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(4/n2)) * 
                prior_normal(x, w, Delta1, Delta2, in1, in2, a, b)
            })
          },  kappa, Inf)$value   
        })
      },  0, Inf)$value)
    }
    if(size == "large"){
      return(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            sapply(y, function(y){
              ( 4 * pmvnorm(lower = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                      qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                      qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                      0), 
                            upper = c(Inf, 
                                      Inf, 
                                      Inf, 
                                      Inf), 
                            mean = c(x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c)), 
                            sigma = SIGMA)  - 
                  3 * pmvnorm(lower = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                        qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                        qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                        qnorm(1 - alpha) + 0.8/sqrt(y^2/c)), 
                              upper = c(Inf, 
                                        Inf, 
                                        Inf, 
                                        Inf), 
                              mean = c(x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c)), 
                              sigma = SIGMA) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(4/n2)) * 
                prior_normal(x, w, Delta1, Delta2, in1, in2, a, b) 
            })
          },  kappa, Inf)$value
        })
      },  0, Inf)$value)
    }
    if(size == "all"){
      return(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            sapply(y, function(y){
              ( 4 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha), 
                                      qnorm(1 - alpha), 
                                      0), 
                            upper = c(Inf, 
                                      Inf, 
                                      Inf, 
                                      Inf), 
                            mean = c(x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c)), 
                            sigma = SIGMA) - 
                  3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                        qnorm(1 - alpha), 
                                        qnorm(1 - alpha), 
                                        qnorm(1 - alpha)), 
                              upper = c(Inf, 
                                        Inf, 
                                        Inf, 
                                        Inf), 
                              mean = c(x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c)), 
                              sigma = SIGMA) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(4/n2)) * 
                prior_normal(x, w, Delta1, Delta2, in1, in2, a, b)
            })
          },  kappa, Inf)$value
        })
      },  0, Inf)$value)
    }
    
  }
  
}

#' @rdname utility_multitrial_normal 
#' @export
utility4_normal <-  function(n2, kappa, w, Delta1, Delta2, in1, in2, a, b,
                             alpha, beta, 
                             c2, c3, c02, c03, 
                             K, N, S,
                             b1, b2, b3,
                             case, fixed){ 
  
  n3  <-  En3_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                     w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                     fixed = fixed)
  
  n3  <- ceiling(n3)
  
  if(round(n3/2) != n3 / 2) {n3 = n3 + 1}
  
  if(n2+4*n3>N){
    
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
    
  }else{
    
    pg    <-  Epgo_normal(kappa = kappa, n2 = n2,
                          w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                          fixed = fixed)
    
    K2    <-  c02 + c2 * n2  #cost phase II
    K3    <-  c03 * pg + c3 * n3  #cost phase III
    
    if(K2+4*K3>K){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
    }else{
      # probability of a successful program; small, medium, large effect size
      prob1 <-  EPsProg4_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                                w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                                case = case, size = "small", fixed = fixed)
      prob3 <-  EPsProg4_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                                w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                                case = case, size = "large", fixed = fixed)
      prob2 <-  EPsProg4_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                                w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                                case = case, size = "all", fixed = fixed) - prob1 - prob3
      
      SP    <-  prob1 + prob2 + prob3
      
        
        if(SP<S){
          
          return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
          
        }else{
          
          G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
          
          EU    <-  - K2 - 4*K3 + G
          
          return(c(EU, 4*n3, SP, pg, K2, 4*K3, prob1, prob2, prob3))
          #output: expected utility Eud, En3, EsP, Epgo, cost phase II and III
        }
      }
    
  }
  
}


#################################
# Two or three phase III trials #
#################################
# Case 2: Strategy 2/2( + 1); at least two trials significant (and the 
# treatment effect of the other one at least showing in the same direction)

#' Expected probability to do third phase III trial
#' 
#' In the setting of Case 2: Strategy 2/2( + 1); at least two trials significant (and the 
#' treatment effect of the other one at least showing in the same direction) this function calculates the probability that a third phase III trial is necessary.
#' @param kappa threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param w weight for mixture prior distribution
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param in1 amount of information for `Delta1` in terms of sample size
#' @param in2 amount of information for `Delta2` in terms of sample size
#' @return The output of the the function `Epgo23_normal()` is the probability to to a third phase III trial.
#' @examples res <- Epgo23_normal(kappa = 0.1, n2 = 50, w = 0.3, alpha = 0.025, beta = 0.1, a = 0.25, b=0.75,
#'                                Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600)
#' @export
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
Epgo23_normal <-  function(kappa, n2, alpha, beta, a, b,  w, Delta1, Delta2, in1, in2){
  
  SIGMA <-  diag(2)
  c     <-  (qnorm(1 - alpha) + qnorm(1 - beta))^2
  
  integrate(function(x){                    
    sapply(x, function(x){
      integrate(function(y){
        sapply(y, function(y){
          2 * (pmvnorm(lower = c(qnorm(1 - alpha), 
                                 0), 
                       upper = c(Inf, 
                                 qnorm(1 - alpha)), 
                       mean = c(x/sqrt(y^2/c), 
                                x/sqrt(y^2/c)), 
                       sigma = SIGMA)) * 
            dnorm(y, 
                  mean = x, 
                  sd = sqrt(4/n2)) * 
            prior_normal(x, w, Delta1, Delta2, in1, in2, a, b)
        })
      },  kappa, Inf)$value
    })
  },  -Inf, Inf)$value
} 

#' Expected probability of a successful program deciding between two or three phase III trials for a normally distributed outcome
#'
#' The function `EPsProg23_normal()` calculates the expected probability of a successful program
#' with a normally distributed outcome. This function follows a special decision rule in order to determine
#' whether two or three phase III trials should be conducted. First, two phase III trials are performed. Depending
#' on their success, the decision for a third phase III trial is made:
#' - If both trials are successful, no third phase III trial will be conducted.
#' - If only one of the two trials is successful and the other trial has a treatment effect that points in the same direction,
#' a third phase III trial will be conducted with a sample size of N3 = N3(ymin), which depends on an assumed minimal clinical relevant effect (`ymin`).
#' The third trial then has to be significant at level `alpha`
#' - If only one of the two trials is successful and the treatment effect of the other points in opposite direction or 
#' if none of the two trials are successful, then no third trial is performed and the drug development development program is not successful. 
#' In the utility function, this will lead to a utility of -9999.
#' 
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
#'                                  Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#'                                  a = 0.25, b = 0.75, 
#'                                  case = 2, size = "small", ymin = 0.5)
#' @export
#' @editor Johannes Cepicka, Lukas D. Sauer
#' @editDate 2022-05-06
EPsProg23_normal <-  function(kappa, n2, alpha, beta, w, Delta1, Delta2, in1, in2, a, b, case, size, ymin){
  # Option 2.1: first two phase III trials are successful: no third phase III trial
  # Option 2.2: one of the two first phase III trials successful, the treatment
  #  effect of the other one points in the same direction: 
  #  conduct third phase III trial with N3 = N3(ymin)
  
  SIGMA <-  diag(2)
  SIGMA3<-  diag(3)
  c     <-  (qnorm(1 - alpha) + qnorm(1 - beta))^2
  
  if(case == 2){ # Option 2.1
    if(size == "small"){
      return(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            sapply(y, function(y){
              ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                  qnorm(1 - alpha)), 
                        upper = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                  qnorm(1 - alpha) + 0.8/sqrt(y^2/c)), 
                        mean = c(x/sqrt(y^2/c), 
                                 x/sqrt(y^2/c)), 
                        sigma = SIGMA)  - 
                  pmvnorm(lower = c(qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                    qnorm(1 - alpha) + 0.5/sqrt(y^2/c)), 
                          upper = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                    qnorm(1 - alpha) + 0.8/sqrt(y^2/c)), 
                          mean = c(x/sqrt(y^2/c), 
                                   x/sqrt(y^2/c)), 
                          sigma = SIGMA)) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(4/n2)) * 
                prior_normal(x, w, Delta1, Delta2, in1, in2, a, b) 
            })
          },  kappa, Inf)$value
        })
      },  0, Inf)$value) 
    }
    if(size == "large"){
      return(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            sapply(y, function(y){
              ( pmvnorm(lower = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                  qnorm(1 - alpha) + 0.8/sqrt(y^2/c)), 
                        upper = c(Inf, Inf), 
                        mean = c(x/sqrt(y^2/c), 
                                 x/sqrt(y^2/c)), 
                        sigma = SIGMA) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(4/n2)) * 
                prior_normal(x, w, Delta1, Delta2, in1, in2, a, b) 
            })
          },  kappa, Inf)$value
        })
      },  0, Inf)$value)    
    }
    if(size == "all"){
      return(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            sapply(y, function(y){
              ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                  qnorm(1 - alpha)), 
                        upper = c(Inf, 
                                  Inf), 
                        mean = c(x/sqrt(y^2/c), 
                                 x/sqrt(y^2/c)), 
                        sigma = SIGMA) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(4/n2)) * 
                prior_normal(x, w, Delta1, Delta2, in1, in2, a, b) 
            })
          },  kappa, Inf)$value
        })
      },  0, Inf)$value)    
    }
  }
  if(case == 3){# Option 2.2
    if(size == "small"){
      return(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            sapply(y, function(y){
              2 * ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                      0, 
                                      qnorm(1 - alpha)), 
                            upper = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                      qnorm(1 - alpha), 
                                      qnorm(1 - alpha) + 0.8/sqrt(ymin^2/c)), 
                            mean = c(x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c), 
                                     x/sqrt(ymin^2/c)), 
                            sigma = SIGMA3)  - 
                      pmvnorm(lower = c(qnorm(1 - alpha) + 0.5/sqrt(y^2/c), 
                                        0, 
                                        qnorm(1 - alpha) + 0.5/sqrt(ymin^2/c)), 
                              upper = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                        qnorm(1 - alpha), 
                                        qnorm(1 - alpha) + 0.8/sqrt(ymin^2/c)), 
                              mean = c(x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c), 
                                       x/sqrt(ymin^2/c)), 
                              sigma = SIGMA3)) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(4/n2)) * 
                prior_normal(x, w, Delta1, Delta2, in1, in2, a, b) 
            })
          },  kappa, Inf)$value
        })
      },  0, Inf)$value) 
    }
    if(size == "large"){
      return(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            sapply(y, function(y){
              2 * ( pmvnorm(lower = c(qnorm(1 - alpha) + 0.8/sqrt(y^2/c), 
                                      0, 
                                      qnorm(1 - alpha) + 0.8/sqrt(ymin^2/c)), 
                            upper = c(Inf, 
                                      qnorm(1 - alpha), 
                                      Inf), 
                            mean = c(x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c), 
                                     x/sqrt(ymin^2/c)), 
                            sigma = SIGMA3) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(4/n2)) * 
                prior_normal(x, w, Delta1, Delta2, in1, in2, a, b) 
            })
          },  kappa, Inf)$value
        })
      },  0, Inf)$value)    
    }
    if(size == "all"){
      return(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            sapply(y, function(y){
              2 * ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                      0, 
                                      qnorm(1 - alpha)), 
                            upper = c(Inf, 
                                      qnorm(1 - alpha), 
                                      Inf), 
                            mean = c(x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c), 
                                     x/sqrt(ymin^2/c)), 
                            sigma = SIGMA3) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(4/n2)) * 
                prior_normal(x, w, Delta1, Delta2, in1, in2, a, b) 
            })
          },  kappa, Inf)$value
        })
      },  0, Inf)$value)    
    }
  }
  
}

#' Utility function for multitrial programs deciding between two or three phase III trials for a normally distributed outcome
#'
#' The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters and the expected probability of a successful program. 
#' The utility is in a further step maximized by the `optimal_multitrial_normal()` function.
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
#' @return The output of the the function utility23_normal() is the expected utility of the program depending on whether two or three phase III trials are performed.
#' @examples #res <- utility23_normal(n2 = 50, kappa = 0.2, w = 0.3, alpha = 0.025, beta = 0.1,
#'       #                           Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
#'       #                           a = 0.25, b = 0.75, 
#'      #                            c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,
#'      #                            b1 = 3000, b2 = 8000, b3 = 10000)
#' @export
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
utility23_normal <-  function(n2, kappa, w, Delta1, Delta2, in1, in2, a, b,
                       alpha, beta, 
                       c2, c3, c02, c03, 
                       b1, b2, b3){ 
  
  pg    <-  Epgo_normal(kappa = kappa, n2 = n2, 
                        w = w , Delta1 = Delta1, Delta2 = Delta2 , 
                        in1 = in1, in2 = in2,
                        a = a, b = b, fixed = FALSE)
  
  
  n3  <-  En3_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                     w = w, Delta1 = Delta1, Delta2 = Delta2, 
                     in1 = in1, in2 = in2, 
                     a = a, b = b, fixed = FALSE)
  
  n3  <- ceiling(n3)
  
  if(round(n3/2) != n3 / 2) {n3 = n3 + 1}
  
  
  # probability of a successful program:
  # small, medium and large effect size
  prob1 <- EPsProg23_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                            w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b, 
                            case = 2, size = "small", ymin = ymin)
  prob3 <- EPsProg23_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                            w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b, 
                            case = 2, size = "large", ymin = ymin)
  prob2 <- EPsProg23_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                            w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b, 
                            case = 2, size = "all", ymin = ymin) - prob1 - prob3
  
  # prob to do third phase III trial
  pg3   <-  Epgo23_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta, 
                          a = a, b = b, w = w, Delta1 = Delta1, Delta2 = Delta2,
                          in1, in2) 
  
  # n3 for third phase III trial
  n33   <-  (4 * (qnorm(1 - alpha) + qnorm(1 - beta))^2)/(ymin^2) 
  
  n33  <- ceiling(n33*pg3)
  
  if(round(n33/2) != n33 / 2) {n33 <- n33 + 1}
  
  
  # probability of a successful program: effect sizes, 
  # for program with third phase III trial
  
  # small 
  prob13   <-  EPsProg23_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                                w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b, 
                                case = 3, size = "small", ymin = ymin) 
  # large
  prob33   <-  EPsProg23_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                                w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b, 
                                case = 3, size = "large", ymin = ymin) 
  # medium
  prob23   <-  EPsProg23_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                                w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b, 
                                case = 3, size = "all", ymin = ymin) - prob13 - prob33
  
  K1    <-  c02 + c2 * n2 # cost phase II 
  
  # cost for one of the first two phase III trials in case of go decision
  K2    <-  c03 * pg + c3 * n3
  
  # cost for the third phase III trial in case of third phase III trial
  K3    <-  pg3 * c03 + c3 * n33
  
  G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 +
    b1 * prob13 + b2 * prob23 + b3 * prob33 # gain
  
  EU    <-  - K1 - 2 * K2 - K3 + G
  SP    <-  prob1 + prob2 + prob3 + 
      prob13 + prob23 + prob33
  
  return(
    c(EU, 2*n3, SP, pg, 2*K2, K3, prob1, prob2, prob3, pg3, n33, prob13, prob23, prob33 )
  )
  
}

