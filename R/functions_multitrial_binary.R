# auxiliary functions
t1 <- function(x, p0){((1-p0)/p0) + ((1-x)/x)}
t2 <- function(x, p0){sqrt(2*(1-((p0 + x)/2))/((p0 + x)/2))}
t3 <- function(x, p0){sqrt(((1-p0)/p0) + ((1-x)/x))}


########################
# Two phase III trials #
########################
# Case 1: Strategy 1/2; at least one trial significant, the treatment effect of the other one at least showing in the same direction 
# Case 2: Strategy 2/2; both trials significant 

 #' Expected probability of a successful program for multitrial programs with binary distributed outcomes
 #' 
 #' These functions calculate the expected probability of a successful program given the parameters. 
 #' Each function represents a specific strategy, e.g. the function EpsProg3_binary() calculates the expected probability if three phase III trials are performed. 
 #' The parameter case specifies how many of the trials have to be successful, i.e. how many trials show a significantly relevant positive treatment effect.
 #' 
 #' The following cases can be investigated by the software:
 #' - Two phase III trials
 #'   -  Case 1: Strategy 1/2; at least one trial significant, the treatment effect of the other one at least showing in the same direction 
 #'   - Case 2: Strategy 2/2; both trials significant
 #' - Three phase III trials 
 #'   - Case 2: Strategy 2/3; at least two trials significant, the treatment effect of the other one at least showing in the same direction
 #'   - Case 3: Strategy 3/3; all trials significant 
 #' - Four phase III trials 
 #'   - Case 3: Strategy 3/4; at least three trials significant, the treatment effect of the other one at least showing in the same direction
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
 #' @return The output of the the function EPsProg2_binary(), EPsProg3_binary() and EPsProg4_binary() is the expected probability of a successful program when performing several phase III trials (2, 3 or 4 respectively)
 #' @examples res <- EPsProg2_binary(RRgo = 0.8, n2 = 50,  alpha = 0.025, beta = 0.1, 
 #'                                  p0 = 0.6,  w = 0.3, p11 =  0.3, p12 = 0.5, 
 #'                                  in1 = 300, in2 = 600, case = 2, size = "small",
 #'                                  fixed = FALSE)
 #'           res <- EPsProg3_binary(RRgo = 0.8, n2 = 50,  alpha = 0.025, beta = 0.1, 
 #'                                  p0 = 0.6,  w = 0.3, p11 =  0.3, p12 = 0.5, 
 #'                                  in1 = 300, in2 = 600, case = 2, size = "small",
 #'                                  fixed = FALSE)
 #'           res <- EPsProg4_binary(RRgo = 0.8, n2 = 50,  alpha = 0.025, beta = 0.1, 
 #'                                  p0 = 0.6,  w = 0.3, p11 =  0.3, p12 = 0.5, 
 #'                                  in1 = 300, in2 = 600, case = 3, size = "small",
 #'                                  fixed = FALSE)
 #' @export                                 
 #' @name EPsProg_multitrial_binary  
 #' @editor Johannes Cepicka
 #' @editDate 2022-04-23
 
EPsProg2_binary <-  function(RRgo, n2, alpha, beta, p0, w, p11, p12, in1, in2, case, size, fixed){
  
  SIGMA <-  diag(2)
  const <- function(x, p0){(qnorm(1-alpha)*t2(x, p0) + qnorm(1-beta)*t3(x, p0))^2}
  
  if(fixed){
    
    rho <- -log(p11/p0)
    
    if(case == 1){
      if(size == "small"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( pmvnorm(lower = c(0, 
                                0), 
                      upper = c(qnorm(1 - alpha) - log(0.95)/sqrt(t1(p11, p0)*y^2/const(p11, p0)),
                                qnorm(1 - alpha) - log(0.95)/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                      mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                               rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                      sigma = SIGMA) - 
                pmvnorm(lower = c(0, 
                                  0), 
                        upper = c(qnorm(1 - alpha), 
                                  qnorm(1 - alpha)), 
                        mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                 rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                        sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = rho, 
                    sd = sqrt(((2/n2))*t1(p11, p0))) 
          })
        },  - log(RRgo), Inf)$value)  
      }
      if(size == "large"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( pmvnorm(lower = c(0, 
                                0), 
                      upper = c(Inf, 
                                Inf), 
                      mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                               rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                      sigma = SIGMA) - 
                pmvnorm(lower = c(0, 
                                  0), 
                        upper = c(qnorm(1 - alpha) - 
                                    log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                  qnorm(1 - alpha) - 
                                    log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                        mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                 rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                        sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = rho, 
                    sd = sqrt(((2/n2))*t1(p11, p0))) 
          })
        },  - log(RRgo), Inf)$value)  
      }
      if(size == "all"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( pmvnorm(lower = c(0, 
                                0), 
                      upper = c(Inf, 
                                Inf), 
                      mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                               rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                      sigma = SIGMA) - 
                pmvnorm(lower = c(0, 
                                  0), 
                        upper = c(qnorm(1 - alpha), 
                                  qnorm(1 - alpha)), 
                        mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                 rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                        sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = rho, 
                    sd = sqrt(((2/n2))*t1(p11, p0)))
          })
        },  - log(RRgo), Inf)$value)     
      }
    }
    if(case == 2){
      if(size == "small"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                qnorm(1 - alpha)), 
                      upper = c(qnorm(1 - alpha) - 
                                  log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                qnorm(1 - alpha) - 
                                  log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                      mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                               rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                      sigma = SIGMA) - 
                pmvnorm(lower = c(qnorm(1 - alpha) - 
                                    log(0.95)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                  qnorm(1 - alpha) - 
                                    log(0.95)/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                        upper = c(qnorm(1 - alpha) - 
                                    log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                  qnorm(1 - alpha) - 
                                    log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                        mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                 rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                        sigma = SIGMA)) * 
              dnorm(y, 
                    mean = rho, 
                    sd = sqrt(((2/n2))*t1(p11, p0))) 
          })
        },  - log(RRgo), Inf)$value) 
      }
      if(size == "large"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( pmvnorm(lower = c(qnorm(1 - alpha) - 
                                  log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                qnorm(1 - alpha) - 
                                  log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                      upper = c(Inf, 
                                Inf), 
                      mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                               rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                      sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = rho, 
                    sd = sqrt(((2/n2))*t1(p11, p0))) 
          })
        },  - log(RRgo), Inf)$value)    
      }
      if(size == "all"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                qnorm(1 - alpha)), 
                      upper = c(Inf, 
                                Inf), 
                      mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                               rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                      sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = rho, 
                    sd = sqrt(((2/n2))*t1(p11, p0))) 
          })
        },  - log(RRgo), Inf)$value)    
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
                          upper = c(qnorm(1 - alpha) - log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0)),
                                    qnorm(1 - alpha) - log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                          mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                   x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                          sigma = SIGMA) - 
                    pmvnorm(lower = c(0, 
                                      0), 
                            upper = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha)), 
                            mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                     x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                            sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*t1(x, p0))) * 
                  prior_binary(x, w, p11, p12, in1, in2)  
              })
            },  - log(RRgo), Inf)$value  
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
                          mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                   x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                          sigma = SIGMA) - 
                    pmvnorm(lower = c(0, 
                                      0), 
                            upper = c(qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                            mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                     x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                            sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*t1(x, p0))) * 
                  prior_binary(x, w, p11, p12, in1, in2)  
              })
            },  - log(RRgo), Inf)$value  
          })
        },  0, 1)$value)  
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
                          mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                   x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                          sigma = SIGMA) - 
                    pmvnorm(lower = c(0, 
                                      0), 
                            upper = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha)), 
                            mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                     x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                            sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*t1(x, p0))) * 
                  prior_binary(x, w, p11, p12, in1, in2)  
              })
            },  - log(RRgo), Inf)$value  
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
                          upper = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                          mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                   x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                          sigma = SIGMA) - 
                    pmvnorm(lower = c(qnorm(1 - alpha) - 
                                        log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                      qnorm(1 - alpha) - 
                                        log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                            upper = c(qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                            mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                     x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                            sigma = SIGMA)) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*t1(x, p0))) * 
                  prior_binary(x, w, p11, p12, in1, in2) 
              })
            },  - log(RRgo), Inf)$value
          })
        },  0, 1)$value) 
      }
      if(size == "large"){
        return(integrate(function(x){
          sapply(x, function(x){
            integrate(function(y){
              sapply(y, function(y){
                ( pmvnorm(lower = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                          upper = c(Inf, 
                                    Inf), 
                          mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                   x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                          sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*t1(x, p0))) * 
                  prior_binary(x, w, p11, p12, in1, in2) 
              })
            },  - log(RRgo), Inf)$value
          })
        },  0, 1)$value)    
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
                          mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                   x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                          sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*t1(x, p0))) * 
                  prior_binary(x, w, p11, p12, in1, in2) 
              })
            },  - log(RRgo), Inf)$value
          })
        },  0, 1)$value)    
      }
    }
    
  }
  
  
  
}


 #' Utility function for multitrial programs with binary distributed outcomes
 #' 
 #' The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters and the expected probability of a successful program. 
 #' The utility is in further step maximized by the `optimal_multitrial_binary()` function.
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
 #' @return The output of the the `functions utility2_binary()`, `utility3_binary()` and `utility4_binary()` is the expected utility of the program when 2, 3 or 4 phase III trials are performed.
 #' @examples res <- utility2_binary(n2 = 50, RRgo = 0.8,  w = 0.3, 
 #'                                  p0 = 0.6, p11 =  0.3, p12 = 0.5, 
 #'                                  in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
 #'                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
 #'                                  K = Inf, N = Inf, S = -Inf,
 #'                                  b1 = 1000, b2 = 2000, b3 = 3000, 
 #'                                  case = 2, fixed = TRUE)
 #'           res <- utility3_binary(n2 = 50, RRgo = 0.8,  w = 0.3, 
 #'                                  p0 = 0.6, p11 =  0.3, p12 = 0.5, 
 #'                                  in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
 #'                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
 #'                                  K = Inf, N = Inf, S = -Inf,
 #'                                  b1 = 1000, b2 = 2000, b3 = 3000, 
 #'                                  case = 2, fixed = TRUE)
 #'          res <- utility4_binary(n2 = 50, RRgo = 0.8,  w = 0.3, 
 #'                                  p0 = 0.6, p11 =  0.3, p12 = 0.5, 
 #'                                  in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
 #'                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
 #'                                  K = Inf, N = Inf, S = -Inf,
 #'                                  b1 = 1000, b2 = 2000, b3 = 3000, 
 #'                                  case = 3, fixed = TRUE)
 #' @name utility_multitrial_binary 
 #' @export                                
 #' @editor Johannes Cepicka
 #' @editDate 2022-04-23
utility2_binary <-  function(n2, RRgo, w, p0, p11, p12, in1, in2,
                            alpha, beta, 
                            c2, c3, c02, c03, 
                            K, N, S,
                            b1, b2, b3,
                            case, fixed){
  
  
  n3  <-  En3_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                     p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, fixed = fixed)
  
  n3  <- ceiling(n3)
  
  if(round(n3/2) != n3 / 2) {n3 = n3 + 1}
  
  if(n2+2*n3>N){
    
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
    
  }else{
    
    pg    <-  Epgo_binary(RRgo = RRgo, n2 = n2, p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, fixed = fixed)
    
    K2    <-  c02 + c2 * n2         # cost phase II
    K3    <-  c03 * pg + c3 * n3    # cost phase III
    
    if(K2+2*K3>K){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
    }else{
      # probability of a successful program; small, medium, large effect size
      prob1 <-  EPsProg2_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                               p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                               case = case, size = "small", fixed = fixed)
      prob3 <-  EPsProg2_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                               p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                               case = case, size = "large", fixed = fixed)
      prob2 <-  EPsProg2_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                               p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                               case = case, size = "all", fixed = fixed) - prob1 - prob3
      
      SP    <-  prob1 + prob2 + prob3
      
      if(SP<S){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        
      }else{
        
        G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
        
        EU    <-  - K2 - 2*K3 + G
        
        return(c(EU, 2*n3, SP, pg, K2, 2*K3, prob1, prob2, prob3))
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

#' @rdname EPsProg_multitrial_binary 
#' @export
EPsProg3_binary <-  function(RRgo, n2, alpha, beta, p0, w, p11, p12, in1, in2, case, size, fixed){
  
  SIGMA <-  diag(3)
  
  const <- function(x, p0){(qnorm(1-alpha)*t2(x, p0) + qnorm(1-beta)*t3(x, p0))^2}
  
  if(fixed){
    
    rho <- -log(p11/p0)
    
    if(case == 2){
      if(size == "small"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( 3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    0), 
                          upper = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                    qnorm(1 - alpha) - 
                                      log(0.95)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                    qnorm(1 - alpha) - 
                                      log(0.95)/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                          mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                   rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                   rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                          sigma = SIGMA) - 
                2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha), 
                                      qnorm(1 - alpha)), 
                            upper = c(qnorm(1 - alpha) - 
                                        log(0.95)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                      qnorm(1 - alpha) - 
                                        log(0.95)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                      qnorm(1 - alpha) - 
                                        log(0.95)/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                            mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                     rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                     rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                            sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = rho, 
                    sd = sqrt(((2/n2))*t1(p11, p0))) 
          })
        },  - log(RRgo), Inf)$value)
      }
      if(size == "large"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( 3 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                    0), 
                          upper = c(Inf, 
                                    Inf, 
                                    Inf), 
                          mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                   rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                   rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                          sigma = SIGMA) - 
                2 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                            upper = c(Inf, 
                                      Inf, 
                                      Inf), 
                            mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                     rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                     rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                            sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = rho, 
                    sd = sqrt(((2/n2))*t1(p11, p0))) 
          })
        },  - log(RRgo), Inf)$value)
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
                          mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                   rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                   rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                          sigma = SIGMA) - 
                2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha), 
                                      qnorm(1 - alpha)), 
                            upper = c(Inf, 
                                      Inf, 
                                      Inf), 
                            mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                     rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                     rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                            sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = rho, 
                    sd = sqrt(((2/n2))*t1(p11, p0)))  
          })
        },  - log(RRgo), Inf)$value)
      }
    }
    if(case == 3){
      if(size == "small"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( 3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    qnorm(1 - alpha)), 
                          upper = c(qnorm(1 - alpha) - 
                                      log(0.95)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                    qnorm(1 - alpha) - 
                                      log(0.95)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                          mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                   rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                   rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                          sigma = SIGMA) - 
                2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha), 
                                      qnorm(1 - alpha)), 
                            upper = c(qnorm(1 - alpha) - 
                                        log(0.95)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                      qnorm(1 - alpha) - 
                                        log(0.95)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                      qnorm(1 - alpha) - 
                                        log(0.95)/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                            mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                     rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                     rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                            sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = rho, 
                    sd = sqrt(((2/n2))*t1(p11, p0)))
          })
        },  - log(RRgo), Inf)$value)
      }
      if(size == "large"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( 3 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                    qnorm(1 - alpha)), 
                          upper = c(Inf, 
                                    Inf, 
                                    Inf), 
                          mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                   rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                   rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                          sigma = SIGMA) - 
                2 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                            upper = c(Inf, 
                                      Inf, 
                                      Inf), 
                            mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                     rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                     rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                            sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = rho, 
                    sd = sqrt(((2/n2))*t1(p11, p0))) 
          })
        },  - log(RRgo), Inf)$value)
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
                      mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                               rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                               rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                      sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = rho, 
                    sd = sqrt(((2/n2))*t1(p11, p0))) 
          })
        },  - log(RRgo), Inf)$value)
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
                              upper = c(qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                        qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                        qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                              mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                       x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                       x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                              sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                          qnorm(1 - alpha), 
                                          qnorm(1 - alpha)), 
                                upper = c(qnorm(1 - alpha) - 
                                            log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                          qnorm(1 - alpha) - 
                                            log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                          qnorm(1 - alpha) - 
                                            log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                                mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                         x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                         x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*t1(x, p0))) * 
                  prior_binary(x, w, p11, p12, in1, in2)
              })
            },  - log(RRgo), Inf)$value   
          })
        },  0, 1)$value)
      }
      if(size == "large"){
        return(integrate(function(x){
          sapply(x, function(x){
            integrate(function(y){
              sapply(y, function(y){
                ( 3 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                        qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                        0), 
                              upper = c(Inf, 
                                        Inf, 
                                        Inf), 
                              mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                       x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                       x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                              sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                            log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                          qnorm(1 - alpha) - 
                                            log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                          qnorm(1 - alpha) - 
                                            log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                                upper = c(Inf, 
                                          Inf, 
                                          Inf), 
                                mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                         x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                         x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*t1(x, p0))) * 
                  prior_binary(x, w, p11, p12, in1, in2) 
              })
            },  - log(RRgo), Inf)$value
          })
        },  0, 1)$value)
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
                              mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                       x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                       x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                              sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                          qnorm(1 - alpha), 
                                          qnorm(1 - alpha)), 
                                upper = c(Inf, 
                                          Inf, 
                                          Inf), 
                                mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                         x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                         x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*t1(x, p0))) * 
                  prior_binary(x, w, p11, p12, in1, in2) 
              })
            },  - log(RRgo), Inf)$value
          })
        },  0, 1)$value)
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
                              upper = c(qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                        qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                        qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                              mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                       x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                       x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                              sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                          qnorm(1 - alpha), 
                                          qnorm(1 - alpha)), 
                                upper = c(qnorm(1 - alpha) - 
                                            log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                          qnorm(1 - alpha) - 
                                            log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                          qnorm(1 - alpha) - 
                                            log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                                mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                         x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                         x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*t1(x, p0))) * 
                  prior_binary(x, w, p11, p12, in1, in2) 
              })
            },  - log(RRgo), Inf)$value
          })
        },  0, 1)$value)
      }
      if(size == "large"){
        return(integrate(function(x){
          sapply(x, function(x){
            integrate(function(y){
              sapply(y, function(y){
                ( 3 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                        qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                        qnorm(1 - alpha)), 
                              upper = c(Inf, 
                                        Inf, 
                                        Inf), 
                              mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                       x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                       x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                              sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                            log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                          qnorm(1 - alpha) - 
                                            log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                          qnorm(1 - alpha) - 
                                            log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                                upper = c(Inf, 
                                          Inf, 
                                          Inf), 
                                mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                         x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                         x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*t1(x, p0))) * 
                  prior_binary(x, w, p11, p12, in1, in2) 
              })
            },  - log(RRgo), Inf)$value
          })
        },  0, 1)$value)
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
                          mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                   x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                   x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                          sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*t1(x, p0))) * 
                  prior_binary(x, w, p11, p12, in1, in2) 
              })
            },  - log(RRgo), Inf)$value
          })
        },  0, 1)$value)
      }
    }  
    
  }
  
}

#' @rdname utility_multitrial_binary 
#' @export
utility3_binary <-  function(n2, RRgo, w, p0, p11, p12, in1, in2,
                            alpha, beta, 
                            c2, c3, c02, c03, 
                            K, N, S,
                            b1, b2, b3,
                            case, fixed){
  
  
  n3  <-  En3_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                     p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, fixed = fixed)
  
  n3  <- ceiling(n3)
  
  if(round(n3/2) != n3 / 2) {n3 = n3 + 1}
  
  if(n2+3*n3>N){
    
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
    
  }else{
    
    pg    <-  Epgo_binary(RRgo = RRgo, n2 = n2, p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, fixed = fixed)
    
    K2    <-  c02 + c2 * n2         # cost phase II
    K3    <-  c03 * pg + c3 * n3    # cost phase III
    
    if(K2+3*K3>K){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
    }else{
      # probability of a successful program; small, medium, large effect size
      prob1 <-  EPsProg3_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                               p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                               case = case, size = "small", fixed = fixed)
      prob3 <-  EPsProg3_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                               p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                               case = case, size = "large", fixed = fixed)
      prob2 <-  EPsProg3_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                               p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                               case = case, size = "all", fixed = fixed) - prob1 - prob3 
      
      SP    <-  prob1 + prob2 + prob3
      
      if(SP<S){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        
      }else{
        
        G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
        
        EU    <-  - K2 - 3*K3 + G
        
        return(c(EU, 3*n3, SP, pg, K2, 3*K3, prob1, prob2, prob3))
      }
    }
  }
}




#########################
# Four phase III trials #
#########################
# Case 3: Strategy 3/4; at least three trials significant, the treatment effect 
# of the other one at least showing in the same direction

#' @rdname EPsProg_multitrial_binary
#' @export
EPsProg4_binary <-  function(RRgo, n2, alpha, beta, p0, w, p11, p12, in1, in2, case, size, fixed){
  
  const <- function(x, p0){(qnorm(1-alpha)*t2(x, p0) + qnorm(1-beta)*t3(x, p0))^2}
  
  SIGMA <-  diag(4)
  
  if(fixed){
    
    rho <- -log(p11/p0)
    
    if(size == "small"){
      return(integrate(function(y){
        sapply(y, function(y){
          ( 4 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                  qnorm(1 - alpha), 
                                  qnorm(1 - alpha), 
                                  0), 
                        upper = c(qnorm(1 - alpha) - log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                  qnorm(1 - alpha) - log(0.95)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                  qnorm(1 - alpha) - log(0.95)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                  qnorm(1 - alpha) - log(0.95)/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                        mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                 rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                 rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                 rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                        sigma = SIGMA)  - 
              3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    qnorm(1 - alpha)), 
                          upper = c(qnorm(1 - alpha) - 
                                      log(0.95)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                    qnorm(1 - alpha) - 
                                      log(0.95)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                    qnorm(1 - alpha) - 
                                      log(0.95)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                    qnorm(1 - alpha) - 
                                      log(0.95)/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                          mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                   rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                   rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                   rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                          sigma = SIGMA) ) * 
            dnorm(y, 
                  mean = rho, 
                  sd = sqrt(((2/n2))*t1(p11, p0))) 
        })
      },  - log(RRgo), Inf)$value)
    }
    if(size == "large"){
      return(integrate(function(y){
        sapply(y, function(y){
          ( 4 * pmvnorm(lower = c(qnorm(1 - alpha) - log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                  qnorm(1 - alpha) - log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                  qnorm(1 - alpha) - log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                  0), 
                        upper = c(Inf, 
                                  Inf, 
                                  Inf, 
                                  Inf), 
                        mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                 rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                 rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                 rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                        sigma = SIGMA)  - 
              3 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                          upper = c(Inf, 
                                    Inf, 
                                    Inf, 
                                    Inf), 
                          mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                   rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                   rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                   rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                          sigma = SIGMA) ) * 
            dnorm(y, 
                  mean = rho, 
                  sd = sqrt(((2/n2))*t1(p11, p0))) 
        })
      },  - log(RRgo), Inf)$value)
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
                        mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                 rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                 rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                 rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                        sigma = SIGMA) - 
              3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    qnorm(1 - alpha)), 
                          upper = c(Inf, 
                                    Inf, 
                                    Inf, 
                                    Inf), 
                          mean = c(rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                   rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                   rho/sqrt(t1(p11, p0)*y^2/const(p11, p0)), 
                                   rho/sqrt(t1(p11, p0)*y^2/const(p11, p0))), 
                          sigma = SIGMA) ) * 
            dnorm(y, 
                  mean = rho, 
                  sd = sqrt(((2/n2))*t1(p11, p0))) 
        })
      },  - log(RRgo), Inf)$value)
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
                            upper = c(qnorm(1 - alpha) - log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                      qnorm(1 - alpha) - log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                      qnorm(1 - alpha) - log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                      qnorm(1 - alpha) - log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                            mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                     x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                     x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                     x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                            sigma = SIGMA)  - 
                  3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                        qnorm(1 - alpha), 
                                        qnorm(1 - alpha), 
                                        qnorm(1 - alpha)), 
                              upper = c(qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                        qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                        qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                        qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                              mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                       x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                       x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                       x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                              sigma = SIGMA) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(((2/n2))*t1(x, p0))) * 
                prior_binary(x, w, p11, p12, in1, in2)
            })
          },  - log(RRgo), Inf)$value   
        })
      },  0, 1)$value)
    }
    if(size == "large"){
      return(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            sapply(y, function(y){
              ( 4 * pmvnorm(lower = c(qnorm(1 - alpha) - log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                      qnorm(1 - alpha) - log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                      qnorm(1 - alpha) - log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                      0), 
                            upper = c(Inf, 
                                      Inf, 
                                      Inf, 
                                      Inf), 
                            mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                     x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                     x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                     x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                            sigma = SIGMA)  - 
                  3 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                        qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                        qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                        qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                              upper = c(Inf, 
                                        Inf, 
                                        Inf, 
                                        Inf), 
                              mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                       x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                       x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                       x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                              sigma = SIGMA) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(((2/n2))*t1(x, p0))) * 
                prior_binary(x, w, p11, p12, in1, in2) 
            })
          },  - log(RRgo), Inf)$value
        })
      },  0, 1)$value)
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
                            mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                     x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                     x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                     x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                            sigma = SIGMA) - 
                  3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                        qnorm(1 - alpha), 
                                        qnorm(1 - alpha), 
                                        qnorm(1 - alpha)), 
                              upper = c(Inf, 
                                        Inf, 
                                        Inf, 
                                        Inf), 
                              mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                       x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                       x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                       x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                              sigma = SIGMA) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(((2/n2))*t1(x, p0))) * 
                prior_binary(x, w, p11, p12, in1, in2) 
            })
          },  - log(RRgo), Inf)$value
        })
      },  0, 1)$value)
    }
    
  }
  
}




#' @rdname utility_multitrial_binary 
#' @export
utility4_binary <-  function(n2, RRgo, w, p0, p11, p12, in1, in2,
                             alpha, beta, 
                             c2, c3, c02, c03, 
                             K, N, S,
                             b1, b2, b3,
                             case, fixed){
  
  
  n3  <-  En3_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                     p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, fixed = fixed)
  
  n3  <- ceiling(n3)
  
  if(round(n3/2) != n3 / 2) {n3 = n3 + 1}
  
  if(n2+4*n3>N){
    
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
    
  }else{
    
    pg    <-  Epgo_binary(RRgo = RRgo, n2 = n2, p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, fixed = fixed)
    
    K2    <-  c02 + c2 * n2         # cost phase II
    K3    <-  c03 * pg + c3 * n3    # cost phase III
    
    if(K2+4*K3>K){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
    }else{
      # probability of a successful program; small, medium, large effect size
      prob1 <-  EPsProg4_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                                p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                                case = case, size = "small", fixed = fixed)
      prob3 <-  EPsProg4_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                                p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                                case = case, size = "large", fixed = fixed)
      prob2 <-  EPsProg4_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                                p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                                case = case, size = "all", fixed = fixed) - prob1 - prob3
      
      SP    <-  prob1 + prob2 + prob3
      
      if(SP<S){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        
      }else{
        
        G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
        
        EU    <-  - K2 - 4*K3 + G
        
        return(c(EU, 4*n3, SP, pg, K2, 4*K3, prob1, prob2, prob3))
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
#' @return The output of the the function `Epgo23_binary()` is the probability to to a third phase III trial.
#' @examples res <- Epgo23_binary(RRgo = 0.8, n2 = 50,  p0 = 0.3, w = 0.3, alpha = 0.025, beta = 0.1,
#'                                p11 =  0.3, p12 = 0.5, in1 = 300, in2 = 600)
#' @editor Johannes Cepicka
#' @editDate 2022-05-09
#' @export
Epgo23_binary <-  function(RRgo, n2, alpha, beta, p0, w, p11, p12, in1, in2){
  
  SIGMA <-  diag(2)
  const     <-  (qnorm(1 - alpha) + qnorm(1 - beta))^2
  
  integrate(function(x){                    
    sapply(x, function(x){
      integrate(function(y){
        sapply(y, function(y){
          2 * (pmvnorm(lower = c(qnorm(1 - alpha), 
                                 0), 
                       upper = c(Inf, 
                                 qnorm(1 - alpha)), 
                       mean = c(x/sqrt(y^2/const), 
                                x/sqrt(y^2/const)), 
                       sigma = SIGMA)) * 
            dnorm(y, 
                  mean = x, 
                  sd = sqrt((2/n2)*t1(p11, p0))) * 
            prior_binary(x, w, p11, p12, in1, in2)
        })
      },  RRgo, Inf)$value
    })
  },  - Inf, Inf)$value
} 

#' Expected probability of a successful program deciding between two or three phase III trials for a binary distributed outcome
#'
#' The function `EPsProg23_binary()` calculates the expected probability of a successful program
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
#' @return The output of the the function `EPsProg23_binary()` is the expected probability of a successful program.
#' @examples res <- EPsProg23_binary(RRgo = 0.8, n2 = 50,  alpha = 0.025, beta = 0.1, 
#'                                  w = 0.6,  p0 = 0.3, p11 =  0.3, p12 = 0.5, 
#'                                  in1 = 300, in2 = 600, case = 2, size = "small",
#'                                  ymin = 0.5)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export
EPsProg23_binary <-  function(RRgo, n2, alpha, beta, w, p0, p11, p12, in1, in2, case, size, ymin){
  # Option 2.1: first two phase III trials are successful: no third phase III trial
  # Option 2.2: one of the two first phase III trials successful, the treatment
  #  effect of the other one points in the same direction: 
  #  conduct third phase III trial with N3 = N3(ymin)
  
  SIGMA <-  diag(2)
  SIGMA3<-  diag(3)
    
  const <- function(x, p0){(qnorm(1-alpha)*t2(x, p0) + qnorm(1-beta)*t3(x, p0))^2}
  
  if(case == 2){ # Option 2.1
    if(size == "small"){
      return(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            sapply(y, function(y){
              ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                  qnorm(1 - alpha)), 
                        upper = c(qnorm(1 - alpha) - log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                  qnorm(1 - alpha) - log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                        mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                 x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                        sigma = SIGMA)  - 
                  pmvnorm(lower = c(qnorm(1 - alpha) - log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                    qnorm(1 - alpha) - log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                          upper = c(qnorm(1 - alpha) - log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                    qnorm(1 - alpha) - log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                          mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                   x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                          sigma = SIGMA)) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(((2/n2))*t1(x, p0))) * 
                prior_binary(x, w, p11, p12, in1, in2) 
            })
          },  - log(RRgo), Inf)$value
        })
      },  0, 1)$value) 
    }
    if(size == "large"){
      return(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            sapply(y, function(y){
              ( pmvnorm(lower = c(qnorm(1 - alpha) - 
                                    log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                  qnorm(1 - alpha) - 
                                    log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                        upper = c(Inf, Inf), 
                        mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                 x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                        sigma = SIGMA) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(((2/n2))*t1(x, p0))) * 
                prior_binary(x, w, p11, p12, in1, in2) 
            })
          },  - log(RRgo), Inf)$value
        })
      },  0, 1)$value)    
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
                        mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                 x/sqrt(t1(x, p0)*y^2/const(x, p0))), 
                        sigma = SIGMA) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(((2/n2))*t1(x, p0))) * 
                prior_binary(x, w, p11, p12, in1, in2) 
            })
          },  - log(RRgo), Inf)$value
        })
      },  0, 1)$value)    
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
                            upper = c(qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                      qnorm(1 - alpha), 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(t1(x, p0)*ymin^2/const(x, p0))), 
                            mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                     x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                     x/sqrt(t1(x, p0)*ymin^2/const(x, p0))), 
                            sigma = SIGMA3)  - 
                      pmvnorm(lower = c(qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                        0, 
                                        qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(t1(x, p0)*ymin^2/const(x, p0))), 
                              upper = c(qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                        qnorm(1 - alpha), 
                                        qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(t1(x, p0)*ymin^2/const(x, p0))), 
                              mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                       x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                       x/sqrt(t1(x, p0)*ymin^2/const(x, p0))), 
                              sigma = SIGMA3)) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(((2/n2))*t1(x, p0))) * 
                prior_binary(x, w, p11, p12, in1, in2) 
            })
          },  - log(RRgo), Inf)$value
        })
      },  0, 1)$value) 
    }
    if(size == "large"){
      return(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            sapply(y, function(y){
              2 * ( pmvnorm(lower = c(qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                      0, 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(t1(x, p0)*ymin^2/const(x, p0))), 
                            upper = c(Inf, 
                                      qnorm(1 - alpha), 
                                      Inf), 
                            mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                     x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                     x/sqrt(t1(x, p0)*ymin^2/const(x, p0))), 
                            sigma = SIGMA3) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(((2/n2))*t1(x, p0))) * 
                prior_binary(x, w, p11, p12, in1, in2) 
            })
          },  - log(RRgo), Inf)$value
        })
      },  0, 1)$value)    
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
                            mean = c(x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                     x/sqrt(t1(x, p0)*y^2/const(x, p0)), 
                                     x/sqrt(t1(x, p0)*ymin^2/const(x, p0))), 
                            sigma = SIGMA3) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(((2/n2))*t1(x, p0))) * 
                prior_binary(x, w, p11, p12, in1, in2) 
            })
          },  - log(RRgo), Inf)$value
        })
      },  0, 1)$value)    
    }
  }
  
}




#' Utility function for multitrial programs deciding between two or three phase III trials for a binary distributed outcome
#'
#' The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters and the expected probability of a successful program. 
#' The utility is in further step maximized by the `optimal_multitrial_binary()` function.
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
#' @return The output of the the function `utility23_binary()` is the expected utility of the program depending on whether two or three phase III trials are performed.
#' @examples #res <- utility23_binary(n2 = 50, RRgo = 0.8,  w = 0.3, 
#'           #                       alpha = 0.05, beta = 0.1,
#'           #                       p0 = 0.6, p11 =  0.3, p12 = 0.5, 
#'           #                       in1 = 300, in2 = 600, alpha = 0.025, beta = 0.1,
#'           #                       c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'           #                       b1 = 1000, b2 = 2000, b3 = 3000)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export
utility23_binary <-  function(n2, RRgo, w, p0, p11, p12, in1, in2,
                             alpha, beta, 
                             c2, c3, c02, c03,
                             b1, b2, b3){
  
  pg    <-  Epgo_binary(RRgo = RRgo, n2 = n2, p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, fixed= FALSE)
  
  
  n3  <-  En3_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                     p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, fixed = fixed)
  
  n3  <- ceiling(n3)
  
  if(round(n3/2) != n3 / 2) {n3 = n3 + 1}
  
  prob1 <-  EPsProg23_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                            p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                            case = 2, size = "small", ymin = ymin)
  prob3 <-  EPsProg23_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                            p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                            case = 2, size = "large", ymin = ymin)
  prob2 <-  EPsProg23_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                            p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                            case = 2, size = "all", ymin = ymin) - prob1 - prob3
  
  
    
  pg3   <-  Epgo23_binary(RRgo = RRgo, alpha = alpha, beta = beta, n2 = n2, 
                          p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2)
    
  n33   <-  (4 * (qnorm(1 - alpha) + qnorm(1 - beta))^2)/(ymin^2) 
  
  n33  <- ceiling(n33*pg3)
  
  if(round(n33/2) != n33 / 2) {n33 <- n33 + 1}
  
  prob13 <-  EPsProg23_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                             p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                             case = 3, size = "small", ymin = ymin)
  prob33 <-  EPsProg23_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                             p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                             case = 3, size = "large", ymin = ymin)
  prob23 <-  EPsProg23_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                             p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                             case = 3, size = "all", ymin = ymin) - prob1 - prob3
  
  
   K2   <-  c02 + c2 * n2         # cost phase II
   K3    <-  c03 * pg + c3 * n3    # cost for one of the first two phase III trials in case of go decision
   K33    <-  pg3 * c03 + c3 * n33  # cost for the third phase III trial in case of third phase III trial
    
   
    
      
      SP    <-  prob1 + prob2 + prob3 + prob13 + prob23 + prob33
      
      
        
      G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 +
        b1 * prob13 + b2 * prob23 + b3 * prob33   # gain
      
      EU    <-  - K2 - 2 * K3 - K33 + G
        
        return(c(EU, 2*n3, SP, pg, 2*K3, K33, prob1, prob2, prob3,n2, 2*n3, pg3, n33, prob13, prob23, prob33))
}
