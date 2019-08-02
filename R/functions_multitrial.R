
########################
# Two phase III trials #
########################
# Case 1: Strategy 1/2; at least one trial significant, the treatment effect of the other one at least showing in the same direction 
# Case 2: Strategy 2/2; both trials significant 

# Expected probability of a successful program
EPsProg2 <-  function(HRgo, d2, alpha, beta, w, hr1, hr2, id1, id2, case, size, fixed){
  
  SIGMA <-  diag(2)
  c     <-  (qnorm(1 - alpha) + qnorm(1 - beta))^2
  
  if(fixed){

    theta <- -log(hr1)
    
    if(case == 1){
      if(size == "small"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( pmvnorm(lower = c(0, 
                                0), 
                      upper = c(qnorm(1 - alpha) - log(0.95)/sqrt(y^2/c),
                                qnorm(1 - alpha) - log(0.95)/sqrt(y^2/c)), 
                      mean = c(theta/sqrt(y^2/c), 
                               theta/sqrt(y^2/c)), 
                      sigma = SIGMA) - 
                pmvnorm(lower = c(0, 
                                  0), 
                        upper = c(qnorm(1 - alpha), 
                                  qnorm(1 - alpha)), 
                        mean = c(theta/sqrt(y^2/c), 
                                 theta/sqrt(y^2/c)), 
                        sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = theta, 
                    sd = sqrt(4/d2))  
          })
        },  - log(HRgo), Inf)$value)  
      }
      if(size == "large"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( pmvnorm(lower = c(0, 
                                0), 
                      upper = c(Inf, 
                                Inf), 
                      mean = c(theta/sqrt(y^2/c), 
                               theta/sqrt(y^2/c)), 
                      sigma = SIGMA) - 
                pmvnorm(lower = c(0, 
                                  0), 
                        upper = c(qnorm(1 - alpha) - 
                                    log(0.85)/sqrt(y^2/c), 
                                  qnorm(1 - alpha) - 
                                    log(0.85)/sqrt(y^2/c)), 
                        mean = c(theta/sqrt(y^2/c), 
                                 theta/sqrt(y^2/c)), 
                        sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = theta, 
                    sd = sqrt(4/d2)) 
          })
        },  - log(HRgo), Inf)$value)  
      }
      if(size == "all"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( pmvnorm(lower = c(0, 
                                0), 
                      upper = c(Inf, 
                                Inf), 
                      mean = c(theta/sqrt(y^2/c), 
                               theta/sqrt(y^2/c)), 
                      sigma = SIGMA) - 
                pmvnorm(lower = c(0, 
                                  0), 
                        upper = c(qnorm(1 - alpha), 
                                  qnorm(1 - alpha)), 
                        mean = c(theta/sqrt(y^2/c), 
                                 theta/sqrt(y^2/c)), 
                        sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = theta, 
                    sd = sqrt(4/d2))
          })
        },  - log(HRgo), Inf)$value)     
      }
    }
    if(case == 2){
      if(size == "small"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                qnorm(1 - alpha)), 
                      upper = c(qnorm(1 - alpha) - 
                                  log(0.85)/sqrt(y^2/c), 
                                qnorm(1 - alpha) - 
                                  log(0.85)/sqrt(y^2/c)), 
                      mean = c(theta/sqrt(y^2/c), 
                               theta/sqrt(y^2/c)), 
                      sigma = SIGMA) - 
                pmvnorm(lower = c(qnorm(1 - alpha) - 
                                    log(0.95)/sqrt(y^2/c), 
                                  qnorm(1 - alpha) - 
                                    log(0.95)/sqrt(y^2/c)), 
                        upper = c(qnorm(1 - alpha) - 
                                    log(0.85)/sqrt(y^2/c), 
                                  qnorm(1 - alpha) - 
                                    log(0.85)/sqrt(y^2/c)), 
                        mean = c(theta/sqrt(y^2/c), 
                                 theta/sqrt(y^2/c)), 
                        sigma = SIGMA)) * 
              dnorm(y, 
                    mean = theta, 
                    sd = sqrt(4/d2)) 
          })
        },  - log(HRgo), Inf)$value) 
      }
      if(size == "large"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( pmvnorm(lower = c(qnorm(1 - alpha) - 
                                  log(0.85)/sqrt(y^2/c), 
                                qnorm(1 - alpha) - 
                                  log(0.85)/sqrt(y^2/c)), 
                      upper = c(Inf, 
                                Inf), 
                      mean = c(theta/sqrt(y^2/c), 
                               theta/sqrt(y^2/c)), 
                      sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = theta, 
                    sd = sqrt(4/d2)) 
          })
        },  - log(HRgo), Inf)$value)    
      }
      if(size == "all"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                qnorm(1 - alpha)), 
                      upper = c(Inf, 
                                Inf), 
                      mean = c(theta/sqrt(y^2/c), 
                               theta/sqrt(y^2/c)), 
                      sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = theta, 
                    sd = sqrt(4/d2)) 
          })
        },  - log(HRgo), Inf)$value)    
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
                          upper = c(qnorm(1 - alpha) - log(0.95)/sqrt(y^2/c),
                                    qnorm(1 - alpha) - log(0.95)/sqrt(y^2/c)), 
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
                        sd = sqrt(4/d2)) * 
                  prior_tte(x, w, hr1, hr2, id1, id2)  
              })
            },  - log(HRgo), Inf)$value  
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
                            upper = c(qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(y^2/c)), 
                            mean = c(x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c)), 
                            sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(4/d2)) * 
                  prior_tte(x, w, hr1, hr2, id1, id2)  
              })
            },  - log(HRgo), Inf)$value  
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
                        sd = sqrt(4/d2)) * 
                  prior_tte(x, w, hr1, hr2, id1, id2)  
              })
            },  - log(HRgo), Inf)$value  
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
                                      log(0.85)/sqrt(y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c)), 
                          mean = c(x/sqrt(y^2/c), 
                                   x/sqrt(y^2/c)), 
                          sigma = SIGMA) - 
                    pmvnorm(lower = c(qnorm(1 - alpha) - 
                                        log(0.95)/sqrt(y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.95)/sqrt(y^2/c)), 
                            upper = c(qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(y^2/c)), 
                            mean = c(x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c)), 
                            sigma = SIGMA)) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(4/d2)) * 
                  prior_tte(x, w, hr1, hr2, id1, id2) 
              })
            },  - log(HRgo), Inf)$value
          })
        },  - Inf, Inf)$value) 
      }
      if(size == "large"){
        return(integrate(function(x){
          sapply(x, function(x){
            integrate(function(y){
              sapply(y, function(y){
                ( pmvnorm(lower = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c)), 
                          upper = c(Inf, 
                                    Inf), 
                          mean = c(x/sqrt(y^2/c), 
                                   x/sqrt(y^2/c)), 
                          sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(4/d2)) * 
                  prior_tte(x, w, hr1, hr2, id1, id2) 
              })
            },  - log(HRgo), Inf)$value
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
                        sd = sqrt(4/d2)) * 
                  prior_tte(x, w, hr1, hr2, id1, id2) 
              })
            },  - log(HRgo), Inf)$value
          })
        },  - Inf, Inf)$value)    
      }
    }
    
  }
  

  
}

# Utility function
utility2 <-  function(d2, HRgo, w, hr1, hr2, id1, id2,
                      alpha, beta, xi2, xi3,
                      c2, c3, c02, c03, 
                      K, N, S,
                      b1, b2, b3,
                      case, fixed){ 
  
  
  d3    <-  Ed3_tte(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                    fixed = fixed)
  
  # sample size is rounded up to next even natural number
  n2  <- ceiling(d2 * (1/xi2))
  if(round(n2/2) != n2 / 2) {n2 <- n2 + 1}
  
  n3  <- ceiling(d3 * (1/xi3))
  if(round(n3/2) != n3 / 2) {n3 <- n3 + 1}
  
  # expected number of events is rounded to natural number for presentation
  d3  <- ceiling(d3)
  
  if(n2+2*n3>N){
    
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
    
  }else{
    
    pg    <-  Epgo_tte(HRgo = HRgo, d2 = d2, w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                       fixed = fixed)
    
    K2    <-  c02 + c2 * n2         # cost phase II
    K3    <-  c03 * pg + c3 * n3    # cost phase III
    
    if(K2+2*K3>K){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
    }else{
      # probability of a successful program; small, medium, large effect size
      prob1 <-  EPsProg2(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta,
                            w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                            case = case, size = "small", fixed = fixed)
      prob3 <-  EPsProg2(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta,
                            w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                            case = case, size = "large", fixed = fixed)
      prob2 <-  EPsProg2(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta,
                            w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                            case = case, size = "all", fixed = fixed) - prob1 - prob3
      
      SP    <-  prob1 + prob2 + prob3
      
      if(SP<S){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        
      }else{
        
        G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
        
        EU    <-  - K2 - 2*K3 + G
        
        return(c(EU, 2*d3, SP, pg, K2, 2*K3, prob1, prob2, prob3, n2, 2*n3))
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

# Expected probability of a successful program
EPsProg3 <-  function(HRgo, d2, alpha, beta, w, hr1, hr2, id1, id2, case, size, fixed){
  
  SIGMA <-  diag(3)
  c     <-  (qnorm(1 - alpha) + qnorm(1 - beta))^2
  
  if(fixed){
    
    theta <- -log(hr1)
    
    if(case == 2){
      if(size == "small"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( 3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    0), 
                          upper = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.95)/sqrt(y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.95)/sqrt(y^2/c)), 
                          mean = c(theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c)), 
                          sigma = SIGMA) - 
                2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha), 
                                      qnorm(1 - alpha)), 
                            upper = c(qnorm(1 - alpha) - 
                                        log(0.95)/sqrt(y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.95)/sqrt(y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.95)/sqrt(y^2/c)), 
                            mean = c(theta/sqrt(y^2/c), 
                                     theta/sqrt(y^2/c), 
                                     theta/sqrt(y^2/c)), 
                            sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = theta, 
                    sd = sqrt(4/d2)) 
          })
        },  - log(HRgo), Inf)$value)
      }
      if(size == "large"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( 3 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c), 
                                    0), 
                          upper = c(Inf, 
                                    Inf, 
                                    Inf), 
                          mean = c(theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c)), 
                          sigma = SIGMA) - 
                2 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(y^2/c)), 
                            upper = c(Inf, 
                                      Inf, 
                                      Inf), 
                            mean = c(theta/sqrt(y^2/c), 
                                     theta/sqrt(y^2/c), 
                                     theta/sqrt(y^2/c)), 
                            sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = theta, 
                    sd = sqrt(4/d2)) 
          })
        },  - log(HRgo), Inf)$value)
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
                          mean = c(theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c)), 
                          sigma = SIGMA) - 
                2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha), 
                                      qnorm(1 - alpha)), 
                            upper = c(Inf, 
                                      Inf, 
                                      Inf), 
                            mean = c(theta/sqrt(y^2/c), 
                                     theta/sqrt(y^2/c), 
                                     theta/sqrt(y^2/c)), 
                            sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = theta, 
                    sd = sqrt(4/d2))  
          })
        },  - log(HRgo), Inf)$value)
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
                                      log(0.95)/sqrt(y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.95)/sqrt(y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c)), 
                          mean = c(theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c)), 
                          sigma = SIGMA) - 
                2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha), 
                                      qnorm(1 - alpha)), 
                            upper = c(qnorm(1 - alpha) - 
                                        log(0.95)/sqrt(y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.95)/sqrt(y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.95)/sqrt(y^2/c)), 
                            mean = c(theta/sqrt(y^2/c), 
                                     theta/sqrt(y^2/c), 
                                     theta/sqrt(y^2/c)), 
                            sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = theta, 
                    sd = sqrt(4/d2))
          })
        },  - log(HRgo), Inf)$value)
      }
      if(size == "large"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( 3 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c), 
                                    qnorm(1 - alpha)), 
                          upper = c(Inf, 
                                    Inf, 
                                    Inf), 
                          mean = c(theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c)), 
                          sigma = SIGMA) - 
                2 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(y^2/c)), 
                            upper = c(Inf, 
                                      Inf, 
                                      Inf), 
                            mean = c(theta/sqrt(y^2/c), 
                                     theta/sqrt(y^2/c), 
                                     theta/sqrt(y^2/c)), 
                            sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = theta, 
                    sd = sqrt(4/d2)) 
          })
        },  - log(HRgo), Inf)$value)
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
                      mean = c(theta/sqrt(y^2/c), 
                               theta/sqrt(y^2/c), 
                               theta/sqrt(y^2/c)), 
                      sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = theta, 
                    sd = sqrt(4/d2)) 
          })
        },  - log(HRgo), Inf)$value)
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
                                          log(0.85)/sqrt(y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(y^2/c)), 
                              mean = c(x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c)), 
                              sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                          qnorm(1 - alpha), 
                                          qnorm(1 - alpha)), 
                                upper = c(qnorm(1 - alpha) - 
                                            log(0.95)/sqrt(y^2/c), 
                                          qnorm(1 - alpha) - 
                                            log(0.95)/sqrt(y^2/c), 
                                          qnorm(1 - alpha) - 
                                            log(0.95)/sqrt(y^2/c)), 
                                mean = c(x/sqrt(y^2/c), 
                                         x/sqrt(y^2/c), 
                                         x/sqrt(y^2/c)), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(4/d2)) * 
                  prior_tte(x, w, hr1, hr2, id1, id2)
              })
            },  - log(HRgo), Inf)$value   
          })
        },  - Inf, Inf)$value)
      }
      if(size == "large"){
        return(integrate(function(x){
          sapply(x, function(x){
            integrate(function(y){
              sapply(y, function(y){
                ( 3 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c), 
                                        0), 
                              upper = c(Inf, 
                                        Inf, 
                                        Inf), 
                              mean = c(x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c)), 
                              sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                            log(0.85)/sqrt(y^2/c), 
                                          qnorm(1 - alpha) - 
                                            log(0.85)/sqrt(y^2/c), 
                                          qnorm(1 - alpha) - 
                                            log(0.85)/sqrt(y^2/c)), 
                                upper = c(Inf, 
                                          Inf, 
                                          Inf), 
                                mean = c(x/sqrt(y^2/c), 
                                         x/sqrt(y^2/c), 
                                         x/sqrt(y^2/c)), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(4/d2)) * 
                  prior_tte(x, w, hr1, hr2, id1, id2) 
              })
            },  - log(HRgo), Inf)$value
          })
        },  - Inf, Inf)$value)
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
                        sd = sqrt(4/d2)) * 
                  prior_tte(x, w, hr1, hr2, id1, id2) 
              })
            },  - log(HRgo), Inf)$value
          })
        },  - Inf, Inf)$value)
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
                                          log(0.95)/sqrt(y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c)), 
                              mean = c(x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c)), 
                              sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                          qnorm(1 - alpha), 
                                          qnorm(1 - alpha)), 
                                upper = c(qnorm(1 - alpha) - 
                                            log(0.95)/sqrt(y^2/c), 
                                          qnorm(1 - alpha) - 
                                            log(0.95)/sqrt(y^2/c), 
                                          qnorm(1 - alpha) - 
                                            log(0.95)/sqrt(y^2/c)), 
                                mean = c(x/sqrt(y^2/c), 
                                         x/sqrt(y^2/c), 
                                         x/sqrt(y^2/c)), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(4/d2)) * 
                  prior_tte(x, w, hr1, hr2, id1, id2) 
              })
            },  - log(HRgo), Inf)$value
          })
        },  - Inf, Inf)$value)
      }
      if(size == "large"){
        return(integrate(function(x){
          sapply(x, function(x){
            integrate(function(y){
              sapply(y, function(y){
                ( 3 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c), 
                                        qnorm(1 - alpha)), 
                              upper = c(Inf, 
                                        Inf, 
                                        Inf), 
                              mean = c(x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c)), 
                              sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                            log(0.85)/sqrt(y^2/c), 
                                          qnorm(1 - alpha) - 
                                            log(0.85)/sqrt(y^2/c), 
                                          qnorm(1 - alpha) - 
                                            log(0.85)/sqrt(y^2/c)), 
                                upper = c(Inf, 
                                          Inf, 
                                          Inf), 
                                mean = c(x/sqrt(y^2/c), 
                                         x/sqrt(y^2/c), 
                                         x/sqrt(y^2/c)), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(4/d2)) * 
                  prior_tte(x, w, hr1, hr2, id1, id2) 
              })
            },  - log(HRgo), Inf)$value
          })
        },  - Inf, Inf)$value)
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
                        sd = sqrt(4/d2)) * 
                  prior_tte(x, w, hr1, hr2, id1, id2) 
              })
            },  - log(HRgo), Inf)$value
          })
        },  - Inf, Inf)$value)
      }
    }  
    
  }

}

# Utility function
utility3 <-  function(d2, HRgo, w, hr1, hr2, id1, id2,
                      alpha, beta, xi2, xi3,
                      c2, c3, c02, c03, 
                      K, N, S,
                      b1, b2, b3,
                      case, fixed){ 
  
  d3    <-  Ed3_tte(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                    fixed = fixed)
  
  # sample size is rounded up to next even natural number
  n2  <- ceiling(d2 * (1/xi2))
  if(round(n2/2) != n2 / 2) {n2 <- n2 + 1}
  
  n3  <- ceiling(d3 * (1/xi3))
  if(round(n3/2) != n3 / 2) {n3 <- n3 + 1}
  
  # expected number of events is rounded to natural number for presentation
  d3  <- ceiling(d3)
  
  if(n2+3*n3>N){
    
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
    
  }else{
    
    pg    <-  Epgo_tte(HRgo = HRgo, d2 = d2, w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                       fixed = fixed)
    
    K2    <-  c02 + c2 * n2         # cost phase II
    K3    <-  c03 * pg + c3 * n3    # cost phase III
    
    if(K2+3*K3>K){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
    }else{
      # probability of a successful program; small, medium, large effect size
      prob1 <-  EPsProg3(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta,
                         w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                         case = case, size = "small", fixed = fixed)
      prob3 <-  EPsProg3(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta,
                         w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                         case = case, size = "large", fixed = fixed)
      prob2 <-  EPsProg3(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta,
                         w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                         case = case, size = "all", fixed = fixed) - prob1 - prob3
      
      SP    <-  prob1 + prob2 + prob3
      
      if(SP<S){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        
      }else{
        
        G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
        
        EU    <-  - K2 - 3*K3 + G
        
        return(c(EU, 3*d3, SP, pg, K2, 3*K3, prob1, prob2, prob3, n2, 3*n3))
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

# Expected probability of a successful program
EPsProg4 <-  function(HRgo, d2, alpha, beta, w, hr1, hr2, id1, id2, case, size, fixed){
  
  SIGMA <-  diag(4)
  c     <-  (qnorm(1 - alpha) + qnorm(1 - beta))^2
  
  if(fixed){
    theta <- -log(hr1)
    
    if(size == "small"){
      return(integrate(function(y){
        sapply(y, function(y){
          ( 4 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                  qnorm(1 - alpha), 
                                  qnorm(1 - alpha), 
                                  0), 
                        upper = c(qnorm(1 - alpha) - log(0.85)/sqrt(y^2/c), 
                                  qnorm(1 - alpha) - log(0.95)/sqrt(y^2/c), 
                                  qnorm(1 - alpha) - log(0.95)/sqrt(y^2/c), 
                                  qnorm(1 - alpha) - log(0.95)/sqrt(y^2/c)), 
                        mean = c(theta/sqrt(y^2/c), 
                                 theta/sqrt(y^2/c), 
                                 theta/sqrt(y^2/c), 
                                 theta/sqrt(y^2/c)), 
                        sigma = SIGMA)  - 
              3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    qnorm(1 - alpha)), 
                          upper = c(qnorm(1 - alpha) - 
                                      log(0.95)/sqrt(y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.95)/sqrt(y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.95)/sqrt(y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.95)/sqrt(y^2/c)), 
                          mean = c(theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c)), 
                          sigma = SIGMA) ) * 
            dnorm(y, 
                  mean = theta, 
                  sd = sqrt(4/d2)) 
        })
      },  - log(HRgo), Inf)$value)
    }
    if(size == "large"){
      return(integrate(function(y){
        sapply(y, function(y){
          ( 4 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                    log(0.85)/sqrt(y^2/c), 
                                  qnorm(1 - alpha) - log(0.85)/sqrt(y^2/c), 
                                  qnorm(1 - alpha) - log(0.85)/sqrt(y^2/c), 
                                  0), 
                        upper = c(Inf, 
                                  Inf, 
                                  Inf, 
                                  Inf), 
                        mean = c(theta/sqrt(y^2/c), 
                                 theta/sqrt(y^2/c), 
                                 theta/sqrt(y^2/c), 
                                 theta/sqrt(y^2/c)), 
                        sigma = SIGMA)  - 
              3 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c)), 
                          upper = c(Inf, 
                                    Inf, 
                                    Inf, 
                                    Inf), 
                          mean = c(theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c)), 
                          sigma = SIGMA) ) * 
            dnorm(y, 
                  mean = theta, 
                  sd = sqrt(4/d2)) 
        })
      },  - log(HRgo), Inf)$value)
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
                        mean = c(theta/sqrt(y^2/c), 
                                 theta/sqrt(y^2/c), 
                                 theta/sqrt(y^2/c), 
                                 theta/sqrt(y^2/c)), 
                        sigma = SIGMA) - 
              3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    qnorm(1 - alpha)), 
                          upper = c(Inf, 
                                    Inf, 
                                    Inf, 
                                    Inf), 
                          mean = c(theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c)), 
                          sigma = SIGMA) ) * 
            dnorm(y, 
                  mean = theta, 
                  sd = sqrt(4/d2)) 
        })
      },  - log(HRgo), Inf)$value)
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
                            upper = c(qnorm(1 - alpha) - log(0.85)/sqrt(y^2/c), 
                                      qnorm(1 - alpha) - log(0.95)/sqrt(y^2/c), 
                                      qnorm(1 - alpha) - log(0.95)/sqrt(y^2/c), 
                                      qnorm(1 - alpha) - log(0.95)/sqrt(y^2/c)), 
                            mean = c(x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c)), 
                            sigma = SIGMA)  - 
                  3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                        qnorm(1 - alpha), 
                                        qnorm(1 - alpha), 
                                        qnorm(1 - alpha)), 
                              upper = c(qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(y^2/c)), 
                              mean = c(x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c)), 
                              sigma = SIGMA) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(4/d2)) * 
                prior_tte(x, w, hr1, hr2, id1, id2)
            })
          },  - log(HRgo), Inf)$value   
        })
      },  - Inf, Inf)$value)
    }
    if(size == "large"){
      return(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            sapply(y, function(y){
              ( 4 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(y^2/c), 
                                      qnorm(1 - alpha) - log(0.85)/sqrt(y^2/c), 
                                      qnorm(1 - alpha) - log(0.85)/sqrt(y^2/c), 
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
                  3 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c)), 
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
                      sd = sqrt(4/d2)) * 
                prior_tte(x, w, hr1, hr2, id1, id2) 
            })
          },  - log(HRgo), Inf)$value
        })
      },  - Inf, Inf)$value)
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
                      sd = sqrt(4/d2)) * 
                prior_tte(x, w, hr1, hr2, id1, id2) 
            })
          },  - log(HRgo), Inf)$value
        })
      },  - Inf, Inf)$value)
    }
    
  }
  
}

# Utility function
utility4 <-  function(d2, HRgo, w, hr1, hr2, id1, id2,
                      alpha, beta, xi2, xi3,
                      c2, c3, c02, c03, 
                      K, N, S,
                      b1, b2, b3,
                      case, fixed){ 
  
  d3    <-  Ed3_tte(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                    fixed = fixed)
  
  # sample size is rounded up to next even natural number
  n2  <- ceiling(d2 * (1/xi2))
  if(round(n2/2) != n2 / 2) {n2 <- n2 + 1}
  
  n3  <- ceiling(d3 * (1/xi3))
  if(round(n3/2) != n3 / 2) {n3 <- n3 + 1}
  
  # expected number of events is rounded to natural number for presentation
  d3  <- ceiling(d3)
  
  if(n2+4*n3>N){
    
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
    
  }else{
    
    pg    <-  Epgo_tte(HRgo = HRgo, d2 = d2, w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                       fixed = fixed)
    
    K2    <-  c02 + c2 * n2         # cost phase II
    K3    <-  c03 * pg + c3 * n3    # cost phase III
    
    if(K2+3*K3>K){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
    }else{
      # probability of a successful program; small, medium, large effect size
      prob1 <-  EPsProg4(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta,
                         w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                         case = case, size = "small", fixed = fixed)
      prob3 <-  EPsProg4(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta,
                         w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                         case = case, size = "large", fixed = fixed)
      prob2 <-  EPsProg4(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta,
                         w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                         case = case, size = "all", fixed = fixed) - prob1 - prob3
      
      SP    <-  prob1 + prob2 + prob3
      
      if(SP<S){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        
      }else{
        
        G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
        
        EU    <-  - K2 - 4*K3 + G
        
        return(c(EU, 4*d3, SP, pg, K2, 4*K3, prob1, prob2, prob3, n2, 4*n3))
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

# Expected probability to do third phase III trial: Epgo3
Epgo23 <-  function(HRgo, d2, alpha, beta, w, hr1, hr2, id1, id2){
  
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
                  sd = sqrt(4/d2)) * 
            prior_tte(x, w, hr1, hr2, id1, id2) 
        })
      },  - log(HRgo), Inf)$value
    })
  },  - Inf, Inf)$value
} 

# Expected number of events for third phase III when going to third phase III
Ed323 <-  function(HRgo, d2, alpha, beta, w, ymin){
  
  ( (4 * (qnorm(1 - alpha) + qnorm(1 - beta))^2)/(ymin^2)) * 
    Epgo23(HRgo, d2, alpha, beta, w, hr1, hr2, id1, id2)
  
}

# Expected probability of a successful program
EPsProg23 <-  function(HRgo, d2, alpha, beta, w, hr1, hr2, id1, id2, case, size, ymin){
  # Option 2.1: first two phase III trials are successful: no third phase III trial
  # Option 2.2: one of the two first phase III trials successful, the treatment
  #  effect of the other one shows a minimal clinically relevant effect: 
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
                        upper = c(qnorm(1 - alpha) - 
                                    log(0.85)/sqrt(y^2/c), 
                                  qnorm(1 - alpha) - 
                                    log(0.85)/sqrt(y^2/c)), 
                        mean = c(x/sqrt(y^2/c), 
                                 x/sqrt(y^2/c)), 
                        sigma = SIGMA)  - 
                  pmvnorm(lower = c(qnorm(1 - alpha) - 
                                      log(0.95)/sqrt(y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.95)/sqrt(y^2/c)), 
                          upper = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c)), 
                          mean = c(x/sqrt(y^2/c), 
                                   x/sqrt(y^2/c)), 
                          sigma = SIGMA)) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(4/d2)) * 
                prior_tte(x, w, hr1, hr2, id1, id2) 
            })
          },  - log(HRgo), Inf)$value
        })
      },  - Inf, Inf)$value) 
    }
    if(size == "large"){
      return(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            sapply(y, function(y){
              ( pmvnorm(lower = c(qnorm(1 - alpha) - 
                                    log(0.85)/sqrt(y^2/c), 
                                  qnorm(1 - alpha) - 
                                    log(0.85)/sqrt(y^2/c)), 
                        upper = c(Inf, Inf), 
                        mean = c(x/sqrt(y^2/c), x/sqrt(y^2/c)), 
                        sigma = SIGMA) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(4/d2)) * 
                prior_tte(x, w, hr1, hr2, id1, id2) 
            })
          },  - log(HRgo), Inf)$value
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
                      sd = sqrt(4/d2)) * 
                prior_tte(x, w, hr1, hr2, id1, id2) 
            })
          },  - log(HRgo), Inf)$value
        })
      },  - Inf, Inf)$value)    
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
                                        log(0.85)/sqrt(y^2/c), 
                                      qnorm(1 - alpha), 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(ymin^2/c)), 
                            mean = c(x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c), 
                                     x/sqrt(ymin^2/c)), 
                            sigma = SIGMA3)  - 
                      pmvnorm(lower = c(qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(y^2/c), 
                                        0, 
                                        qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(ymin^2/c)), 
                              upper = c(qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c), 
                                        qnorm(1 - alpha), 
                                        qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(ymin^2/c)), 
                              mean = c(x/sqrt(y^2/c), 
                                       x/sqrt(y^2/c), 
                                       x/sqrt(ymin^2/c)), 
                              sigma = SIGMA3)) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(4/d2)) * 
                prior_tte(x, w, hr1, hr2, id1, id2) 
            })
          },  - log(HRgo), Inf)$value
        })
      },  - Inf, Inf)$value) 
    }
    if(size == "large"){
      return(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            sapply(y, function(y){
              2 * ( pmvnorm(lower = c(qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(y^2/c), 
                                      0, 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt(ymin^2/c)), 
                            upper = c(Inf, 
                                      qnorm(1 - alpha), 
                                      Inf), 
                            mean = c(x/sqrt(y^2/c), 
                                     x/sqrt(y^2/c), 
                                     x/sqrt(ymin^2/c)), 
                            sigma = SIGMA3) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(4/d2)) * 
                prior_tte(x, w, hr1, hr2, id1, id2) 
            })
          },  - log(HRgo), Inf)$value
        })
      },  - Inf, Inf)$value)    
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
                      sd = sqrt(4/d2)) * 
                prior_tte(x, w, hr1, hr2, id1, id2) 
            })
          },  - log(HRgo), Inf)$value
        })
      },  - Inf, Inf)$value)    
    }
  }
  
}
# Utility function
utility23 <-  function(d2, HRgo, w, hr1, hr2, id1, id2,
                       alpha, beta, xi2, xi3,
                       c2, c3, c02, c03, 
                       b1, b2, b3){ 
  
  pg    <-  Epgo_tte(HRgo = HRgo, d2 = d2, 
                     w = w, hr1 = hr1, hr2 = hr2, 
                     id1 = id1, id2 = id2,
                     fixed = fixed)
  
  d3  <- Ed3_tte(HRgo = HRgo, d2 = d2, alpha = alpha, 
                 beta = beta, w = w, hr1 = hr1, hr2 = hr2,
                 id1 = id1, id2 = id2, fixed = fixed)

  # sample size is rounded up to next even natural number
  n2  <- ceiling(d2*(1/xi2))
  #if(round(n2/2) != n2 / 2) {n2 <- n2 + 1}
  
  n3  <- ceiling(d3 * (1/xi3))
  #if(round(n3/2) != n3 / 2) {n3 <- n3 + 1}
  
  # expected number of events is rounded to natural number
  d3  <- ceiling(d3)
  
    
  # probability of a successful program:
  # small, medium and large effect size
  prob1 <- EPsProg23(HRgo = HRgo, d2 = d2,
                       alpha = alpha, beta = beta,
                       w = w, hr1 = hr1, hr2 = hr2, 
                       id1 = id1, id2 = id2, 
                     case = 2, size = "small", ymin = ymin)
  prob3 <- EPsProg23(HRgo = HRgo, d2 = d2, 
                       alpha = alpha, beta = beta,
                       w = w, hr1 = hr1, hr2 = hr2, 
                       id1 = id1, id2 = id2, 
                     case = 2, size = "large", ymin = ymin)
  prob2 <- EPsProg23(HRgo = HRgo, d2 = d2, 
                       alpha = alpha, beta = beta,
                       w = w, hr1 = hr1, hr2 = hr2, 
                       id1 = id1, id2 = id2, 
                     case = 2, size = "all", ymin = ymin) - prob1 -prob3
      
        
  # prob to do third phase III trial
  pg3   <-  Epgo23(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, 
                   w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2) 
  
  # n3 for third phase III trial
  d33   <-  (4 * (qnorm(1 - alpha) + qnorm(1 - beta))^2)/(ymin^2) 
  
  n33  <- ceiling(d33 * (1/xi3))
  #if(round(n33/2) != n33 / 2) {n33 <- n33 + 1}
  
  d33  <- ceiling(d33)
  
  # probability of a successful program: effect sizes, 
  # for program with third phase III trial
  
  # small 
  prob13   <-  EPsProg23(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, 
                         w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                         case = 3, size = "small", ymin = ymin) 
  # large
  prob33   <-  EPsProg23(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta,
                         w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                         case = 3, size = "large", ymin = ymin) 
  # medium
  prob23   <-  EPsProg23(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, 
                         w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                         case = 3, size = "all", ymin = ymin)  -  prob1  -  prob3
  
  K1    <-  c02 + c2 * n2 # cost phase II 
  
  # cost for one of the first two phase III trials in case of go decision
  K2    <-  c03 * pg + c3 * n3
  
  # cost for the third phase III trial in case of third phase III trial
  K3    <-  pg3 * (c03 + c3 * n33)
  
  G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 +
    b1 * prob13 + b2 * prob23 + b3 * prob33 # gain
  
  EU    <-  - K1 - 2 * K2 - K3 + G
  SP    <-  prob1 + prob2 + prob3 +  
    prob13 + prob23 + prob33
  
  return(
    c(EU, d3, SP, pg, K2, K3, prob1, prob2, prob3, n2, n3, pg3, d33, n33, prob13, prob23, prob33 )
  )
  
}





