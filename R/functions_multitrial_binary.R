########################
# Two phase III trials #
########################
# Case 1: Strategy 1/2; at least one trial significant, the treatment effect of the other one at least showing in the same direction 
# Case 2: Strategy 2/2; both trials significant 

# Expected probability of a successful program
EPsProg2_binary <-  function(RRgo, n2, alpha, beta, step1, step2, p0, w, p11, p12, in1, in2,case, size, fixed){
  
  SIGMA <-  diag(2)
  c     <-  (qnorm(1-alpha)*sqrt(2*(1-((p0 + p11)/2))/((p0 + p11)/2)) + qnorm(1-beta)*sqrt(((1-p0)/p0) + ((1-p11)/p11)))^2
  
  if(fixed){
    
    if(case == 1){
      if(size == "small"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( pmvnorm(lower = c(0, 
                                0), 
                      upper = c(qnorm(1 - alpha) - log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c),
                                qnorm(1 - alpha) - log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                      mean = c((-log((p11)/p0))/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                               -log((p11)/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                      sigma = SIGMA) - 
                pmvnorm(lower = c(0, 
                                  0), 
                        upper = c(qnorm(1 - alpha), 
                                  qnorm(1 - alpha)), 
                        mean = c((-log((p11)/p0))/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                 -log((p11)/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                        sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = -log(p11/p0), 
                    sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11)))) 
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
                      mean = c((-log((p11)/p0))/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                               (-log((p11)/p0))/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                      sigma = SIGMA) - 
                pmvnorm(lower = c(0, 
                                  0), 
                        upper = c(qnorm(1 - alpha) - 
                                    log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                  qnorm(1 - alpha) - 
                                    log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                        mean = c((-log((p11)/p0))/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                 (-log((p11)/p0))/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                        sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = -log(p11/p0), 
                    sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11)))) 
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
                      mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                               -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                      sigma = SIGMA) - 
                pmvnorm(lower = c(0, 
                                  0), 
                        upper = c(qnorm(1 - alpha), 
                                  qnorm(1 - alpha)), 
                        mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                 -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                        sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = -log(p11/p0), 
                    sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11))))
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
                                  log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                qnorm(1 - alpha) - 
                                  log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                      mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                               -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                      sigma = SIGMA) - 
                pmvnorm(lower = c(qnorm(1 - alpha) - 
                                    log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                  qnorm(1 - alpha) - 
                                    log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                        upper = c(qnorm(1 - alpha) - 
                                    log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                  qnorm(1 - alpha) - 
                                    log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                        mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                 -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                        sigma = SIGMA)) * 
              dnorm(y, 
                    mean = -log(p11/p0), 
                    sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11)))) 
          })
        },  - log(RRgo), Inf)$value) 
      }
      if(size == "large"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( pmvnorm(lower = c(qnorm(1 - alpha) - 
                                  log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                qnorm(1 - alpha) - 
                                  log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                      upper = c(Inf, 
                                Inf), 
                      mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                               -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                      sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = -log(p11/p0), 
                    sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11)))) 
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
                      mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                               -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                      sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = -log(p11/p0), 
                    sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11)))) 
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
                          upper = c(qnorm(1 - alpha) - log(0.95)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c),
                                    qnorm(1 - alpha) - log(0.95)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                          mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                   x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                          sigma = SIGMA) - 
                    pmvnorm(lower = c(0, 
                                      0), 
                            upper = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha)), 
                            mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                     x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                            sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-x)/x)))) * 
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
                          mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                   x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                          sigma = SIGMA) - 
                    pmvnorm(lower = c(0, 
                                      0), 
                            upper = c(qnorm(1 - alpha) - 
                                        log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                            mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                     x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                            sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-x)/x)))) * 
                  prior_binary(x, w, p11, p12, in1, in2)  
              })
            },  - log(RRgo), Inf)$value  
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
                          mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                   x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                          sigma = SIGMA) - 
                    pmvnorm(lower = c(0, 
                                      0), 
                            upper = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha)), 
                            mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                     x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                            sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-x)/x)))) * 
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
                                      log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                          mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                   x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                          sigma = SIGMA) - 
                    pmvnorm(lower = c(qnorm(1 - alpha) - 
                                        log(0.95)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.95)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                            upper = c(qnorm(1 - alpha) - 
                                        log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                            mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                     x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                            sigma = SIGMA)) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-x)/x)))) * 
                  prior_binary(x, w, p11, p12, in1, in2) 
              })
            },  - log(RRgo), Inf)$value
          })
        },  - Inf, Inf)$value) 
      }
      if(size == "large"){
        return(integrate(function(x){
          sapply(x, function(x){
            integrate(function(y){
              sapply(y, function(y){
                ( pmvnorm(lower = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                          upper = c(Inf, 
                                    Inf), 
                          mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                   x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                          sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-x)/x)))) * 
                  prior_binary(x, w, p11, p12, in1, in2) 
              })
            },  - log(RRgo), Inf)$value
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
                          mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                   x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                          sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-x)/x)))) * 
                  prior_binary(x, w, p11, p12, in1, in2) 
              })
            },  - log(RRgo), Inf)$value
          })
        },  - Inf, Inf)$value)    
      }
    }
    
  }
  
  
  
}

# Utility function
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
                               step1 = steps1, step2 =  steps2,
                               p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                               case = case, size = "small", fixed = fixed)
      prob3 <-  EPsProg2_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                               step1 =  stepm1, step2 =  stepm2,
                               p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                               case = case, size = "large", fixed = fixed)
      prob2 <-  EPsProg2_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                               step1 =  stepl1, step2 = stepl2,
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

# Expected probability of a successful program
EPsProg3_binary <-  function(RRgo, n2, alpha, beta, step1, step2, p0, w, p11, p12, in1, in2,case, size, fixed){
  
  SIGMA <-  diag(3)
  c     <-  (qnorm(1-alpha)*sqrt(2*(1-((p0 + p11)/2))/((p0 + p11)/2)) + qnorm(1-beta)*sqrt(((1-p0)/p0) + ((1-p11)/p11)))^2
  
  if(fixed){
    
    if(case == 2){
      if(size == "small"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( 3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    0), 
                          upper = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                          mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                   -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                   -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                          sigma = SIGMA) - 
                2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha), 
                                      qnorm(1 - alpha)), 
                            upper = c(qnorm(1 - alpha) - 
                                        log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                            mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                     -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                     -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                            sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = -log(p11/p0), 
                    sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11)))) 
          })
        },  - log(RRgo), Inf)$value)
      }
      if(size == "large"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( 3 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                    0), 
                          upper = c(Inf, 
                                    Inf, 
                                    Inf), 
                          mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                   -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                   -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                          sigma = SIGMA) - 
                2 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                        log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                            upper = c(Inf, 
                                      Inf, 
                                      Inf), 
                            mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                     -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                     -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                            sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = -log(p11/p0), 
                    sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11)))) 
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
                          mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                   -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                   -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                          sigma = SIGMA) - 
                2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha), 
                                      qnorm(1 - alpha)), 
                            upper = c(Inf, 
                                      Inf, 
                                      Inf), 
                            mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                     -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                     -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                            sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = -log(p11/p0), 
                    sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11))))  
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
                                      log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                          mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                   -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                   -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                          sigma = SIGMA) - 
                2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha), 
                                      qnorm(1 - alpha)), 
                            upper = c(qnorm(1 - alpha) - 
                                        log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                            mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                     -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                     -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                            sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = -log(p11/p0), 
                    sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11))))
          })
        },  - log(RRgo), Inf)$value)
      }
      if(size == "large"){
        return(integrate(function(y){
          sapply(y, function(y){
            ( 3 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                    qnorm(1 - alpha)), 
                          upper = c(Inf, 
                                    Inf, 
                                    Inf), 
                          mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                   -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                   -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                          sigma = SIGMA) - 
                2 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                        log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                            upper = c(Inf, 
                                      Inf, 
                                      Inf), 
                            mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                     -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                     -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                            sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = -log(p11/p0), 
                    sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11)))) 
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
                      mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                               -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                               -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                      sigma = SIGMA) ) * 
              dnorm(y, 
                    mean = -log(p11/p0), 
                    sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11)))) 
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
                                          log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.95)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.95)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                              mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                       x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                       x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                              sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                          qnorm(1 - alpha), 
                                          qnorm(1 - alpha)), 
                                upper = c(qnorm(1 - alpha) - 
                                            log(0.95)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                          qnorm(1 - alpha) - 
                                            log(0.95)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                          qnorm(1 - alpha) - 
                                            log(0.95)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                                mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                         x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                         x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11)))) * 
                  prior_binary(x, w, p11, p12, in1, in2)
              })
            },  - log(RRgo), Inf)$value   
          })
        },  - Inf, Inf)$value)
      }
      if(size == "large"){
        return(integrate(function(x){
          sapply(x, function(x){
            integrate(function(y){
              sapply(y, function(y){
                ( 3 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                          log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                        0), 
                              upper = c(Inf, 
                                        Inf, 
                                        Inf), 
                              mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                       x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                       x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                              sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                            log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                          qnorm(1 - alpha) - 
                                            log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                          qnorm(1 - alpha) - 
                                            log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                                upper = c(Inf, 
                                          Inf, 
                                          Inf), 
                                mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                         x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                         x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-x)/x)))) * 
                  prior_binary(x, w, p11, p12, in1, in2) 
              })
            },  - log(RRgo), Inf)$value
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
                              mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                       x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                       x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                              sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                          qnorm(1 - alpha), 
                                          qnorm(1 - alpha)), 
                                upper = c(Inf, 
                                          Inf, 
                                          Inf), 
                                mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                         x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                         x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-x)/x)))) * 
                  prior_binary(x, w, p11, p12, in1, in2) 
              })
            },  - log(RRgo), Inf)$value
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
                                          log(0.95)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.95)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                              mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                       x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                       x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                              sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                          qnorm(1 - alpha), 
                                          qnorm(1 - alpha)), 
                                upper = c(qnorm(1 - alpha) - 
                                            log(0.95)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                          qnorm(1 - alpha) - 
                                            log(0.95)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                          qnorm(1 - alpha) - 
                                            log(0.95)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                                mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                         x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                         x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-x)/x)))) * 
                  prior_binary(x, w, p11, p12, in1, in2) 
              })
            },  - log(RRgo), Inf)$value
          })
        },  - Inf, Inf)$value)
      }
      if(size == "large"){
        return(integrate(function(x){
          sapply(x, function(x){
            integrate(function(y){
              sapply(y, function(y){
                ( 3 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                          log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                        qnorm(1 - alpha)), 
                              upper = c(Inf, 
                                        Inf, 
                                        Inf), 
                              mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                       x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                       x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                              sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                            log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                          qnorm(1 - alpha) - 
                                            log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                          qnorm(1 - alpha) - 
                                            log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                                upper = c(Inf, 
                                          Inf, 
                                          Inf), 
                                mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                         x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                         x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-x)/x)))) * 
                  prior_binary(x, w, p11, p12, in1, in2) 
              })
            },  - log(RRgo), Inf)$value
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
                          mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                   x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                   x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                          sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = x, 
                        sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-x)/x)))) * 
                  prior_binary(x, w, p11, p12, in1, in2) 
              })
            },  - log(RRgo), Inf)$value
          })
        },  - Inf, Inf)$value)
      }
    }  
    
  }
  
}


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
                               step1 = steps1, step2 =  steps2,
                               p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                               case = case, size = "small", fixed = fixed)
      prob3 <-  EPsProg3_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                               step1 =  stepm1, step2 =  stepm2,
                               p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                               case = case, size = "large", fixed = fixed)
      prob2 <-  EPsProg3_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                               step1 =  stepl1, step2 = stepl2,
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

# Expected probability of a successful program
EPsProg4_binary <-  function(RRgo, n2, alpha, beta, step1, step2, p0, w, p11, p12, in1, in2,case, size, fixed){
  
  SIGMA <-  diag(4)
  c     <-  (qnorm(1-alpha)*sqrt(2*(1-((p0 + p11)/2))/((p0 + p11)/2)) + qnorm(1-beta)*sqrt(((1-p0)/p0) + ((1-p11)/p11)))^2
  
  if(fixed){
    
    if(size == "small"){
      return(integrate(function(y){
        sapply(y, function(y){
          ( 4 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                  qnorm(1 - alpha), 
                                  qnorm(1 - alpha), 
                                  0), 
                        upper = c(qnorm(1 - alpha) - log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                  qnorm(1 - alpha) - log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                  qnorm(1 - alpha) - log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                  qnorm(1 - alpha) - log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                        mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                 -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                 -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                 -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                        sigma = SIGMA)  - 
              3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    qnorm(1 - alpha)), 
                          upper = c(qnorm(1 - alpha) - 
                                      log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                          mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                   -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                   -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                   -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                          sigma = SIGMA) ) * 
            dnorm(y, 
                  mean = -log(p11/p0), 
                  sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11)))) 
        })
      },  - log(RRgo), Inf)$value)
    }
    if(size == "large"){
      return(integrate(function(y){
        sapply(y, function(y){
          ( 4 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                    log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                  qnorm(1 - alpha) - log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                  qnorm(1 - alpha) - log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                  0), 
                        upper = c(Inf, 
                                  Inf, 
                                  Inf, 
                                  Inf), 
                        mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                 -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                 -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                 -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                        sigma = SIGMA)  - 
              3 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                          upper = c(Inf, 
                                    Inf, 
                                    Inf, 
                                    Inf), 
                          mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                   -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                   -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                   -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                          sigma = SIGMA) ) * 
            dnorm(y, 
                  mean = -log(p11/p0), 
                  sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11)))) 
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
                        mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                 -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                 -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                 -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                        sigma = SIGMA) - 
              3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    qnorm(1 - alpha)), 
                          upper = c(Inf, 
                                    Inf, 
                                    Inf, 
                                    Inf), 
                          mean = c(-log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                   -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                   -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                   -log(p11/p0)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                          sigma = SIGMA) ) * 
            dnorm(y, 
                  mean = -log(p11/p0), 
                  sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11)))) 
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
                            upper = c(qnorm(1 - alpha) - log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                      qnorm(1 - alpha) - log(0.95)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                      qnorm(1 - alpha) - log(0.95)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                      qnorm(1 - alpha) - log(0.95)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                            mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                     x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                     x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                     x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                            sigma = SIGMA)  - 
                  3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                        qnorm(1 - alpha), 
                                        qnorm(1 - alpha), 
                                        qnorm(1 - alpha)), 
                              upper = c(qnorm(1 - alpha) - 
                                          log(0.95)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.95)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.95)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.95)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                              mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                       x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                       x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                       x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                              sigma = SIGMA) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-x)/x)))) * 
                prior_binary(x, w, p11, p12, in1, in2)
            })
          },  - log(RRgo), Inf)$value   
        })
      },  - Inf, Inf)$value)
    }
    if(size == "large"){
      return(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            sapply(y, function(y){
              ( 4 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                        log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                      qnorm(1 - alpha) - log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                      qnorm(1 - alpha) - log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                      0), 
                            upper = c(Inf, 
                                      Inf, 
                                      Inf, 
                                      Inf), 
                            mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                     x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                     x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                     x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                            sigma = SIGMA)  - 
                  3 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                          log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                        qnorm(1 - alpha) - 
                                          log(0.85)/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                              upper = c(Inf, 
                                        Inf, 
                                        Inf, 
                                        Inf), 
                              mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                       x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                       x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                       x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                              sigma = SIGMA) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-x)/x)))) * 
                prior_binary(x, w, p11, p12, in1, in2) 
            })
          },  - log(RRgo), Inf)$value
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
                            mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                     x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                     x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                     x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                            sigma = SIGMA) - 
                  3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                        qnorm(1 - alpha), 
                                        qnorm(1 - alpha), 
                                        qnorm(1 - alpha)), 
                              upper = c(Inf, 
                                        Inf, 
                                        Inf, 
                                        Inf), 
                              mean = c(x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                       x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                       x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c), 
                                       x/sqrt((((1-p0)/p0) + ((1-x)/x))*y^2/c)), 
                              sigma = SIGMA) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-x)/x)))) * 
                prior_binary(x, w, p11, p12, in1, in2) 
            })
          },  - log(RRgo), Inf)$value
        })
      },  - Inf, Inf)$value)
    }
    
  }
  
}




# Utility function
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
                                step1 = steps1, step2 =  steps2,
                                p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                                case = case, size = "small", fixed = fixed)
      prob3 <-  EPsProg4_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                                step1 =  stepm1, step2 =  stepm2,
                                p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                                case = case, size = "large", fixed = fixed)
      prob2 <-  EPsProg4_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                                step1 =  stepl1, step2 = stepl2,
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



EPsProg23_binary <-  function(HRgo, d2, alpha, beta, w, hr1, hr2, id1, id2, case, size, ymin){
  # Option 2.1: first two phase III trials are successful: no third phase III trial
  # Option 2.2: one of the two first phase III trials successful, the treatment
  #  effect of the other one points in the same direction: 
  #  conduct third phase III trial with N3 = N3(ymin)
  
  SIGMA <-  diag(2)
  SIGMA3<-  diag(3)
  c  <- (qnorm(1-alpha)*sqrt(2*(1-((p0 + p11)/2))/((p0 + p11)/2)) + qnorm(1-beta)*sqrt(((1-p0)/p0) + ((1-p11)/p11)))^2
  
  if(case == 2){ # Option 2.1
    if(size == "small"){
      return(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            sapply(y, function(y){
              ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                  qnorm(1 - alpha)), 
                        upper = c(qnorm(1 - alpha) - 
                                    log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                  qnorm(1 - alpha) - 
                                    log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                        mean = c(x/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                 x/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                        sigma = SIGMA)  - 
                  pmvnorm(lower = c(qnorm(1 - alpha) - 
                                      log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                          upper = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                    qnorm(1 - alpha) - 
                                      log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                          mean = c(x/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                   x/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                          sigma = SIGMA)) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11)))) * 
                prior_binary(x, w, p11, p12, in1, in2) 
            })
          },  - log(RRgo), Inf)$value
        })
      },  - Inf, Inf)$value) 
    }
    if(size == "large"){
      return(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            sapply(y, function(y){
              ( pmvnorm(lower = c(qnorm(1 - alpha) - 
                                    log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                  qnorm(1 - alpha) - 
                                    log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                        upper = c(Inf, Inf), 
                        mean = c(x/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                 x/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                        sigma = SIGMA) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11)))) * 
                prior_binary(x, w, p11, p12, in1, in2) 
            })
          },  - log(RRgo), Inf)$value
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
                        mean = c(x/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                 x/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c)), 
                        sigma = SIGMA) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11)))) * 
                prior_binary(x, w, p11, p12, in1, in2) 
            })
          },  - log(RRgo), Inf)$value
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
                                        log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                      qnorm(1 - alpha), 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*ymin^2/c)), 
                            mean = c(x/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                     x/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                     x/sqrt((((1-p0)/p0) + ((1-p11)/p11))*ymin^2/c)), 
                            sigma = SIGMA3)  - 
                      pmvnorm(lower = c(qnorm(1 - alpha) - 
                                          log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                        0, 
                                        qnorm(1 - alpha) - 
                                          log(0.95)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*ymin^2/c)), 
                              upper = c(qnorm(1 - alpha) - 
                                          log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                        qnorm(1 - alpha), 
                                        qnorm(1 - alpha) - 
                                          log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*ymin^2/c)), 
                              mean = c(x/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                       x/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                       x/sqrt((((1-p0)/p0) + ((1-p11)/p11))*ymin^2/c)), 
                              sigma = SIGMA3)) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11)))) * 
                prior_binary(x, w, p11, p12, in1, in2) 
            })
          },  - log(RRgo), Inf)$value
        })
      },  - Inf, Inf)$value) 
    }
    if(size == "large"){
      return(integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            sapply(y, function(y){
              2 * ( pmvnorm(lower = c(qnorm(1 - alpha) - 
                                        log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                      0, 
                                      qnorm(1 - alpha) - 
                                        log(0.85)/sqrt((((1-p0)/p0) + ((1-p11)/p11))*ymin^2/c)), 
                            upper = c(Inf, 
                                      qnorm(1 - alpha), 
                                      Inf), 
                            mean = c(x/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                     x/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                     x/sqrt((((1-p0)/p0) + ((1-p11)/p11))*ymin^2/c)), 
                            sigma = SIGMA3) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11)))) * 
                prior_binary(x, w, p11, p12, in1, in2) 
            })
          },  - log(RRgo), Inf)$value
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
                            mean = c(x/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                     x/sqrt((((1-p0)/p0) + ((1-p11)/p11))*y^2/c), 
                                     x/sqrt((((1-p0)/p0) + ((1-p11)/p11))*ymin^2/c)), 
                            sigma = SIGMA3) ) * 
                dnorm(y, 
                      mean = x, 
                      sd = sqrt(((2/n2))*(((1-p0)/p0) + ((1-p11)/p11)))) * 
                prior_binary(x, w, p11, p12, in1, in2) 
            })
          },  - log(RRgo), Inf)$value
        })
      },  - Inf, Inf)$value)    
    }
  }
  
}




# Utility function
utility23_binary <-  function(n2, RRgo, w, p0, p11, p12, in1, in2,
                             alpha, beta, 
                             c2, c3, c02, c03,
                             b1, b2, b3){
  
  pg    <-  Epgo_binary(RRgo = RRgo, n2 = n2, p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, fixed=fixed)
  
  
  n3  <-  En3_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                     p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, fixed = fixed)
  
  n3  <- ceiling(n3)
  
  if(round(n3/2) != n3 / 2) {n3 = n3 + 1}
  
  prob1 <-  EPsProg23_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                            step1 = steps1, step2 =  steps2,
                            p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                            case = 2, size = "small", ymin = ymin)
  prob3 <-  EPsProg23_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                            step1 =  stepm1, step2 =  stepm2,
                            p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                            case = 2, size = "large", ymin = ymin)
  prob2 <-  EPsProg23_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                            step1 =  stepl1, step2 = stepl2,
                            p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                            case = 2, size = "all", ymin = ymin) - prob1 - prob3
  
  
    
  pg3   <-  Epgo_binary(RRgo = RRgo, n2 = n2, p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2)
    
  n33   <-  (4 * (qnorm(1 - alpha) + qnorm(1 - beta))^2)/(ymin^2) 
  
  n33  <- ceiling(n33*pg3)
  
  if(round(n33/2) != n33 / 2) {n33 <- n33 + 1}
  
  prob13 <-  EPsProg23_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                             step1 = steps1, step2 =  steps2,
                             p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                             case = 3, size = "small", ymin = ymin)
  prob33 <-  EPsProg23_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                             step1 =  stepm1, step2 =  stepm2,
                             p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                             case = 3, size = "large", ymin = ymin)
  prob23 <-  EPsProg23_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                             step1 =  stepl1, step2 = stepl2,
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
