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

# Expected sample size for phase III when going to phase III: En3
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

# Expected probability of a successful program: EsP
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

# Utility function
utility_binary_L <-  function(n2, RRgo, Adj, w, p0, p11, p12, in1, in2,
                              alpha, beta, 
                              c2, c3, c02, c03, 
                              K, N, S,
                              steps1, stepm1, stepl1,
                              b1, b2, b3,
                              fixed){
  
  
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

# Expected probability to go to phase III: Epgo
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

# Expected sample size for phase III when going to phase III: En3
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

# Expected probability of a successful program: EsP
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

# Utility function
utility_binary_L2 <-  function(n2, RRgo, Adj, w, p0, p11, p12, in1, in2,
                               alpha, beta, 
                               c2, c3, c02, c03, 
                               K, N, S,
                               steps1, stepm1, stepl1,
                               b1, b2, b3,
                               fixed){
  
  
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

# Expected sample size for phase III when going to phase III: En3
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

# Expected probability of a successful program: EsP
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

# Utility function
utility_binary_R <-  function(n2, RRgo, Adj, w, p0, p11, p12, in1, in2,
                              alpha, beta, 
                              c2, c3, c02, c03, 
                              K, N, S,
                              steps1, stepm1, stepl1,
                              b1, b2, b3,
                              fixed){
  
 
  
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

# Expected probability to go to phase III: Epgo
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

# Expected sample size for phase III when going to phase III: En3
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

# Expected probability of a successful program: EsP
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

# Utility function
utility_binary_R2 <-  function(n2, RRgo, Adj, w, p0, p11, p12, in1, in2,
                               alpha, beta, 
                               c2, c3, c02, c03, 
                               K, N, S,
                               steps1, stepm1, stepl1,
                               b1, b2, b3,
                               fixed){
  
 
  
  
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
