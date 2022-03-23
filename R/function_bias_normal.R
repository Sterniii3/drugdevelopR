# 1.1. conservative sample size calculation: use lower bound of one-sided confidence intervall
##############################################################################################

# prior distribution
# as above

# expected probability to go to phase III
# as above

# Expected sample size for phase III when going to phase III: En3
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
# Expected probability of a successful program: EsP
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
                    mean = (x)/sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c),
                    sd = 1) -
                pnorm(qnorm(1 - alpha) + step1/sqrt((y-qnorm(1-Adj)*sqrt(4/n2))^2/c),
                      mean = (x)/sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c),
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

# Utility function
utility_normal_L <-  function(n2, kappa, Adj, w, Delta1, Delta2, in1, in2, a, b,
                            alpha, beta, 
                            c2, c3, c02, c03, 
                            K, N, S,
                            steps1, stepm1, stepl1,
                            b1, b2, b3,
                            fixed){
  
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

# Expected probability to go to phase III: Epgo
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


# Expected sample size for phase III when going to phase III: En3
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
    return(
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

# Expected probability of a successful program: EsP
EPsProg_normal_L2 <-  function(kappa, n2, alpha, beta, step1, step2, w, Delta1, Delta2, in1, in2, a, b, fixed){
  
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

# Utility function
utility_normal_L2 <-  function(n2, kappa, Adj, w, Delta1, Delta2, in1, in2, a, b,
                            alpha, beta, 
                            c2, c3, c02, c03, 
                            K, N, S,
                            steps1, stepm1, stepl1,
                            b1, b2, b3,
                            fixed){
  
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
    
    pg    <-  Epgo_normal_L2(kappa = kappa, n2 = n2,
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

# 2.1. conservative sample size calculation: use estimate with retetion factor
##############################################################################################

# prior distribution
# as above

# expected probability to go to phase III
# as above

# Expected sample size for phase III when going to phase III: En3
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

# Expected probability of a successful program: EsP
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
                    mean = (x)/sqrt(y^2/c),
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

# Utility function
utility_normal_R <-  function(n2, kappa, Adj, w, Delta1, Delta2, in1, in2, a, b,
                            alpha, beta, 
                            c2, c3, c02, c03, 
                            K, N, S,
                            steps1, stepm1, stepl1,
                            b1, b2, b3,
                            fixed){
  
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
# use estimate with retetion factor
##############################################################################################

# prior distribution
# as above

# Expected probability to go to phase III: Epgo
Epgo_normal <-  function(kappa, n2, Adj, w, Delta1, Delta2, in1, in2, a, b, fixed){
  
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

# Expected sample size for phase III when going to phase III: En3
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

# Expected probability of a successful program: EsP
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
                    mean = (x)/sqrt(y^2/c),
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


# Utility function
utility_normal_R2 <-  function(n2, kappa, Adj,  w, Delta1, Delta2, in1, in2, a, b,
                            alpha, beta, 
                            c2, c3, c02, c03, 
                            K, N, S,
                            steps1, stepm1, stepl1,
                            b1, b2, b3,
                            fixed){
  
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

