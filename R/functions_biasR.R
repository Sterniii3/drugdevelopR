# Prior distribution for theta
prior_tte<-function(x, w, hr1, hr2, id1, id2){
    w * dnorm(x, -log(hr1), sqrt(4/id1)) + (1 - w) * dnorm(x, -log(hr2), sqrt(4/id2))
}

# Expected probability to go to phase III: Epgo
Epgo_tte <-  function(HRgo, d2, w, hr1, hr2, id1, id2, fixed){
  
  if(!fixed){
    
    return(  
      integrate(function(x){
        sapply(x, function(x){
          pnorm((log(HRgo) + x)/sqrt(4/d2)) *
            prior_tte(x, w, hr1, hr2, id1, id2)
        })
      },  - Inf, Inf)$value
    ) 
    
  }else{
    
    return(
      pnorm((log(HRgo) - log(hr1))/sqrt(4/d2))
      )
  }

}

# 1.1. conservative sample size calculation: use lower bound of one-sided confidence intervall
##############################################################################################

# prior distribution
# as above

# expected probability to go to phase III
# as above

# expected number of events for phase III when going to phase III
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
    return(ceiling(integrated))
    
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
    return(ceiling(integrated))
  }

} 

# expected probability of a successful program
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

# Utility function
utility_L <-  function(d2, HRgo, Adj, w, hr1, hr2, id1, id2,
                         alpha, beta, xi2, xi3,
                         c2, c3, c02, c03, 
                         K, N, S,
                         steps1, stepm1, stepl1,
                         b1, b2, b3,
                         fixed){

  d3    <-  Ed3_L(HRgo = HRgo, d2 = d2, Adj = Adj,
                  alpha = alpha, beta = beta, w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                  fixed = fixed)

  if(is.na(d3)){
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))  
  }else{
    # round up to next even natural number
    n2 = ceiling(d2 * (1/xi2))
    if(round(n2/2) != n2 / 2) {n2 = n2 + 1}
    
    n3 = ceiling(d3 * (1/xi3))
    if(round(n3/2) != n3 / 2) {n3 = n3 + 1}
    
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

# expected probability to go to phase III
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
EPsProg_L2<-function(HRgo, d2, Adj, alpha, beta, step1, step2, w, hr1, hr2, id1, id2, fixed){
  
  c=(qnorm(1-alpha)+qnorm(1-beta))^2
  
  if(fixed){
    return(

          integrate(function(y){ 
            
            ( pnorm(qnorm(1-alpha)+step2/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
                    mean=-log(hr1)/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
                    sd=1)-
                pnorm(qnorm(1-alpha)+step1/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
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

              ( pnorm(qnorm(1-alpha)+step2/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
                      mean=x/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
                      sd=1)-
                  pnorm(qnorm(1-alpha)+step1/(sqrt((y-qnorm(1-Adj)*sqrt(4/d2))^2/c)),
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
utility_L2 <-  function(d2, HRgo, Adj, w, hr1, hr2, id1, id2,
                       alpha, beta, xi2, xi3,
                       c2, c3, c02, c03, 
                       K, N, S,
                       steps1, stepm1, stepl1,
                       b1, b2, b3,
                       fixed){
  
  d3    <-  Ed3_L2(HRgo = HRgo, d2 = d2, Adj = Adj,
                  alpha = alpha, beta = beta, w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                  fixed = fixed)
  if(is.na(d3)){
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))  
  }else{
    # round up to next even natural number
    n2 = ceiling(d2 * (1/xi2))
    if(round(n2/2) != n2 / 2) {n2 = n2 + 1}
    
    n3 = ceiling(d3 * (1/xi3))
    if(round(n3/2) != n3 / 2) {n3 = n3 + 1}
    
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
    return(ceiling(integrated))
    
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
    return(ceiling(integrated))
  }
  
} 

# expected probability of a successful program
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
utility_R <-  function(d2, HRgo, Adj, w, hr1, hr2, id1, id2,
                       alpha, beta, xi2, xi3,
                       c2, c3, c02, c03, 
                       K, N, S,
                       steps1, stepm1, stepl1,
                       b1, b2, b3,
                       fixed){
  
  d3    <-  Ed3_R(HRgo = HRgo, d2 = d2, Adj = Adj,
                  alpha = alpha, beta = beta, w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                  fixed = fixed)
  
  if(is.na(d3)){
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))  
  }else{
  
    # round up to next even natural number
    n2 = ceiling(d2 * (1/xi2))
    if(round(n2/2) != n2 / 2) {n2 = n2 + 1}
    
    n3 = ceiling(d3 * (1/xi3))
    if(round(n3/2) != n3 / 2) {n3 = n3 + 1}
    
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
EPsProg_R2<-function(HRgo, d2, Adj, alpha, beta, step1, step2, w, hr1, hr2, id1, id2, fixed){
  
  c=(qnorm(1-alpha)+qnorm(1-beta))^2
  
  if(fixed){
      return(
        integrate(function(y){ 
          ( pnorm(qnorm(1-alpha)+step2/(sqrt((y*Adj)^2/c)),
                  mean=-log(hr1)/(sqrt((y*Adj)^2/c)),
                  sd=1)-
              pnorm(qnorm(1-alpha)+step1/(sqrt((y*Adj)^2/c)),
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
              ( pnorm(qnorm(1-alpha)+step2/(sqrt((y*Adj)^2/c)),
                      mean=x/(sqrt((y*Adj)^2/c)),
                      sd=1)-
                  pnorm(qnorm(1-alpha)+step1/(sqrt((y*Adj)^2/c)),
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
utility_R2 <-  function(d2, HRgo, Adj, w, hr1, hr2, id1, id2,
                       alpha, beta, xi2, xi3,
                       c2, c3, c02, c03, 
                       K, N, S,
                       steps1, stepm1, stepl1,
                       b1, b2, b3,
                       fixed){
  
  d3    <-  Ed3_R2(HRgo = HRgo, d2 = d2, Adj = Adj,
                  alpha = alpha, beta = beta, w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                  fixed = fixed)
  
  if(is.na(d3)){
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))  
  }else{
  
    # round up to next even natural number
    n2 = ceiling(d2 * (1/xi2))
    if(round(n2/2) != n2 / 2) {n2 = n2 + 1}
    
    n3 = ceiling(d3 * (1/xi3))
    if(round(n3/2) != n3 / 2) {n3 = n3 + 1}
    
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
