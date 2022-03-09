# Prior distribution for p1
prior_binary<-function(x, w, p11, p12, in1, in2){
  w * dnorm(x, p11, sqrt(p11*(1-p11)/in1)) + (1-w) * dnorm(x, p12, sqrt(p12*(1-p12)/in2))
}

# 10000 realizations of the prior distribution
box_binary<-function(w, p11, p12, in1, in2){
  w*rnorm(1000000,p11,sqrt(p11*(1-p11)/in1))+(1-w)*rnorm(1000000,p12,sqrt(p12*(1-p12)/in2))
}

# auxiliary functions
t1 <- function(x, p0){((1-p0)/p0) + ((1-x)/x)}
t2 <- function(x, p0){sqrt(2*(1-((p0 + x)/2))/((p0 + x)/2))}
t3 <- function(x, p0){sqrt(((1-p0)/p0) + ((1-x)/x))}

# Expected probability to go to phase III: Epgo
Epgo_binary <-  function(RRgo, n2, p0, w, p11, p12, in1, in2, fixed){
  if(fixed){
    return(
      pnorm((-log(p11/p0) + log(RRgo))/sqrt((2/n2)*t1(p11, p0))) 
    )
  }else{
    return(
      integrate(function(x){
        sapply(x, function(x){
          pnorm((-log(x/p0) + log(RRgo))/sqrt((2/n2)*t1(x, p0)))  *
            prior_binary(x, w, p11, p12, in1, in2)
        })
      },  0, 1)$value
    )
  }
}

# Expected sample size for phase III when going to phase III: En3
En3_binary <-  function(RRgo, n2, alpha, beta, p0, w, p11, p12, in1, in2, fixed){
  if(fixed){
    return(
      integrate(function(y){
        ((2*(qnorm(1-alpha)*t2(p11, p0)+qnorm(1-beta)*t3(p11, p0))^2)/y^2) *
          dnorm(y,
                mean = -log(p11/p0),
                sd = sqrt((2/n2)*t1(p11, p0)))
      }, - log(RRgo), Inf)$value
    )
  }else{
    return(
      integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            ((2*(qnorm(1-alpha)*t2(x, p0)+qnorm(1-beta)*t3(x, p0))^2)/y^2) *
              dnorm(y,
                    mean = -log(x/p0),
                    sd = sqrt((2/n2)*t1(x, p0)))*
              prior_binary(x, w, p11, p12, in1, in2)
          }, - log(RRgo), Inf)$value
        })
      }, 0, 1)$value
    )
  }
}

# Expected probability of a successful program: EsP
EPsProg_binary <-  function(RRgo, n2, alpha, beta, step1, step2, p0, w, p11, p12, in1, in2, gamma, fixed){

  if(fixed){
    return(
      integrate(function(y){
        ( pnorm(qnorm(1 - alpha) -
                  log(step2)/sqrt((t1(p11, p0)*y^2)/(qnorm(1-alpha)*t2(p11, p0) +
                                                     qnorm(1-beta)*t3(p11, p0))^2),
                mean = -log((p11+gamma)/p0)/sqrt((t1(p11, p0)*y^2)/(qnorm(1-alpha)*t2(p11, p0) +
                                                                  qnorm(1-beta)*t3(p11, p0))^2),
                sd = 1) -
            pnorm(qnorm(1 - alpha) -
                    log(step1)/sqrt((t1(p11, p0)*y^2)/(qnorm(1-alpha)*t2(p11, p0) +
                                                       qnorm(1-beta)*t3(p11, p0))^2),
                  mean = -log((p11+gamma)/p0)/sqrt((t1(p11, p0)*y^2)/(qnorm(1-alpha)*t2(p11, p0) +
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
                      log(step2)/sqrt((t1(x, p0)*y^2)/(qnorm(1-alpha)*t2(x, p0) +
                                                         qnorm(1-beta)*t3(x, p0))^2),
                    mean = -log((x+gamma)/p0)/sqrt((t1(x, p0)*y^2)/(qnorm(1-alpha)*t2(x, p0) +
                                                                      qnorm(1-beta)*t3(x, p0))^2),
                    sd = 1) -
                pnorm(qnorm(1 - alpha) -
                        log(step1)/sqrt((t1(x, p0)*y^2)/(qnorm(1-alpha)*t2(x, p0) +
                                                           qnorm(1-beta)*t3(x, p0))^2),
                      mean = -log((x+gamma)/p0)/sqrt((t1(x, p0)*y^2)/(qnorm(1-alpha)*t2(x, p0) +
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
utility_binary <-  function(n2, RRgo, w, p0, p11, p12, in1, in2,
                            alpha, beta, 
                            c2, c3, c02, c03, 
                            K, N, S,
                            steps1, stepm1, stepl1,
                            b1, b2, b3,
                            gamma, fixed){

  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 0

  n3  <-  En3_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                       p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, fixed = fixed)

  n3  <- ceiling(n3)
  
  if(round(n3/2) != n3 / 2) {n3 = n3 + 1}

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
      prob1 <-  EPsProg_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                               step1 = steps1, step2 =  steps2,
                               p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                               gamma = gamma, fixed = fixed)
      prob2 <-  EPsProg_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                               step1 =  stepm1, step2 =  stepm2,
                               p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                               gamma = gamma, fixed = fixed)
      prob3 <-  EPsProg_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta,
                               step1 =  stepl1, step2 = stepl2,
                               p0 = p0, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2, 
                               gamma = gamma, fixed = fixed)
  
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

#################
# skip phase II #
#################

# number of events for phase III based on median_prior
n3_skipII_binary <-function(alpha, beta, p0, median_prior){

  median_RR = -log(median_prior/p0)

  return(((2*(qnorm(1-alpha)*sqrt(2*(1-((p0 + median_prior)/2))/((p0 + median_prior)/2))+qnorm(1-beta)*sqrt((1-p0)/p0+(1-median_prior)/median_prior))^2)/median_RR^2))

}

# expected probability of a successful program based on median_prior
EPsProg_skipII_binary <-function(alpha, beta, step1, step2, p0, median_prior, w, p11, p12, in1, in2, gamma, fixed){

  c=(qnorm(1-alpha)+qnorm(1-beta))^2

  median_RR = -log(median_prior/p0)

  if(fixed){
    return(
      pnorm(qnorm(1-alpha) - log(step2)/(sqrt(median_RR^2/c)),
              mean=-log((p11+gamma)/p0)/(sqrt(median_RR^2/c)),
              sd=1)-
          pnorm(qnorm(1-alpha) - log(step1)/(sqrt(median_RR^2/c)),
                mean=-log((p11+gamma)/p0)/(sqrt(median_RR^2/c)),
                sd=1)
    )
  }else{
    return(
      integrate(function(x){
        sapply(x,function(x){
          ( pnorm(qnorm(1-alpha) - log(step2)/(sqrt(median_RR^2/c)),
                  mean=-log((x+gamma)/p0)/(sqrt(median_RR^2/c)),
                  sd=1)-
              pnorm(qnorm(1-alpha) - log(step1)/(sqrt(median_RR^2/c)),
                    mean=-log((x+gamma)/p0)/(sqrt(median_RR^2/c)),
                    sd=1) )*
            prior_binary(x, w, p11, p12, in1, in2)
        })
      }, 0, 1)$value 
    )
  }
}

#utility function
utility_skipII_binary <-function(alpha, beta, c03, c3, b1, b2, b3, p0, median_prior, 
                                K, N, S,
                                steps1, stepm1,  stepl1, 
                                w, p11, p12, in1, in2, gamma, fixed){

  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 0
  
  n3  <- n3_skipII_binary(alpha = alpha, beta = beta, p0 = p0, median_prior = median_prior)
  
  n3  <- ceiling(n3)
  
  if(round(n3/2) != n3 / 2) {n3 = n3 + 1}

  if(n3>N){
    
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
    
  }else{
  
  K2  <- 0 #cost phase II
  K3  <- c03 + c3 * n3 #cost phase III

    if(K2+K3>K){
  
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999))
  
    }else{
  
      # probability of a successful program; small, medium, large effect size
      prob1 <- EPsProg_skipII_binary(alpha = alpha, beta = beta, step1 = steps1, step2 = steps2,
                              p0 = p0, median_prior = median_prior, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2,
                              gamma = gamma, fixed = fixed)
      prob2 <- EPsProg_skipII_binary(alpha = alpha, beta = beta, step1 = stepm1, step2 = stepm2,
                              p0 = p0, median_prior = median_prior, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2,
                              gamma = gamma, fixed = fixed)
      prob3 <- EPsProg_skipII_binary(alpha = alpha, beta = beta, step1 = stepl1, step2 = stepl2,
                              p0 = p0, median_prior = median_prior, w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2,
                              gamma = gamma, fixed = fixed)
  
      SP    <-  prob1 + prob2 + prob3
      
        if(SP<S){
          
          return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
          
        }else{
        
        G     <- b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
    
        EU    <-  - K2 - K3 + G
    
    
        return(c(EU, n3, SP, K3, prob1, prob2, prob3))
      }
    }
  }
}



