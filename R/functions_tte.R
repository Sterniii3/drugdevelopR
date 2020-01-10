# prior distribution for theta
prior_tte<-function(x, w, hr1, hr2, id1, id2){
    w * dnorm(x, -log(hr1), sqrt(4/id1)) + 
    (1 - w) * dnorm(x, -log(hr2), sqrt(4/id2))
}

# 10000 realizations of the prior distribution
box_tte<-function(w, hr1, hr2, id1, id2){
  w * rnorm(1000000, -log(hr1),sqrt(4/id1)) + 
    (1 - w) * rnorm(1000000, -log(hr2), sqrt(4/id2))
}

# expected probability to go to phase III
Epgo_tte <-  function(HRgo, d2, w, hr1, hr2, id1, id2, fixed){
  if(!fixed){
    return(  
      integrate(function(x){
        sapply(x, function(x){
          pnorm((log(HRgo) + x)/sqrt(4/d2))*
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

# expected number of events for phase III 
# in before phase II perspective
Ed3_tte <-  function(HRgo, d2, alpha, beta, 
                     w, hr1, hr2, id1, id2, fixed){
  if(!fixed){
    return(  
      integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            ((4*(qnorm(1-alpha)+qnorm(1-beta))^2)/(y^2))*
              dnorm(y,
                  mean = x,
                  sd = sqrt(4/d2))*
            prior_tte(x, w, hr1, hr2, id1, id2)
          }, -log(HRgo), Inf)$value
        })
      },  - Inf, Inf)$value
    )
  }else{
    return(  
      integrate(function(y){
        ((4*(qnorm(1-alpha)+qnorm(1-beta))^2)/(y^2))*
            dnorm(y,
                  mean = -log(hr1),
                  sd = sqrt(4/d2)) 
      }, -log(HRgo), Inf)$value
    )
  }
}

# expected probability of a successful program
EPsProg_tte <-  function(HRgo, d2, alpha, beta, 
                         step1, step2, 
                         w, hr1, hr2, id1, id2, 
                         gamma, fixed){

  c = (qnorm(1 - alpha) + qnorm(1 - beta))^2

  if(!fixed){
    return(
      integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            (pnorm(qnorm(1-alpha)-log(step2)/(sqrt(y^2/c)),
                    mean = (x+gamma)/(sqrt(y^2/c)),
                    sd = 1) -
                pnorm(qnorm(1-alpha)-log(step1)/(sqrt(y^2/c)),
                      mean = (x+gamma)/(sqrt(y^2/c)),
                      sd = 1) )*
                  dnorm(y,
                      mean = x,
                      sd = sqrt(4/d2))*
                    prior_tte(x, w, hr1, hr2, id1, id2)
          }, -log(HRgo), Inf)$value
        })
      },  - Inf, Inf)$value
    )
  }else{
    return(
      integrate(function(y){
        (pnorm(qnorm(1-alpha)-log(step2)/(sqrt(y^2/c)),
                mean = (-log(hr1)+gamma)/(sqrt(y^2/c)),
                sd = 1) -
            pnorm(qnorm(1-alpha)-log(step1)/(sqrt(y^2/c)),
                  mean = (-log(hr1)+gamma)/(sqrt(y^2/c)),
                  sd = 1))*
              dnorm(y,
                mean = -log(hr1),
                sd = sqrt(4/d2)) 
      },  - log(HRgo), Inf)$value
    )
  }
}

# utility function
utility_tte <-  function(d2, HRgo, w, hr1, hr2, id1, id2,
                         alpha, beta, xi2, xi3,
                         c2, c3, c02, c03, 
                         K, N, S,
                         steps1, stepm1, stepl1,
                         b1, b2, b3,
                         gamma, fixed){

  d3  <- Ed3_tte(HRgo = HRgo, d2 = d2, alpha = alpha, 
                 beta = beta, w = w, hr1 = hr1, hr2 = hr2,
                 id1 = id1, id2 = id2, fixed = fixed)

  # sample size is rounded up to next even natural number
  n2  <- ceiling(d2*(1/xi2))
  if(round(n2/2) != n2 / 2) {n2 <- n2 + 1}

  n3  <- ceiling(d3 * (1/xi3))
  if(round(n3/2) != n3 / 2) {n3 <- n3 + 1}

  # expected number of events is rounded to natural number
  d3  <- ceiling(d3)
  
  if(n2+n3>N){
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, 
             -9999, -9999, -9999, -9999, -9999))
  }else{
    pg    <-  Epgo_tte(HRgo = HRgo, d2 = d2, 
                       w = w, hr1 = hr1, hr2 = hr2, 
                       id1 = id1, id2 = id2,
                       fixed = fixed)
    
    K2    <-  c02 + c2 * n2         # cost phase II
    K3    <-  c03 * pg + c3 * n3    # cost phase III
    
    if(K2+K3>K){
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, 
               -9999, -9999, -9999, -9999, -9999))
    }else{
      # probability of a successful program:
      # small, medium and large effect size
      prob1 <- EPsProg_tte(HRgo = HRgo, d2 = d2,
                           alpha = alpha, beta = beta,
                           step1 = steps1, step2 =  steps2,
                           w = w, hr1 = hr1, hr2 = hr2, 
                           id1 = id1, id2 = id2, 
                           gamma = gamma, fixed = fixed)
      prob2 <- EPsProg_tte(HRgo = HRgo, d2 = d2, 
                           alpha = alpha, beta = beta,
                           step1 =  stepm1, step2 =  stepm2,
                           w = w, hr1 = hr1, hr2 = hr2, 
                           id1 = id1, id2 = id2, 
                           gamma = gamma, fixed = fixed)
      prob3 <- EPsProg_tte(HRgo = HRgo, d2 = d2, 
                           alpha = alpha, beta = beta,
                           step1 =  stepl1, step2 = stepl2,
                           w = w, hr1 = hr1, hr2 = hr2, 
                           id1 = id1, id2 = id2, 
                           gamma = gamma, fixed = fixed)
      
      SP    <-  prob1 + prob2 + prob3
      
      if(SP<S){
        return(c(-9999, -9999, -9999, -9999, -9999, 
                 -9999, -9999, -9999, -9999, -9999, -9999))
      }else{
        G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 
        EU    <-  - K2 - K3 + G

        return(
          c(EU, d3, SP, pg, K2, K3, prob1, prob2, prob3, n2, n3)
        )
      } 
    } 
  }
}

#################
# skip phase II #
#################

# number of events for phase III based on median_prior
d3_skipII_tte <-function(alpha, beta, median_prior){
  return(
    (4*(qnorm(1-alpha)+qnorm(1-beta))^2)/(median_prior^2)
  )
}

# expected probability of a successful program 
# based on median_prior
EPsProg_skipII_tte <-function(alpha, beta, step1, step2, 
                              median_prior, w, hr1, hr2, 
                              id1, id2, gamma, fixed){

  c=(qnorm(1-alpha)+qnorm(1-beta))^2

  if(!fixed){
    return(
      integrate(function(x){
        sapply(x,function(x){
          (pnorm(qnorm(1-alpha)-
                   log(step2)/(sqrt(median_prior^2/c)),
                  mean=(x+gamma)/(sqrt(median_prior^2/c)),
                  sd=1)-
              pnorm(qnorm(1-alpha)-
                      log(step1)/(sqrt(median_prior^2/c)),
                    mean=(x+gamma)/(sqrt(median_prior^2/c)),
                    sd=1))*
                prior_tte(x, w, hr1, hr2, id1, id2)
        })
      }, -Inf, Inf)$value
    )  
  }else{
    return(
      pnorm(qnorm(1-alpha)- 
              log(step2)/(sqrt(median_prior^2/c)),
            mean=(-log(hr1)+gamma)/(sqrt(median_prior^2/c)),
            sd=1)-
        pnorm(qnorm(1-alpha)- 
                log(step1)/(sqrt(median_prior^2/c)),
              mean=(-log(hr1)+gamma)/(sqrt(median_prior^2/c)),
              sd=1) 
    )
  }
}

#utility function
utility_skipII_tte <-function(alpha, beta, xi3, c03, c3, 
                              b1, b2, b3, median_prior,
                              K, N, S, 
                              steps1, steps2, 
                              stepm1, stepm2, 
                              stepl1, stepl2,
                              w, hr1, hr2, id1, id2, 
                              gamma, fixed){

  d3  <- d3_skipII_tte(alpha = alpha, beta = beta, 
                       median_prior = median_prior)

  n3  <- ceiling(d3*(1/xi3))
  if(round(n3/2) != n3 / 2) {n3 = n3 + 1}

  d3  <- ceiling(d3)
  
  if(n3>N){
    return(c(-9999, -9999, -9999, -9999,
             -9999, -9999, -9999, -9999))
  }else{
    K2  <- 0 
    K3  <- c03 + c3*n3
    if(K2+K3>K){
      return(c(-9999, -9999, -9999, -9999,
               -9999, -9999, -9999, -9999))
    }else{

      # probability of a successful program:
      # small, medium, large effect size
      prob1 <- EPsProg_skipII_tte(alpha = alpha, beta = beta, 
                                  step1 = steps1, 
                                  step2 = steps2,
                                  median_prior = median_prior, 
                                  w = w, hr1 = hr1, hr2 = hr2, 
                                  id1 = id1, id2 = id2,
                                  gamma = gamma, fixed = fixed)
      prob2 <- EPsProg_skipII_tte(alpha = alpha, beta = beta, 
                                  step1 = stepm1, 
                                  step2 = stepm2,
                                  median_prior = median_prior, 
                                  w = w, hr1 = hr1, hr2 = hr2, 
                                  id1 = id1, id2 = id2,
                                  gamma = gamma, fixed = fixed)
      prob3 <- EPsProg_skipII_tte(alpha = alpha, beta = beta, 
                                  step1 = stepl1, 
                                  step2 = stepl2,
                                  median_prior = median_prior, 
                                  w = w, hr1 = hr1, hr2 = hr2, 
                                  id1 = id1, id2 = id2,
                                  gamma = gamma, fixed = fixed)
      
      SP    <-  prob1 + prob2 + prob3
      
      if(SP<S){
        return(c(-9999, -9999, -9999, -9999,
                 -9999, -9999, -9999, -9999))
        
      }else{
        G     <- b1 * prob1 + b2 * prob2 + b3 * prob3 
        EU    <-  - K2 - K3 + G
        
        return(
          c(EU, d3, n3, SP, K3, prob1, prob2, prob3)
        )
      }
    }
  }
}

