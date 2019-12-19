# Prior distribution for Delta
prior_normal <-function(x, w, Delta1, Delta2, in1, in2, a, b){
  w * dtnorm(x, Delta2, sqrt(4/in1), lower = a, upper = b) +
    (1-w) * dtnorm(x, Delta1, sqrt(4/in2), lower = a, upper = b)
}

# 10000 realizations of the prior distribution
box_normal <-function(w, Delta1, Delta2, in1, in2, a, b){
  w*rtnorm(1000000,Delta2,sqrt(4/in1), lower = a, upper = b) +
    (1-w)*rtnorm(1000000,Delta1,sqrt(4/in2), lower = a, upper = b)
}

# Expected probability to go to phase III: Epgo
Epgo_normal <-  function(kappa, n2, w, Delta1, Delta2, in1, in2, a, b, fixed){
  
  if(fixed){
    return(
      pnorm((Delta1-kappa)/sqrt(4/n2))  
    )
  }else{
    return(
      integrate(function(x){
        sapply(x, function(x){
          pnorm((x-kappa)/sqrt(4/n2)) *
            prior_normal(x, w, Delta1, Delta2, in1, in2, a, b)
        })
      },  - Inf, Inf)$value
    )
  }
}

# Expected sample size for phase III when going to phase III: En3
En3_normal <-  function(kappa, n2, alpha, beta, w, Delta1, Delta2, in1, in2, a, b, taus, fixed){

  if(fixed){
    
    int   = try(integrate(function(y){
      sapply(y,function(y){
        ( (4*(qnorm(1-alpha)+qnorm(1-beta))^2)/((y-taus)^2))*
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
    return(
      integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            ((4*(qnorm(1-alpha)+qnorm(1-beta))^2)/(y-taus)^2) *
              dnorm(y,
                    mean = x,
                    sd = sqrt(4/n2))*
              prior_normal(x, w, Delta1, Delta2, in1, in2, a, b)
          }, kappa, Inf)$value
        })
      },  - Inf, Inf)$value
    )
  }
}

# Expected probability of a successful program: EsP
EPsProg_normal <-  function(kappa, n2, alpha, beta, step1, step2, w, Delta1, Delta2, in1, in2, a, b, taus, gamma, fixed){

  c = (qnorm(1 - alpha) + qnorm(1 - beta))^2

  if(fixed){
    return(
          integrate(function(y){
            ( pnorm(qnorm(1 - alpha) + step2/sqrt((y-taus)^2/c),
                    mean = (Delta1+gamma)/sqrt((y-taus)^2/c),
                    sd = 1) -
                pnorm(qnorm(1 - alpha) + step1/sqrt((y-taus)^2/c),
                      mean = (Delta1+gamma)/sqrt((y-taus)^2/c),
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
            ( pnorm(qnorm(1 - alpha) + step2/sqrt((y-taus)^2/c),
                    mean = (x+gamma)/sqrt((y-taus)^2/c),
                    sd = 1) -
                pnorm(qnorm(1 - alpha) + step1/sqrt((y-taus)^2/c),
                      mean = (x+gamma)/sqrt((y-taus)^2/c),
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
utility_normal <-  function(n2, kappa, w, Delta1, Delta2, in1, in2, a, b,
                            alpha, beta, 
                            c2, c3, c02, c03, 
                            K, N, S,
                            steps1, stepm1, stepl1,
                            b1, b2, b3, taus,
                            gamma, fixed, case){

  n3  <-  En3_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                      w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b, taus = taus,
                      fixed = fixed)
  
  n3 <- ifelse(is.na(n3),99999,n3)
  
  n3  <- ceiling(n3)
  
  if(round(n3/2) != n3 / 2) {n3 = n3 + 1}

  if(n3==100000||n2+n3>N){
    
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
      prob1 <-  EPsProg_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                             step1 = steps1, step2 =  steps2,
                             w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                             taus = taus, gamma = gamma, fixed = fixed)
      prob2 <-  EPsProg_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                             step1 =  stepm1, step2 =  stepm2,
                             w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                             taus = taus, gamma = gamma, fixed = fixed)
      prob3 <-  EPsProg_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                             step1 =  stepl1, step2 = stepl2,
                             w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                             taus = taus, gamma = gamma, fixed = fixed)

      SP    <-  prob1 + prob2 + prob3
      
      if(SP<S){
      
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
      }else{
      
      if(case =="2b"){
        
        prob0 <-  EPsProg_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta,
                                 step1 = taus, step2 =  steps1,
                                 w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                                 taus = taus, gamma = gamma, fixed = fixed)
        
        SP0    <- prob0 + prob1 + prob2 + prob3
        
        G0     <- b1/2 *prob0 + b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
        
        EU0    <-  - K2 - K3 + G0
        
        return(c(EU0, n3, SP0, pg, K2, K3, prob1, prob2, prob3))
       
      }else{
        
        G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
        
        EU    <-  - K2 - K3 + G
        
        return(c(EU, n3, SP, pg, K2, K3, prob1, prob2, prob3))  
       
      }
      
      #output: expected utility Eud, En3, EsP, Epgo, cost phase II and III
      }
    }
  }
}

