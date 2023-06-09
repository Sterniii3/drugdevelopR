###############################################
# Potentially more than two arms in phase III #
###############################################
# Delta1; fixed true effects for treatment 1 and 2
# y: hat_Delta1_2; estimator in phase II
# z: hat_T_3; normalized estimator in phase III
# n2: total sample size in phase II 
# n3: total sample size in phase III

# 1. Strategy: Only best promising treatment goes to phase III
# -> Phase III is always 2 arm trial (1:1 sample size allocation)
# 2. Strategy: All promising treatments go to phase III
# -> Phase III is 2 or 3 arm trial (1:1 or 1:1:1 sample size allocation)


#' Probability to go to phase III for multiarm programs with normally distributed outcomes 
#' 
#' Given our parameters this function calculates the probability to go to phase III after the second phase was conducted. The considered strategies are as follows:
#' - 1. Strategy: Only best promising treatment goes to phase III
#    -> Phase III is always two-arm trial (1:1 sample size allocation)
#  - 2. Strategy: All promising treatments go to phase III
#    -> Phase III is two- or three-arm trial (1:1 or 1:1:1 sample size allocation)
#' @param kappa threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be an even number
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising")
#' @param case different cases: 1 ("nogo"), 21 (treatment 1 is promising, treatment 2 is not), 22 (treatment 2 is promising, treatment 1 is not), 31 (both treatments are promising, treatment 1 is better), 32 (both treatments are promising, treatment 2 is better)
#' @return The function pgo_normal() returns the probability to go to phase III for multiarm programs with normally distributed outcomes 
#' @examples res <- pgo_normal(kappa = 0.1, n2 = 50, Delta1 = 0.375, Delta2 = 0.625, strategy = 2, case = 31)
#' @keywords internal
#' @export
pgo_normal<-function(kappa,n2,Delta1,Delta2,strategy,case){
  
  # distribution of y, yk~N(thetak,sigmak^2) and correlation rho = 1/2 (equal sample size allocation)
  MEANY    = c(Delta1,Delta2)
  sigma1   = sqrt((6/n2))   # sd of y1 (equal sample size allocation)
  sigma2   = sqrt((6/n2))   # sd of y2 (equal sample size allocation)
  SIGMAY   = matrix(c(sigma1^2,1/2*sigma1*sigma2,1/2*sigma1*sigma2,sigma2^2), nrow = 2, ncol = 2)
  
  if(case==1){# no go

    return(mvtnorm::pmvnorm(lower = c(-Inf,-Inf),
                   upper = c(kappa,kappa),
                   mean  = MEANY,
                   sigma = SIGMAY))
  }
  if(strategy==1){# best promising
    if(case==21){# treatment 1 is promising and better than treatment 2
      
      return(integrate(function(y1){
        sapply(y1,function(y1){ 
          integrate(function(y2){
            mvtnorm::dmvnorm(cbind(y1,y2),
                    mean  = MEANY,
                    sigma = SIGMAY)
          }, -Inf, y1)$value   
        })
      }, kappa, Inf)$value)  
      
    }
    if(case==22){# treatment 2 is promising and better than treatment 1
      
      return(integrate(function(y2){
        sapply(y2,function(y2){ 
          integrate(function(y1){
            mvtnorm::dmvnorm(cbind(y1,y2),
                    mean  = MEANY,
                    sigma = SIGMAY)
          }, -Inf, y2)$value   
        })
      }, kappa, Inf)$value)   
      
    }
  }
  if(strategy==2){# all promising
    if(case==21){# treatment 1 is promising, treatment 2 is not
      
      return(mvtnorm::pmvnorm(lower=c(kappa,-Inf),
                     upper=c(Inf,kappa),
                     mean=MEANY,
                     sigma=SIGMAY))   
      
    }
    if(case==22){# treatment 2 is promising, treatment 1 is not
      return(mvtnorm::pmvnorm(lower=c(-Inf,kappa),
                     upper=c(kappa,Inf),
                     mean=MEANY,
                     sigma=SIGMAY))         
    }
    if(case==31){# both treatments are promising, treatment 1 is better
      
      return(integrate(function(y1){
        sapply(y1,function(y1){ 
          integrate(function(y2){
            mvtnorm::dmvnorm(cbind(y1,y2),
                    mean  = MEANY,
                    sigma = SIGMAY)
          }, kappa, y1)$value   
        })
      }, kappa, Inf)$value)  
      
    }
    if(case==32){# both treatments are promising, treatment 2 is better
      
      return(integrate(function(y2){
        sapply(y2,function(y2){ 
          integrate(function(y1){
            mvtnorm::dmvnorm(cbind(y1,y2),
                    mean  = MEANY,
                    sigma = SIGMAY)
          }, kappa,y2)$value   
        })
      }, kappa, Inf)$value) 
      
    }
  }
  
}

#' Total sample size for phase III trial with l treatments and equal allocation ratio for normally distributed outcomes
#'
#'  Depending on the results of phase II and our strategy ,i.e. whether we proceed only with the best promising treatment (l = 1) or with all promising treatments (l = 2), this program calculates the number of participants in phase III.
#'  
#'  l=1: according to Schoenfeld to guarantee power for the log rank test to detect treatment effect of phase II; 
#'  
#'  l=2: according to Dunnett to guarantee y any-pair power (Horn & Vollandt)
#' @param alpha significance level
#' @param beta  1-'beta' power for calculation of sample size for phase III
#' @param y  y_hat_theta_2; estimator in phase II
#' @param l number of treatments in phase III:
#' - l=1: according to Schoenfeld to guarantee power for the log rank test to detect treatment effect of phase II;  
#' - l=2: according to Dunnett to guarantee y any-pair power (Horn & Vollandt)
#' @return the function ss_normal() returns the total sample size for phase III trial with l treatments and equal allocation ratio
#' @examples res <- ss_normal(alpha = 0.05, beta = 0.1, y = 0.5, l = 1)
#' @keywords internal
#' @export
ss_normal<-function(alpha,beta,y,l){
  
  if(l==1){calpha = qnorm(1-alpha)}
  if(l==2){calpha = as.numeric(mvtnorm::qmvnorm(1-alpha, mean=c(0,0), sigma=matrix(c(1,1/2,1/2,1), nrow=2, ncol=2))[1])}
  
  return(((l+1)*(calpha+qnorm(1-beta))^2)/(y^2)*2)
}

#' Expected sample size for phase III for multiarm programs with normally distributed outcomes
#' 
#' Given phase II results are promising enough to get the "go"-decision to go to phase III this function now calculates the expected sample size for phase III given the cases and strategies listed below.
#' The results of this function are necessary for calculating the utility of the program, which is then in a further step maximized by the `optimal_multiarm_normal()` function 
#' 
#' @param kappa threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta  1-beta power for calculation of sample size for phase III
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)
#' @param case different cases: 1 ("nogo"), 21 (treatment 1 is promising, treatment 2 is not), 22 (treatment 2 is promising, treatment 1 is not), 31 (both treatments are promising, treatment 1 is better), 32 (both treatments are promising, treatment 2 is better)
#' @return The function Ess_normal() returns the expected sample size for phase III when going to phase III when outcomes are normally distributed and we consider multiarm programs, i.e. several phase III trials with different doses or different treatments are performed
#' @examples res <- Ess_normal(kappa = 0.1 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
#'                             Delta1 = 0.375, Delta2 = 0.625, strategy = 3, case = 31)
#' @keywords internal
#' @export
Ess_normal<-function(kappa,n2,alpha,beta,Delta1,Delta2,strategy,case){
  
  
  # distribution of y, yk~N(thetak,sigmak^2) and correlation rho = 1/2 (equal sample size allocation)
  MEANY    = c(Delta1,Delta2)
  sigma1   = sqrt((6/n2))   # sd of y1 (equal sample size allocation)
  sigma2   = sqrt((6/n2))   # sd of y2 (equal sample size allocation)
  SIGMAY   = matrix(c(sigma1^2,1/2*sigma1*sigma2,1/2*sigma1*sigma2,sigma2^2), nrow = 2, ncol = 2)
  
  if(case==1){# no go
    
    return(0)
    
  }
  if(strategy==1){# best promising
    
    if(case==21){# treatment 1 is promising and better than treatment 2
      
      return(integrate(function(y1){
        sapply(y1,function(y1){ 
          integrate(function(y2){
            ss_normal(alpha,beta,y1,1)*
              mvtnorm::dmvnorm(cbind(y1,y2),
                      mean  = MEANY,
                      sigma = SIGMAY)
          }, -Inf, y1)$value   
        })
      }, kappa, Inf)$value)  
      
    }
    if(case==22){# treatment 2 is promising and better than treatment 1
      
      return(integrate(function(y2){
        sapply(y2,function(y2){ 
          integrate(function(y1){
            ss_normal(alpha,beta,y2,1)*
              mvtnorm::dmvnorm(cbind(y1,y2),
                      mean  = MEANY,
                      sigma = SIGMAY)
          },-Inf, y2)$value   
        })
      }, kappa, Inf)$value)   
      
    }
  }
  if(strategy==2){# all promising
    
    if(case==21){# treatment 1 is promising, treatment 2 is not
      
      f <- function(y){ 
        ss_normal(alpha,beta,y[1],1)*mvtnorm::dmvnorm(c(y[1],y[2]), mean  = MEANY, sigma = SIGMAY)
      }
      
      return(cubature::adaptIntegrate(f, lowerLimit = c(kappa, -Inf), upperLimit = c(Inf, kappa))$integral)
      
    }
    if(case==22){# treatment 2 is promising, treatment 1 is not
      
      f <- function(y){ 
        ss_normal(alpha,beta,y[2],1)*mvtnorm::dmvnorm(c(y[1],y[2]), mean  = MEANY, sigma = SIGMAY)
      }
      
      return(cubature::adaptIntegrate(f, lowerLimit = c(-Inf, kappa), upperLimit = c(kappa, Inf))$integral)
      
    }
    if(case==31){# both treatments are promising, treatment 1 is better
      
      return(integrate(function(y1){
        sapply(y1,function(y1){ 
          integrate(function(y2){
            ss_normal(alpha,beta,y2,2)*
              mvtnorm::dmvnorm(cbind(y1,y2),
                      mean  = MEANY,
                      sigma = SIGMAY)
          }, kappa, y1)$value   
        })
      }, kappa, Inf)$value)  
      
    }
    if(case==32){# both treatments are promising, treatment 2 is better
      
      return(integrate(function(y2){
        sapply(y2,function(y2){ 
          integrate(function(y1){
            ss_normal(alpha,beta,y1,2)*
              mvtnorm::dmvnorm(cbind(y1,y2),
                      mean  = MEANY,
                      sigma = SIGMAY)
          }, kappa, y2)$value   
        })
      }, kappa, Inf)$value)  
      
    }
  }  
  
} 

#' Probability of a successful program for multiarm programs with normally distributed outcomes 
#' 
#' Given we get the "go"-decision in phase II, this functions now calculates the probability that the results of the confirmatory trial (phase III) are significant, i.e. we have a statistically relevant positive effect of the treatment.
#' @param kappa threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta  1-beta power for calculation of sample size for phase III
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param step1 lower boundary for effect size
#' @param step2 upper boundary for effect size
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)
#' @param case different cases: 1 ("nogo"), 21 (treatment 1 is promising, treatment 2 is not), 22 (treatment 2 is promising, treatment 1 is not), 31 (both treatments are promising, treatment 1 is better), 32 (both treatments are promising, treatment 2 is better)
#' @return The function PsProg_normal() returns the probability of a successful program.
#' @examples res <- PsProg_normal(kappa = 0.1 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
#'                             Delta1 = 0.375, Delta2 = 0.625,  step1 = 0, step2 = 0.5,
#'                             strategy = 3, case = 31)
#' @keywords internal
#' @export
PsProg_normal<-function(kappa,n2,alpha,beta,Delta1,Delta2,step1,step2,strategy,case){
  

  
  # distribution of y, yk~N(thetak,sigmak^2) and correlation rho = 1/2 (equal sample size allocation)
  MEANY    = c(Delta1,Delta2)
  sigma1   = sqrt((6/n2))   # sd of y1 (equal sample size allocation)
  sigma2   = sqrt((6/n2))   # sd of y2 (equal sample size allocation)
  SIGMAY   = matrix(c(sigma1^2,1/2*sigma1*sigma2,1/2*sigma1*sigma2,sigma2^2), nrow = 2, ncol = 2)
  
  if(case==1){# no go
    
    return(0)
  }
  if(strategy==1){# best promising
    if(case==21){# treatment 1 is promising and better than treatment 2
      
      c     = (qnorm(1-alpha)+qnorm(1-beta))^2 
      
      return(integrate(function(y1){ 
        sapply(y1,function(y1){
          integrate(function(y2){ 
            ( pnorm(qnorm(1-alpha)+step2/(sqrt(y1^2/c)),
                    mean=Delta1/(sqrt(y1^2/c)),
                    sd=1) -
                pnorm(qnorm(1-alpha)+step1/(sqrt(y1^2/c)),
                      mean=Delta1/(sqrt(y1^2/c)),
                      sd=1) )*
              mvtnorm::dmvnorm(cbind(y1,y2),
                      mean=MEANY,
                      sigma=SIGMAY) 
          }, -Inf,y1)$value
        })
      }, kappa,Inf)$value)
      
    }
    if(case==22){# treatment 2 is promising and better than treatment 1
      
      c     = (qnorm(1-alpha)+qnorm(1-beta))^2 
      
      return(integrate(function(y2){ 
        sapply(y2,function(y2){
          integrate(function(y1){ 
            ( pnorm(qnorm(1-alpha)+step2/(sqrt(y2^2/c)),
                    mean=Delta2/(sqrt(y2^2/c)),
                    sd=1) -
                pnorm(qnorm(1-alpha)+step1/(sqrt(y2^2/c)),
                      mean=Delta2/(sqrt(y2^2/c)),
                      sd=1) )*
              mvtnorm::dmvnorm(cbind(y1,y2),
                      mean=MEANY,
                      sigma=SIGMAY) 
          }, -Inf,y2)$value
        })
      }, kappa,Inf)$value)
      
    }
  }
  if(strategy==2){# all promising
    
    if(case==21){# treatment 1 is promising, treatment 2 is not
      
      c     = (qnorm(1-alpha)+qnorm(1-beta))^2 
      
      f <- function(y){ 
        ( pnorm(qnorm(1-alpha)+step2/(sqrt(y[1]^2/c)),
                mean=Delta1/(sqrt(y[1]^2/c)),
                sd=1)-
            pnorm(qnorm(1-alpha)+step1/(sqrt(y[1]^2/c)),
                  mean=Delta1/(sqrt(y[1]^2/c)),
                  sd=1) )*
          mvtnorm::dmvnorm(c(y[1],y[2]),
                  mean=MEANY,
                  sigma=SIGMAY)
      }
      
      return(cubature::adaptIntegrate(f, lowerLimit = c(kappa, -Inf), upperLimit = c(Inf, kappa))$integral)
      
    }
    if(case==22){# treatment 2 is promising, treatment 1 is not 
      
      c     = (qnorm(1-alpha)+qnorm(1-beta))^2 
      
      f <- function(y){ 
        ( pnorm(qnorm(1-alpha)+step2/(sqrt(y[2]^2/c)),
                mean=Delta2/(sqrt(y[2]^2/c)),
                sd=1)-
            pnorm(qnorm(1-alpha)+step1/(sqrt(y[2]^2/c)),
                  mean=Delta2/(sqrt(y[2]^2/c)),
                  sd=1) )*
          mvtnorm::dmvnorm(c(y[1],y[2]),
                  mean=MEANY,
                  sigma=SIGMAY)
      }
      
      return(cubature::adaptIntegrate(f, lowerLimit = c(-Inf, kappa), upperLimit = c(kappa, Inf))$integral)
      
    }
    if(case==31){# both treatments are promising, treatment 1 is better
      
      SIGMAZ   = matrix(c(1,1/2,1/2,1), nrow = 2, ncol = 2)
      calpha   = as.numeric(mvtnorm::qmvnorm(1-alpha, mean=c(0,0), sigma= SIGMAZ)[1])
      c        = (calpha+qnorm(1-beta))^2 
      
      return(integrate(function(y1){ 
        sapply(y1,function(y1){
          integrate(function(y2){ 
            sapply(y2,function(y2){ # How to erase??
              ( mvtnorm::pmvnorm(lower=c(-Inf,-Inf),
                        upper=c(calpha+step2/(sqrt(y2^2/c)),
                                calpha+step2/(sqrt(y2^2/c))),
                        mean=c(Delta1/(sqrt(y2^2/c)),Delta2/(sqrt(y2^2/c))),
                        sigma=SIGMAZ)-
                  mvtnorm::pmvnorm(lower=c(-Inf,-Inf),
                          upper=c(calpha+step1/(sqrt(y2^2/c)),
                                  calpha+step1/(sqrt(y2^2/c))),
                          mean=c(Delta1/(sqrt(y2^2/c)),Delta2/(sqrt(y2^2/c))),
                          sigma=SIGMAZ) )*
                mvtnorm::dmvnorm(c(y1,y2),
                        mean=MEANY,
                        sigma=SIGMAY) 
            })
          }, kappa,y1)$value
        })
      }, kappa,Inf)$value)
      
      
    }
    if(case==32){# both treatments are promising, treatment 2 is better
      
      SIGMAZ   = matrix(c(1,1/2,1/2,1), nrow = 2, ncol = 2) 
      calpha   = as.numeric(mvtnorm::qmvnorm(1-alpha, mean=c(0,0), sigma= SIGMAZ)[1])
      c        = (calpha+qnorm(1-beta))^2 
      
      return(integrate(function(y2){ 
        sapply(y2,function(y2){
          integrate(function(y1){ 
            sapply(y1,function(y1){ # How to erase??
              ( mvtnorm::pmvnorm(lower=c(-Inf,-Inf),
                        upper=c(calpha+step2/(sqrt(y1^2/c)),
                                calpha+step2/(sqrt(y1^2/c))),
                        mean=c(Delta1/(sqrt(y1^2/c)),Delta2/(sqrt(y1^2/c))),
                        sigma=SIGMAZ)-
                  mvtnorm::pmvnorm(lower=c(-Inf,-Inf),
                          upper=c(calpha+step1/(sqrt(y1^2/c)),
                                  calpha+step1/(sqrt(y1^2/c))),
                          mean=c(Delta1/(sqrt(y1^2/c)),Delta2/(sqrt(y1^2/c))),
                          sigma=SIGMAZ) )*
                mvtnorm::dmvnorm(c(y1,y2),
                        mean=MEANY,
                        sigma=SIGMAY) 
            })
          }, kappa,y2)$value
        })
      }, kappa,Inf)$value)
      
    }
  }
  
} 

#' Utility function for multiarm programs with normally distributed outcomes
#' 
#' The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters as on the sample size and expected probability of a successful program. 
#' The utility is in further step maximized by the `optimal_multiarm_normal()` function.
#' 
#' @param kappa threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta  1-beta power for calculation of sample size for phase III
#' @param Delta1 assumed true treatment effect for standardized difference in means
#' @param Delta2 assumed true treatment effect for standardized difference in means
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising") 
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param steps1 lower boundary for effect size category "small", default: 0
#' @param stepm1 lower boundary for effect size category "medium" = upper boundary for effect size category "small" default: 0.5
#' @param stepl1 lower boundary for effect size category "large" = upper boundary for effect size category "medium", default: 0.8
#' @param b1 expected gain for effect size category "small"
#' @param b2 expected gain for effect size category "medium"
#' @param b3 expected gain for effect size category "large"
#' @return The output of the function `utility_multiarm_normal()` is the expected utility of the program.
#' @examples res <- utility_multiarm_normal(n2 = 50, kappa = 0.8, alpha = 0.05, beta = 0.1,
#'                             Delta1 = 0.375, Delta2 = 0.625, strategy = 1,
#'                             c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                             K = Inf, N = Inf, S = -Inf,  
#'                             steps1 = 0, stepm1 = 0.5,   stepl1 = 0.8,
#'                             b1 = 1000, b2 = 2000, b3 = 3000)
#' @keywords internal
#' @export
utility_multiarm_normal<-function(n2,kappa,alpha,beta,
                                  Delta1,Delta2,strategy,
                                  c2,c02,c3,c03,K,N,S,
                                  steps1, stepm1, stepl1,b1, b2, b3){ 
  
  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- Inf
  
  if(strategy==1){
    
    n321    = Ess_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                  strategy=strategy,case=21)
    
    n322    = Ess_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                  strategy=strategy,case=22)
    
    n3      = ceiling(n321+n322)           # total expected sample size for phase III
    
    if(n2+n3>N){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
    }else{
      
      pnogo   = pgo_normal(kappa=kappa,n2=n2,Delta1=Delta1,Delta2=Delta2,strategy=strategy,case=1)
      
      K2    <-  c02 + c2 * n2  #cost phase II
      K3    <-  c03 * (1-pnogo) + c3 * n3  #cost phase III
      
      if(K2+K3>K){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        #output: expected utility Eud, En3, EsP, Epgo
        
      }else{
        # probability of a successful program; small, medium, large effect size
        prob121 = PsProg_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                         step1=steps1,step2=steps2,strategy=strategy,case=21)
        prob221 = PsProg_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                         step1=stepm1,step2=stepm2,strategy=strategy,case=21)
        prob321 = PsProg_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                         step1=stepl1,step2=stepl2,strategy=strategy,case=21)
        
        prob122 = PsProg_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                         step1=steps1,step2=steps2,strategy=strategy,case=22)
        prob222 = PsProg_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                         step1=stepm1,step2=stepm2,strategy=strategy,case=22)
        prob322 = PsProg_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                         step1=stepl1,step2=stepl2,strategy=strategy,case=22)
        
        SP    = prob121+prob221+prob321 +                           # probability of a successful program
          prob122+prob222+prob322 
        
        if(SP<S){
          
          return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
          
        }else{
          
          G     = b1*prob121+b2*prob221+b3*prob321 +                  # gain
            b1*prob122+b2*prob222+b3*prob322                   
          EU    = -K2-K3+G                                            # total expected utility
          
          SP2 = SP
          SP3 = 0
          
          return(c(EU, n3, SP, 1-pnogo, SP2, SP3, K2, K3))
          
        }
      }
    }
    
  }
  
  if(strategy==2){
    
    n321    = Ess_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                  strategy=strategy,case=21)
    n322    = Ess_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                  strategy=strategy,case=22)
    n331    = Ess_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                  strategy=strategy,case=31)
    n332    = Ess_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                  strategy=strategy,case=32)
    n3      = ceiling(n321+n322+n331+n332)   # total expected sample size for phase III
    
    if(n2+n3>N){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
    }else{
      
      pnogo   = pgo_normal(kappa=kappa,n2=n2,Delta1=Delta1,Delta2=Delta2,strategy=strategy,case=1)
      
      K2    <-  c02 + c2 * n2  #cost phase II
      K3    <-  c03 * (1-pnogo) + c3 * n3  #cost phase III
      
      if(K2+K3>K){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        #output: expected utility Eud, En3, EsP, Epgo
        
      }else{
        
        # probability of a successful program; small, medium, large effect size
        prob121 = PsProg_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                         step1=steps1,step2=steps2,strategy=strategy,case=21)
        prob221 = PsProg_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                         step1=stepm1,step2=stepm2,strategy=strategy,case=21)
        prob321 = PsProg_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                         step1=stepl1,step2=stepl2,strategy=strategy,case=21)
        
        prob122 = PsProg_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                         step1=steps1,step2=steps2,strategy=strategy,case=22)
        prob222 = PsProg_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                         step1=stepm1,step2=stepm2,strategy=strategy,case=22)
        prob322 = PsProg_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                         step1=stepl1,step2=stepl2,strategy=strategy,case=22)
        
        prob131 = PsProg_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                         step1=steps1,step2=steps2,strategy=strategy,case=31)
        prob231 = PsProg_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                         step1=stepm1,step2=stepm2,strategy=strategy,case=31)
        prob331 = PsProg_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                         step1=stepl1,step2=stepl2,strategy=strategy,case=31)
        
        prob132 = PsProg_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                         step1=steps1,step2=steps2,strategy=strategy,case=32)
        prob232 = PsProg_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                         step1=stepm1,step2=stepm2,strategy=strategy,case=32)
        prob332 = PsProg_normal(kappa=kappa,n2=n2,alpha=alpha,beta=beta,Delta1=Delta1,Delta2=Delta2,
                         step1=stepl1,step2=stepl2,strategy=strategy,case=32)
        
        SP2     = prob121+prob221+prob321 +                   # probability of a successful program with 
          prob122+prob222+prob322                            # two arms phase III trial 
        SP3     = prob131+prob231+prob331 +                   # with three arms phase III trial
          prob132+prob232+prob332
        SP      =  SP2 + SP3                                  # probability of a successful program
        
        if(SP<S){
          
          return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
          
        }else{
          
          G       = b1*prob121+b2*prob221+b3*prob321 +                # gain
            b1*prob122+b2*prob222+b3*prob322 +
            b1*prob131+b2*prob231+b3*prob331 + 
            b1*prob132+b2*prob232+b3*prob332                   
          EU    = -K2-K3+G                                            # total expected utility
          
          
          return(c(EU, n3, SP, 1-pnogo, SP2, SP3, K2, K3))  
          
        }
      }
    }
  }
}
