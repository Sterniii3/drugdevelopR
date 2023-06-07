###############################################
# Potentially more than two arms in phase III #
###############################################
# theta; fixed true effects for treatment 1 and 2
# y: hat_theta_2; estimator in phase II
# z: hat_T_3; normalized estimator in phase III
# n2: total sample size in phase II 
# n3: total sample size in phase III

# 1. Strategy: Only best promising treatment goes to phase III
# -> Phase III is always 2 arm trial (1:1 sample size allocatiob)
# 2. Strategy: All promising treatments go to phase III
# -> Phase III is 2 or 3 arm trial (1:1 or 1:1:1 sample size allocation)


#' Probability to go to phase III for multiarm programs with time-to-event outcomes 
#' 
#' Given our parameters this function calculates the probability to go to phase III after the second phase was conducted. The considered strategies are as follows:
#' - 1. Strategy: Only best promising treatment goes to phase III
#    -> Phase III is always two-arm trial (1:1 sample size allocation)
#  - 2. Strategy: All promising treatments go to phase III
#    -> Phase III is two- or three-arm trial (1:1 or 1:1:1 sample size allocation)
#' @param HRgo threshold value for the go/no-go decision rule
#' @param n2 total sample size in phase II, must be divisible by 3
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param ec control arm event rate for phase II and III
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising") 
#' @param case different cases: 1 ("nogo"), 21 (treatment 1 is promising, treatment 2 is not), 22 (treatment 2 is promising, treatment 1 is not), 31 (both treatments are promising, treatment 1 is better), 32 (both treatments are promising, treatment 2 is better)
#' @return The function pgo_tte() returns the probability to go to phase III.
#' @examples res <- pgo_tte(HRgo = 0.8, n2 = 48 , ec = 0.6, hr1 = 0.7, hr2 = 0.8, strategy = 2, case = 31)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @keywords internal
#' @export
pgo_tte<-function(HRgo,n2,ec,hr1,hr2,strategy,case){
  
  et1      = 1 - (1-ec)^hr1     # event rate in arm 1
  et2      = 1 - (1-ec)^hr2     # event rate in arm 2
  
  # distribution of y, yk~N(thetak,sigmak^2) and correlation rho = 1/2 (equal sample size allocation)
  MEANY    = -log(c(hr1,hr2))
  sigma1   = sqrt((3/n2)*((1/ec)+(1/et1)))   # sd of y1 (equal sample size allocation)
  sigma2   = sqrt((3/n2)*((1/ec)+(1/et2)))   # sd of y2 (equal sample size allocation)
  SIGMAY   = matrix(c(sigma1^2,1/2*sigma1*sigma2,1/2*sigma1*sigma2,sigma2^2), nrow = 2, ncol = 2)
  
  if(case==1){# no go
    
    return(mvtnorm::pmvnorm(lower = c(-Inf,-Inf),
                   upper = c(-log(HRgo),-log(HRgo)),
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
      }, -log(HRgo), Inf)$value)  
      
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
      }, -log(HRgo), Inf)$value)   
      
    }
  }
  if(strategy==2){# all promising
    if(case==21){# treatment 1 is promising, treatment 2 is not
      
      return(mvtnorm::pmvnorm(lower=c(-log(HRgo),-Inf),
                     upper=c(Inf,-log(HRgo)),
                     mean=MEANY,
                     sigma=SIGMAY))   
      
    }
    if(case==22){# treatment 2 is promising, treatment 1 is not
      return(mvtnorm::pmvnorm(lower=c(-Inf,-log(HRgo)),
                     upper=c(-log(HRgo),Inf),
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
          }, -log(HRgo), y1)$value   
        })
      }, -log(HRgo), Inf)$value)  
      
    }
    if(case==32){# both treatments are promising, treatment 2 is better
      
      return(integrate(function(y2){
        sapply(y2,function(y2){ 
          integrate(function(y1){
            mvtnorm::dmvnorm(cbind(y1,y2),
                    mean  = MEANY,
                    sigma = SIGMAY)
          }, -log(HRgo),y2)$value   
        })
      }, -log(HRgo), Inf)$value) 
      
    }
  }
  
}

#' Total sample size for phase III trial with l treatments and equal allocation ratio for time-to-event outcomes 
#'
#'  Depending on the results of phase II and our strategy ,i.e. whether we proceed only with the best promising treatment (l = 1) or with all promising treatments (l = 2), this program calculates the number of participants in phase III.
#'  
#'  l=1: according to Schoenfeld to guarantee power for the log rank test to detect treatment effect of phase II; 
#'  
#'  l=2: according to Dunnett to guarantee y any-pair power (Horn & Vollandt)
#' @param alpha significance level
#' @param beta  1-beta power for calculation of sample size for phase III
#' @param ec control arm event rate for phase II and III
#' @param ek event rate of arm k (either arm 1 or arm 2)
#' @param y hat_theta_2; estimator in phase II
#' @param l number of treatments in phase III:
#' - l=1: according to Schoenfeld to guarantee power for the log rank test to detect treatment effect of phase II;  
#' - l=2: according to Dunnett to guarantee y any-pair power (Horn & Vollandt)
#' @return the function ss_tte() returns the total sample size for phase III trial with l treatments and equal allocation ratio
#' @examples res <- ss_tte(alpha = 0.05, beta = 0.1, ec = 0.6, ek = 0.8, y = 0.5, l=1)
#' @editor Johannes Cepicka
#' @editDate 2022-05-08
#' @keywords internal
#' @export
ss_tte<-function(alpha,beta,ec,ek,y,l){

  if(l==1){calpha = qnorm(1-alpha)}
  if(l==2){calpha = as.numeric(mvtnorm::qmvnorm(1-alpha, mean=c(0,0), sigma=matrix(c(1,1/2,1/2,1), nrow=2, ncol=2))[1])}
  
  return(((l+1)*(calpha+qnorm(1-beta))^2)/(y^2)*((1/ec)+(1/ek)))
}

#' Expected sample size for phase III for multiarm programs with time-to-event outcomes
#' 
#' Given phase II results are promising enough to get the "go"-decision to go to phase III this function now calculates the expected sample size for phase III given the cases and strategies listed below.
#' The results of this function are necessary for calculating the utility of the program, which is then in a further step maximized by the `optimal_multiarm()` function 
#' 
#' @param HRgo threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta  1-beta power for calculation of sample size for phase III
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param ec control arm event rate for phase II and III
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)
#' @param case different cases: 1 ("nogo"), 21 (treatment 1 is promising, treatment 2 is not), 22 (treatment 2 is promising, treatment 1 is not), 31 (both treatments are promising, treatment 1 is better), 32 (both treatments are promising, treatment 2 is better)
#' @return the function Ess_tte() returns the expected sample size for phase III when going to phase III
#' @examples res <- Ess_tte(HRgo = 0.8 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
#'                             ec = 0.6, hr1 = 0.7, hr2 = 0.8, strategy = 2, case = 21)
#' @editor Johannes Cepicka
#' @editDate 2022-05-08
#' @keywords internal
#' @export
Ess_tte<-function(HRgo,n2,alpha,beta,ec,hr1,hr2,strategy,case){
  
  et1      = 1 - (1-ec)^hr1     # event rate in arm 1
  et2      = 1 - (1-ec)^hr2     # event rate in arm 2
  
  # distribution of y, yk~N(thetak,sigmak^2) and correlation rho = 1/2 (equal sample size allocation)
  MEANY    = -log(c(hr1,hr2))
  sigma1   = sqrt((3/n2)*((1/ec)+(1/et1)))   # sd of y1 (equal sample size allocation)
  sigma2   = sqrt((3/n2)*((1/ec)+(1/et2)))   # sd of y2 (equal sample size allocation)
  SIGMAY   = matrix(c(sigma1^2,1/2*sigma1*sigma2,1/2*sigma1*sigma2,sigma2^2), nrow = 2, ncol = 2)
  
  if(case==1){# no go
    
    return(0)
    
  }
  if(strategy==1){# best promising
    
    if(case==21){# treatment 1 is promising and better than treatment 2
      
      return(integrate(function(y1){
        sapply(y1,function(y1){ 
          integrate(function(y2){
            ss_tte(alpha,beta,ec,et1,y1,1)*
              mvtnorm::dmvnorm(cbind(y1,y2),
                      mean  = MEANY,
                      sigma = SIGMAY)
          }, -Inf, y1)$value   
        })
      }, -log(HRgo), Inf)$value)  
      
    }
    if(case==22){# treatment 2 is promising and better than treatment 1
      
      return(integrate(function(y2){
        sapply(y2,function(y2){ 
          integrate(function(y1){
            ss_tte(alpha,beta,ec,et2,y2,1)*
              mvtnorm::dmvnorm(cbind(y1,y2),
                      mean  = MEANY,
                      sigma = SIGMAY)
          },-Inf, y2)$value   
        })
      }, -log(HRgo), Inf)$value)   
      
    }
  }
  if(strategy==2){# all promising
    
    if(case==21){# treatment 1 is promising, treatment 2 is not
      
      f <- function(y){ 
        ss_tte(alpha,beta,ec,et1,y[1],1)*mvtnorm::dmvnorm(c(y[1],y[2]), mean  = MEANY, sigma = SIGMAY)
      }
      
      return(cubature::adaptIntegrate(f, lowerLimit = c(-log(HRgo), -Inf), upperLimit = c(Inf, -log(HRgo)))$integral)
      
    }
    if(case==22){# treatment 2 is promising, treatment 1 is not
      
      f <- function(y){ 
        ss_tte(alpha,beta,ec,et2,y[2],1)*mvtnorm::dmvnorm(c(y[1],y[2]), mean  = MEANY, sigma = SIGMAY)
      }
      
      return(cubature::adaptIntegrate(f, lowerLimit = c(-Inf, -log(HRgo)), upperLimit = c(-log(HRgo), Inf))$integral)
      
    }
    if(case==31){# both treatments are promising, treatment 1 is better
      
      return(integrate(function(y1){
        sapply(y1,function(y1){ 
          integrate(function(y2){
            ss_tte(alpha,beta,ec,et2,y2,2)*
              mvtnorm::dmvnorm(cbind(y1,y2),
                      mean  = MEANY,
                      sigma = SIGMAY)
          }, -log(HRgo), y1)$value   
        })
      }, -log(HRgo), Inf)$value)  
      
    }
    if(case==32){# both treatments are promising, treatment 2 is better
      
      return(integrate(function(y2){
        sapply(y2,function(y2){ 
          integrate(function(y1){
            ss_tte(alpha,beta,ec,et1,y1,2)*
              mvtnorm::dmvnorm(cbind(y1,y2),
                      mean  = MEANY,
                      sigma = SIGMAY)
          }, -log(HRgo), y2)$value   
        })
      }, -log(HRgo), Inf)$value)  
      
    }
  }  
  
} 

#' Probability of a successful program for multiarm programs with time-to-event outcomes 
#' 
#' Given we get the "go"-decision in phase II, this functions now calculates the probability that the results of the confirmatory trial (phase III) are significant, i.e. we have a statistically relevant positive effect of the treatment.
#' @param HRgo threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be divisible by three
#' @param alpha significance level
#' @param beta  1-beta power for calculation of sample size for phase III
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param ec control arm event rate for phase II and III
#' @param step1 lower boundary for effect size
#' @param step2 upper boundary for effect size
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)
#' @param case different cases: 1 ("nogo"), 21 (treatment 1 is promising, treatment 2 is not), 22 (treatment 2 is promising, treatment 1 is not), 31 (both treatments are promising, treatment 1 is better), 32 (both treatments are promising, treatment 2 is better)
#' @return The function PsProg_tte() returns the probability of a successful program
#' @examples res <- PsProg_tte(HRgo = 0.8 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
#'                             ec = 0.6, hr1 = 0.7, hr2 = 0.8, step1 = 1, step2 = 0.95,
#'                             strategy = 2, case = 21)
#' @editor Johannes Cepicka
#' @editDate 2022-05-08
#' @keywords internal
#' @export
PsProg_tte<-function(HRgo,n2,alpha,beta,ec,hr1,hr2,step1,step2,strategy,case){
  
  et1      = 1 - (1-ec)^hr1    # event rate in arm 1
  et2      = 1 - (1-ec)^hr2     # event rate in arm 2
  
  # distribution of y, yk~N(thetak,sigmak^2) and correlation rho = 1/2 (equal sample size allocation)
  MEANY    = -log(c(hr1,hr2))
  sigma1   = sqrt((3/n2)*((1/ec)+(1/et1)))   # sd of y1 (equal sample size allocation)
  sigma2   = sqrt((3/n2)*((1/ec)+(1/et2)))   # sd of y2 (equal sample size allocation)
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
            ( pnorm(qnorm(1-alpha)-log(step2)/(sqrt(y1^2/c)),
                    mean=-log(hr1)/(sqrt(y1^2/c)),
                    sd=1) -
                pnorm(qnorm(1-alpha)-log(step1)/(sqrt(y1^2/c)),
                      mean=-log(hr1)/(sqrt(y1^2/c)),
                      sd=1) )*
              mvtnorm::dmvnorm(cbind(y1,y2),
                      mean=MEANY,
                      sigma=SIGMAY) 
          }, -Inf,y1)$value
        })
      }, -log(HRgo),Inf)$value)
      
    }
    if(case==22){# treatment 2 is promising and better than treatment 1
      
      c     = (qnorm(1-alpha)+qnorm(1-beta))^2 
      
      return(integrate(function(y2){ 
        sapply(y2,function(y2){
          integrate(function(y1){ 
            ( pnorm(qnorm(1-alpha)-log(step2)/(sqrt(y2^2/c)),
                    mean=-log(hr2)/(sqrt(y2^2/c)),
                    sd=1) -
                pnorm(qnorm(1-alpha)-log(step1)/(sqrt(y2^2/c)),
                      mean=-log(hr2)/(sqrt(y2^2/c)),
                      sd=1) )*
              mvtnorm::dmvnorm(cbind(y1,y2),
                      mean=MEANY,
                      sigma=SIGMAY) 
          }, -Inf,y2)$value
        })
      }, -log(HRgo),Inf)$value)
      
    }
  }
  if(strategy==2){# all promising
    
    if(case==21){# treatment 1 is promising, treatment 2 is not
      
      c     = (qnorm(1-alpha)+qnorm(1-beta))^2 
      
      f <- function(y){ 
        ( pnorm(qnorm(1-alpha)-log(step2)/(sqrt(y[1]^2/c)),
                mean=-log(hr1)/(sqrt(y[1]^2/c)),
                sd=1)-
            pnorm(qnorm(1-alpha)-log(step1)/(sqrt(y[1]^2/c)),
                  mean=-log(hr1)/(sqrt(y[1]^2/c)),
                  sd=1) )*
          mvtnorm::dmvnorm(c(y[1],y[2]),
                  mean=MEANY,
                  sigma=SIGMAY)
      }
      
      return(cubature::adaptIntegrate(f, lowerLimit = c(-log(HRgo), -Inf), upperLimit = c(Inf, -log(HRgo)))$integral)
      
    }
    if(case==22){# treatment 2 is promising, treatment 1 is not 
      
      c     = (qnorm(1-alpha)+qnorm(1-beta))^2 
      
      f <- function(y){ 
        ( pnorm(qnorm(1-alpha)-log(step2)/(sqrt(y[2]^2/c)),
                mean=-log(hr2)/(sqrt(y[2]^2/c)),
                sd=1)-
            pnorm(qnorm(1-alpha)-log(step1)/(sqrt(y[2]^2/c)),
                  mean=-log(hr2)/(sqrt(y[2]^2/c)),
                  sd=1) )*
          mvtnorm::dmvnorm(c(y[1],y[2]),
                  mean=MEANY,
                  sigma=SIGMAY)
      }
      
      return(cubature::adaptIntegrate(f, lowerLimit = c(-Inf, -log(HRgo)), upperLimit = c(-log(HRgo), Inf))$integral)
      
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
                        upper=c(calpha-log(step2)/(sqrt(y2^2/c)),
                                calpha-log(step2)/(sqrt(y2^2/c))),
                        mean=c(-log(hr1)/(sqrt(y2^2/c)),-log(hr2)/(sqrt(y2^2/c))),
                        sigma=SIGMAZ)-
                  mvtnorm::pmvnorm(lower=c(-Inf,-Inf),
                          upper=c(calpha-log(step1)/(sqrt(y2^2/c)),
                                  calpha-log(step1)/(sqrt(y2^2/c))),
                          mean=c(-log(hr1)/(sqrt(y2^2/c)),-log(hr2)/(sqrt(y2^2/c))),
                          sigma=SIGMAZ) )*
                mvtnorm::dmvnorm(c(y1,y2),
                        mean=MEANY,
                        sigma=SIGMAY) 
            })
          }, -log(HRgo),y1)$value
        })
      }, -log(HRgo),Inf)$value)
      
      
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
                        upper=c(calpha-log(step2)/(sqrt(y1^2/c)),
                                calpha-log(step2)/(sqrt(y1^2/c))),
                        mean=c(-log(hr1)/(sqrt(y1^2/c)),-log(hr2)/(sqrt(y1^2/c))),
                        sigma=SIGMAZ)-
                  mvtnorm::pmvnorm(lower=c(-Inf,-Inf),
                          upper=c(calpha-log(step1)/(sqrt(y1^2/c)),
                                  calpha-log(step1)/(sqrt(y1^2/c))),
                          mean=c(-log(hr1)/(sqrt(y1^2/c)),-log(hr2)/(sqrt(y1^2/c))),
                          sigma=SIGMAZ) )*
                mvtnorm::dmvnorm(c(y1,y2),
                        mean=MEANY,
                        sigma=SIGMAY) 
            })
          }, -log(HRgo),y2)$value
        })
      }, -log(HRgo),Inf)$value)
      
    }
  }
  
} 

#' Utility function for multiarm programs with time-to-event outcomes
#' 
#' The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters as on the sample size and expected probability of a successful program. 
#' The utility is in further step maximized by the `optimal_multiarm()` function.
#' 
#' @param HRgo threshold value for the go/no-go decision rule
#' @param n2 total sample size for phase II; must be divisible by three
#' @param alpha significance level
#' @param beta  1-beta power for calculation of sample size for phase III
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param ec control arm event rate for phase II and III
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising")
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param steps1 lower boundary for effect size category `"small"` in HR scale, default: 1
#' @param stepm1 lower boundary for effect size category `"medium"` in HR scale = upper boundary for effect size category "small" in RR scale, default: 0.95
#' @param stepl1 lower boundary for effect size category `"large"` in HR scale = upper boundary for effect size category "medium" in RR scale, default: 0.85
#' @param b1 expected gain for effect size category `"small"`
#' @param b2 expected gain for effect size category `"medium"`
#' @param b3 expected gain for effect size category `"large"`
#' @return The output of the function `utility_multiarm()` is the expected utility of the program
#' @examples \dontrun{utility_multiarm(n2 = 50, HRgo = 0.8, alpha = 0.05, beta = 0.1,
#'                             hr1 = 0.7, hr2 = 0.8, strategy = 2, ec = 0.6,
#'                             c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                             K = Inf, N = Inf, S = -Inf,  
#'                             steps1 = 1, stepm1 = 0.95,  stepl1 = 0.85,
#'                             b1 = 1000, b2 = 2000, b3 = 3000)}
#' @editor Johannes Cepicka
#' @editDate 2022-05-08
#' @keywords internal
#' @export
utility_multiarm<-function(n2,HRgo,alpha,beta,hr1,hr2,strategy,ec,c2,c02,c3,c03,K,N,S,steps1, stepm1, stepl1,b1, b2, b3){ 
  
  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 0
  
  if(strategy==1){
    
    n321    = Ess_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                  strategy=strategy,case=21)
    
    n322    = Ess_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                  strategy=strategy,case=22)
    
    n3      = ceiling(n321+n322)           # total expected sample size for phase III
  
    if(n2+n3>N){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
    }else{
      
      pnogo   = pgo_tte(HRgo=HRgo,n2=n2,ec=ec,hr1=hr1,hr2=hr2,strategy=strategy,case=1)
      
      K2    <-  c02 + c2 * n2  #cost phase II
      K3    <-  c03 * (1-pnogo) + c3 * n3  #cost phase III
      
      if(K2+K3>K){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        #output: expected utility Eud, En3, EsP, Epgo
        
      }else{
        # probability of a successful program; small, medium, large effect size
        prob121 = PsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                         step1=steps1,step2=steps2,strategy=strategy,case=21)
        prob221 = PsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                         step1=stepm1,step2=stepm2,strategy=strategy,case=21)
        prob321 = PsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                         step1=stepl1,step2=stepl2,strategy=strategy,case=21)
        
        prob122 = PsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                         step1=steps1,step2=steps2,strategy=strategy,case=22)
        prob222 = PsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                         step1=stepm1,step2=stepm2,strategy=strategy,case=22)
        prob322 = PsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
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
    
    n321    = Ess_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                  strategy=strategy,case=21)
    n322    = Ess_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                  strategy=strategy,case=22)
    n331    = Ess_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                  strategy=strategy,case=31)
    n332    = Ess_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                  strategy=strategy,case=32)
    n3      = ceiling(n321+n322+n331+n332)   # total expected sample size for phase III
    
    if(n2+n3>N){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      
    }else{
      
      pnogo   = pgo_tte(HRgo=HRgo,n2=n2,ec=ec,hr1=hr1,hr2=hr2,strategy=strategy,case=1)
      
      K2    <-  c02 + c2 * n2  #cost phase II
      K3    <-  c03 * (1-pnogo) + c3 * n3  #cost phase III
      
      if(K2+K3>K){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        #output: expected utility Eud, En3, EsP, Epgo
        
      }else{

        # probability of a successful program; small, medium, large effect size
        prob121 = PsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                         step1=steps1,step2=steps2,strategy=strategy,case=21)
        prob221 = PsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                         step1=stepm1,step2=stepm2,strategy=strategy,case=21)
        prob321 = PsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                         step1=stepl1,step2=stepl2,strategy=strategy,case=21)
        
        prob122 = PsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                         step1=steps1,step2=steps2,strategy=strategy,case=22)
        prob222 = PsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                         step1=stepm1,step2=stepm2,strategy=strategy,case=22)
        prob322 = PsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                         step1=stepl1,step2=stepl2,strategy=strategy,case=22)
        
        prob131 = PsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                         step1=steps1,step2=steps2,strategy=strategy,case=31)
        prob231 = PsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                         step1=stepm1,step2=stepm2,strategy=strategy,case=31)
        prob331 = PsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                         step1=stepl1,step2=stepl2,strategy=strategy,case=31)
        
        prob132 = PsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                         step1=steps1,step2=steps2,strategy=strategy,case=32)
        prob232 = PsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
                         step1=stepm1,step2=stepm2,strategy=strategy,case=32)
        prob332 = PsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,
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
