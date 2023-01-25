################################################################################
## R-Code for maximizing expected utility
# Utility-based optimization of phase II/III programs with two time-to-event endpoints
# endpoints: overall survival (OS) and progression free survival (PFS)
#
# Author: Marietta Kirchner
# Date: 25.01.2017
################################################################################

#load functions

#' Density for the maximum of two normally distributed random variables
#' 
#' The function `fmax()` will return the value of f(z), which is the value of the density function of the
#' maximum of two normally distributed random variables.
#'
#'  Z = max(X,Y) with X ~ N(mu1,sigma1^2), Y ~ N(mu2,sigma2^2)
#' 
#' f(z)=f1(-z)+f2(-z)
#'@param z integral variable
#'@param mu1 mean of second endpoint 
#'@param mu2 mean of first endpoint
#'@param sigma1 standard deviation of first endpoint
#'@param sigma2 standard deviation of second endpoint
#'@param rho correlation between the two endpoints
#'@return The function `fmax()` will return the value of f(z), which is the value of the density function of the
#'maximum of two normally distributed random variables.
#'@examples res <- fmax(z = 0.5, mu1 = 0.375, mu2 = 0.25, sigma1 = 8, sigma2 = 12, rho = 0.4 )
#'@editor Johannes Cepicka
#'@editDate 2022-04-23
#' @export
fmax<-function (z,mu1,mu2,sigma1,sigma2,rho){ 
  t1<-dnorm(-z,mean=-mu1,sd=sigma1)
  tt<-rho*(mu1-z)/(sigma1*sqrt(1-rho*rho))
  tt<-tt-(mu2-z)/(sigma2*sqrt(1-rho*rho))
  t1<-t1*pnorm(tt)
  t2<-dnorm(-z,mean=-mu2,sd=sigma2)
  tt<-rho*(mu2-z)/(sigma2*sqrt(1-rho*rho))
  tt<-tt-(mu1-z)/(sigma1*sqrt(1-rho*rho))
  t2<-t2*pnorm(tt)
  return(t1+t2)
}

#'@rdname dbivanorm
#'@export
dbivanorm <- function(x,y, mu1,mu2,sigma1,sigma2,rho){ 
  covariancemat <- matrix(c(sigma1,rho*sqrt(sigma1)*sqrt(sigma2), rho*sqrt(sigma1)*sqrt(sigma2),sigma2),ncol=2)
  ff <- dmvnorm(cbind(x,y), mean=c(mu1,mu2),sigma=covariancemat)
  return(ff)
}


#' Probability to go to phase III for multiple endpoints in the time-to-event setting
#' 
#' This function calculates the probability that we go to phase III, i.e. that results of phase II are promising enough to
#' get a successful drug development program. Successful means that at least one endpoint shows a statistically significant positive treatment effect in phase III. 
#' @param HRgo threshold value for the go/no-go decision rule; 
#' @param n2 total sample size for phase II; must be even number
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param id1 amount of information for `hr1` in terms of number of events
#' @param id2 amount of information for `hr2` in terms of number of events
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `hr1` is used as fixed effect
#' @param rho correlation between the two endpoints
#' @return The output of the the function `pgo_multiple_tte()` is the probability to go to phase III.
#' @examples res <- pgo_multiple_tte(HRgo = 0.8, n2 = 50,
#'                                hr1 = 0.75, hr2 = 0.80, id1 = 300, id2 = 600, 
#'                                fixed = TRUE, rho = 0.3)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export
pgo_multiple_tte<-function(HRgo,n2,hr1,hr2,id1,id2,fixed,rho){
  
  e21<-hr1*n2 # number of events phase II for endpoint 1 (PFS) = event rate * sample size
  e22<-hr2*n2 # number of events phase II for endpoint 2 (OS)
  hr <- c(hr1,hr2)
  var1<-4/e21
  var2<-4/e22
  vartrue1 <- sqrt(4/id1)
  vartrue2 <- sqrt(4/id2)
  covmat<-matrix(c(var1, rho*sqrt(var1)*sqrt(var2), rho*sqrt(var1)*sqrt(var2), var2), ncol=2) #covariance-Matrix of c(hr1,hr2)
  
 if(fixed) {
      return(integrate(function(x){
       sapply(x,function(x)
        fmax(x,hr1,hr2,sqrt(var1),sqrt(var2),rho))
        },-log(HRgo),Inf)$value)
    }
  
  else  {
      return(integrate(function(u){
        sapply(u,function(u){
         integrate(function(v){
            sapply(v,function(v){
              integrate(function(x){sapply(x,function(x){
              (fmax(x,u,v,sqrt(var1),sqrt(var2),rho))*(dbivanorm(u,v,-log(hr1),-log(hr2),vartrue1,vartrue2,rho))
            })
          },-log(HRgo),Inf)$value
        })
      },-Inf,Inf )$value
    })
   },-Inf,Inf)$value)
 
 }

  
}


#' Expected sample size for phase III for multiple endpoints with normally distributed outcomes
#' 
#' Given phase II results are promising enough to get the "go"-decision to go to phase III this function now calculates the expected sample size for phase III.
#' The results of this function are necessary for calculating the utility of the program, which is then in a further step maximized by the `optimal_multiple_tte()` function 
#' @param HRgo threshold value for the go/no-go decision rule; 
#' @param n2 total sample size for phase II; must be even number
#' @param beta `1-beta` power for calculation of the number of events for phase III by Schoenfeld (1981) formula
#' @param alpha one- sided significance level
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param id1 amount of information for `hr1` in terms of number of events
#' @param id2 amount of information for `hr2` in terms of number of events
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `hr1` is used as fixed effect
#' @param rho correlation between the two endpoints
#' @return the output of the the function `Ess_multiple_tte()` is the expected number of participants in phase III
#' @examples res <- Ess_multiple_tte(HRgo = 0.8, n2 = 50, alpha = 0.05, beta = 0.1,
#'                                hr1 = 0.75, hr2 = 0.80, 
#'                                id1 = 300, id2 = 600, 
#'                                fixed = TRUE, rho = 0.3)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export
Ess_multiple_tte<-function(HRgo,n2,alpha,beta,hr1,hr2,id1,id2,fixed,rho){

  e21<-hr1*n2 # number of events phase II for endpoint 1 (PFS) = event rate * sample size
  e22<-hr2*n2 # number of events phase II for endpoint 2 (OS)
  hr <- c(hr1,hr2)
  var1<-4/e21
  var2<-4/e22
  vartrue1 <- sqrt(4/id1)
  vartrue2 <- sqrt(4/id2)
  covmat<-matrix(c(var1, rho*sqrt(var1)*sqrt(var2), rho*sqrt(var1)*sqrt(var2), var2), ncol=2) #covariance-Matrix of c(true1,true2)
  
 if(fixed)  {
         return(integrate(function(x){
          sapply(x,function(x){
            ((4*(qnorm(1 - alpha) + qnorm(1 - beta))^2)/x^2)*fmax(x,-log(hr1),-log(hr2),sqrt(var1),sqrt(var2),rho)
            })
        },-log(HRgo),Inf)$value)
      }

  else   {
          return(integrate(function(u){
            sapply(u,function(u){
              integrate(function(v){
                sapply(v,function(v){
                  integrate(function(x){
                    sapply(x,function(x){
                      (((4*(qnorm(1 - alpha) + qnorm(1 - beta))^2)/x^2)*fmax(x,u,v,sqrt(var1),sqrt(var2),rho))*(dbivanorm(u,v,-log(hr1),-log(hr2),vartrue1,vartrue2,rho))
                    })
                  },-log(HRgo),Inf)$value #divide value by pgo to get E(e3|GO)
                })
              },-Inf,Inf)$value
            })
          },-Inf,Inf)$value)
    
  }
}

#' Probabilty that effect in endpoint one is larger than in endpoint two
#' 
#' This function calculated the probability that the treatment effect in endpoint one (or endpoint x) is larger than in endpoint two (or endpoint y), i.e. P(x>y) = P(x-y>0)
#' 
#' Z=X-Y is normally distributed with expectation mu_x - mu_y and variance sigma_x + sigma_y- 2 rho sdx sdy
#' @param n2 total sample size for phase II; must be even number
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param id1 amount of information for `hr1` in terms of number of events
#' @param id2 amount of information for `hr2` in terms of number of events
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `hr1` is used as fixed effect
#' @param rho correlation between the two endpoints
#' @return The output of the the function `pw()` is the probability that endpoint one has a better result than endpoint two
#' @examples res <- pw(n2 = 50,hr1 = 0.75, hr2 = 0.80, id1 = 300, id2 = 600, 
#'                     fixed = FALSE, rho = 0.3)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export
pw <- function(n2,hr1,hr2,id1,id2,fixed,rho){
  
  e21<-hr1*n2 # number of events phase II for endpoint 1 (PFS) = event rate * sample size
  e22<-hr2*n2 # number of events phase II for endpoint 2 (OS)
  hr <- c(hr1,hr2)
  var1<-4/e21
  var2<-4/e22
  vartrue1 <- sqrt(4/id1)
  vartrue2 <- sqrt(4/id2)
  covmat<-matrix(c(var1, rho*sqrt(var1)*sqrt(var2), rho*sqrt(var1)*sqrt(var2), var2), ncol=2) #covariance-Matrix of c(true1,true2)

if(fixed) {
  return(pnorm(0,mean=(-log(hr1)+log(hr2)),sd=sqrt(var1+var2-2*rho*sqrt(var1)*sqrt(var2)),lower.tail=FALSE))
   }

else {
   return (integrate(function(u){
      sapply(u,function(u){
        integrate(function(v){
          sapply(v,function(v){
            pnorm(0,mean=(u-v),sd=sqrt(var1+var2-2*rho*sqrt(var1)*sqrt(var2)),lower.tail=FALSE)*(dbivanorm(u,v,-log(hr1),-log(hr2),vartrue1,vartrue2,rho))
          })
        },-Inf,Inf)$value
      })
   },-Inf,Inf)$value)
  }

}

#E(n3|GO)
expn3go_tte<-function(HRgo,n2,alpha,beta,hr1,hr2,id1,id2,fixed,rho){
  
  expe3go_tte<-Ess_multiple_tte(HRgo,n2,alpha,beta,hr1,hr2,id1,id2,fixed,rho)/pgo_multiple_tte(HRgo,n2,hr1,hr2,id1,id2,fixed,rho)
  
  return(expe3go_tte/hr[1])*pw(n2,hr1,hr2,id1,id2,fixed,rho)+(expe3go_tte/hr[2])*(1-pw(n2,hr1,hr2,id1,id2,fixed,rho))
  }
  
  
#' Expected probability of a successful program for multiple endpoints in a time-to-event setting 
#' 
#' This function calculates the probability that our drug development program is successful.
#' Successful is defined as at least one endpoint showing a statistically significant positive treatment effect in phase III. 
#' @param HRgo threshold value for the go/no-go decision rule; 
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param ec control arm event rate for phase II and III
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param id1 amount of information for `hr1` in terms of sample size
#' @param id2 amount of information for `hr2` in terms of sample size
#' @param step1 lower boundary for effect size
#' @param step2 upper boundary for effect size
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `hr1` is used as fixed effect
#' @param rho correlation between the two endpoints
#' @return The output of the the function `EPsProg_multiple_tte()` is the expected probability of a successful program, when going to phase III.
#' @examples res <- EPsProg_multiple_tte(HRgo = 0.8, n2 = 50, alpha = 0.025, beta = 0.1,
#'                                ec = 1, hr1 = 0.75, hr2 = 0.80,
#'                                id1 = 300, id2 = 600, 
#'                                step1 = 1, step2 = 0.95,
#'                                fixed = TRUE, rho = 0.3)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export
EPsProg_multiple_tte<-function(HRgo,n2,alpha,beta,ec,hr1,hr2,id1,id2,step1,step2,fixed,rho){
 
  e21<-hr1*n2 # number of events phase II for endpoint 1 (PFS) = event rate * sample size
  e22<-hr2*n2 # number of events phase II for endpoint 2 (OS)
  hr <- c(hr1,hr2)
  var1 <- 4/e21
  var2 <- 4/e22
  vartrue1 <- sqrt(4/id1)
  vartrue2 <- sqrt(4/id2)
  covmat<-matrix(c(var1, rho*sqrt(var1)*sqrt(var2), rho*sqrt(var1)*sqrt(var2), var2), ncol=2) #covariance-Matrix of c(true1,true2) 
  
  c     = (qnorm(1-alpha)+qnorm(1-beta))^2 
  
   if(fixed) {return(integrate(function(x){ 
          sapply(x,function(x){ 
            integrate(function(y){ 
              sapply(y,function(y){
                fmax(y,-log(hr1)/sqrt((x^2/c)*(ec/hr[1])),-log(hr2)/sqrt((x^2/c)*(ec/hr[2])),1,1,rho)*fmax(x,-log(hr1),-log(hr2),sqrt(var1),sqrt(var2),rho)
              })
            },qnorm(1-alpha)-log(step1)/sqrt((x^2/c)),qnorm(1-alpha)-log(step2)/sqrt((x^2/c)))$value
          })
        }, -log(HRgo),Inf)$value)
   }
  
  
   else {

     return(integrate(function(u){
        sapply(u,function(u){
          integrate(function(v){
            sapply(v,function(v){
              integrate(function(x){ 
                sapply(x,function(x){ 
                  integrate(function(y){ 
                    sapply(y,function(y){
                      (fmax(y,-log(hr1)/sqrt((x^2/c)*(hr[ec]/hr[1])),-log(hr2)/sqrt((x^2/c)*(hr[ec]/hr[2])),1,1,rho)
                       *fmax(x,-log(hr1),-log(hr2),sqrt(var1),sqrt(var2),rho)*dbivanorm(u,v,hr1,hr2,vartrue1,vartrue2,rho))
                    })
                  },qnorm(1-alpha)-log(step1)/sqrt((x^2/c)),qnorm(1-alpha)-log(step2)/sqrt((x^2/c)))$value
                })
              },-log(HRgo),Inf)$value
            })
          },-1000,1000)$value
        })
      },-1000,1000)$value)
     
   }
}


#' Probability that endpoint OS significant
#' 
#' This function calculate the probability that the endpoint OS is statistically significant.
#' In the context of cancer research OS stands for overall survival, a positive treatment effect in this endpoints is thus sufficient for a successful program.
#' @param HRgo threshold value for the go/no-go decision rule; 
#' @param n2 total sample size for phase II; must be even number
#' @param beta 1-beta power for calculation of the number of events for phase III by Schoenfeld (1981) formula
#' @param alpha one- sided significance level
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param id1 amount of information for `hr1` in terms of number of events
#' @param id2 amount of information for `hr2` in terms of number of events
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `hr1` is used as fixed effect
#' @param rho correlation between the two endpoints
#' @return The output of the the function `os_tte()` is the probability that endpoint OS significant.
#' @examples res <- os_tte(HRgo = 0.8, n2 = 50, alpha = 0.05, beta = 0.1,
#'                                hr1 = 0.75, hr2 = 0.80, 
#'                                id1 = 300, id2 = 600, 
#'                                fixed = TRUE, rho = 0.3)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export
os_tte<-function(HRgo, n2, alpha, beta, hr1, hr2, id1, id2, fixed, rho){
  
  e21<-hr1*n2 # number of events phase II for endpoint 1 (PFS) = event rate * sample size
  e22<-hr2*n2 # number of events phase II for endpoint 2 (OS)
  hr <- c(hr1,hr2)
  var1<-4/e21
  var2<-4/e22
  vartrue1 <- sqrt(4/id1)
  vartrue2 <- sqrt(4/id2)
  covmat<-matrix(c(var1, rho*sqrt(var1)*sqrt(var2), rho*sqrt(var1)*sqrt(var2), var2), ncol=2) #covariance-Matrix of c(true1,true2) 
  
  
  c     = (qnorm(1-alpha)+qnorm(1-beta))^2

      if(fixed) {
        
        os1_tte<-integrate(function(x){ 
        sapply(x,function(x){ 
          pnorm(-qnorm(1-alpha)-log(hr2)/sqrt((x^2/c)*(hr[1]/hr[2])))*fmax(x,-log(hr1),-log(hr2),sqrt(var1),sqrt(var2),rho)
        })
      },-Inf,Inf)$value

      os2_tte<-integrate(function(x){ 
        sapply(x,function(x){ 
          pnorm(-qnorm(1-alpha)-log(hr2)/sqrt((x^2/c)*(hr[2]/hr[2])))*fmax(x,-log(hr1),-log(hr2),sqrt(var1),sqrt(var2),rho)
        })
      },-Inf,Inf)$value}

  else {
  
      os1_tte<-integrate(function(u){
      sapply(u,function(u){
        integrate(function(v){
          sapply(v,function(v){
            integrate(function(x){ 
              sapply(x,function(x){ 
                pnorm(-qnorm(1-alpha)-log(hr2)/sqrt((x^2/c)*(hr[1]/hr[2])))*fmax(x,u,v,sqrt(var1),sqrt(var2),rho)*dbivanorm(u,v,-log(hr1),-log(hr2),vartrue1,vartrue2,rho)
              })
            },-Inf,Inf)$value
          })
        },-1000,1000)$value
      })
    },-1000,1000)$value

    os2_tte<-integrate(function(u){
      sapply(u,function(u){
        integrate(function(v){
          sapply(v,function(v){
            integrate(function(x){ 
              sapply(x,function(x){ 
                pnorm(-qnorm(1-alpha)-log(hr2)/sqrt((x^2/c)*(hr[2]/hr[2])))*fmax(x,u,v,sqrt(var1),sqrt(var2),rho)*dbivanorm(u,v,-log(hr1),-log(hr2),vartrue1,vartrue2,rho)
              })
            },-Inf,Inf)$value
          })
        },-1000,1000)$value
      })
    },-1000,1000)$value}

    return(os_tte <- os1_tte*pw(n2,hr1,hr2,id1,id2,fixed,rho) + os2_tte*(1-pw(n2,hr1,hr2,id1,id2,fixed,rho)))
}



#' Utility function for multiple endpoints in a time-to-event-setting
#'
#' The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters and the expected probability of a successful program. 
#' The utility is in a further step maximized by the `optimal_multiple_tte()` function. 
#' Note, that for calculating the utility of the program, two different benefit triples are necessary: 
#'  - one triple for the case that the more important endpoint overall survival (OS) shows a significant positive treatment effect 
#'  - one triple when only the endpoint progression-free survival (PFS) shows a significant positive treatment effect
#' @param HRgo threshold value for the go/no-go decision rule; 
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param id1 amount of information for `hr1` in terms of sample size
#' @param id2 amount of information for `hr2` in terms of sample size
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param steps1 lower boundary for effect size category `"small"` in HR scale, default: 1
#' @param stepm1 lower boundary for effect size category `"medium"` in HR scale = upper boundary for effect size category `"small"` in HR scale, default: 0.95
#' @param stepl1 lower boundary for effect size category `"large"` in HR scale = upper boundary for effect size category `"medium"` in HR scale, default: 0.85
#' @param b11 expected gain for effect size category `"small"` if endpoint OS is significant
#' @param b21 expected gain for effect size category `"medium"`if endpoint OS is significant
#' @param b31 expected gain for effect size category `"large"` if endpoint OS is significant 
#' @param b12 expected gain for effect size category `"small"` if endpoint OS is not significant
#' @param b22 expected gain for effect size category `"medium"`if endpoint OS is not significant
#' @param b32 expected gain for effect size category `"large"` if endpoint OS is not significant
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `hr1` is used as fixed effect
#' @param rho correlation between the two endpoints
#' @return The output of the the function `utility_multiple_tte()` is the expected utility of the program.
#' @examples res <- utility_multiple_tte(n2 = 50, HRgo = 0.8, alpha = 0.025, beta = 0.1,
#'                                hr1 = 0.75, hr2 = 0.80,
#'                                id1 = 300, id2 = 600,
#'                                c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,               
#'                                K = Inf, N = Inf, S = -Inf, 
#'                                steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
#'                                b11 = 1000, b21 = 2000, b31 = 3000,
#'                                b12 = 1000, b22 = 1500, b32 = 2000, 
#'                                fixed = TRUE, rho = 0.3)
#' @editor Johannes Cepicka
#' @editDate 2022-04-23  
#' @export
utility_multiple_tte<-function(n2, HRgo, alpha, beta, hr1, hr2, id1, id2,
                               c2, c02, c3, c03, K, N, S,
                               steps1, stepm1, stepl1, 
                               b11, b21, b31, b12, b22, b32, fixed, rho){ 
  
  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 0

   n3 <- Ess_multiple_tte(HRgo=HRgo,n2=n2,
                          alpha=alpha,beta=beta,
                          hr1=hr1,hr2=hr2,id1=id1,id2=id2,
                          fixed=fixed,rho=rho)
   
   OS <- os_tte(HRgo = HRgo, n2 = n2,
                alpha = alpha, beta = beta,
                hr1 = hr1, hr2 = hr2,
                id1 = id1, id2 = id2,
                fixed = fixed, rho = rho)
  
   pw <- pw(n2=n2,hr1=hr1,hr2=hr2,id1=id1,id2=id2,fixed=fixed,rho=rho)
   
   if(n2+n3>N){
     
     return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
     
   }else{
     
     pnogo   = pgo_multiple_tte(HRgo=HRgo,n2=n2,hr1=hr1,hr2=hr2,id1=id1,id2=id2,fixed=fixed,rho=rho)
     
     K2    <-  c02 + c2 * n2  #cost phase II
     K3    <-  c03 * (1-pnogo) + c3 * n3  #cost phase III
     
     if(K2+K3>K){
       
       return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
       #output: expected utility Eud, En3, EsP, Epgo
       
     }else{
       # probability of a successful program; small, medium, large effect size
       prob11 = EPsProg_multiple_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=1,hr1=hr1,hr2=hr2,id1=id1,id2=id2,
                                     step1=steps1,step2=steps2,fixed=fixed,rho=rho)
       prob21 = EPsProg_multiple_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=1,hr1=hr1,hr2=hr2,id1=id1,id2=id2,
                                     step1=stepm1,step2=stepm2,fixed=fixed,rho=rho)
       prob31 = EPsProg_multiple_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=1,hr1=hr1,hr2=hr2,id1=id1,id2=id2,
                                     step1=stepl1,step2=stepl2,fixed=fixed,rho=rho)
       
       prob12 = EPsProg_multiple_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=2,hr1=hr1,hr2=hr2,id1=id1,id2=id2,
                                     step1=steps1,step2=steps2,fixed=fixed,rho=rho)
       prob22 = EPsProg_multiple_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=2,hr1=hr1,hr2=hr2,id1=id1,id2=id2,
                                     step1=stepm1,step2=stepm2,fixed=fixed,rho=rho)
       prob32 = EPsProg_multiple_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=2,hr1=hr1,hr2=hr2,id1=id1,id2=id2,
                                     step1=stepl1,step2=stepl2,fixed=fixed,rho=rho)
       
       
        prob1 <- prob11*pw + prob12*(1-pw)
        prob2 <- prob21*pw + prob22*(1-pw)
        prob3 <- prob31*pw + prob32*(1-pw)
       
        SP    = (prob11+prob21+prob31)*pw + (prob12+prob22+prob32)*(1-pw)                         # probability of a successful program
          
       
       if(SP<S){
         
         return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
         
       }else{
         
         G     = (b11*prob1+b21*prob2+b31*prob3)*OS + (b12*prob1+b22*prob2+b32*prob3)*(1-OS)                 # gain
                             
         EU    = -K2-K3+G                                            # total expected utility
         
         SP2 = SP
         SP3 = 0
         
         return(c(EU, n3, SP, 1-pnogo, SP2, SP3, K2, K3, OS))
         
       }
     }
   }
   
}

