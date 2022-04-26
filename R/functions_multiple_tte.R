################################################################################
## R-Code for maximizing expected utility
# Utility-based optimization of phase II/III programs with two time-to-event endpoints
# endpoints: overall survival (OS) and progression free survival (PFS)
#
# Author: Marietta Kirchner
# Date: 25.01.2017
################################################################################

#load packages
library(mvtnorm)
library(MASS)

#load functions

#density for maximum of two normally distributed, random variables
#Z=max(X,Y) with X~N(mu1,sigma1^2),Y~N(mu2,sigma2^2)
#f(z)=f1(-z)+f2(-z)
#function fmax will return the value of f(z)
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

#density of bivariate normal distribution
dbivanorm <- function(x,y, mu1,mu2,sigma1,sigma2,rho){ 
  covariancemat <- matrix(c(sigma1,rho*sqrt(sigma1)*sqrt(sigma2), rho*sqrt(sigma1)*sqrt(sigma2),sigma2),ncol=2)
  ff <- dmvnorm(cbind(x,y), mean=c(mu1,mu2),sigma=covariancemat)
  return(ff)
}


# probability to go to phase III
pgo_tte<-function(HRgo,n2,ec,hr1,hr2,id1,id2,fixed,rho){
  
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
          },kappa,Inf)$value
        })
      },-Inf,Inf )$value
    })
   },-Inf,Inf)$value)
 
 }

  
}


# expected sample size for phase III when going to phase III
Ess_tte<-function(HRgo,n2,alpha,beta,ec,hr1,hr2,id1,id2,fixed,rho){

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
          sapply(x,function(x)
            ((4*(qnorm(1 - alpha) + qnorm(1 - beta))^2)/x^2)*fmax(x,-log(hr1),-log(hr2),sqrt(var1),sqrt(var2),rho))
        },-log(HRgo),Inf,abs.tol=1e-2)$value)
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
                  },-log(HRgo),Inf)$value #devide value by pgo to get E(e3|GO)
                })
              },-Inf,Inf)$value
            })
          },-Inf,Inf)$value)
    
  }
}

pw <- function(n2,ec,hr1,hr2,id1,id2,fixed,rho){
  
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
expn3go_tte<-function(HRgo,n2,alpha,beta,ec,hr1,hr2,id1,id2,fixed,rho){
  
  expe3go_tte<-Ess_tte(HRgo,n2,alpha,beta,ec,hr1,hr2,id1,id2,fixed,rho)/pgo_tte(HRgo,n2,ec,hr1,hr2,id1,id2,fixed,rho)
  
  return(expe3go_tte/hr[1])*pw(n2,ec,hr1,hr2,id1,id2,fixed,rho)+(expe3go_tte/hr[2])*(1-pw(n2,ec,hr1,hr2,id1,id2,fixed,rho))
  }
  
  
  
# probability of a successful program

EPsProg_tte<-function(HRgo,n2,alpha,beta,ec,hr1,hr2,id1,id2,step1,step2,fixed,rho){
 
  e21<-hr1*n2 # number of events phase II for endpoint 1 (PFS) = event rate * sample size
  e22<-hr2*n2 # number of events phase II for endpoint 2 (OS)
  hr <- c(hr1,hr2)
  var1<-4/e21
  var2<-4/e22
  vartrue1 <- sqrt(4/id1)
  vartrue2 <- sqrt(4/id2)
  covmat<-matrix(c(var1, rho*sqrt(var1)*sqrt(var2), rho*sqrt(var1)*sqrt(var2), var2), ncol=2) #covariance-Matrix of c(true1,true2) 
  
  c     = (qnorm(1-alpha)+qnorm(1-beta))^2 
  
   if(fixed) {return(  integrate(function(x){ 
          sapply(x,function(x){ 
            integrate(function(y){ 
              sapply(y,function(y){
                fmax(y,hr1/sqrt((x^2/c)*(hr[ec]/hr[1])),hr2/sqrt((x^2/c)*(hr[ec]/hr[2])),1,1,rho)*fmax(x,-log(hr1),-log(hr2),sqrt(var1),sqrt(var2),rho)
              })
            },qnorm(1-alpha)+step1/sqrt((x^2/c)),qnorm(1-alpha)+step2/sqrt((x^2/c)))$value
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
                      (fmax(y,-log(hr1)/sqrt((x^2/c)*(hr[ec]/hr[1])),-log(hr2)/sqrt((x^2/c)*(hr[ec]/hr[2])),1,1,rho)*fmax(x,-log(hr1),-log(hr2),sqrt(var1),sqrt(var2),rho)*dbivanorm(u,v,hr1,hr2,vartrue1,vartrue2,rho))
                    })
                  },qnorm(1-alpha)+step1/sqrt((x^2/c)),qnorm(1-alpha)+step2/sqrt((x^2/c)))$value
                })
              },kappa,Inf)$value
            })
          },-Inf,Inf)$value
        })
      },-Inf,Inf)$value)
     
   }
}


#propability endpoint OS significant

os_tte<-function(HRgo,n2,alpha,beta,hr1,hr2,id1,id2,fixed,rho){
  
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
        },-Inf,Inf)$value
      })
    },-Inf,Inf)$value

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
        },-Inf,Inf)$value
      })
    },-Inf,Inf)$value}

    return(os_tte <- os1_tte*pw(n2,ec,hr1,hr2,id1,id2,fixed,rho) + os2_tte*(1-pw(n2,ec,hr1,hr2,id1,id2,fixed,rho)))
}


  

utility_multiple_tte<-function(n2,HRgo,alpha,beta,hr1,hr2,id1,id2,ec,
                               c2,c02,c3,c03,K,N,S,
                               steps1, stepm1, stepl1,b1, b2, b3,fixed,rho){ 

   n3 <- Ess_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,hr1=hr1,hr2=hr2,id1=id1,id2=id2,fixed=fixed,rho=rho)
   OS <- os_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,hr1=hr1,hr2=hr2,id1=id1,id2=id2,fixed=fixed,rho=rho)
   pw <- pw(n2=n2,ec=ec,hr1=hr1,hr2=hr2,id1=id1,id2=id2,fixed=fixed,rho=rho)
   
   if(n2+n3>N){
     
     return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
     
   }else{
     
     pnogo   = pgo_tte(HRgo=HRgo,n2=n2,ec=ec,hr1=hr1,hr2=hr2,id1=id1,id2=id2,fixed=fixed,rho=rho)
     
     K2    <-  c02 + c2 * n2  #cost phase II
     K3    <-  c03 * (1-pnogo) + c3 * n3  #cost phase III
     
     if(K2+K3>K){
       
       return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
       #output: expected utility Eud, En3, EsP, Epgo
       
     }else{
       # probability of a successful program; small, medium, large effect size
       prob11 = EPsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=1,hr1=hr1,hr2=hr2,id1=id1,id2=id2,
                        step1=steps1,step2=steps2,fixed=fixed,rho=rho)
       prob21 = EPsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=1,hr1=hr1,hr2=hr2,id1=id1,id2=id2,
                        step1=stepm1,step2=stepm2,fixed=fixed,rho=rho)
       prob31 = EPsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=1,hr1=hr1,hr2=hr2,id1=id1,id2=id2,
                        step1=stepl1,step2=stepl2,fixed=fixed,rho=rho)
       
       prob12 = EPsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=2,hr1=hr1,hr2=hr2,id1=id1,id2=id2,
                           step1=steps1,step2=steps2,fixed=fixed,rho=rho)
       prob22 = EPsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=2,hr1=hr1,hr2=hr2,id1=id1,id2=id2,
                           step1=stepm1,step2=stepm2,fixed=fixed,rho=rho)
       prob32 = EPsProg_tte(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=2,hr1=hr1,hr2=hr2,id1=id1,id2=id2,
                           step1=stepl1,step2=stepl2,fixed=fixed,rho=rho)
       
       
        prob1 <- prob11*pw + prob12(1-pw)
        prob2 <- prob21*pw + prob22(1-pw)
        prob3 <- prob31*pw + prob32(1-pw)
       
        SP    = (prob11+prob21+prob31)*pw + (prob12+prob22+prob32)*(1-pw)                         # probability of a successful program
          
       
       if(SP<S){
         
         return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
         
       }else{
         
         G     = b1*prob1+b2*prob2+b3*prob3                   # gain
                             
         EU    = -K2-K3+G                                            # total expected utility
         
         SP2 = SP
         SP3 = 0
         
         return(c(EU, n3, SP, 1-pnogo, SP2, SP3, K2, K3, OS))
         
       }
     }
   }
   
}
   
