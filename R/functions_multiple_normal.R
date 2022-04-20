############################################################################################
## R-Code for maximizing expected utility
# Utility-based optimization of phase II/III programs with two endpoints (continuous, normally distributed)
# Success: both endpoints significant in Phase III
# Two normally distributed endpoints
# endpoints: ADAS-cog, ADCS-ADL (Alzheimer's Disease)
#
# Author: Marietta Kirchner
# Date: 13.12.2016
############################################################################################

#load packages
library(mvtnorm)
library(MASS)
library(doParallel)
library(parallel)

#load functions
#density for minimum of two normally distributed, random variables
#Z=min(X,Y) with X~N(mu1,sigma1^2),Y~N(mu2,sigma2^2)
#f(z)=f1(z)+f2(z)
#function fmin will return the value of f(z)
fmin<-function (y,mu1,mu2,sigma1,sigma2,rho)
{t1<-dnorm(y,mean=mu1,sd=sigma1)
tt<-rho*(y-mu1)/(sigma1*sqrt(1-rho*rho))
tt<-tt-(y-mu2)/(sigma2*sqrt(1-rho*rho))
t1<-t1*pnorm(tt)
t2<-dnorm(y,mean=mu2,sd=sigma2)
tt<-rho*(y-mu2)/(sigma2*sqrt(1-rho*rho))
tt<-tt-(y-mu1)/(sigma1*sqrt(1-rho*rho))
t2<-t2*pnorm(tt)
return(t1+t2)}

#density of bivariate normal distribution
dbivanorm <- function(x,y, mu1,mu2,sigma1,sigma2,rho){ 
  covariancemat <- matrix(c(sigma1,rho*sqrt(sigma1)*sqrt(sigma2), rho*sqrt(sigma1)*sqrt(sigma2),sigma2),ncol=2)
  ff <- dmvnorm(cbind(x,y), mean=c(mu1,mu2),sigma=covariancemat)
  return(ff)
}

pgo_normal<-function(kappa, n2, Delta1, Delta2, in1, in2, sigma1, sigma2, fixed, rho){
  
  Sigma <- c(sigma1,sigma2)
  r<-c(4*Sigma[1]^2,4*Sigma[2]^2) #(r1,r2) known constant for endpoint i
  var1<-r[1]/n2 #variance of effect for endpoint 1
  var2<-r[2]/n2 #variance of effect for endpoint 2
  covmat<-matrix(c(var1, rho*sqrt(var1)*sqrt(var2), rho*sqrt(var1)*sqrt(var2), var2), ncol=2) #covariance-Matrix of c(true1,true2)
  
  
   if(fixed) {
    return(pmvnorm(lower=kappa, upper=c(Inf,Inf), mean=c(Delta1, Delta2),sigma=covmat)[1])
  }
  
  else  {
    return(integrate(function(u){
      sapply(u,function(u){
        integrate(function(v){
          sapply(v,function(v){
            (pmvnorm(lower=kappa, upper=c(Inf,Inf), mean=c(u,v),sigma=covmat)[1])*dbivanorm(u,v,Delta1,Delta2, 4/in1, 4/in2, rho)
          })
        },-Inf,Inf )$value
      })
    },-Inf,Inf )$value)
    
     }
  
}


Ess_normal<-function(kappa, n2, alpha, beta, Delta1, Delta2, in1, in2, sigma1, sigma2, fixed, rho){
  
  Sigma <- c(sigma1,sigma2)
  r<-c(4*Sigma[1]^2,4*Sigma[2]^2) #(r1,r2) known constant for endpoint i
  var1<-r[1]/n2 #variance of effect for endpoint 1
  var2<-r[2]/n2 #variance of effect for endpoint 2
  covmat<-matrix(c(var1, rho*sqrt(var1)*sqrt(var2), rho*sqrt(var1)*sqrt(var2), var2), ncol=2) #covariance-Matrix of c(true1,true2)
  
  c <- qnorm(1-alpha)+qnorm(1-beta)
 
  
   if(fixed)  {
    return(integrate(function(x){ sapply(x,function(x){
      (4*(qnorm(1-alpha)+qnorm(1-beta))^2/x^2)*fmin(x,Delta1,Delta2,var1,var2,rho)
    })
    },kappa,Inf)$value)
  }
  
  else   {
    return(integrate(function(u){ sapply(u,function(u){
      integrate(function(v){sapply(v,function(v){
        integrate(function(y){ sapply(y,function(y){ 
          integrate(function(x){ sapply(x,function(x){
            max((r[1]*(qnorm(1-alpha)+qnorm(1-beta))^2/x^2),(r[2]*(qnorm(1-alpha)+qnorm(1-beta))^2/y^2))*dbivanorm(x,y,u,v,var1,var2,rho)*dbivanorm(u,v,Delta1,Delta2,4/in1, 4/in2,rho)
          })
          },kappa[1],Inf)$value
        })
        },kappa[2],Inf)$value
      })
      },-Inf,Inf )$value
    })
    },-Inf,Inf )$value)
    
  }
}

posp_normal <- function(kappa, n2, alpha, beta, Delta1,Delta2, sigma1, sigma2, in1, in2, fixed, rho){
  
  Sigma <- c(sigma1,sigma2)
  r<-c(4*Sigma[1]^2,4*Sigma[2]^2) #(r1,r2) known constant for endpoint i
  var1<-r[1]/n2 #variance of effect for endpoint 1
  var2<-r[2]/n2 #variance of effect for endpoint 2
  covmat<-matrix(c(var1, rho*sqrt(var1)*sqrt(var2), rho*sqrt(var1)*sqrt(var2), var2), ncol=2) #covariance-Matrix of c(true1,true2)
  covmatt3<-matrix(c(1, rho, rho, 1), ncol=2)
  
  c <- qnorm(1-alpha)+qnorm(1-beta)
  
  if(fixed) {
    return(integrate(function(y){ 
      sapply(y,function(y){ 
        integrate(function(x){ 
          sapply(x,function(x){
            pmvnorm(lower=c(qnorm(1-alpha),qnorm(1-alpha)), 
                    upper=c(Inf,Inf), 
                    mean=c(Delta1/sqrt(r[1]/max((r[1]*c/x^2),(r[2]*c/y^2))),Delta2/sqrt(r[2]/max((r[1]*c/x^2),(r[2]*c/y^2)))),
                    sigma=covmatt3)[1]*dbivanorm(x,y,Delta1,Delta2,var1,var2,rho)
          })
        }, kappa[1],Inf)$value #x
      })
    }, kappa[2],Inf)$value)
  }
  
  else {
    return (integrate(function(u){ sapply(u,function(u){
      integrate(function(v){sapply(v,function(v){
        integrate(function(y){ sapply(y,function(y){ 
          integrate(function(x){ sapply(x,function(x){
            pmvnorm(lower=c(qnorm(1-alpha),qnorm(1-alpha)),
                    upper=c(Inf,Inf),
                    mean=c(u/sqrt(r[1]/max((r[1]*c/x^2),(r[2]*c/y^2))),v/sqrt(r[2]/max((r[1]*c/x^2),(r[2]*c/y^2)))),
                    sigma=covmatt3)[1]*dbivanorm(x,y,u,v,var1,var2,rho)*dbivanorm(u,v,Delta1,Delta2,4/in1,4/in2,rho)
          })
          }, kappa[1],Inf)$value #x
        })
        }, kappa[2],Inf)$value #y
      })
      },-Inf,Inf )$value
    })
    },-Inf,Inf )$value)
  }
  
}

#E(n3|GO)
  
  expn3go_normal<-Ess_normal(kappa, n2, alpha, beta, Delta1, Delta2, in1, in2, sigma1, sigma2, fixed, rho)/pgo_normal(kappa, n2, Delta1, Delta2, in1, in2, sigma1, sigma2, fixed, rho)
  




# probability of a successful program

EPsProg_normal<-function(kappa,n2,alpha,beta,Delta1,Delta2, sigma1, sigma2,
                      step11, step12, step21, step22, 
                      in1, in2, fixed,rho){
  
  Sigma <- c(sigma1,sigma2)
  r<-c(4*Sigma[1]^2,4*Sigma[2]^2) #(r1,r2) known constant for endpoint i
  var1<-r[1]/n2 #variance of effect for endpoint 1
  var2<-r[2]/n2 #variance of effect for endpoint 2
  covmat<-matrix(c(var1, rho*sqrt(var1)*sqrt(var2), rho*sqrt(var1)*sqrt(var2), var2), ncol=2) #covariance-Matrix of c(true1,true2)
  covmatt3<-matrix(c(1, rho, rho, 1), ncol=2)
  
  
  c     = (qnorm(1-alpha)+qnorm(1-beta))^2 
  
  if(fixed) {
    return(integrate(function(y){
      sapply(y,function(y){
        integrate(function(x){ 
          sapply(x,function(x){ 
            pmvnorm(lower=c(qnorm(1-alpha)+step11*sqrt(max((r[1]*c/x^2),(r[2]*c/y^2)))/2,
                            qnorm(1-alpha)+step12*sqrt(max((r[1]*c/x^2),(r[2]*c/y^2)))/2), 
                    upper=c(qnorm(1-alpha)+step21*sqrt(max((r[1]*c/x^2),(r[2]*c/y^2)))/2,
                            qnorm(1-alpha)+step22*sqrt(max((r[1]*c/x^2),(r[2]*c/y^2)))/2), 
                    mean=c(Delta1/sqrt(r[1]/max((r[1]*c/x^2),(r[2]*c/y^2))),Delta2/sqrt(r[2]/max((r[1]*c/x^2),(r[2]*c/y^2)))),
                    sigma=covmatt3)[1]*dbivanorm(x,y,Delta1,Delta2,var1,var2,rho)
          })
        },kappa[1],Inf)$value
      })
    },kappa[2],Inf)$value)
  }
  
  
  else {
    
    return(integrate(function(u){ sapply(u,function(u){
      integrate(function(v){sapply(v,function(v){
        integrate(function(y){ sapply(y,function(y){ 
          integrate(function(x){ sapply(x,function(x){
            pmvnorm(lower=c(qnorm(1-alpha)+step11*sqrt(max((r[1]*c/x^2),(r[2]*c/y^2)))/2,
                            qnorm(1-alpha)+step12*sqrt(max((r[1]*c/x^2),(r[2]*c/y^2)))/2), 
                    upper=c(qnorm(1-alpha)+step21*sqrt(max((r[1]*c/x^2),(r[2]*c/y^2)))/2,
                            qnorm(1-alpha)+step22*sqrt(max((r[1]*c/x^2),(r[2]*c/y^2)))/2), 
                    mean=c(u/sqrt(r[1]/max((r[1]*c/x^2),(r[2]*c/y^2))),v/sqrt(r[2]/max((r[1]*c/x^2),(r[2]*c/y^2)))),
                    sigma=covmatt3)[1]*dbivanorm(x,y,u,v,var1,var2,rho)*dbivanorm(u,v,Delta1,Delta2,4/in1,4/in2,rho)
          })
          },kappa[1],Inf)$value
        })
        },kappa[2],Inf)$value
      })
      },-Inf,Inf )$value
    })
    },-Inf,Inf )$value)
    
  }
}



utility_multiple_normal<-function(n2,kappa,alpha,beta,Delta1,Delta2, in1, in2, sigma1, sigma2,
                               c2,c02,c3,c03,K,N,S,
                               steps1, stepm1, stepl1,b1, b2, b3,fixed,rho,relaxed){ 
  
  n3 = Ess_normal(kappa = kappa, n2 = n2, alpha = alpha , beta = beta, 
                  Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, 
                  sigma1 = sigma1, sigma2=sigma2, fixed = fixed, rho = rho)
  
  POSP = posp_normal(kappa = kappa, n2=n2, alpha= alpha, beta=beta, 
                     Delta1 = Delta1, Delta2 = Delta2, in1, in2,
                     sigma1=sigma1, sigma2=sigma2, fixed=fixed, rho=rho)
 
  
  if(n2+n3>N){
    
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
    
  }else{
    
    pnogo   = pgo_normal(kappa=kappa, n2=n2, Delta1=Delta1, Delta2=Delta2, in1=in1, in2=in2, sigma1 = sigma1, sigma2 = sigma2, fixed=fixed, rho=rho)
    
    K2    <-  c02 + c2 * n2  #cost phase II
    K3    <-  c03 * (1-pnogo) + c3 * n3  #cost phase III
    
    if(K2+K3>K){
      
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
      #output: expected utility Eud, En3, EsP, Epgo
      
    }else{
      # probability of a successful program; small, medium, large effect size
      prob11 = EPsProg_normal(kappa = kappa ,n2 = n2, alpha=alpha, beta=beta, Delta1 = Delta1, Delta2=Delta2,
                           sigma1 = sigma1, sigma2 = sigma2,
                           step11 = steps1 , step12 = stepm1  , step21 = stepm1, step22 = Inf, 
                           in1 = in1, in2 = in2, fixed = fixed ,rho = rho)
      prob21 = EPsProg_normal(kappa = kappa ,n2 = n2, alpha=alpha, beta=beta, Delta1 = Delta1, Delta2=Delta2,
                              sigma1 = sigma1, sigma2 = sigma2,
                              step11 = stepm1, step12 = steps1 , step21 = Inf , step22 = stepm1,
                              in1 = in1, in2 = in2, fixed = fixed ,rho = rho)
      prob31 = EPsProg_normal(kappa = kappa ,n2 = n2, alpha=alpha, beta=beta, Delta1 = Delta1, Delta2=Delta2,
                              sigma1 = sigma1, sigma2 = sigma2,
                              step11 = steps1, step12 = steps1, step21 = stepm1, step22 = stepm1,
                              in1 = in1, in2 = in2, fixed = fixed ,rho = rho)
      
      proba1 <- prob11 + prob21 + prob31
      proba3 = EPsProg_normal(kappa = kappa ,n2 = n2, alpha=alpha, beta=beta, Delta1 = Delta1, Delta2=Delta2,
                              sigma1 = sigma1, sigma2 = sigma2,
                              step11 = stepl1, step12 = stepl1, step21 = Inf, step22 = Inf,
                              in1 = in1, in2 = in2, fixed = fixed ,rho = rho)
      proba2 <- POSP - proba1 - proba3
      
      prob13 = EPsProg_normal(kappa = kappa ,n2 = n2, alpha=alpha, beta=beta, Delta1 = Delta1, Delta2=Delta2,
                              sigma1 = sigma1, sigma2 = sigma2,
                              step11 = stepl1, step12 = steps1, step21 = Inf, step22 = stepl1,
                              in1 = in1, in2 = in2, fixed = fixed ,rho = rho)
      prob23 = EPsProg_normal(kappa = kappa ,n2 = n2, alpha=alpha, beta=beta, Delta1 = Delta1, Delta2=Delta2,
                              sigma1 = sigma1, sigma2 = sigma2,
                              step11 = steps1, step12 = stepl1, step21 = stepl1, step22 = Inf,
                              in1 = in1, in2 = in2, fixed = fixed ,rho = rho)
      prob33 = EPsProg_normal(kappa = kappa ,n2 = n2, alpha=alpha, beta=beta, Delta1 = Delta1, Delta2=Delta2,
                              sigma1 = sigma1, sigma2 = sigma2,
                              step11 = stepl1, step12 = stepl1, step21 = Inf, step22 = Inf,
                              in1 = in1, in2 = in2, fixed = fixed ,rho = rho)
      
      probb3 <- prob13 + prob23 + prob33
      probb1 = EPsProg_normal(kappa = kappa ,n2 = n2, alpha=alpha, beta=beta, Delta1 = Delta1, Delta2=Delta2,
                              sigma1 = sigma1, sigma2 = sigma2,
                              step11 = steps1, step12 = steps1, step21 = stepm1, step22 = stepm1,
                              in1 = in1, in2 = in2, fixed = fixed ,rho = rho)
      probb2 <- POSP - probb1 - probb3
      
      if (relaxed == "TRUE"){
        prob1 <- probb1
        prob2 <- probb2
        prob3 <- probb3
      } else {
        prob1 <- proba1
        prob2 <- proba2
        prob3 <- proba3
      }
      
      
      
       SP    = prob1 + prob2 + prob3                         # probability of a successful program
      
      
      if(SP<S){
        
        return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))
        
      }else{
        
        G     = b1*prob1+b2*prob2+b3*prob3                      # gain
        
        EU    = -K2-K3+G                                            # total expected utility
        
        SP2 = SP
        SP3 = 0
        
        return(c(EU, n3, SP, 1-pnogo, SP2, SP3, K2, K3))
        
      }
    }
  }
  
}

