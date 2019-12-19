library("shiny")
library("ggplot2")
library("gridExtra")
library("mvtnorm")
library("doParallel")
library("parallel")
library("cubature")
library("plotly")

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
# -> Phase III is 2 or 3 arm trial (1:1 or 1:1:1 sample size allocatiob)

mainPath <- "/opt/shiny-server/samplesizr/multiarm/"

#mainPath <- ""


# probability to go to phase III
pgo<-function(HRgo,n2,ec,theta1,theta2,strategy,case){
   
   et1      = 1 - (1-ec)^exp(-theta1)     # event rate in arm 1
   et2      = 1 - (1-ec)^exp(-theta2)     # event rate in arm 2
   
   # distribution of y, yk~N(thetak,sigmak^2) and correlation rho = 1/2 (equal sample size allocation)
   MEANY    = c(theta1,theta2)
   sigma1   = sqrt((3/n2)*((1/ec)+(1/et1)))   # sd of y1 (equal sample size allocation)
   sigma2   = sqrt((3/n2)*((1/ec)+(1/et2)))   # sd of y2 (equal sample size allocation)
   SIGMAY   = matrix(c(sigma1^2,1/2*sigma1*sigma2,1/2*sigma1*sigma2,sigma2^2), nrow = 2, ncol = 2)
   
   if(case==1){# no go
      
      return(pmvnorm(lower = c(-Inf,-Inf),
                     upper = c(-log(HRgo),-log(HRgo)),
                     mean  = MEANY,
                     sigma = SIGMAY))
   }
   if(strategy==1){# best promising
      if(case==21){# treatment 1 is promising and better than treatment 2
         
         return(integrate(function(y1){
            sapply(y1,function(y1){ 
               integrate(function(y2){
                  dmvnorm(cbind(y1,y2),
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
                  dmvnorm(cbind(y1,y2),
                          mean  = MEANY,
                          sigma = SIGMAY)
               }, -Inf, y2)$value   
            })
         }, -log(HRgo), Inf)$value)   
         
      }
   }
   if(strategy==2){# all promising
      if(case==21){# treatment 1 is promising, treatment 2 is not
         
         return(pmvnorm(lower=c(-log(HRgo),-Inf),
                        upper=c(Inf,-log(HRgo)),
                        mean=MEANY,
                        sigma=SIGMAY))   
         
      }
      if(case==22){# treatment 2 is promising, treatment 1 is not
         return(pmvnorm(lower=c(-Inf,-log(HRgo)),
                        upper=c(-log(HRgo),Inf),
                        mean=MEANY,
                        sigma=SIGMAY))         
      }
      if(case==31){# both treatments are promising, treatment 1 is better
         
         return(integrate(function(y1){
            sapply(y1,function(y1){ 
               integrate(function(y2){
                  dmvnorm(cbind(y1,y2),
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
                  dmvnorm(cbind(y1,y2),
                          mean  = MEANY,
                          sigma = SIGMAY)
               }, -log(HRgo),y2)$value   
            })
         }, -log(HRgo), Inf)$value) 
         
      }
   }
   
}

# total sample size for phase III trial with l treatments and equal allocation ratio
# l=1: according to Schoenfeld to guarantee power for the log rank test to detect treatment effect of phase II; 
# l=2: according to Dunnett to guarantee y any-pair power (Horn & Vollandt)
ss<-function(alpha,beta,ec,ek,y,l){
   
   if(l==1){calpha = qnorm(1-alpha)}
   if(l==2){calpha = as.numeric(qmvnorm(1-alpha, mean=c(0,0), sigma=matrix(c(1,1/2,1/2,1), nrow=2, ncol=2))[1])}
   
   return(((l+1)*(calpha+qnorm(1-beta))^2)/(y^2)*((1/ec)+(1/ek)))
}

# Expected sample size for phase III when going to phase III
Ess<-function(HRgo,n2,alpha,beta,ec,theta1,theta2,strategy,case){
   
   et1      = 1 - (1-ec)^exp(-theta1)     # event rate in arm 1
   et2      = 1 - (1-ec)^exp(-theta2)     # event rate in arm 2
   
   # distribution of y, yk~N(thetak,sigmak^2) and correlation rho = 1/2 (equal sample size allocation)
   MEANY    = c(theta1,theta2)
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
                  ss(alpha,beta,ec,et1,y1,1)*
                     dmvnorm(cbind(y1,y2),
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
                  ss(alpha,beta,ec,et2,y2,1)*
                     dmvnorm(cbind(y1,y2),
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
            ss(alpha,beta,ec,et1,y[1],1)*dmvnorm(c(y[1],y[2]), mean  = MEANY, sigma = SIGMAY)
         }
         
         return(adaptIntegrate(f, lowerLimit = c(-log(HRgo), -Inf), upperLimit = c(Inf, -log(HRgo)))$integral)
         
      }
      if(case==22){# treatment 2 is promising, treatment 1 is not
         
         f <- function(y){ 
            ss(alpha,beta,ec,et2,y[2],1)*dmvnorm(c(y[1],y[2]), mean  = MEANY, sigma = SIGMAY)
         }
         
         return(adaptIntegrate(f, lowerLimit = c(-Inf, -log(HRgo)), upperLimit = c(-log(HRgo), Inf))$integral)
         
      }
      if(case==31){# both treatments are promising, treatment 1 is better
         
         return(integrate(function(y1){
            sapply(y1,function(y1){ 
               integrate(function(y2){
                  ss(alpha,beta,ec,et2,y2,2)*
                     dmvnorm(cbind(y1,y2),
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
                  ss(alpha,beta,ec,et1,y1,2)*
                     dmvnorm(cbind(y1,y2),
                             mean  = MEANY,
                             sigma = SIGMAY)
               }, -log(HRgo), y2)$value   
            })
         }, -log(HRgo), Inf)$value)  
         
      }
   }  
   
} 

# Probability of a successful program
PsProg<-function(HRgo,n2,alpha,beta,ec,theta1,theta2,step1,step2,strategy,case){
   
   et1      = 1 - (1-ec)^exp(-theta1)     # event rate in arm 1
   et2      = 1 - (1-ec)^exp(-theta2)     # event rate in arm 2
   
   # distribution of y, yk~N(thetak,sigmak^2) and correlation rho = 1/2 (equal sample size allocation)
   MEANY    = c(theta1,theta2)
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
                  ( pnorm(qnorm(1-alpha)+step2/(sqrt(y1^2/c)),
                          mean=theta1/(sqrt(y1^2/c)),
                          sd=1) -
                       pnorm(qnorm(1-alpha)+step1/(sqrt(y1^2/c)),
                             mean=theta1/(sqrt(y1^2/c)),
                             sd=1) )*
                     dmvnorm(cbind(y1,y2),
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
                  ( pnorm(qnorm(1-alpha)+step2/(sqrt(y2^2/c)),
                          mean=theta2/(sqrt(y2^2/c)),
                          sd=1) -
                       pnorm(qnorm(1-alpha)+step1/(sqrt(y2^2/c)),
                             mean=theta2/(sqrt(y2^2/c)),
                             sd=1) )*
                     dmvnorm(cbind(y1,y2),
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
            ( pnorm(qnorm(1-alpha)+step2/(sqrt(y[1]^2/c)),
                    mean=theta1/(sqrt(y[1]^2/c)),
                    sd=1)-
                 pnorm(qnorm(1-alpha)+step1/(sqrt(y[1]^2/c)),
                       mean=theta1/(sqrt(y[1]^2/c)),
                       sd=1) )*
               dmvnorm(c(y[1],y[2]),
                       mean=MEANY,
                       sigma=SIGMAY)
         }
         
         return(adaptIntegrate(f, lowerLimit = c(-log(HRgo), -Inf), upperLimit = c(Inf, -log(HRgo)))$integral)
         
      }
      if(case==22){# treatment 2 is promising, treatment 1 is not 
         
         c     = (qnorm(1-alpha)+qnorm(1-beta))^2 
         
         f <- function(y){ 
            ( pnorm(qnorm(1-alpha)+step2/(sqrt(y[2]^2/c)),
                    mean=theta2/(sqrt(y[2]^2/c)),
                    sd=1)-
                 pnorm(qnorm(1-alpha)+step1/(sqrt(y[2]^2/c)),
                       mean=theta2/(sqrt(y[2]^2/c)),
                       sd=1) )*
               dmvnorm(c(y[1],y[2]),
                       mean=MEANY,
                       sigma=SIGMAY)
         }
         
         return(adaptIntegrate(f, lowerLimit = c(-Inf, -log(HRgo)), upperLimit = c(-log(HRgo), Inf))$integral)
         
      }
      if(case==31){# both treatments are promising, treatment 1 is better
         
         SIGMAZ   = matrix(c(1,1/2,1/2,1), nrow = 2, ncol = 2)
         calpha   = as.numeric(qmvnorm(1-alpha, mean=c(0,0), sigma= SIGMAZ)[1])
         c        = (calpha+qnorm(1-beta))^2 
         
         return(integrate(function(y1){ 
            sapply(y1,function(y1){
               integrate(function(y2){ 
                  sapply(y2,function(y2){ # How to erase??
                     ( pmvnorm(lower=c(-Inf,-Inf),
                               upper=c(calpha+step2/(sqrt(y2^2/c)),
                                       calpha+step2/(sqrt(y2^2/c))),
                               mean=c(theta1/(sqrt(y2^2/c)),theta2/(sqrt(y2^2/c))),
                               sigma=SIGMAZ)-
                          pmvnorm(lower=c(-Inf,-Inf),
                                  upper=c(calpha+step1/(sqrt(y2^2/c)),
                                          calpha+step1/(sqrt(y2^2/c))),
                                  mean=c(theta1/(sqrt(y2^2/c)),theta2/(sqrt(y2^2/c))),
                                  sigma=SIGMAZ) )*
                        dmvnorm(c(y1,y2),
                                mean=MEANY,
                                sigma=SIGMAY) 
                  })
               }, -log(HRgo),y1)$value
            })
         }, -log(HRgo),Inf)$value)
         
         
      }
      if(case==32){# both treatments are promising, treatment 2 is better
         
         SIGMAZ   = matrix(c(1,1/2,1/2,1), nrow = 2, ncol = 2) 
         calpha   = as.numeric(qmvnorm(1-alpha, mean=c(0,0), sigma= SIGMAZ)[1])
         c        = (calpha+qnorm(1-beta))^2 
         
         return(integrate(function(y2){ 
            sapply(y2,function(y2){
               integrate(function(y1){ 
                  sapply(y1,function(y1){ # How to erase??
                     ( pmvnorm(lower=c(-Inf,-Inf),
                               upper=c(calpha+step2/(sqrt(y1^2/c)),
                                       calpha+step2/(sqrt(y1^2/c))),
                               mean=c(theta1/(sqrt(y1^2/c)),theta2/(sqrt(y1^2/c))),
                               sigma=SIGMAZ)-
                          pmvnorm(lower=c(-Inf,-Inf),
                                  upper=c(calpha+step1/(sqrt(y1^2/c)),
                                          calpha+step1/(sqrt(y1^2/c))),
                                  mean=c(theta1/(sqrt(y1^2/c)),theta2/(sqrt(y1^2/c))),
                                  sigma=SIGMAZ) )*
                        dmvnorm(c(y1,y2),
                                mean=MEANY,
                                sigma=SIGMAY) 
                  })
               }, -log(HRgo),y2)$value
            })
         }, -log(HRgo),Inf)$value)
         
      }
   }
   
} 

#utility function
utility<-function(HRgo,n2,alpha,beta,theta1,theta2,strategy,ec,c2,c02,c3,c03,b1,b2,b3,Nc){ 
   
   if(strategy==1){
      
      K2    = c02+c2*n2             # fixed and variable per-patient cost phase II 
      
      n321    = Ess(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                    strategy=strategy,case=21)
      n322    = Ess(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                    strategy=strategy,case=22)
      
      n3    <- ceiling(n321+n322)           # total expected sample size for phase III
    
      if(n2+n3>Nc){
         return(c(-9999, n2, n3, -9999, -9999, -9999, -9999, K2, -9999))
      }else{
         
         pnogo   = pgo(HRgo=HRgo,n2=n2,ec=ec,theta1=theta1,theta2=theta2,strategy=strategy,case=1)
         
         K3    = c03*(1-pnogo)+c3*n3   # fixed and variable per-patient cost phase III
         
         # probability of a successful program; small, medium, large effect size
         prob121 = PsProg(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                          step1=0,step2=-log(0.95),strategy=strategy,case=21)
         prob221 = PsProg(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                          step1=-log(0.95),step2=-log(0.85),strategy=strategy,case=21)
         prob321 = PsProg(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                          step1=-log(0.85),step2=Inf,strategy=strategy,case=21)
         
         prob122 = PsProg(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                          step1=0,step2=-log(0.95),strategy=strategy,case=22)
         prob222 = PsProg(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                          step1=-log(0.95),step2=-log(0.85),strategy=strategy,case=22)
         prob322 = PsProg(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                          step1=-log(0.85),step2=Inf,strategy=strategy,case=22)
         
         G     = b1*prob121+b2*prob221+b3*prob321 +                  # gain
            b1*prob122+b2*prob222+b3*prob322                   
         EU    = -K2-K3+G                                            # total expected utility
         SP    = prob121+prob221+prob321 +                           # probability of a successful program
            prob122+prob222+prob322 
         
         SP2 = SP
         SP3 = 0
         return(c(EU, n2, n3, SP, 1-pnogo, SP2, SP3, K2, K3))  
         # output: total expected utility, 
         # total sample size for phase II,
         # total expected sample size for phase III,
         # probability of a successful program,
         # probability to go to phase III  
         # probability of a successful two arm phase III trial
         # probability of a successful three arm phase III trial 
         # cost phase II
         # cost phase III
      }
   }
   if(strategy==2){
      
      K2      = c02+c2*n2             # fixed and variable per-patient cost phase II 
      
      n321    = Ess(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                    strategy=strategy,case=21)
      n322    = Ess(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                    strategy=strategy,case=22)
      n331    = Ess(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                    strategy=strategy,case=31)
      n332    = Ess(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                    strategy=strategy,case=32)
      n3      = ceiling(n321+n322+n331+n332)   # total expected sample size for phase III
      
      if(n2+n3>Nc){
         return(c(-9999, n2, n3, -9999, -9999, -9999, -9999, K2, -9999))
      }else{
         
         pnogo   = pgo(HRgo=HRgo,n2=n2,ec=ec,theta1=theta1,theta2=theta2,strategy=strategy,case=1)
         
         K3    = c03*(1-pnogo)+c3*n3   # fixed and variable per-patient cost phase III
         
         # probability of a successful program; small, medium, large effect size
         prob121 = PsProg(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                          step1=0,step2=-log(0.95),strategy=strategy,case=21)
         prob221 = PsProg(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                          step1=-log(0.95),step2=-log(0.85),strategy=strategy,case=21)
         prob321 = PsProg(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                          step1=-log(0.85),step2=Inf,strategy=strategy,case=21)
         
         prob122 = PsProg(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                          step1=0,step2=-log(0.95),strategy=strategy,case=22)
         prob222 = PsProg(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                          step1=-log(0.95),step2=-log(0.85),strategy=strategy,case=22)
         prob322 = PsProg(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                          step1=-log(0.85),step2=Inf,strategy=strategy,case=22)
         
         prob131 = PsProg(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                          step1=0,step2=-log(0.95),strategy=strategy,case=31)
         prob231 = PsProg(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                          step1=-log(0.95),step2=-log(0.85),strategy=strategy,case=31)
         prob331 = PsProg(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                          step1=-log(0.85),step2=Inf,strategy=strategy,case=31)
         
         prob132 = PsProg(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                          step1=0,step2=-log(0.95),strategy=strategy,case=32)
         prob232 = PsProg(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                          step1=-log(0.95),step2=-log(0.85),strategy=strategy,case=32)
         prob332 = PsProg(HRgo=HRgo,n2=n2,alpha=alpha,beta=beta,ec=ec,theta1=theta1,theta2=theta2,
                          step1=-log(0.85),step2=Inf,strategy=strategy,case=32)
         
         G       = b1*prob121+b2*prob221+b3*prob321 +                # gain
            b1*prob122+b2*prob222+b3*prob322 +
            b1*prob131+b2*prob231+b3*prob331 + 
            b1*prob132+b2*prob232+b3*prob332                   
         EU    = -K2-K3+G                                            # total expected utility
         SP2     = prob121+prob221+prob321 +                   # probability of a successful program with 
            prob122+prob222+prob322                            # two arms phase III trial 
         SP3     = prob131+prob231+prob331 +                   # with three arms phase III trial
            prob132+prob232+prob332
         SP      =  SP2 + SP3                                  # probability of a successful program
         
         return(c(EU, n2, n3, SP, 1-pnogo, SP2, SP3, K2, K3))  
         # output: total expected utility, 
         # total sample size for phase II,
         # total expected sample size for phase III,
         # probability of a successful program,
         # probability to go to phase III  
         # probability of a successful two arm phase III trial
         # probability of a successful three arm phase III trial 
         # cost phase II
         # cost phase III
      }
   }
}


shinyServer(function(input, output,session) {
   
     output$table <- renderTable({
        
         input$go
         Select = isolate(input$Select)
         input$go
         N2 = isolate(input$N2)
         input$go
         HRgo = isolate(input$HRgo)
         input$go
         stepN2 = isolate(input$stepN2)
         input$go
         stepHRgo = isolate(input$stepHRgo)
         input$go
         alpha = isolate(input$alpha)
         input$go
         beta = isolate(input$beta)
         input$go
         ec = isolate(input$ec)
         input$go
         R = isolate(input$R)
         input$go
         c02 = isolate(input$c02)
         input$go
         c03 = isolate(input$c03)
         input$go
         c2 = isolate(input$c2)
         input$go
         c3 = isolate(input$c3)
         input$go
         b1 = isolate(input$b1)
         input$go
         b2 = isolate(input$b2)
         input$go
         b3 = isolate(input$b3)
         input$go
         HR1 = isolate(input$HR1)
         input$go
         HR2 = isolate(input$HR2)
         input$go
         refresh = isolate(input$refresh)
         input$go
         Nc = isolate(input$Nc)
         
         theta1   = -log(HR1)
         theta2   = -log(HR2)
         
         N2seq = seq(from = N2[1], to = N2[2], by = stepN2)
         hrgo = seq(from = HRgo[1], to = HRgo[2], by = stepHRgo)
         
         if(Select==1){STRATEGY = 1}
         if(Select==2){STRATEGY = 2}
         if(Select==3){STRATEGY = c(1,2)}
         
         y = function(x){
            result = utility(HRgo,x,alpha,beta,theta1,theta2,strategy,ec,c2,c02,c3,c03,b1,b2,b3,Nc)
            return(result)
         }
         
         for(strategy in STRATEGY){
         
         ufkt <- n2fkt <- HRgofkt <- n3fkt<- spfkt <- pgofkt  <- sp2fkt <- sp3fkt <- 
            nfkt <- K2fkt <- K3fkt<- matrix(0,length(N2seq),length(hrgo))
         
         progress <- Progress$new(session, min=1, max=(length(hrgo)))
         on.exit(progress$close())
         
         progress$set(message = 'Optimization progress',
                      detail=paste("for Strategy",strategy))
      
         for(j in c(1:(length(hrgo)))){
            
            HRgo = hrgo[j]

            result <-sapply(N2seq,y)
            
            progress$set(value = j)
            
            #c(EU, n2, n3, SP, 1-pnogo, SP2, SP3, K2, K3)
            
            ufkt[,j]       <- result[1,]
            
            HRgofkt[,j]    <- rep(HRgo,length(N2seq))
            n2fkt[,j]      <- result[2,]
            n3fkt[,j]      <- result[3,]
            spfkt[,j]      <- result[4,]
            pgofkt[,j]     <- result[5,]
            sp2fkt[,j]     <- result[6,]
            sp3fkt[,j]     <- result[7,]
            K2fkt[,j]      <- result[8,]
            K3fkt[,j]      <- result[9,]

            nfkt[,j]       <- ceiling(n3fkt[,j])+N2seq
         }

         if(strategy==1){

            save(ufkt,file=paste0(mainPath, "ufkt1.RData"))
            save(n2fkt,file=paste0(mainPath, "n2fkt1.RData"))
            save(HRgofkt,file=paste0(mainPath, "HRgofkt1.RData"))   
         }
         if(strategy==2){

            save(ufkt,file=paste0(mainPath, "ufkt2.RData"))
            save(n2fkt,file=paste0(mainPath, "n2fkt2.RData"))
            save(HRgofkt,file=paste0(mainPath, "HRgofkt2.RData"))   
         }
         
         ind   <- which(ufkt == max(ufkt), arr.ind <- TRUE)
         
         I <- as.vector(ind[1,1])
         J <- as.vector(ind[1,2])
         
         Eud   <- ufkt[I,J]
         n2tot    <- n2fkt[I,J]
         n3tot    <- ceiling(n3fkt[I,J])
         ntot     <- nfkt[I,J]
         prob  <- spfkt[I,J]
         pg    <- pgofkt[I,J]
         prob2 <- sp2fkt[I,J]
         prob3 <- sp3fkt[I,J]
         k2    <- K2fkt[I,J]
         k3    <- K3fkt[I,J]
         
         
         if(refresh==0||(Select==3&strategy==2)){

            load(file=paste0(mainPath, "optimizationresults.RData"))
            
            }else{DF = NULL}
         
         DF <- rbind(DF,data.frame(Strategy=format(strategy,digits=0),u=Eud,HR1=exp(-theta1),HR2=exp(-theta2),ec=ec,
                                   HRgo=hrgo[J],n2=format(n2tot,digits=0),
                                   n3=format(n3tot,digits=0),n=format(ntot,digits=0),N=format(Nc,digits=0),
                                   pgo=pg,sProg=prob,Sp2=prob2,Sp3=prob3,K2=k2,K3=k3,
                                   alpha=format(alpha,digits=3),beta=beta,
                                   c02=c02,c03=c03,c2=c2,c3=c3,b1=b1,b2=b2,b3=b3))
         
         save(DF,file=paste0(mainPath, "optimizationresults.RData"))
         }
         return(DF)
               })

     output$plot <- renderPlotly({
           
           input$go
           Select = isolate(input$Select)
           input$go
           Plot = isolate(input$Plot)
           if(Plot==1){
              if(Select==1){

                bests = 1
                 load(file=paste0(mainPath, "ufkt1.RData"))
                 load(file=paste0(mainPath, "n2fkt1.RData"))
                 load(file=paste0(mainPath, "HRgofkt1.RData")) 
              }
              if(Select==2){

                bests = 2
                 load(file=paste0(mainPath, "ufkt2.RData"))
                 load(file=paste0(mainPath, "n2fkt2.RData"))
                 load(file=paste0(mainPath, "HRgofkt2.RData"))   
              }
              if(Select==3){
                 
                 load(file=paste0(mainPath, "ufkt1.RData"))
                 ufkt1  = ufkt
                 load(file=paste0(mainPath, "ufkt2.RData"))
                 ufkt2  = ufkt
                 
                 if(max(ufkt1)>max(ufkt2)){
                    
                   bests= 1
                    load(file=paste0(mainPath, "ufkt1.RData"))
                    load(file=paste0(mainPath, "n2fkt1.RData"))
                    load(file=paste0(mainPath, "HRgofkt1.RData"))  
                 }else{
                    
                   bests=2
                    load(file=paste0(mainPath, "ufkt2.RData"))
                    load(file=paste0(mainPath, "n2fkt2.RData"))
                    load(file=paste0(mainPath, "HRgofkt2.RData"))   
                 }
                 
              }
              plot_ly(x=HRgofkt,y=n2fkt,z=ufkt, type="surface")   %>% layout(title=
                      paste0("Optimization region of Strategy ", bests),
                       scene = list(
                         xaxis = list(title = "HRgo"),
                         yaxis = list(title = "n2"),
                         zaxis = list(title = "expected utility")))
              
              
           }

           
        })   
        
     
})
