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
# theta; fithetaed true effects for treatment 1 and 2
# y: hat_theta_2; estimator in phase II
# z: hat_T_3; normalized estimator in phase III
# d2: total sample size in phase II 
# d3: total sample size in phase III

# 1. Strategy: Only best promising treatment goes to phase III
# -> Phase III is always 2 arm trial (1:1 sample size allocatiob)
# 2. Strategy: All promising treatments go to phase III
# -> Phase III is 2 or 3 arm trial (1:1 or 1:1:1 sample size allocatiob)

mainPath <- "/opt/shiny-server/samplesizr/multitrial/"

#mainPath <- ""


######################################
# One phase III trial - fixed effect #
######################################
# Ursprung: basic R Shiny App Code
# theta: true underlying effect
# y: hat_theta_2; estimator in phase II
# z: T_3; normalized estimator in phase III

# Probability to go to phase III: pgo
pgof <-  function(HRgo, d2, theta){
   
   pnorm((log(HRgo) + theta)/sqrt(4/d2))
}

# Ethetapected number of events for phase III when going to phase III: Ed3
Ed3f <-  function(HRgo, d2, alpha, beta, theta){
   
   ceiling(integrate(function(y){
      ( (4 * (qnorm(1 - alpha) + qnorm(1 - beta))^2)/(y^2)) * 
         dnorm(y, 
               mean = theta, 
               sd = sqrt(4/d2))  
   },  - log(HRgo), Inf)$value)
} 

# Ethetapected probability of a successful program: EsP
EPsProgf <-  function(HRgo, d2, alpha, beta, theta, step1, step2){
   
   c = (qnorm(1 - alpha) + qnorm(1 - beta))^2
   
   integrate(function(y){ 
      ( pnorm(qnorm(1 - alpha) + step2/(sqrt(y^2/c)), 
              mean = theta/(sqrt(y^2/c)), 
              sd = 1) - 
           pnorm(qnorm(1 - alpha) + step1/(sqrt(y^2/c)), 
                 mean = theta/(sqrt(y^2/c)), 
                 sd = 1) ) * 
         dnorm(y, 
               mean = theta, 
               sd = sqrt(4/d2)) 
   },  - log(HRgo), Inf)$value
   
   
}

# Utility function
utilityf <-  function(HRgo,d2,alpha,beta,theta,p2,p3,c2,c02,c3,c03,b1,b2,b3){ 
   
   pg    <-  pgof(HRgo = HRgo, d2 = d2, theta = theta)
   d3    <-  Ed3f(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta)
   
   K2    <-  c02 + c2 * d2 * (1/p2) #cost phase II 
   K3    <-  c03 * pg + c3 * d3 * (1/p3) #cost phase III
      
      steps1   <- 0
      steps2   <- - log(0.95)
      stepm1   <- - log(0.95)
      stepm2   <- - log(0.85)
      stepl1   <- - log(0.85)
      stepl2   <- Inf
      
      # probability of a successful program; small, medium, large effect size
      prob1 <-  EPsProgf(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta,
                         step1 = steps1, step2 =  steps2)
      prob2 <-  EPsProgf(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta,
                         step1 =  stepm1, step2 =  stepm2)
      prob3 <-  EPsProgf(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta,
                         step1 =  stepl1, step2 = stepl2)
      
      G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
      
      EU    <-  - K2 - K3 + G
      SP    <-  prob1 + prob2 + prob3
      
      return(c(EU, d2, d3, SP, pg, K2, K3, prob1, prob2, prob3)) 
      #output: ethetapected utility Eud, Ed3, EsP, Epgo, cost phase II and III   
   
}

# Ursprung: U:\Dissi\BiometricalJournal\austausch_preussler\Code zum Einreichen\Code\functions 

########################################
# Two phase III trials - fithetaed effects #
########################################
# Case 1: Strategy 1/2; at least one trial significant, the treatment effect of the other one at least showing in the same direction 
# Case 2: Strategy 2/2; both trials significant 

# Ethetapected probability of a successful program
EPsProgf2 <-  function(HRgo, d2, alpha, beta, theta, case, size){
   
   SIGMA <-  diag(2)
   c     <-  (qnorm(1 - alpha) + qnorm(1 - beta))^2
   
   if(case == 1){
      if(size == "small"){
         return(integrate(function(y){
            sapply(y, function(y){
               ( pmvnorm(lower = c(0, 
                                   0), 
                         upper = c(qnorm(1 - alpha) - log(0.95)/sqrt(y^2/c),
                                   qnorm(1 - alpha) - log(0.95)/sqrt(y^2/c)), 
                         mean = c(theta/sqrt(y^2/c), 
                                  theta/sqrt(y^2/c)), 
                         sigma = SIGMA) - 
                    pmvnorm(lower = c(0, 
                                      0), 
                            upper = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha)), 
                            mean = c(theta/sqrt(y^2/c), 
                                     theta/sqrt(y^2/c)), 
                            sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = theta, 
                        sd = sqrt(4/d2))  
            })
         },  - log(HRgo), Inf)$value)  
      }
      if(size == "large"){
         return(integrate(function(y){
            sapply(y, function(y){
               ( pmvnorm(lower = c(0, 
                                   0), 
                         upper = c(Inf, 
                                   Inf), 
                         mean = c(theta/sqrt(y^2/c), 
                                  theta/sqrt(y^2/c)), 
                         sigma = SIGMA) - 
                    pmvnorm(lower = c(0, 
                                      0), 
                            upper = c(qnorm(1 - alpha) - 
                                         log(0.85)/sqrt(y^2/c), 
                                      qnorm(1 - alpha) - 
                                         log(0.85)/sqrt(y^2/c)), 
                            mean = c(theta/sqrt(y^2/c), 
                                     theta/sqrt(y^2/c)), 
                            sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = theta, 
                        sd = sqrt(4/d2)) 
            })
         },  - log(HRgo), Inf)$value)  
      }
      if(size == "all"){
         return(integrate(function(y){
            sapply(y, function(y){
               ( pmvnorm(lower = c(0, 
                                   0), 
                         upper = c(Inf, 
                                   Inf), 
                         mean = c(theta/sqrt(y^2/c), 
                                  theta/sqrt(y^2/c)), 
                         sigma = SIGMA) - 
                    pmvnorm(lower = c(0, 
                                      0), 
                            upper = c(qnorm(1 - alpha), 
                                      qnorm(1 - alpha)), 
                            mean = c(theta/sqrt(y^2/c), 
                                     theta/sqrt(y^2/c)), 
                            sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = theta, 
                        sd = sqrt(4/d2))
            })
         },  - log(HRgo), Inf)$value)     
      }
   }
   if(case == 2){
      if(size == "small"){
         return(integrate(function(y){
            sapply(y, function(y){
               ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                   qnorm(1 - alpha)), 
                         upper = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c), 
                                   qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c)), 
                         mean = c(theta/sqrt(y^2/c), 
                                  theta/sqrt(y^2/c)), 
                         sigma = SIGMA) - 
                    pmvnorm(lower = c(qnorm(1 - alpha) - 
                                         log(0.95)/sqrt(y^2/c), 
                                      qnorm(1 - alpha) - 
                                         log(0.95)/sqrt(y^2/c)), 
                            upper = c(qnorm(1 - alpha) - 
                                         log(0.85)/sqrt(y^2/c), 
                                      qnorm(1 - alpha) - 
                                         log(0.85)/sqrt(y^2/c)), 
                            mean = c(theta/sqrt(y^2/c), 
                                     theta/sqrt(y^2/c)), 
                            sigma = SIGMA)) * 
                  dnorm(y, 
                        mean = theta, 
                        sd = sqrt(4/d2)) 
            })
         },  - log(HRgo), Inf)$value) 
      }
      if(size == "large"){
         return(integrate(function(y){
            sapply(y, function(y){
               ( pmvnorm(lower = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c), 
                                   qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c)), 
                         upper = c(Inf, 
                                   Inf), 
                         mean = c(theta/sqrt(y^2/c), 
                                  theta/sqrt(y^2/c)), 
                         sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = theta, 
                        sd = sqrt(4/d2)) 
            })
         },  - log(HRgo), Inf)$value)    
      }
      if(size == "all"){
         return(integrate(function(y){
            sapply(y, function(y){
               ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                   qnorm(1 - alpha)), 
                         upper = c(Inf, 
                                   Inf), 
                         mean = c(theta/sqrt(y^2/c), 
                                  theta/sqrt(y^2/c)), 
                         sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = theta, 
                        sd = sqrt(4/d2)) 
            })
         },  - log(HRgo), Inf)$value)    
      }
   }
   
}

# Utility function
utilityf2 <-  function(HRgo,d2,alpha,beta,theta,p2,p3,c2,c02,c3,c03,b1,b2,b3,case){ 
   
   pg  <-  pgof(HRgo = HRgo, d2 = d2, theta = theta)
   d3  <-  Ed3f(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta)
   
   # probability of a successful program: effect sizes
   
   # small 
   prob1 <-  EPsProgf2(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta,
                       case = case, size = "small") 
   # large
   prob3 <-  EPsProgf2(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta,
                       case = case, size = "large") 
   # medium
   prob2 <-  EPsProgf2(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta,
                       case = case, size = "all")  -  prob1  -  prob3
   
   K2    <-  c02 + c2 * d2 * (1/p2) # cost phase II 
   K3    <-  c03 * pg + c3 * d3 * (1/p3) # cost for one phase III
   G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 # gain
   
   EU    <-  - K2 - 2 * K3 + G
   SP    <-  prob1 + prob2 + prob3
   
   return(c(EU, d2, d3, SP, pg, K2, K3, prob1, prob2, prob3)) # output: expected utility Eud, Ed3, EsP, Epgo
   
}

##########################
# Three phase III trials #
##########################
# Case 2: Strategy 2/3; at least two trials significant, the treatment effect 
# of the other one at least showing in the same direction
# Case 3: Strategy 3/3; all trials significant

# Expected probability of a successful program
EPsProgf3 <-  function(HRgo, d2, alpha, beta, theta, case, size){
   
   SIGMA <-  diag(3)
   c     <-  (qnorm(1 - alpha) + qnorm(1 - beta))^2
   
   if(case == 2){
      if(size == "small"){
         return(integrate(function(y){
            sapply(y, function(y){
               ( 3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                       qnorm(1 - alpha), 
                                       0), 
                             upper = c(qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c), 
                                       qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(y^2/c), 
                                       qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(y^2/c)), 
                             mean = c(theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c)), 
                             sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                          qnorm(1 - alpha), 
                                          qnorm(1 - alpha)), 
                                upper = c(qnorm(1 - alpha) - 
                                             log(0.95)/sqrt(y^2/c), 
                                          qnorm(1 - alpha) - 
                                             log(0.95)/sqrt(y^2/c), 
                                          qnorm(1 - alpha) - 
                                             log(0.95)/sqrt(y^2/c)), 
                                mean = c(theta/sqrt(y^2/c), 
                                         theta/sqrt(y^2/c), 
                                         theta/sqrt(y^2/c)), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = theta, 
                        sd = sqrt(4/d2)) 
            })
         },  - log(HRgo), Inf)$value)
      }
      if(size == "large"){
         return(integrate(function(y){
            sapply(y, function(y){
               ( 3 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c), 
                                       qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c), 
                                       0), 
                             upper = c(Inf, 
                                       Inf, 
                                       Inf), 
                             mean = c(theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c)), 
                             sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                             log(0.85)/sqrt(y^2/c), 
                                          qnorm(1 - alpha) - 
                                             log(0.85)/sqrt(y^2/c), 
                                          qnorm(1 - alpha) - 
                                             log(0.85)/sqrt(y^2/c)), 
                                upper = c(Inf, 
                                          Inf, 
                                          Inf), 
                                mean = c(theta/sqrt(y^2/c), 
                                         theta/sqrt(y^2/c), 
                                         theta/sqrt(y^2/c)), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = theta, 
                        sd = sqrt(4/d2)) 
            })
         },  - log(HRgo), Inf)$value)
      }
      if(size == "all"){
         return(integrate(function(y){
            sapply(y, function(y){
               ( 3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                       qnorm(1 - alpha), 
                                       0), 
                             upper = c(Inf, 
                                       Inf, 
                                       Inf), 
                             mean = c(theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c)), 
                             sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                          qnorm(1 - alpha), 
                                          qnorm(1 - alpha)), 
                                upper = c(Inf, 
                                          Inf, 
                                          Inf), 
                                mean = c(theta/sqrt(y^2/c), 
                                         theta/sqrt(y^2/c), 
                                         theta/sqrt(y^2/c)), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = theta, 
                        sd = sqrt(4/d2))  
            })
         },  - log(HRgo), Inf)$value)
      }
   }
   if(case == 3){
      if(size == "small"){
         return(integrate(function(y){
            sapply(y, function(y){
               ( 3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                       qnorm(1 - alpha), 
                                       qnorm(1 - alpha)), 
                             upper = c(qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(y^2/c), 
                                       qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(y^2/c), 
                                       qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c)), 
                             mean = c(theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c)), 
                             sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                          qnorm(1 - alpha), 
                                          qnorm(1 - alpha)), 
                                upper = c(qnorm(1 - alpha) - 
                                             log(0.95)/sqrt(y^2/c), 
                                          qnorm(1 - alpha) - 
                                             log(0.95)/sqrt(y^2/c), 
                                          qnorm(1 - alpha) - 
                                             log(0.95)/sqrt(y^2/c)), 
                                mean = c(theta/sqrt(y^2/c), 
                                         theta/sqrt(y^2/c), 
                                         theta/sqrt(y^2/c)), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = theta, 
                        sd = sqrt(4/d2))
            })
         },  - log(HRgo), Inf)$value)
      }
      if(size == "large"){
         return(integrate(function(y){
            sapply(y, function(y){
               ( 3 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c), 
                                       qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c), 
                                       qnorm(1 - alpha)), 
                             upper = c(Inf, 
                                       Inf, 
                                       Inf), 
                             mean = c(theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c)), 
                             sigma = SIGMA) - 
                    2 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                             log(0.85)/sqrt(y^2/c), 
                                          qnorm(1 - alpha) - 
                                             log(0.85)/sqrt(y^2/c), 
                                          qnorm(1 - alpha) - 
                                             log(0.85)/sqrt(y^2/c)), 
                                upper = c(Inf, 
                                          Inf, 
                                          Inf), 
                                mean = c(theta/sqrt(y^2/c), 
                                         theta/sqrt(y^2/c), 
                                         theta/sqrt(y^2/c)), 
                                sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = theta, 
                        sd = sqrt(4/d2)) 
            })
         },  - log(HRgo), Inf)$value)
      }
      if(size == "all"){
         return(integrate(function(y){
            sapply(y, function(y){
               ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                   qnorm(1 - alpha), 
                                   qnorm(1 - alpha)), 
                         upper = c(Inf, 
                                   Inf, 
                                   Inf), 
                         mean = c(theta/sqrt(y^2/c), 
                                  theta/sqrt(y^2/c), 
                                  theta/sqrt(y^2/c)), 
                         sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = theta, 
                        sd = sqrt(4/d2)) 
            })
         },  - log(HRgo), Inf)$value)
      }
   }
   
}

# Utility function
utilityf3 <-  function(HRgo,d2,alpha,beta,theta,p2,p3,c2,c02,c3,c03,b1,b2,b3,case){ 
   
   pg  <-  pgof(HRgo = HRgo, d2 = d2, theta = theta)
   d3  <-  Ed3f(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta)
   
   # probability of a successful program: effect sizes
   
   # small 
   prob1 <-  EPsProgf3(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta,
                       case = case, size = "small") 
   # large
   prob3 <-  EPsProgf3(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta,
                       case = case, size = "large") 
   # medium
   prob2 <-  EPsProgf3(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta,
                       case = case, size = "all")  -  prob1  - prob3
   
   K2    <-  c02 + c2 * d2 * (1/p2) # cost phase II 
   K3    <-  c03 * pg + c3 * d3 * (1/p3) # cost for one phase III
   G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 # gain
   
   EU    <-  - K2 - 3 * K3 + G
   SP    <-  prob1 + prob2 + prob3
   
   return(c(EU, d2, d3, SP, pg, K2, K3, prob1, prob2, prob3)) # output: expected utility Eud, Ed3, EsP, Epgo 
   
}

#########################
# Four phase III trials #
#########################
# Case 3: Strategy 3/4; at least three trials significant, the treatment effect 
# of the other one at least showing in the same direction

# Expected probability of a successful program
EPsProgf4 <-  function(HRgo, d2, alpha, beta, theta, size){
   
   SIGMA <-  diag(4)
   c     <-  (qnorm(1 - alpha) + qnorm(1 - beta))^2
   
   if(size == "small"){
      return(integrate(function(y){
         sapply(y, function(y){
            ( 4 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    0), 
                          upper = c(qnorm(1 - alpha) - log(0.85)/sqrt(y^2/c), 
                                    qnorm(1 - alpha) - log(0.95)/sqrt(y^2/c), 
                                    qnorm(1 - alpha) - log(0.95)/sqrt(y^2/c), 
                                    qnorm(1 - alpha) - log(0.95)/sqrt(y^2/c)), 
                          mean = c(theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c)), 
                          sigma = SIGMA)  - 
                 3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                       qnorm(1 - alpha), 
                                       qnorm(1 - alpha), 
                                       qnorm(1 - alpha)), 
                             upper = c(qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(y^2/c), 
                                       qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(y^2/c), 
                                       qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(y^2/c), 
                                       qnorm(1 - alpha) - 
                                          log(0.95)/sqrt(y^2/c)), 
                             mean = c(theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c)), 
                             sigma = SIGMA) ) * 
               dnorm(y, 
                     mean = theta, 
                     sd = sqrt(4/d2)) 
         })
      },  - log(HRgo), Inf)$value)
   }
   if(size == "large"){
      return(integrate(function(y){
         sapply(y, function(y){
            ( 4 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                       log(0.85)/sqrt(y^2/c), 
                                    qnorm(1 - alpha) - log(0.85)/sqrt(y^2/c), 
                                    qnorm(1 - alpha) - log(0.85)/sqrt(y^2/c), 
                                    0), 
                          upper = c(Inf, 
                                    Inf, 
                                    Inf, 
                                    Inf), 
                          mean = c(theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c)), 
                          sigma = SIGMA)  - 
                 3 * pmvnorm(lower = c(qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c), 
                                       qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c), 
                                       qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c), 
                                       qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c)), 
                             upper = c(Inf, 
                                       Inf, 
                                       Inf, 
                                       Inf), 
                             mean = c(theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c)), 
                             sigma = SIGMA) ) * 
               dnorm(y, 
                     mean = theta, 
                     sd = sqrt(4/d2)) 
         })
      },  - log(HRgo), Inf)$value)
   }
   if(size == "all"){
      return(integrate(function(y){
         sapply(y, function(y){
            ( 4 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    qnorm(1 - alpha), 
                                    0), 
                          upper = c(Inf, 
                                    Inf, 
                                    Inf, 
                                    Inf), 
                          mean = c(theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c), 
                                   theta/sqrt(y^2/c)), 
                          sigma = SIGMA) - 
                 3 * pmvnorm(lower = c(qnorm(1 - alpha), 
                                       qnorm(1 - alpha), 
                                       qnorm(1 - alpha), 
                                       qnorm(1 - alpha)), 
                             upper = c(Inf, 
                                       Inf, 
                                       Inf, 
                                       Inf), 
                             mean = c(theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c)), 
                             sigma = SIGMA) ) * 
               dnorm(y, 
                     mean = theta, 
                     sd = sqrt(4/d2)) 
         })
      },  - log(HRgo), Inf)$value)
   }
   
}

# Utility function
utilityf4 <-  function(HRgo,d2,alpha,beta,theta,p2,p3,c2,c02,c3,c03,b1,b2,b3,case){ 
   
   pg  <-  pgof(HRgo = HRgo, d2 = d2, theta = theta)
   d3  <-  Ed3f(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta)
   
   # probability of a successful program: effect sizes
   
   # small 
   prob1 <-  EPsProgf4(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta,
                       size = "small") 
   # large
   prob3 <-  EPsProgf4(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta,
                       size = "large") 
   # medium
   prob2 <-  EPsProgf4(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta,
                       size = "all")  -  prob1  -  prob3
   
   K2    <-  c02 + c2 * d2 * (1/p2) # cost phase II 
   K3    <-  c03 * pg + c3 * d3 * (1/p3) # cost for one phase III
   G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 # gain
   
   EU    <-   - K2 - 4 * K3 + G
   SP    <-  prob1 + prob2 + prob3
   
   return(c(EU, d2, d3, SP, pg, K2, K3, prob1, prob2, prob3)) # output: expected utility Eud, Ed3, EsP, Epgo
   
}

#################################
# Two or three phase III trials #
#################################
# Case 2: Strategy 2/2( + 1); at least two trials significant (and the 
# treatment effect of the other one at least showing in the same direction)

# Expected probability to do third phase III trial: Epgo3
pgof23 <-  function(HRgo, d2, alpha, beta, theta){
   
   SIGMA <-  diag(2)
   c     <-  (qnorm(1 - alpha) + qnorm(1 - beta))^2
   
   return(integrate(function(y){
      sapply(y, function(y){
         2 * (pmvnorm(lower = c(qnorm(1 - alpha), 
                                0), 
                      upper = c(Inf, 
                                qnorm(1 - alpha)), 
                      mean = c(theta/sqrt(y^2/c), 
                               theta/sqrt(y^2/c)), 
                      sigma = SIGMA)) * 
            dnorm(y, 
                  mean = theta, 
                  sd = sqrt(4/d2)) 
      })
   },  - log(HRgo), Inf)$value)
} 

# Expected number of events for third phase III when going to third phase III: Ed33
Ed3f23 <-  function(HRgo, d2, alpha, beta, theta, ymin){
   
   ( (4 * (qnorm(1 - alpha) + qnorm(1 - beta))^2)/(ymin^2)) * 
      pgo23f(HRgo, d2, alpha, beta, w)
   
}

# Ethetapected probability of a successful program
EPsProgf23 <-  function(HRgo, d2, alpha, beta, theta, case, size, ymin){
   # Option 2.1: first two phase III trials are successful: no third phase III trial
   # Option 2.2: one of the two first phase III trials successful, the treatment
   #  effect of the other one shows a minimal clinically relevant effect: 
   #  conduct third phase III trial with d3 = d3(ymin)
   
   SIGMA <-  diag(2)
   SIGMA3<-  diag(3)
   c     <-  (qnorm(1 - alpha) + qnorm(1 - beta))^2
   
   if(case == 2){ # Option 2.1
      if(size == "small"){
         return(integrate(function(y){
            sapply(y, function(y){
               ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                   qnorm(1 - alpha)), 
                         upper = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c), 
                                   qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c)), 
                         mean = c(theta/sqrt(y^2/c), 
                                  theta/sqrt(y^2/c)), 
                         sigma = SIGMA)  - 
                    pmvnorm(lower = c(qnorm(1 - alpha) - 
                                         log(0.95)/sqrt(y^2/c), 
                                      qnorm(1 - alpha) - 
                                         log(0.95)/sqrt(y^2/c)), 
                            upper = c(qnorm(1 - alpha) - 
                                         log(0.85)/sqrt(y^2/c), 
                                      qnorm(1 - alpha) - 
                                         log(0.85)/sqrt(y^2/c)), 
                            mean = c(theta/sqrt(y^2/c), 
                                     theta/sqrt(y^2/c)), 
                            sigma = SIGMA)) * 
                  dnorm(y, 
                        mean = theta, 
                        sd = sqrt(4/d2)) 
            })
         },  - log(HRgo), Inf)$value) 
      }
      if(size == "large"){
         return(integrate(function(y){
            sapply(y, function(y){
               ( pmvnorm(lower = c(qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c), 
                                   qnorm(1 - alpha) - 
                                      log(0.85)/sqrt(y^2/c)), 
                         upper = c(Inf, Inf), 
                         mean = c(theta/sqrt(y^2/c), theta/sqrt(y^2/c)), 
                         sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = theta, 
                        sd = sqrt(4/d2)) 
            })
         },  - log(HRgo), Inf)$value)    
      }
      if(size == "all"){
         return(integrate(function(y){
            sapply(y, function(y){
               ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                   qnorm(1 - alpha)), 
                         upper = c(Inf, 
                                   Inf), 
                         mean = c(theta/sqrt(y^2/c), 
                                  theta/sqrt(y^2/c)), 
                         sigma = SIGMA) ) * 
                  dnorm(y, 
                        mean = theta, 
                        sd = sqrt(4/d2))
            })
         },  - log(HRgo), Inf)$value)    
      }
   }
   if(case == 3){# Option 2.2
      if(size == "small"){
         return(integrate(function(y){
            sapply(y, function(y){
               2 * ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                       0, 
                                       qnorm(1 - alpha)), 
                             upper = c(qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c), 
                                       qnorm(1 - alpha), 
                                       qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(ymin^2/c)), 
                             mean = c(theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c), 
                                      theta/sqrt(ymin^2/c)), 
                             sigma = SIGMA3)  - 
                        pmvnorm(lower = c(qnorm(1 - alpha) - 
                                             log(0.95)/sqrt(y^2/c), 
                                          0, 
                                          qnorm(1 - alpha) - 
                                             log(0.95)/sqrt(ymin^2/c)), 
                                upper = c(qnorm(1 - alpha) - 
                                             log(0.85)/sqrt(y^2/c), 
                                          qnorm(1 - alpha), 
                                          qnorm(1 - alpha) - 
                                             log(0.85)/sqrt(ymin^2/c)), 
                                mean = c(theta/sqrt(y^2/c), 
                                         theta/sqrt(y^2/c), 
                                         theta/sqrt(ymin^2/c)), 
                                sigma = SIGMA3)) * 
                  dnorm(y, 
                        mean = theta, 
                        sd = sqrt(4/d2))  
            })
         },  - log(HRgo), Inf)$value) 
      }
      if(size == "large"){
         return(integrate(function(y){
            sapply(y, function(y){
               2 * ( pmvnorm(lower = c(qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(y^2/c), 
                                       0, 
                                       qnorm(1 - alpha) - 
                                          log(0.85)/sqrt(ymin^2/c)), 
                             upper = c(Inf, 
                                       qnorm(1 - alpha), 
                                       Inf), 
                             mean = c(theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c), 
                                      theta/sqrt(ymin^2/c)), 
                             sigma = SIGMA3) ) * 
                  dnorm(y, 
                        mean = theta, 
                        sd = sqrt(4/d2)) 
            })
         },  - log(HRgo), Inf)$value)    
      }
      if(size == "all"){
         return(integrate(function(y){
            sapply(y, function(y){
               2 * ( pmvnorm(lower = c(qnorm(1 - alpha), 
                                       0, 
                                       qnorm(1 - alpha)), 
                             upper = c(Inf, 
                                       qnorm(1 - alpha), 
                                       Inf), 
                             mean = c(theta/sqrt(y^2/c), 
                                      theta/sqrt(y^2/c), 
                                      theta/sqrt(ymin^2/c)), 
                             sigma = SIGMA3) ) * 
                  dnorm(y, 
                        mean = theta, 
                        sd = sqrt(4/d2)) 
            })
         },  - log(HRgo), Inf)$value)    
      }
   }
   
}

# Utility function
utilityf23 <-  function(HRgo,d2,alpha,beta,theta,p2,p3,c2,c02,c3,c03,b1,b2,b3,case){ 
   
   pg  <-  pgof23(HRgo = HRgo, d2 = d2, theta = theta)
   d3  <-  Ed3f23(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta)
   
   # probability of a successful program: effect sizes, 
   # for program with two phase III trials
   
   # small 
   prob1    <-  EPsProgf23(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta,
                           case = 2, size = "small", ymin = ymin) 
   # large
   prob3    <-  EPsProgf23(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta,
                           case = 2, size = "large", ymin = ymin) 
   # medium
   prob2    <-  EPsProgf23(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta,
                           case = 2, size = "all", ymin = ymin)  - prob1  -  prob3
   
   # prob to do third phase III trial
   pg3   <-  pgof23(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta) 
   
   # d3 for third phase III trial
   d33   <-  (4 * (qnorm(1 - alpha) + qnorm(1 - beta))^2)/(ymin^2) 
   
   
   # probability of a successful program: effect sizes, 
   # for program with third phase III trial
   
   # small 
   prob13   <-  EPsProgf23(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta,
                           case = 3, size = "small", ymin = ymin) 
   # large
   prob33   <-  EPsProgf23(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta,
                           case = 3, size = "large", ymin = ymin) 
   # medium
   prob23   <-  EPsProgf23(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta,
                           case = 3, size = "all", ymin = ymin)  -  prob1  -  prob3
   
   K1    <-  c02 + c2 * d2 * (1/p2) # cost phase II 
   
   # cost for one of the first two phase III trials in case of go decision
   K2    <-  c03 * pg + c3 * d3 * (1/P3[2]) 
   
   # cost for the third phase III trial in case of third phase III trial
   K3    <-  pg3 * (c03 + c3 * d33 * (1/P3[1])) 
   
   G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 
   +  b1 * prob13 + b2 * prob23 + b3 * prob33 # gain
   
   EU    <-  - K1 - 2 * K2 - K3 + G
   SP    <-  prob1 + prob2 + prob3 +  
      prob13 + prob23 + prob33
   
   # output: expected utility Eud, Ed3, EsP, Epgo, Epgo3, Ed33    
   return(c(EU, d3, SP, pg, pg3, d33 * pg3)) 
   
}

shinyServer(function(input, output,session) {
   
     output$table <- renderTable({
        
         input$go
         case = isolate(input$Case)
         input$go
         d2 = isolate(input$D2)
         input$go
         HRgo = isolate(input$HRgo)
         input$go
         stepd2 = isolate(input$stepD2)
         input$go
         stepHRgo = isolate(input$stepHRgo)

         input$go
         beta = isolate(input$beta)
         input$go
         p2 = isolate(input$p2)
         input$go
         p3 = isolate(input$p3)
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
         HR = isolate(input$HR)
         input$go
         refresh = isolate(input$refresh)

         
         theta   = -log(HR)

         
         d2seq = seq(from = d2[1], to = d2[2], by = stepd2)
         hrgo = seq(from = HRgo[1], to = HRgo[2], by = stepHRgo)
         
         
         if(case==1){
            # Strategy 1alpha vs. Strategy 1/2,
            STRATEGY = c(1, 2)
            }
         if(case==2){
            # Strategy 1alpha^2 vs. Strategy 2/2 vs. Strategy 2/3 vs. Strategy 2/2( + 1)
            STRATEGY = c(1, 2, 3)
            }
         if(case==3){
            # Strategy 1alpha^3 vs. Strategy 3/3 vs. Strategy 3/4
            STRATEGY = c(1, 3, 4)
            }
         
         
         for(strategy in STRATEGY){
         
            
         if(strategy==1){
            
            if(case==1){
               input$go
               alpha = isolate(input$alpha)
            }
            if(case==2){
               input$go
               alpha = isolate(input$alpha)^2
            }
            if(case==3){
               input$go
               alpha = isolate(input$alpha)^3
            }
            
            y = function(x){
               result = utilityf(HRgo,x,alpha,beta,theta,p2,p3,c2,c02,c3,c03,b1,b2,b3)
               return(result)
            }
         }  
         if(strategy==2){
            
            input$go
            alpha = isolate(input$alpha)
            
            y = function(x){
               result = utilityf2(HRgo,x,alpha,beta,theta,p2,p3,c2,c02,c3,c03,b1,b2,b3,case)
               return(result)
            }
         }  
         if(strategy==3){
            
            input$go
            alpha = isolate(input$alpha)
            
            y = function(x){
               result = utilityf3(HRgo,x,alpha,beta,theta,p2,p3,c2,c02,c3,c03,b1,b2,b3,case)
               return(result)
            }
         }  
         if(strategy==4){
            
            input$go
            alpha = isolate(input$alpha)
            
            y = function(x){
               result = utilityf4(HRgo,x,alpha,beta,theta,p2,p3,c2,c02,c3,c03,b1,b2,b3,case)
               return(result)
            }
         }     
            
         ufkt <- d2fkt <- HRgofkt <- d3fkt<- spfkt <- pgofkt <- sp1fkt <- sp2fkt <- sp3fkt <- 
            dfkt <- K2fkt <- K3fkt<- matrix(0,length(d2seq),length(hrgo))
         
         progress <- Progress$new(session, min=1, max=(length(hrgo)))
         on.exit(progress$close())
         
         progress$set(message = 'Optimization progress',
                      detail=paste("for Strategy",strategy))
      
         for(j in c(1:(length(hrgo)))){
            
            HRgo = hrgo[j]

            result <-sapply(d2seq,y)

            progress$set(value = j)
            
            #c(EU, d2, d3, SP, pg, K2, K3, prob1, prob2, prob3)
            
            ufkt[,j]       <- result[1,]
            
            HRgofkt[,j]    <- rep(HRgo,length(d2seq))
            d2fkt[,j]      <- result[2,]
            d3fkt[,j]      <- result[3,]
            spfkt[,j]      <- result[4,]
            pgofkt[,j]     <- result[5,]
            K2fkt[,j]      <- result[6,]
            K3fkt[,j]      <- result[7,]
            sp1fkt[,j]      <- result[8,]
            sp2fkt[,j]      <- result[9,]
            sp3fkt[,j]      <- result[10,]
            
            dfkt[,j]       <- d3fkt[,j]+d2seq
         }

         if(strategy==1){

            save(ufkt,file=paste0(mainPath, "ufkt1.RData"))
            save(d2fkt,file=paste0(mainPath, "d2fkt1.RData"))
            save(HRgofkt,file=paste0(mainPath, "HRgofkt1.RData"))   
         }
         if(strategy==2){

            save(ufkt,file=paste0(mainPath, "ufkt2.RData"))
            save(d2fkt,file=paste0(mainPath, "d2fkt2.RData"))
            save(HRgofkt,file=paste0(mainPath, "HRgofkt2.RData"))   
         }
         if(strategy==3){
            
            save(ufkt,file=paste0(mainPath, "ufkt3.RData"))
            save(d2fkt,file=paste0(mainPath, "d2fkt3.RData"))
            save(HRgofkt,file=paste0(mainPath, "HRgofkt3.RData"))   
         }
         if(strategy==4){
            
            save(ufkt,file=paste0(mainPath, "ufkt4.RData"))
            save(d2fkt,file=paste0(mainPath, "d2fkt4.RData"))
            save(HRgofkt,file=paste0(mainPath, "HRgofkt4.RData"))   
         }
         
         ind   <- which(ufkt == max(ufkt), arr.ind <- TRUE)
         
         I <- as.vector(ind[1,1])
         J <- as.vector(ind[1,2])
         
         Eud   <- ufkt[I,J]
         d2tot <- d2fkt[I,J]
         d3tot <- ceiling(d3fkt[I,J])
         dtot  <- dfkt[I,J]
         prob  <- spfkt[I,J]
         pg    <- pgofkt[I,J]
         prob1 <- sp1fkt[I,J]
         prob2 <- sp2fkt[I,J]
         prob3 <- sp3fkt[I,J]
         k2    <- K2fkt[I,J]
         k3    <- K3fkt[I,J]
         
         
         if(refresh==1&strategy==1){
            DF = NULL
            }else{load(file=paste0(mainPath, "optimizationresults.RData"))}
         
      
         
         DF <- rbind(DF,data.frame(Strategy=format(strategy,digits=0),Case=case,u=Eud,HR=exp(-theta),
                                   HRgo=hrgo[J],d2=format(d2tot,digits=0),
                                   d3=format(d3tot,digits=0),d=format(d2tot+d3tot*as.numeric(strategy),digits=0),
                                   n2= format(d2tot/p2,digits=0),n3= format(d3tot/p3,digits=0),
                                   n=format(d2tot/p2+d3tot/p3*as.numeric(strategy),digits=0),
                                   pgo=pg,sProg=prob,sProg1=prob1,sProg2=prob2,sProg3=prob3,K2=k2,K3=k3,
                                   xi2=p2,xi3=p3,alpha=format(alpha,digits=3),beta=beta,
                                   c02=c02,c03=c03,c2=c2,c3=c3,b1=b1,b2=b2,b3=b3))
         
         save(DF,file=paste0(mainPath, "optimizationresults.RData"))
         }
         return(DF)
               })

     output$plot <- renderPlotly({
           
           input$go
           case = isolate(input$Case)
           input$go
           Plot = isolate(input$Plot)
           if(Plot==1){
    
              if(case==1){
                 
                 load(file=paste0(mainPath, "ufkt1.RData"))
                 ufkt1  = ufkt
                 load(file=paste0(mainPath, "ufkt2.RData"))
                 ufkt2  = ufkt
                 
                 if(max(ufkt1)>max(ufkt2)){
                    
                    load(file=paste0(mainPath, "ufkt1.RData"))
                    load(file=paste0(mainPath, "d2fkt1.RData"))
                    load(file=paste0(mainPath, "HRgofkt1.RData"))  
                 }else{
                    
                    load(file=paste0(mainPath, "ufkt2.RData"))
                    load(file=paste0(mainPath, "d2fkt2.RData"))
                    load(file=paste0(mainPath, "HRgofkt2.RData"))   
                 }
                 
              }
              
              if(case==2){
                 
                 load(file=paste0(mainPath, "ufkt1.RData"))
                 ufkt1  = ufkt
                 load(file=paste0(mainPath, "ufkt2.RData"))
                 ufkt2  = ufkt
                 load(file=paste0(mainPath, "ufkt3.RData"))
                 ufkt3  = ufkt
                 
                 u1 = max(ufkt1)
                 u2 = max(ufkt2)
                 u3 = max(ufkt3)
                 
                 u = c(u1,u2,u3)
                 
                 if(max(u)==max(ufkt1)){
                    load(file=paste0(mainPath, "ufkt1.RData"))
                    load(file=paste0(mainPath, "d2fkt1.RData"))
                    load(file=paste0(mainPath, "HRgofkt1.RData"))  
                 }
                 if(max(u)==max(ufkt2)){
                    load(file=paste0(mainPath, "ufkt2.RData"))
                    load(file=paste0(mainPath, "d2fkt2.RData"))
                    load(file=paste0(mainPath, "HRgofkt2.RData"))   
                 }
                 if(max(u)==max(ufkt3)){
                    load(file=paste0(mainPath, "ufkt3.RData"))
                    load(file=paste0(mainPath, "d2fkt3.RData"))
                    load(file=paste0(mainPath, "HRgofkt3.RData"))   
                 }
                 
              }
              
              if(case==3){
                 
                 load(file=paste0(mainPath, "ufkt1.RData"))
                 ufkt1  = ufkt
                 load(file=paste0(mainPath, "ufkt3.RData"))
                 ufkt3  = ufkt
                 load(file=paste0(mainPath, "ufkt4.RData"))
                 ufkt4  = ufkt
                 
                 u1 = max(ufkt1)
                 u3 = max(ufkt3)
                 u4 = max(ufkt4)
                 
                 u = c(u1,u3,u4)
                 
                 if(max(u)==max(ufkt1)){
                    load(file=paste0(mainPath, "ufkt1.RData"))
                    load(file=paste0(mainPath, "d2fkt1.RData"))
                    load(file=paste0(mainPath, "HRgofkt1.RData"))  
                 }
                 if(max(u)==max(ufkt3)){
                    load(file=paste0(mainPath, "ufkt3.RData"))
                    load(file=paste0(mainPath, "d2fkt3.RData"))
                    load(file=paste0(mainPath, "HRgofkt3.RData"))   
                 }
                 if(max(u)==max(ufkt4)){
                    load(file=paste0(mainPath, "ufkt4.RData"))
                    load(file=paste0(mainPath, "d2fkt4.RData"))
                    load(file=paste0(mainPath, "HRgofkt4.RData"))   
                 }
                 
              }
              
              plot_ly(x=HRgofkt,y=d2fkt,z=ufkt, type="surface")   
           }

           
        })   
        
     
})
