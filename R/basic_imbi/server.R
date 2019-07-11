library("shiny")
library("ggplot2")
library("gridExtra")
library("mvtnorm")
library("doParallel")
library("parallel")
library("cubature")
library("plotly")

mainPath <- "/opt/shiny-server/samplesizr/basic/"

#mainPath <- ""

######################################
# One phase III trial - fixed effect #
######################################
# theta: true underlying effect
# y: hat_theta_2; estimator in phase II
# z: T_3; normalized estimator in phase III

# Probability to go to phase III: pgo
pgof_time <-  function(HRgo, d2, theta){

   pnorm((log(HRgo) + theta)/sqrt(4/d2))
}

# Expected number of events for phase III when going to phase III: Ed3
Ed3f_time <-  function(HRgo, d2, alpha, beta, theta){

   ceiling(integrate(function(y){
      ( (4 * (qnorm(1 - alpha/2) + qnorm(1 - beta))^2)/(y^2)) *
         dnorm(y,
               mean = theta,
               sd = sqrt(4/d2))
   },  - log(HRgo), Inf)$value)
}

# Expected probability of a successful program: EsP
EPsProgf_time <-  function(HRgo, d2, alpha, beta, theta, step1, step2){

   c = (qnorm(1 - alpha/2) + qnorm(1 - beta))^2

   integrate(function(y){
      ( pnorm(qnorm(1 - alpha/2) - log(step2)/(sqrt(y^2/c)),
              mean = theta/(sqrt(y^2/c)),
              sd = 1) -
           pnorm(qnorm(1 - alpha/2) - log(step1)/(sqrt(y^2/c)),
                 mean = theta/(sqrt(y^2/c)),
                 sd = 1) ) *
         dnorm(y,
               mean = theta,
               sd = sqrt(4/d2))
   },  - log(HRgo), Inf)$value


}

# Utility function
utilityf_time <-  function(HRgo,d2,alpha,beta,theta,p2,p3,c2,c02,c3,c03,b1,b2,b3,K,
                           steps1,steps2,stepm1,stepm2,stepl1,stepl2){

   pg    <-  pgof_time(HRgo = HRgo, d2 = d2, theta = theta)
   d3    <-  Ed3f_time(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta)

   # round up to next even natural number
   n2 = ceiling(d2 * (1/p2))
   if(round(n2/2) != n2 / 2) {n2 = n2 + 1}
   
   n3 = ceiling(d3 * (1/p3))
   if(round(n3/2) != n3 / 2) {n3 = n3 + 1}
  
   K2    <-  c02 + c2 * n2 #cost phase II
   K3    <-  c03 * pg + c3 * n3 #cost phase III

   if(K2+K3>K){

      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999))


   }else{

      # probability of a successful program; small, medium, large effect size
      prob1 <-  EPsProgf_time(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta,
                         step1 = steps1, step2 =  steps2)
      prob2 <-  EPsProgf_time(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta,
                         step1 =  stepm1, step2 =  stepm2)
      prob3 <-  EPsProgf_time(HRgo = HRgo, d2 = d2, alpha = alpha, beta = beta, theta = theta,
                         step1 =  stepl1, step2 = stepl2)

      G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain

      EU    <-  - K2 - K3 + G
      SP    <-  prob1 + prob2 + prob3

      return(c(EU, d2, d3, SP, pg, K2, K3, prob1, prob2, prob3, n2, n3))
      #output: expected utility Eud, En3, EsP, Epgo, cost phase II and III
   }

}



# Probability to go to phase III: pgo
pgof_binary <-  function(RRgo, n2, p0, p1){

   rho = -log(p1/p0)
   Var   = (2/n2)*((1-p0)/p0+(1-p1)/p1)
   pnorm((rho + log(RRgo))/sqrt(Var))
}

# Expected number of events for phase III when going to phase III: En3
En3f_binary <-  function(RRgo, n2, alpha, beta, p0, p1){

   rho   = -log(p1/p0)
   Var   = (2/n2)*((1-p0)/p0+(1-p1)/p1)
   P     = (p0 + p1)/2
   c     = (qnorm(1-alpha/2)*sqrt(2*(1-P)/P)+qnorm(1-beta)*sqrt((1-p0)/p0+(1-p1)/p1))^2

   ceiling(integrate(function(y){
      ((2*c)/y^2) *
         dnorm(y,
               mean = rho,
               sd = sqrt(Var))
   }, - log(RRgo), Inf)$value)
}

# Expected probability of a successful program: EsP
EPsProgf_binary <-  function(RRgo, n2, alpha, beta, p0, p1, step1, step2){

   # Phase II
   rho   = -log(p1/p0)
   Var   = (2/n2)*((1-p0)/p0+(1-p1)/p1)

   # Phase III
   V     = ((1-p0)/p0+(1-p1)/p1)
   P     = (p0 + p1)/2
   c     = (qnorm(1-alpha/2)*sqrt(2*(1-P)/P)+qnorm(1-beta)*sqrt((1-p0)/p0+(1-p1)/p1))^2

   integrate(function(y){
      ( pnorm(qnorm(1 - alpha/2) - log(step2)/(sqrt(V*y^2/c)),
              mean = rho/(sqrt(V*y^2/c)),
              sd = 1) -
           pnorm(qnorm(1 - alpha/2) - log(step1)/(sqrt(V*y^2/c)),
                 mean = rho/(sqrt(V*y^2/c)),
                 sd = 1) ) *
         dnorm(y,
               mean = rho,
               sd = sqrt(Var))
   },  - log(RRgo), Inf)$value

}

# Utility function
utilityf_binary <-  function(RRgo,n2,alpha,beta,p0,p1,c2,c02,c3,c03,b1,b2,b3,K,
                             steps1,steps2,stepm1,stepm2,stepl1,stepl2){

   pg    <-  pgof_binary(RRgo = RRgo, n2 = n2, p0 = p0, p1 = p1)
   n3    <-  En3f_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta, p0 = p0, p1 = p1)

   # round up to next even natural number
   if(round(n3/2) != n3 / 2) {n3 = n3 + 1}
   
   K2    <-  c02 + c2 * n2       #cost phase II
   K3    <-  c03 * pg + c3 * n3  #cost phase III

   if(K2+K3>K){

      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999)) #output: expected utility Eud, En3, EsP, Epgo

   }else{

      # probability of a successful program; small, medium, large effect size
      prob1 <-  EPsProgf_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta, p0 = p0, p1 = p1,
                         step1 = steps1, step2 =  steps2)
      prob2 <-  EPsProgf_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta, p0 = p0, p1 = p1,
                         step1 =  stepm1, step2 =  stepm2)
      prob3 <-  EPsProgf_binary(RRgo = RRgo, n2 = n2, alpha = alpha, beta = beta, p0 = p0, p1 = p1,
                         step1 =  stepl1, step2 = stepl2)

      G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain

      EU    <-  - K2 - K3 + G
      SP    <-  prob1 + prob2 + prob3

      return(c(EU, n2, n3, SP, pg, K2, K3, prob1, prob2, prob3, n2, n3))
      #output: expected utility Eud, En3, EsP, Epgo, cost phase II and III
   }

}


# x: mu; true underlying standardized difference in means: Delta := (mu0 - mu1)/sigma; H0: Delta <= 0, H1: Delta>0 with pooled variance sigma^2
# y: hat_Delta_2; estimator in phase II
# z: hat_Delta_3; normalized estimator in phase III

# Probability to go to phase III: pgo
pgof_normal <-  function(kappa, n2, Delta){
   pnorm((Delta-kappa)/sqrt(4/n2))
}

# Expected number of events for phase III when going to phase III: En3
En3f_normal <-  function(kappa, n2, alpha, beta, Delta){
   ceiling(integrate(function(y){
      ((4*(qnorm(1-alpha/2)+qnorm(1-beta))^2)/y^2) *
         dnorm(y,
               mean = Delta,
               sd = sqrt(4/n2))
   }, kappa, Inf)$value)
}

# Expected probability of a successful program: EsP
EPsProgf_normal <-  function(kappa, n2, alpha, beta, Delta, step1, step2){

   c     <- (qnorm(1-alpha/2)+qnorm(1-beta))^2

   integrate(function(y){
      ( pnorm(qnorm(1 - alpha/2) + step2/sqrt(y^2/c),
              mean = Delta/sqrt(y^2/c),
              sd = 1) -
           pnorm(qnorm(1 - alpha/2) + step1/sqrt(y^2/c),
                 mean = Delta/sqrt(y^2/c),
                 sd = 1) ) *
         dnorm(y,
               mean = Delta,
               sd = sqrt(4/n2))
   },  kappa, Inf)$value

}

# Utility function
utilityf_normal <-  function(kappa,n2,alpha,beta,Delta,c2,c02,c3,c03,b1,b2,b3,K,
                             steps1,steps2,stepm1,stepm2,stepl1,stepl2){

   pg    <-  pgof_normal(kappa = kappa, n2 = n2, Delta = Delta)
   n3    <-  En3f_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta, Delta = Delta)

   if(round(n3/2) != n3 / 2) {n3 = n3 + 1}
   
   K2    <-  c02 + c2 * n2       #cost phase II
   K3    <-  c03 * pg + c3 * n3  #cost phase III

   if(K2+K3>K){

      return(c(-9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999, -9999)) #output: expected utility Eud, En3, EsP, Epgo

   }else{

      # probability of a successful program; small, medium, large effect size
      prob1 <-  EPsProgf_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta, Delta = Delta,
                         step1 = steps1, step2 =  steps2)
      prob2 <-  EPsProgf_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta, Delta = Delta,
                         step1 =  stepm1, step2 =  stepm2)
      prob3 <-  EPsProgf_normal(kappa = kappa, n2 = n2, alpha = alpha, beta = beta, Delta = Delta,
                         step1 =  stepl1, step2 = stepl2)

      G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain

      EU    <-  - K2 - K3 + G
      SP    <-  prob1 + prob2 + prob3

      return(c(EU, n2, n3, SP, pg, K2, K3, prob1, prob2, prob3, n2, n3))
      #output: expected utility Eud, En3, EsP, Epgo, cost phase II and III
   }

}





shinyServer(function(input, output,session) {

output$table <- renderTable({


  Select = input$Select

   if(Select==1){
     input$go1
     Num1 = isolate(input$Num1)
     input$go1
     HRgo = isolate(input$HRgo)
     input$go1
     stepNum1 = isolate(input$stepNum1)
     input$go1
     stepHRgo = isolate(input$stepHRgo)
     input$go1
     HR = isolate(input$HR)

     input$go1
     alpha = isolate(input$alpha1)
     input$go1
     beta = isolate(input$beta1)
     input$go1
     p2 = isolate(input$p2)
     input$go1
     p3 = isolate(input$p3)
     input$go1
     c02 = isolate(input$c021)
     input$go1
     c03 = isolate(input$c031)
     input$go1
     c2 = isolate(input$c21)
     input$go1
     c3 = isolate(input$c31)
     input$go1
     steps1 = isolate(input$small1)
     input$go1
     stepm1 = isolate(input$medium1)
     input$go1
     stepl1 = isolate(input$large1)
     input$go1
     b1 = isolate(input$b11)
     input$go1
     b2 = isolate(input$b21)
     input$go1
     b3 = isolate(input$b31)
     input$go1
     refresh = isolate(input$refresh1)
     input$go1
     K = isolate(input$K1)


  }

   if(Select==2){

     input$go2
     p0 = isolate(input$p0)
     input$go2
     p1 = isolate(input$p1)
     input$go2
     Num2 = isolate(input$Num2)
     input$go2
     RRgo = isolate(input$RRgo)
     input$go2
     stepNum2 = isolate(input$stepNum2)
     input$go2
     stepRRgo = isolate(input$stepRRgo)


     input$go2
     alpha = isolate(input$alpha2)
     input$go2
     beta = isolate(input$beta2)
     input$go2
     p2 = isolate(input$p22)
     input$go2
     p3 = isolate(input$p32)
     input$go2
     c02 = isolate(input$c022)
     input$go2
     c03 = isolate(input$c032)
     input$go2
     c2 = isolate(input$c22)
     input$go2
     c3 = isolate(input$c32)
     input$go2
     steps1 = isolate(input$small2)
     input$go2
     stepm1 = isolate(input$medium2)
     input$go2
     stepl1 = isolate(input$large2)
     input$go2
     b1 = isolate(input$b12)
     input$go2
     b2 = isolate(input$b22)
     input$go2
     b3 = isolate(input$b32)
     input$go2
     refresh = isolate(input$refresh2)
     input$go2
     K = isolate(input$K2)

  }

   if(Select==3){

     input$go3
     Delta = isolate(input$Delta)
     input$go3
     Num3 = isolate(input$Num3)
     input$go3
     kappa = isolate(input$kappa)
     input$go3
     stepNum3 = isolate(input$stepNum3)
     input$go3
     stepkappa = isolate(input$stepkappa)



     input$go3
     alpha = isolate(input$alpha3)
     input$go3
     beta = isolate(input$beta3)
     input$go3
     p2 = isolate(input$p23)
     input$go3
     p3 = isolate(input$p33)
     input$go3
     c02 = isolate(input$c023)
     input$go3
     c03 = isolate(input$c033)
     input$go3
     c2 = isolate(input$c23)
     input$go3
     c3 = isolate(input$c33)
     input$go3
     steps1 = isolate(input$small3)
     input$go3
     stepm1 = isolate(input$medium3)
     input$go3
     stepl1 = isolate(input$large3)
     input$go3
     b1 = isolate(input$b13)
     input$go3
     b2 = isolate(input$b23)
     input$go3
     b3 = isolate(input$b33)
     input$go3
     refresh = isolate(input$refresh3)
     input$go3
     K = isolate(input$K3)

  }

   if(Select==1){

     theta   = -log(HR)

     NUM = seq(from = Num1[1], to = Num1[2], by = stepNum1)
     TRES = seq(from = HRgo[1], to = HRgo[2], by = stepHRgo)

     steps2 <- stepm1
     stepm2 <- stepl1
     stepl2 <- 0


     y = function(x){
        result = utilityf_time(tres,x,alpha,beta,theta,p2,p3,c2,c02,c3,c03,b1,b2,b3,
                               K,steps1,steps2,stepm1,stepm2,stepl1,stepl2)
        return(result)
     }
  }

   if(Select==2){

      NUM = seq(from = Num2[1], to = Num2[2], by = stepNum2)
      TRES = seq(from = RRgo[1], to = RRgo[2], by = stepRRgo)

      steps2 <- stepm1
      stepm2 <- stepl1
      stepl2 <- 0

      y = function(x){
         result = utilityf_binary(tres,x,alpha,beta,p0,p1,c2,c02,c3,c03,b1,b2,b3,
                                  K,steps1,steps2,stepm1,stepm2,stepl1,stepl2)
         return(result)
      }
   }

   if(Select==3){

     NUM = seq(from = Num3[1], to = Num3[2], by = stepNum3)
     TRES = seq(from = kappa[1], to = kappa[2], by = stepkappa)

     steps2 <- stepm1
     stepm2 <- stepl1
     stepl2 <- Inf

     y = function(x){
        result = utilityf_normal(tres,x,alpha,beta,Delta,c2,c02,c3,c03,b1,b2,b3,
                                 K,steps1,steps2,stepm1,stepm2,stepl1,stepl2)
        return(result)
     }
  }

   ufkt <- numfkt <- tresfkt <- d3fkt<- spfkt <- pgofkt  <- sp1fkt <-sp2fkt <- sp3fkt <-
      dfkt <- K2fkt <- K3fkt<- n2fkt<- n3fkt<- matrix(0,length(NUM),length(TRES))

   progress <- Progress$new(session, min=1, max=(length(TRES)))
   on.exit(progress$close())

   progress$set(message = 'Optimization progress')

   for(j in c(1:(length(TRES)))){

      tres = TRES[j]

      result <-sapply(NUM,y)

      progress$set(value = j)

      #      return(c(EU, d2, d3, SP, pg, K2, K3, prob1, prob2, prob3, n2, n3))

      ufkt[,j]       <- result[1,]

      tresfkt[,j]    <- rep(tres,length(NUM))
      numfkt[,j]     <- result[2,]
      d3fkt[,j]      <- result[3,]
      spfkt[,j]      <- result[4,]
      pgofkt[,j]     <- result[5,]
      K2fkt[,j]      <- result[6,]
      K3fkt[,j]      <- result[7,]
      sp1fkt[,j]     <- result[8,]
      sp2fkt[,j]     <- result[9,]
      sp3fkt[,j]     <- result[10,]
      n2fkt[,j]      <- result[11,]
      n3fkt[,j]      <- result[12,]

      dfkt[,j]       <- d3fkt[,j]+NUM
   }

   save(ufkt,file=paste0(mainPath, "ufkt.RData"))
   save(numfkt,file=paste0(mainPath, "numfkt.RData"))
   save(tresfkt,file=paste0(mainPath, "tresfkt.RData"))

   ind   <- which(ufkt == max(ufkt), arr.ind <- TRUE)

   I <- as.vector(ind[1,1])
   J <- as.vector(ind[1,2])

   Eud   <- ufkt[I,J]
   num2tot <- numfkt[I,J]
   num3tot <- d3fkt[I,J]
   numtot  <- dfkt[I,J]
   prob  <- spfkt[I,J]
   pg    <- pgofkt[I,J]
   prob1 <- sp1fkt[I,J]
   prob2 <- sp2fkt[I,J]
   prob3 <- sp3fkt[I,J]
   k2    <- K2fkt[I,J]
   k3    <- K3fkt[I,J]
   n2tot <- n2fkt[I,J]
   n3tot <- n3fkt[I,J]

   if(refresh==0){
      load(file=paste0(mainPath, "optimizationresults.RData"))
      }else{DF = NULL}

   if(Select==1){
      DF <- rbind(DF,data.frame(u=Eud,HR=HR,
                                HRgo=TRES[J],d2=format(num2tot,digits=0),
                                d3=format(num3tot,digits=0),d=format(numtot,digits=0),
                                n2=format(n2tot,digits=0),n3=format(n3tot,digits=0),
                                n=format(n2tot+n3tot,digits=0),
                                K=format(K,digits=0),
                                pgo=pg,sProg=prob,sProg1=prob1,sProg2=prob2,sProg3=prob3,K2=k2,K3=k3,
                                xi2=p2,xi3=p3,alpha=format(alpha,digits=3),beta=beta,
                                c02=c02,c03=c03,c2=c2,c3=c3,b1=b1,b2=b2,b3=b3,
                                steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2)))

      save(DF,file=paste0(mainPath, "optimizationresults.RData"))
   }

   if(Select==2){
      DF <- rbind(DF,data.frame(u=Eud,P0=p0,P1=p1,RR=p1/p0,
                                RRgo=TRES[J],
                                n2=format(num2tot,digits=0),n3=format(num3tot,digits=0),
                                n=format(num2tot+num3tot,digits=0),
                                K=format(K,digits=0),
                                pgo=pg,sProg=prob,sProg1=prob1,sProg2=prob2,sProg3=prob3,K2=k2,K3=k3,
                                alpha=format(alpha,digits=3),beta=beta,
                                c02=c02,c03=c03,c2=c2,c3=c3,b1=b1,b2=b2,b3=b3,
                                steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2)))

      save(DF,file=paste0(mainPath, "optimizationresults.RData"))
   }

   if(Select==3){
      DF <- rbind(DF,data.frame(u=Eud,Delta=Delta,
                                Kappa=TRES[J],
                                n2=format(num2tot,digits=0),n3=format(num3tot,digits=0),
                                n=format(num2tot+num3tot,digits=0),
                                K=format(K,digits=0),
                                pgo=pg,sProg=prob,sProg1=prob1,sProg2=prob2,sProg3=prob3,K2=k2,K3=k3,
                                alpha=format(alpha,digits=3),beta=beta,
                                c02=c02,c03=c03,c2=c2,c3=c3,b1=b1,b2=b2,b3=b3,
                                steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2)))

      save(DF,file=paste0(mainPath, "optimizationresults.RData"))
   }

   return(DF)
         })

     output$plot <- renderPlotly({


        Select = input$Select

        if(Select==1){
           input$go1
           Plot = input$Plot1
        }
        if(Select==2){
           input$go2
           Plot = input$Plot2
        }
        if(Select==3){
           input$go3
           Plot = input$Plot3
        }

           if(Plot==1){

              load(file=paste0(mainPath, "ufkt.RData"))
              load(file=paste0(mainPath, "numfkt.RData"))
              load(file=paste0(mainPath, "tresfkt.RData"))

              plot_ly(x=tresfkt,y=numfkt,z=ufkt, type="surface")
           }


        })


})
