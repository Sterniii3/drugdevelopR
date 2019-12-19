library("shiny")
library("ggplot2")
library("gridExtra")
library("mvtnorm")
library("doParallel")
library("parallel")
library("cubature")
library("plotly")



mainPath <- "/opt/shiny-server/samplesizr/bias/"

#mainPath <- ""

# w; weight for true underlying fixed effect
# y: hat_theta_2; estimator in phase II
# z: hat_t_3; normalized estimator in phase III

# 0. unadjusted setting with option to skip phase II
######################################################

#  probability to go to phase III
pgof<-function(HRgo,d2,theta){
   return(pnorm((theta+log(HRgo))/sqrt(4/d2)))
}

# expected number of events for phase III when going to phase III
Ed3f<-function(HRgo,d2,alpha,beta,theta){
   ceiling(integrate(function(y){
            sapply(y,function(y){
               ( (4*(qnorm(1-alpha)+qnorm(1-beta))^2)/(y^2))*
                  dnorm(y,
                        mean=theta,
                        sd=sqrt(4/d2)) 
            })
         }, -log(HRgo),Inf)$value)
} 

# expected probability of a successful program
EPsProgf<-function(HRgo,d2,alpha,beta,theta,step1,step2){
   
   c=(qnorm(1-alpha)+qnorm(1-beta))^2
   
   integrate(integrate(function(y){ 
            sapply(y,function(y){
               ( pnorm(qnorm(1-alpha)+step2/(sqrt(y^2/c)),
                       mean=theta/(sqrt(y^2/c)),
                       sd=1)-
                    pnorm(qnorm(1-alpha)+step1/(sqrt(y^2/c)),
                          mean=theta/(sqrt(y^2/c)),
                          sd=1) )*
                  dnorm(y,
                        mean=theta,
                        sd=sqrt(4/d2))
            })
         }, -log(HRgo),Inf)$value)$value
   
}

#utility function
utilityf<-function(HRgo,d2,alpha,beta,theta,p2,p3,c2,c02,c3,c03,b1,b2,b3){ 
   
   pg    = pgof(HRgo=HRgo,d2=d2,theta=theta)
   d3    = Ed3f(HRgo=HRgo,d2=d2,alpha=alpha,beta=beta,theta=theta)
   
   K2    = c02+c2*d2*(1/p2) #cost phase II 
   K3    = c03*pg+c3*d3*(1/p3) #cost phase III
   
   # probability of a successful program; small, medium, large effect size
   prob1 = EPsProgf(HRgo=HRgo,d2=d2,alpha=alpha,beta=beta,theta=theta,step1=stepls,step2=stepus)
   prob2 = EPsProgf(HRgo=HRgo,d2=d2,alpha=alpha,beta=beta,theta=theta,step1=steplm,step2=stepum)
   prob3 = EPsProgf(HRgo=HRgo,d2=d2,alpha=alpha,beta=beta,theta=theta,step1=stepll,step2=stepul)
   
   G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
   
   EU    <-  - K2 - K3 + G
   SP    <-  prob1 + prob2 + prob3

   return(c(EU, d2, d3, SP, pg, K2, K3, prob1, prob2, prob3)) #output: pgo, n3, sProg, expected utility, cost phase II and III, gain
   
}


# 1.1. conservative sample size calculation: use lower bound of one-sided confidence intervall
##############################################################################################

# prior distribution
# as above

# expected probability to go to phase III
# as above

# expected number of events for phase III when going to phase III
Ed3fL<-function(HRgo,d2,alpha,beta,theta,alphaL){
   
   int   = try(integrate(function(y){
            sapply(y,function(y){
               ( (4*(qnorm(1-alpha)+qnorm(1-beta))^2)/((y-qnorm(1-alphaL)*sqrt(4/d2))^2))*
                  dnorm(y,
                        mean=theta,
                        sd=sqrt(4/d2))
            })
         }, -log(HRgo),Inf),silent=TRUE)
   if(inherits(int ,'try-error')){
      warning(as.vector(int))
      integrated <- NA_real_
   } else {
      integrated <- int$value
   }
   return(ceiling(integrated))
} 

# expected probability of a successful program
EPsProgfL<-function(HRgo,d2,alpha,beta,theta,step1,step2,alphaL){
   
   c=(qnorm(1-alpha)+qnorm(1-beta))^2
   
         integrate(function(y){ 
            sapply(y,function(y){
               ( pnorm(qnorm(1-alpha)+step2/(sqrt((y-qnorm(1-alphaL)*sqrt(4/d2))^2/c)),
                       mean=theta/(sqrt((y-qnorm(1-alphaL)*sqrt(4/d2))^2/c)),
                       sd=1)-
                    pnorm(qnorm(1-alpha)+step1/(sqrt((y-qnorm(1-alphaL)*sqrt(4/d2))^2/c)),
                          mean=theta/(sqrt((y-qnorm(1-alphaL)*sqrt(4/d2))^2/c)),
                          sd=1) )*
                  dnorm(y,
                        mean=theta,
                        sd=sqrt(4/d2))
            })
         }, -log(HRgo),Inf)$value
   
}

# utility function
utilityfL<-function(HRgo,d2,alpha,beta,theta,p2,p3,c2,c02,c3,c03,b1,b2,b3,adj){ 
   
   alphaL = adj
   
   d3    = Ed3fL(HRgo=HRgo,d2=d2,alpha=alpha,beta=beta,theta=theta,alphaL=alphaL)
   
   if(is.na(d3)==FALSE){
      
      pg    = pgof(HRgo=HRgo,d2=d2,theta=theta)
      
      K2    = c02+c2*d2*(1/p2) #cost phase II 
      K3    = c03*pg+c3*d3*(1/p3) #cost phase III
      
      # probability of a successful program; small, medium, large effect size
      prob1 = EPsProgfL(HRgo=HRgo,d2=d2,alpha=alpha,beta=beta,theta=theta,step1=stepls,step2=stepus,alphaL=alphaL)
      prob2 = EPsProgfL(HRgo=HRgo,d2=d2,alpha=alpha,beta=beta,theta=theta,step1=steplm,step2=stepum,alphaL=alphaL)
      prob3 = EPsProgfL(HRgo=HRgo,d2=d2,alpha=alpha,beta=beta,theta=theta,step1=stepll,step2=stepul,alphaL=alphaL)
      
      G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
      
      EU    <-  - K2 - K3 + G
      SP    <-  prob1 + prob2 + prob3
      
      return(c(EU, d2, d3, SP, pg, K2, K3, prob1, prob2, prob3)) 
      #output: pgo, n3, sProg, expected utility, cost phase II and III, gain
      
   }else{
      
      return(c(-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999))
      
   }
   
   
   
}


# 2. conservative sample size calculation: use estimate with retetion factor
##############################################################################################

# prior distribution
# as above

# expected probability to go to phase III
# as above

# expected number of events for phase III when going to phase III
Ed3fR<-function(HRgo,d2,alpha,beta,theta,r){
   
   int   = try(integrate(function(y){
            sapply(y,function(y){
               ( (4*(qnorm(1-alpha)+qnorm(1-beta))^2)/((y*r)^2))*
                  dnorm(y,
                        mean=theta,
                        sd=sqrt(4/d2)) 
            })
         }, -log(HRgo),Inf),silent=TRUE)
   if(inherits(int ,'try-error')){
      warning(as.vector(int))
      integrated <- NA_real_
   } else {
      integrated <- int$value
   }
   return(ceiling(integrated))
} 

# expected probability of a successful program
EPsProgfR<-function(HRgo,d2,alpha,beta,theta,step1,step2,r){
   
   c=(qnorm(1-alpha)+qnorm(1-beta))^2
   
      integrate(function(y){ 
            sapply(y,function(y){
               ( pnorm(qnorm(1-alpha)+step2/(sqrt((y*r)^2/c)),
                       mean=theta/(sqrt((y*r)^2/c)),
                       sd=1)-
                    pnorm(qnorm(1-alpha)+step1/(sqrt((y*r)^2/c)),
                          mean=theta/(sqrt((y*r)^2/c)),
                          sd=1) )*
                  dnorm(y,
                        mean=theta,
                        sd=sqrt(4/d2))
            })
         }, -log(HRgo),Inf)$value
}

# utility function
utilityfR<-function(HRgo,d2,alpha,beta,theta,p2,p3,c2,c02,c3,c03,b1,b2,b3,adj){ 
   
   r = adj
   
   d3    = Ed3fR(HRgo=HRgo,d2=d2,alpha=alpha,beta=beta,theta=theta,r=r)
   
   if(is.na(d3)==FALSE){
      
      pg    = pgof(HRgo=HRgo,d2=d2,theta=theta)
      
      K2    = c02+c2*d2*(1/p2) #cost phase II in 10^5
      K3    = c03*pg+c3*d3*(1/p3) #cost phase III in 10^5
   
      
      # probability of a successful program; small, medium, large effect size
      prob1 = EPsProgfR(HRgo=HRgo,d2=d2,alpha=alpha,beta=beta,theta=theta,step1=stepls,step2=stepus,r=r)
      prob2 = EPsProgfR(HRgo=HRgo,d2=d2,alpha=alpha,beta=beta,theta=theta,step1=steplm,step2=stepum,r=r)
      prob3 = EPsProgfR(HRgo=HRgo,d2=d2,alpha=alpha,beta=beta,theta=theta,step1=stepll,step2=stepul,r=r)

      G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 #gain
      
      EU    <-  - K2 - K3 + G
      SP    <-  prob1 + prob2 + prob3
      
      return(c(EU, d2, d3, SP, pg, K2, K3, prob1, prob2, prob3)) 
      #output: pgo, n3, sProg, expected utility, cost phase II and III, gain
      
   }else{
      
      return(c(-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999))
      
   }
   
   
   
}

# Threshold values for effect sizes
stepls = 0
stepus = -log(0.95)
steplm = -log(0.95)
stepum = -log(0.85)
stepll = -log(0.85)
stepul = Inf   



shinyServer(function(input, output,session) {
   
     output$table <- renderTable({
        
         input$go
         Select = isolate(input$Select)
         input$go
         D2 = isolate(input$D2)
         input$go
         R = isolate(input$R)
         input$go
         alphaL = isolate(input$alphaL)
         input$go
         HRgo = isolate(input$HRgo)
         input$go
         stepD2 = isolate(input$stepD2)
         input$go
         stepHRgo = isolate(input$stepHRgo)
         input$go
         stepR = isolate(input$stepR)
         input$go
         stepalphaL = isolate(input$stepalphaL)
         input$go
         alpha = isolate(input$alpha)
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
         
         d2seq = seq(from = D2[1], to = D2[2], by = stepD2)
         hrgo = seq(from = HRgo[1], to = HRgo[2], by = stepHRgo)
         
         if(Select==1){STRATEGY = 1}
         if(Select==2){STRATEGY = 2}
         if(Select==3){STRATEGY = c(2,1)}
         
         
         if(refresh==1){
            
            DF = NULL   
            
         }else{
            
            load(file=paste0(mainPath, "optimizationresults.RData"))
         }
         
         
         for(strategy in STRATEGY){
         
         if(strategy==1){
            
            ADJ = seq(from = R[1], to = R[2], by = stepR)
            
            meth = "multipl." 
            
            y = function(x){
               result = utilityfR(HRgo,x,alpha,beta,theta,p2,p3,c2,c02,c3,c03,b1,b2,b3,adj)
               return(result)
            }
         }
         
         if(strategy==2){
            
            ADJ = seq(from = alphaL[1], to = alphaL[2], by = stepalphaL)
            
            meth = "additive"
            
            y = function(x){
               result = utilityfL(HRgo,x,alpha,beta,theta,p2,p3,c2,c02,c3,c03,b1,b2,b3,adj)
               return(result)
            }
         }  
         
         DF0 = NULL
         
         progress <- Progress$new(session, min=1, max=(length(ADJ)))
         on.exit(progress$close())
         
         progress$set(message = 'Optimization progress',
                      detail=paste("for", meth, "method"))
           
         for(i in c(1:length(ADJ))){
             
            adj = ADJ[i]
            progress$set(value = i)
         
         ufkt <- d2fkt <- HRgofkt <- d3fkt<- spfkt <- pgofkt  <- sp1fkt<- sp2fkt <- sp3fkt <- 
            dfkt <- K2fkt <- K3fkt<- matrix(0,length(d2seq),length(hrgo))
         
      
         for(j in c(1:(length(hrgo)))){
            
            HRgo = hrgo[j]

            result <-sapply(d2seq,y)
            

            #c(EU, d2, d3, SP, pg, K2, K3, prob1, prob2, prob3)
            
            ufkt[,j]       <- result[1,]
            
            HRgofkt[,j]    <- rep(HRgo,length(d2seq))
            d2fkt[,j]      <- result[2,]
            d3fkt[,j]      <- result[3,]
            spfkt[,j]      <- result[4,]
            pgofkt[,j]     <- result[5,]
            K2fkt[,j]      <- result[6,]
            K3fkt[,j]      <- result[7,]
            sp1fkt[,j]     <- result[8,]
            sp2fkt[,j]     <- result[9,]
            sp3fkt[,j]     <- result[10,]


            dfkt[,j]       <- d3fkt[,j]+d2seq
         }
         
         ind   <- which(ufkt == max(ufkt), arr.ind <- TRUE)
         
         I <- as.vector(ind[1,1])
         J <- as.vector(ind[1,2])
         
         Eud   <- ufkt[I,J]
         d2tot    <- d2fkt[I,J]
         d3tot    <- d3fkt[I,J]
         dtot     <- dfkt[I,J]
         prob  <- spfkt[I,J]
         pg    <- pgofkt[I,J]
         prob1 <- sp1fkt[I,J]
         prob2 <- sp2fkt[I,J]
         prob3 <- sp3fkt[I,J]
         k2    <- K2fkt[I,J]
         k3    <- K3fkt[I,J]

         DF0 <- rbind(DF0,data.frame(Method=meth, u=Eud, HR=exp(-theta),
                                   Adj = adj,HRgo=hrgo[J],d2=format(d2tot,digits=0),
                                   d3=format(d3tot,digits=0),d=format(dtot,digits=0),
                                   n2=format(d2tot/p2,digits=0),
                                   n3=format(d3tot/p3,digits=0),n=format(d2tot/p2+d3tot/p3,digits=0),
                                   pgo=pg,sProg=prob,Sprog1=prob1,Sprog2=prob2,Sprog3=prob3,K2=k2,K3=k3,
                                   alpha=format(alpha,digits=3),beta=beta,p2=p2,p3=p3,
                                   c02=c02,c03=c03,c2=c2,c3=c3,b1=b1,b2=b2,b3=b3))
         
         }

         if(strategy==3){
            index = which(subset(DF0,Method=="multipl.")$u == 
                             max(subset(DF0,Method=="multipl.")$u)|subset(DF0,Method=="additive")$u == 
                             max(subset(DF0,Method=="additive")$u))   
         }else{
            
            index = which(subset(DF0,Method==meth)$u == 
                             max(subset(DF0,Method==meth)$u))  
         }

         

         DF = rbind(DF,subset(DF0,Method==meth)[index,])  
            
         
         
         save(DF,file=paste0(mainPath, "optimizationresults.RData"))

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

                 load(file=paste0(mainPath, "ufkt1.RData"))
                 load(file=paste0(mainPath, "d2fkt1.RData"))
                 load(file=paste0(mainPath, "HRgofkt1.RData")) 
                 showplot = "multiplicatively"
              }
              if(Select==2){

                 load(file=paste0(mainPath, "ufkt2.RData"))
                 load(file=paste0(mainPath, "d2fkt2.RData"))
                 load(file=paste0(mainPath, "HRgofkt2.RData")) 
                 showplot = "additively"
              }
              if(Select==3){
                 
                 load(file=paste0(mainPath, "ufkt1.RData"))
                 ufkt1  = ufkt
                 load(file=paste0(mainPath, "ufkt2.RData"))
                 ufkt2  = ufkt
                 
                 if(max(ufkt1)>max(ufkt2)){
                    
                    load(file=paste0(mainPath, "ufkt1.RData"))
                    load(file=paste0(mainPath, "d2fkt1.RData"))
                    load(file=paste0(mainPath, "HRgofkt1.RData")) 
                    showplot = "multiplicatively"
                 }else{
                    
                    load(file=paste0(mainPath, "ufkt2.RData"))
                    load(file=paste0(mainPath, "d2fkt2.RData"))
                    load(file=paste0(mainPath, "HRgofkt2.RData"))   
                    showplot = "additively"
                 }
                 
              }
              plot_ly(x=HRgofkt,y=d2fkt,z=ufkt, type="surface")  %>% 
                layout(title=
                         paste0("Optimization region of ", showplot, " adjusted setting"),
                       scene = list(
                         xaxis = list(title = "HRgo"),
                         yaxis = list(title = "d2"),
                         zaxis = list(title = "expected utility")))  
           }

           
        })   
        
     
})
