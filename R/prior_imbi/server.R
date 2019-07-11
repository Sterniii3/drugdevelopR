library("shiny")
library("ggplot2")
library("gridExtra")
library("mvtnorm")
library("doParallel")
library("parallel")
library("cubature")
library("plotly")
library("viridis")


mainPath <- "/opt/shiny-server/samplesizr/bias/"

#mainPath <- ""

# Prior distribution for theta
prior_tte<-function(x, w, hr1, hr2, id1, id2){
  w * dnorm(x, -log(hr1), sqrt(4/id1)) + (1 - w) * dnorm(x, -log(hr2), sqrt(4/id2))
}

# 10000 realizations of the prior distribution
box_tte<-function(w, hr1, hr2, id1, id2){
  w * rnorm(100000, -log(hr1),sqrt(4/id1)) + (1 - w) * rnorm(100000, -log(hr2), sqrt(4/id2))
}


# Prior distribution for p1
prior_binary<-function(x, w, p11, p12, in1, in2){
  w * dnorm(x, p11, sqrt(p11*(1-p11)/in1)) + (1-w) * dnorm(x, p12, sqrt(p12*(1-p12)/in2))
}

# 10000 realizations of the prior distribution
box_binary<-function(w, p11, p12, in1, in2){
  w*rnorm(100000,p11,sqrt(p11*(1-p11)/in1))+(1-w)*rnorm(100000,p12,sqrt(p12*(1-p12)/in2))
}


# Prior distribution for Delta
prior_normal <-function(x, w, Delta1, Delta2, in1, in2, a, b){
  w * dtnorm(x, Delta1, sqrt(4/in1), lower = a, upper = b) +
    (1-w) * dtnorm(x, Delta2, sqrt(4/in2), lower = a, upper = b)
}

# 10000 realizations of the prior distribution
box_normal <-function(w, Delta1, Delta2, in1, in2, a, b){
  w*rtnorm(100000,Delta1,sqrt(4/in1), lower=a, upper=b) +
    (1-w)*rtnorm(100000,Delta2,sqrt(4/in2), lower=a, upper=b)
}

col =  viridis(2)
col1 = viridis(2,alpha = 0.5)

shinyServer(function(input, output,session) {


     output$plot1 <- renderPlot({


        Select = input$Select

        if(Select==1){
           HR = input$HR

           n = input$n1

           w = input$w1

           gamma = input$gamma1
           
           x <- seq(-log(1.5), -log(0.1), length=100)
           hx <- prior_tte(x,w[1],HR[1],HR[2],n[1],n[2])
           
           if(gamma==0){
             
             hx1 <- prior_tte(x,w[2],HR[1],HR[2],n[1],n[2])
             
             y=x[round(hx+hx1,3)!=0]
             
             hy = hx[round(hx+hx1,3)!=0]
             hy1 = hx1[round(hx+hx1,3)!=0]
           
           }else{
             
           y=x[round(hx,3)!=0]
           hy = hx[round(hx,3)!=0]
}
           if(gamma==0){
             plot(exp(-y), hy, type="l", lty=2, xlab="HR",
                  ylab="Density", col = col[1], lwd=2, ylim=c(0,max(hx,hx1)))
             lines(exp(-y),hy1,col=col[2], lwd=2)
            legend(x="topleft",legend=w, col=col, lwd=2, lty=c(2,1),title="w", x.intersp=1, y.intersp=1, bty="n")
           
             }else{
               plot(exp(-y), hy, type="l",  xlab="HR",
                    ylab="Density", col = col[1], lwd=2, ylim=c(0,max(hx)))
               lines(exp(-y-gamma),hy,col=col1[1], lwd=2, lty = 3)
               legend(x="topleft",legend=c(0,gamma), col=c(col[1],col1[1]), lwd=2, lty=c(1,3),title=paste0("w = ",w[1]), x.intersp=1, y.intersp=1, bty="n")
             
           }
             
           
           
        }

        if(Select==2){

           p0 = input$p0
           p1 = input$p1

           n = input$n2

           w = input$w2
           
           gamma = input$gamma2

           x <- seq(0, 1, length=100)
           hx <- prior_binary(x,w[1],p1[1],p1[2],n[1],n[2])
           
           if(gamma==0){
             hx1 <- prior_binary(x,w[2],p1[1],p1[2],n[1],n[2])
             
             y=x[round(hx+hx1,3)!=0]
             hy = hx[round(hx+hx1,3)!=0]
             hy1 = hx1[round(hx+hx1,3)!=0]
             
             plot(y/p0, hy, type="l", lty=2, xlab="RR",
                  ylab="Density", col = col[1], lwd=2, ylim=c(0,max(hx,hx1)))
             lines(y/p0,hy1,col=col[2], lwd=2)
             legend(x="topleft",legend=w, col=col, lwd=2, lty=c(2,1),title="w", x.intersp=1, y.intersp=1, bty="n") 
           
             if(p0==1){
               plot(y, hy, type="l", lty=2, xlab=expression("p"[1]),
                    ylab="Density", col = col[1], lwd=2, ylim=c(0,max(hx,hx1)))
               lines(y,hy1,col=col[2], lwd=2)
               legend(x="topleft",legend=w, col=col, lwd=2, lty=c(2,1),title="w", x.intersp=1, y.intersp=1, bty="n")
             }
           }else{
               
             y=x[round(hx,3)!=0]
             hy = hx[round(hx,3)!=0]
             
             plot(y/p0, hy, type="l",  xlab="RR",
                  ylab="Density", col = col[1], lwd=2, ylim=c(0,max(hx)),lty = 1)
             
             lines((y-gamma)/p0,hy,col=col1[1], lwd=2, lty=3)
             
             legend(x="topleft",legend=c(0,gamma), col=c(col[1],col1[1]), lwd=2, lty=c(1,3),title=paste0("w = ",w[1]), x.intersp=1, y.intersp=1, bty="n")             
             
             if(p0==1){
               plot(y, hy, type="l", lty=1, xlab=expression("p"[1]),
                    ylab="Density", col = col[1], lwd=2, ylim=c(0,max(hx)))
               lines(y-gamma,hy,col=col1[1], lwd=2,lty=3)
               legend(x="topleft",legend=c(0,gamma), col=c(col[1],col1[1]), lwd=2, lty=c(1,3),title=paste0("w = ",w[1]), x.intersp=1, y.intersp=1, bty="n")             
               }

             
             }

        }


        if(Select==3){

           Delta = input$delta

           n = input$n3

           w = input$w3
           
           gamma = input$gamma3

           B = input$b

           x <- seq(0, 1, length=100)

           hx1 <- prior_normal(x,w[1],Delta[2],Delta[1],n[1],n[2],B[1],B[2])
           
           
          if(gamma==0){
            hx <- prior_normal(x,w[2],Delta[2],Delta[1],n[1],n[2],B[1],B[2])
            
            plot(x, hx, type="l", lty=1, xlab=expression(Delta),
                 ylab="Density", col = col[2], lwd=2, ylim=c(0,max(hx,hx1)), xlim=c(B[1]-0.05,B[2]+0.05))
            lines(x,hx1,col=col[1], lwd=2,lty=2)
            legend(x="topleft",legend=w, col=col, lwd=2, lty=c(2,1),title="w", x.intersp=1, y.intersp=1, bty="n")
            
          }else{
            plot(x, hx1, type="l", lty=1, xlab=expression(Delta),
                 ylab="Density", col = col[1], lwd=2, ylim=c(0,max(hx1)), xlim=c(B[1]-0.05,B[2]+0.05))
            lines(x+gamma,hx1,col=col1[1], lwd=2,lty=3)
            legend(x="topleft",legend=c(0,gamma), col=c(col[1],col1[1]), lwd=2, lty=c(1,3),title=paste0("w = ",w[1]), x.intersp=1, y.intersp=1, bty="n")
            
            
            
          }

           
        }


        }, height = 400, width = 600 )



     output$plot <- renderPlot({


        Select = input$Select

        if(Select==1){
           HR = input$HR

           n = input$n1

           w = input$w1
           
           gamma = input$gamma1

           x <- seq(-log(1.5), -log(0.1), length=100)
           px <- prior_tte(x,w[1],HR[1],HR[2],n[1],n[2])
           
           hx <- box_tte(w[1],HR[1],HR[2],n[1],n[2])
           
           if(gamma==0){
             
             px1 <- prior_tte(x,w[2],HR[1],HR[2],n[1],n[2])
             y=x[round(px+px1,3)!=0]
             
             hx1 <- box_tte(w[2],HR[1],HR[2],n[1],n[2])

             boxplot(exp(-hx), exp(-hx1),
                     names = c(paste0(round(median(exp(-hx)),2)), paste0(round(median(exp(-hx1)),2))),
                     xlab="HR", ylab="median HR", col = col,   ylim=c(min(exp(-y)),max(exp(-y))), horizontal = TRUE)
             
           }else{
             y=x[round(px,3)!=0]
             boxplot( exp(-hx), exp(-hx-gamma),
                     names = c( paste0(round(median(exp(-hx)),2)),paste0(round(median(exp(-hx-gamma)),2))),
                     xlab="HR", ylab="median HR", col = c(col[1],col1[1]),   ylim=c(min(exp(-y)),max(exp(-y))), 
                     horizontal = TRUE)
             
           }

        }

        if(Select==2){

           p0 = input$p0
           p1 = input$p1

           n = input$n2

           w = input$w2

           gamma = input$gamma2
           
           x <- seq(0, 1, length=100)
           px <- prior_binary(x,w[1],p1[1],p1[2],n[1],n[2])
           px1 <- prior_binary(x,w[2],p1[1],p1[2],n[1],n[2])
           y=x[round(px+px1,3)!=0]

           
           hx <- box_binary(w[1],p1[1],p1[2],n[1],n[2])
           
           if(gamma==0){

           hx1 <- box_binary(w[2],p1[1],p1[2],n[1],n[2])

           boxplot(hx/p0,hx1/p0, names=c( paste0(round(median(hx/p0),2)),paste0(round(median(hx1/p0),2))), 
                   horizontal = TRUE,ylab = "median RR",
                xlab="RR", col = col,  ylim=c(min(y/p0),max(y/p0)))

           if(p0==1){
              boxplot(hx,hx1, names=c( paste0(round(median(hx),2)),paste0(round(median(hx1),2))), horizontal = TRUE,
                      ylab= expression("median p"[1]),
                      xlab=expression("p"[1]), col = col,  ylim=c(min(y),max(y)))
           }
           }else{
          
             
             y=x[round(px,3)!=0]
             boxplot( hx/p0,(hx-gamma)/p0,
                      names = c( paste0(round(median(hx/p0),2)),paste0(round(median((hx-gamma)/p0),2))),
                      xlab="RR", ylab="median RR", col = c(col[1],col1[1]),  ylim=c(min(y/p0),max(y/p0)), 
                      horizontal = TRUE)
            
             
             if(p0==1){
               boxplot( hx/p0,(hx-gamma)/p0,
                        names = c( paste0(round(median(hx/p0),2)),paste0(round(median((hx-gamma)/p0),2))),
                        ylab= expression("median p"[1]),
                        xlab=expression("p"[1]),
                         col = c(col[1],col1[1]),  ylim=c(min(y/p0),max(y/p0)), 
                        horizontal = TRUE)
             } 
             
        }

        }


        if(Select==3){

           Delta = input$delta

           n = input$n3

           w = input$w3

           gamma = input$gamma3
           
           B = input$b

           x <- seq(0, 1, length=100)
           hx1 <- box_normal(w[1],Delta[2],Delta[1],n[1],n[2],B[1],B[2])  
           
           if(gamma==0){
             hx <- box_normal(w[2],Delta[2],Delta[1],n[1],n[2],B[1],B[2])

             boxplot(hx1, hx, 
                     names = c( paste0(round(median(hx1),2)),paste0(round(median(hx),2))),
                     ylab= expression("median "~Delta),
                     xlab=expression(Delta),
                     horizontal = TRUE,
                      col = c(col[1],col[2]), ylim=c(B[1]-0.05,B[2]+0.05))
             
           }else{
             
      
             
             boxplot(hx1, hx1+gamma, 
                     names = c( paste0(round(median(hx1),2)),paste0(round(median(hx1+gamma),2))),
                     ylab= expression("median "~Delta),
                     xlab=expression(Delta),
                     horizontal = TRUE,
                      col = c(col[1],col1[1]), ylim=c(B[1]-0.05,B[2]+0.05))  
             
           }
           






        }


     }, height = 400, width = 600 )

})
