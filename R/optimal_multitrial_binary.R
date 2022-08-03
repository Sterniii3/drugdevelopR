#' Optimal phase II/III drug development planning where several phase III trials are performed
#'
#' The drugdevelopR package enables planning of phase II/III drug development programs with optimal sample size allocation and go/no-go decision rules. For binary endpoints the treatment effect is measured by the risk ratio (RR).The assumed true treatment effects can be assumed fixed or modelled by a prior distribution. The R Shiny application \href{https://web.imbi.uni-heidelberg.de/prior/}{prior} visualizes the prior distributions used in this package. Fast coputing is enabled by parallel programming.
#' 
#' @name optimal_multitrial_binary
#' @param w weight for mixture prior distribution
#' @param p0 assumed true rate of control group
#' @param p11 assumed true rate of treatment group
#' @param p12 assumed true rate of treatment group
#' @param in1 amount of information for p11 in terms of sample size
#' @param in2 amount of information for p12 in terms of sample size
#' @param n2min minimal total sample size for phase II; must be even number
#' @param n2max maximal total sample size for phase II, must be even number
#' @param stepn2 stepsize for the optimization over n2; must be even number
#' @param rrgomin minimal threshold value for the go/no-go decision rule
#' @param rrgomax maximal threshold value for the go/no-go decision rule
#' @param steprrgo stepsize for the optimization over RRgo
#' @param beta 1-beta power for calculation of sample size for phase III
#' @param alpha significance level
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param b1 expected gain for effect size category "small"
#' @param b2 expected gain for effect size category "medium"
#' @param b3 expected gain for effect size category "large"
#' @param case choose case: "at least 1, 2 or 3 significant trials needed for approval"
#' @param strategy choose strategy: "conduct 1, 2, 3 or 4 trial"; TRUE calculates all strategies of the selected Case
#' @param fixed choose if true treatment effects are fixed or random, if TRUE p11 is used as fixed effect for p1
#' @param num_cl number of clusters used for parallel computing, default: 1
#' @return
#' The output of the function \code{\link{optimal_binary}} is a data.frame containing the optimization results:
#' \describe{
#'   \item{Case}{Case: "number of significant trials needed"}
#'   \item{Strategy}{Strategy: "number of conducted trials"}
#'   \item{u}{maximal expected utility}
#'   \item{RRgo}{optimal threshold value for the decision rule to go to phase III}
#'   \item{n2}{total sample size for phase II}
#'   \item{n3}{total sample size for phase III; rounded to the next even natural number}
#'   \item{n}{total sample size in the program; n = n2 + n3}
#'   \item{K}{maximal costs of the program}
#'   \item{pgo}{probability to go to phase III}
#'   \item{sProg}{probability of a successful program}
#'   \item{sProg1}{probability of a successful program with "small" treatment effect in Phase III}
#'   \item{sProg2}{probability of a successful program with "medium" treatment effect in Phase III}
#'   \item{sProg3}{probability of a successful program with "large" treatment effect in Phase III }
#'   \item{K2}{expected costs for phase II}
#'   \item{K3}{expected costs for phase III}
#' }
#' and further input parameters.
#'
#' Taking cat(comment()) of the data.frame object lists the used optimization sequences, start and finish date of the optimization procedure.
#' @examples
#' res <- optimal_multitrial_binary(w = 0.3,                # define parameters for prior
#'   p0 = 0.6, p11 =  0.3, p12 = 0.5, in1 = 30, in2 = 60,   # (https://web.imbi.uni-heidelberg.de/prior/)
#'   n2min = 20, n2max = 100, stepn2 = 4,                   # define optimization set for n2
#'   rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.05,         # define optimization set for RRgo
#'   alpha = 0.05, beta = 0.1,                              # drug development planning parameters
#'   c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,               # define fixed and variable costs for phase II and III,
#'   K = Inf, N = Inf, S = -Inf,                            # set constraints
#'   b1 = 1000, b2 = 2000, b3 = 3000,                       # define expected benefit for a "small", "medium" and "large" treatment effect
#    case = 1, strategy = TRUE,                             # chose Case and Strategy                                   
#'   fixed = FALSE,                                         # choose if true treatment effects are fixed or random
#'   num_cl = 1)                                            # set number of cores used for parallelized computing 
#' res
#' cat(comment(res))                                        # displays the optimization sequence, start/finish date of procedure
#' @section drugdevelopR functions:
#' The drugdevelopR package provides the functions
#' \itemize{
#'   \item \code{\link{optimal_tte}},
#'   \item \code{\link{optimal_binary}} or
#'   \item \code{\link{optimal_normal}}
#' }
#' to plan optimal phase II/III drug development programs with
#' \itemize{
#'   \item time-to-event (treatment effect measured by hazard ratio (HR)),
#'   \item binary (treatment effect measured by risk ratio (RR)) and
#'   \item normally distributed (treatment effect measured by standardized difference in means (Delta))
#' }
#' endpoint, where the treatment effect is modelled by a \href{https://web.imbi.uni-heidelberg.de/prior/}{prior}. Optimal phase II/III drug development planning with fixed treatment effects can be done with the help of the R Shiny application \href{https://web.imbi.uni-heidelberg.de/basic/}{basic}. Extensions are 
#' \itemize{
#'   \item optimal planning of programs including methods for discounting of phase II results (function: \code{\link{optimal_bias}}, App: \href{https://web.imbi.uni-heidelberg.de/bias/}{bias}),
#'   \item optimal planning of programs with several phase III trials (function: \code{\link{optimal_multitrial}}, App: \href{https://web.imbi.uni-heidelberg.de/multitrial/}{multitrial}) and
#'   \item optimal planning of programs with multiple arms (function: \code{\link{optimal_multiarm}}, App: \href{https://web.imbi.uni-heidelberg.de/multiarm/}{multiarm}).
#' }
#' @references
#' IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at \href{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}, assessed last 15.05.19.
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export

optimal_multitrial_binary <- function(w, p0, p11, p12, in1, in2,
                               n2min, n2max, stepn2,
                               rrgomin, rrgomax, steprrgo,
                               alpha, beta, 
                               c2, c3, c02, c03, 
                               K = Inf, N = Inf, S = -Inf,
                               b1, b2, b3,
                               case, strategy = TRUE,
                               fixed = FALSE,  num_cl = 1){
  
  result <- result23 <- NULL
  
  # spezifications for one phase III trial
  steps1 = 1; stepm1 = 0.95; stepl1 = 0.85
  steps2  <- stepm1
  stepm2  <- stepl1
  stepl2  <- 0
  gamma   <- 0
  ymin <- -log(0.8)
  
  
  alpha_in <- alpha
  
  date <- Sys.time()
  
  RRGO <- seq(rrgomin, rrgomax, steprrgo)
  N2   <- seq(n2min, n2max, stepn2)
  
  if(!is.numeric(strategy)){
    if(case==1){
      # Strategy 1alpha vs. Strategy 1/2,
      STRATEGY = c(1, 2)
    }
    if(case==2){
      # Strategy 1alpha^2 vs. Strategy 2/2 vs. Strategy 2/3 vs. Strategy 2/2( + 1)
      STRATEGY = c(1, 2, 3, 23)
    }
    if(case==3){
      # Strategy 1alpha^3 vs. Strategy 3/3 vs. Strategy 3/4
      STRATEGY = c(1, 3, 4)
    }  
  }else{
    STRATEGY = strategy
  }
  
  for(Strategy in STRATEGY){
    
    ufkt <- spfkt <- pgofkt <- K2fkt <- K3fkt <-
      sp1fkt <- sp2fkt <- sp3fkt <- n2fkt <- n3fkt <- pgo3fkt <- 
      n33fkt<-  sp13fkt   <- sp23fkt   <- sp33fkt <- matrix(0, length(N2), length(RRGO))
    
    cat("", fill = TRUE)
    cat(paste("Case ", case,": Optimization progess for Strategy ", Strategy), fill = TRUE)
    cat("", fill = TRUE)
    pb <- txtProgressBar(min = 0, max = length(RRGO), style = 3, label = "Optimization progess")
    
    for(j in 1:length(RRGO)){
      
      RRgo <- RRGO[j]
      
      cl <-  parallel::makeCluster(getOption("cl.cores", num_cl)) #define cluster
      
      ###################
      # Strategy 1alpha #
      ###################
      if(Strategy == 1){
        if(case==1){
          alpha <- alpha_in
        }
        if(case==2){
          alpha <- alpha_in^2
        }
        if(case==3){
          alpha <- alpha_in^3
        }  
      }else{
        alpha <- alpha_in
      }
      
      parallel::clusterExport(cl, c("pmvnorm", "dmvnorm", "prior_binary", "Epgo_binary", "Epgo23_binary", "En3_binary",
                          "EPsProg_binary", "EPsProg2_binary", "EPsProg3_binary", "EPsProg4_binary", "EPsProg23_binary",
                          "alpha", "beta",
                          "steps1", "steps2", "stepm1", "stepm2", "stepl1", "stepl2",
                          "K", "N", "S", "fixed",
                          "c2", "c3", "c02", "c03",
                          "b1", "b2", "b3", "w", "RRgo",
                          "p0", "p11", "p12", "in1", "in2"), envir = environment())
      
      if(Strategy==1){
        res <- parallel::parSapply(cl, N2, utility_binary, RRgo, w, p0, p11, p12, in1, in2,
                         alpha, beta, 
                         c2, c3, c02, c03, K, N, S,
                         steps1, stepm1, stepl1,
                         b1, b2, b3,
                         fixed)  
      }
      if(Strategy==2){
        res <- parallel::parSapply(cl, N2, utility2_binary, RRgo, w, p0, p11, p12, in1, in2,
                         alpha, beta, 
                         c2, c3, c02, c03, 
                         K, N, S,
                         b1, b2, b3,
                         case, fixed)  
      }
      if(Strategy==3){
        res <- parallel::parSapply(cl, N2, utility3_binary, RRgo, w, p0, p11, p12, in1, in2,
                         alpha, beta, 
                         c2, c3, c02, c03,
                         K, N, S,
                         b1, b2, b3,
                         case, fixed)  
      }
      if(Strategy==23){
        res <- parallel::parSapply(cl, N2, utility23_binary, RRgo, w, p0, p11, p12, in1, in2,
                         alpha, beta, 
                         c2, c3, c02, c03, 
                         b1, b2, b3)  
      }
      if(Strategy==4){
        res <- parallel::parSapply(cl, N2, utility4_binary, RRgo, w, p0, p11, p12, in1, in2,
                         alpha, beta, 
                         c2, c3, c02, c03, 
                         K, N, S,
                         b1, b2, b3,
                         case, fixed)  
      }
      
      setTxtProgressBar(title= "i", pb, j)
      parallel::stopCluster(cl)
      
      ufkt[, j]      <-  res[1, ]
      n3fkt[, j]     <-  res[2, ]
      spfkt[, j]     <-  res[3, ]
      pgofkt[, j]    <-  res[4, ]
      K2fkt[, j]     <-  res[5, ]
      K3fkt[, j]     <-  res[6, ]
      sp1fkt[, j]    <-  res[7, ]
      sp2fkt[, j]    <-  res[8, ]
      sp3fkt[, j]    <-  res[9, ]
      
      if(Strategy==23){
        pgo3fkt[, j]    <-  res[10, ]
        n33fkt[, j]     <-  res[11, ]
        sp13fkt[, j]    <-  res[12, ]
        sp23fkt[, j]    <-  res[13, ]
        sp33fkt[, j]    <-  res[14, ]
      }
      
    }
    
    ind   <-  which(ufkt  ==  max(ufkt), arr.ind <-  TRUE)
    
    I <-  as.vector(ind[1, 1])
    J <-  as.vector(ind[1, 2])
    
    Eud   <- ufkt[I, J]
    n3    <- n3fkt[I, J]
    prob  <- spfkt[I, J]
    pg    <- pgofkt[I, J]
    k2    <- K2fkt[I, J]
    k3    <- K3fkt[I, J]
    prob1 <- sp1fkt[I, J]
    prob2 <- sp2fkt[I, J]
    prob3 <- sp3fkt[I, J]
    
    if(Strategy==23){
      n33    <- n33fkt[I, J]
      pg3    <- pgo3fkt[I, J]
      prob13 <- sp13fkt[I, J]
      prob23 <- sp23fkt[I, J]
      prob33 <- sp33fkt[I, J]
      
    }else{
      n33    <- 0
      pg3    <- 0
      prob13 <- 0
      prob23 <- 0
      prob33 <- 0
    }
    
    
    close(pb)
    
    if(!fixed){
      
      result <-  rbind(result, data.frame(Case = case, Strategy = Strategy, 
                                          u = round(Eud,2), RRgo = RRGO[J], n2 = N2[I],
                                          n3 = n3, n = N2[I] + n3,
                                          pgo = round(pg,2), sProg = round(prob,2),
                                          w = w, p0 = p0, p11 = p11, p12 = p12, in1 = in1, in2 = in2,
                                          K = K, N = N, S = S, K2 = round(k2), K3 = round(k3),
                                          sProg1 = round(prob1,2), sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                                          steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                                          pgo3 = round(pg3,2), d33= d33, n33 = n33,
                                          sProg13 = round(prob13,2), sProg23 = round(prob23,2), sProg33 = round(prob33,2),
                                          alpha = alpha, beta = beta, c02 = c02,
                                          c03 = c03, c2 = c2, c3 = c3, b1 = b1, b2 = b2, b3 = b3))
    }else{
      
      result <-  rbind(result, data.frame(Case = case, Strategy = Strategy, 
                                          u = round(Eud,2), RRgo = RRGO[J], n2 = N2[I],
                                          n3 = n3, n = N2[I] + n3,
                                          pgo = round(pg,2), sProg = round(prob,2),
                                          p0 = p0, p1 = p11,
                                          K = K, N = N, S = S, K2 = round(k2), K3 = round(k3),
                                          sProg1 = round(prob1,2), sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                                          steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                                          alpha = alpha, beta = beta, c02 = c02,
                                          c03 = c03, c2 = c2, c3 = c3, b1 = b1, b2 = b2, b3 = b3))
    }
    
    
    
    comment(result) <-   c("\noptimization sequence RRgo:", RRGO,
                           "\noptimization sequence n2:", N2,
                           "\nset on date:", as.character(date),
                           "\nfinish date:", as.character(Sys.time()))
    
  }
  cat("", fill = TRUE)
  cat("", fill = TRUE)
  cat("Optimization result:", fill = TRUE)
  cat("", fill = TRUE)
  print(result)
  cat("", fill = TRUE)
  
  return(result)
  
}
