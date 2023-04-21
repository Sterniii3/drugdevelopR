#' Optimal phase II/III drug development planning when discounting phase II results with binary endpoint
#'
#' The function \code{\link{optimal_bias_binary}} of the drugdevelopR package enables planning of phase II/III drug development programs with optimal sample size allocation and go/no-go decision rules including methods for discounting of phase II results for binary endpoints (Preussler et. al, 2020). 
#' The discounting may be necessary as programs that proceed to phase III can be overoptimistic about the treatment effect (i.e. they are biased).
#' The assumed true treatment effects can be assumed fixed or modelled by a prior distribution.
#' The R Shiny application \href{https://web.imbi.uni-heidelberg.de/prior/}{prior} visualizes the prior distributions used in this package. 
#' Fast computing is enabled by parallel programming.
#' 
#' @name optimal_bias_binary
#' @inheritParams optimal_binary_generic
#' @inheritParams optimal_bias_generic
#' 
#' @return
#' `r optimal_return_doc(type = "binary", setting = "bias")`
#'
#' @examples
#' res <- optimal_bias_binary(w = 0.3,                   # define parameters for prior
#'   p0 = 0.6, p11 =  0.3, p12 = 0.5,
#'    in1 = 30, in2 = 60,                                # (https://web.imbi.uni-heidelberg.de/prior/)
#'   n2min = 20, n2max = 100, stepn2 = 10,               # define optimization set for n2
#'   rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.05,      # define optimization set for RRgo
#'   adj = "both",                                       # choose type of adjustment
#'   alpha = 0.05, beta = 0.1,                           # drug development planning parameters
#'   lambdamin = 0.2, lambdamax = 1, steplambda = 0.05,  # define optimization set for lambda
#'   alphaCImin = 0.025, alphaCImax = 0.5,
#'   stepalphaCI = 0.025,                                # define optimization set for alphaCI
#'   c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,            # fixed and variable costs for phase II/III
#'   K = Inf, N = Inf, S = -Inf,                         # set constraints
#'   steps1 = 1,                                         # define lower boundary for "small"
#'   stepm1 = 0.95,                                      # "medium"
#'   stepl1 = 0.85,                                      # and "large" effect size categories
#'   b1 = 1000, b2 = 2000, b3 = 3000,                    # define expected benefits
#'   fixed = TRUE,                                       # true treatment effects are fixed/random
#'   num_cl = 1)                                         # number of cores for parallelized computing
#' res
#' cat(comment(res))                                       
#' 
#' @references
#' IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at \href{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}, assessed last 15.05.19.
#' @editor Johannes Cepicka
#' @editDate 2022-04-23
#' @export

optimal_bias_binary <- function(w, p0, p11, p12, in1, in2,
                           n2min, n2max, stepn2,
                           rrgomin, rrgomax, steprrgo,
                           adj = "both",
                           lambdamin = NULL, lambdamax = NULL, steplambda = NULL,
                           alphaCImin = NULL, alphaCImax = NULL, stepalphaCI = NULL,
                           alpha, beta, 
                           c2, c3, c02, c03, 
                           K = Inf, N = Inf, S = -Inf,
                           steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                           b1, b2, b3,
                           fixed = FALSE, num_cl = 1){
  
  date <- Sys.time()
  
  result <- NULL
  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 0
  
  
  RRGO <- seq(rrgomin, rrgomax, steprrgo)
  N2   <- seq(n2min, n2max, stepn2)
  
  if(adj=="both"){
    STRATEGY = c(1,2)
  }
  if(adj=="multiplicative"){
    STRATEGY = 1
  }  
  if(adj=="additive"){
    STRATEGY = 2
  } 
  if(adj=="all"){
    STRATEGY = c(1,2,3,4)
  }
  
  
  for (strategy in STRATEGY){
    
    calresults <- NULL
    
    if(strategy == 1|strategy==3){
      proz <- "multiplicative"
      ADJ <- seq(lambdamin, lambdamax, steplambda)
    }
    if(strategy == 2|strategy==4){
      proz <- "additive"
      ADJ <- seq(alphaCImin, alphaCImax, stepalphaCI)
    }
    
    cat("", fill = TRUE)
    cat(paste("Optimization progess for adjustment type ", proz), fill = TRUE)
    cat("", fill = TRUE)
    pb <- txtProgressBar(min = 0, max = length(ADJ), style = 3, label = "Optimization progess") 
    
    for(a in 1:length(ADJ)){
      
      Adj <- ADJ[a]
  
      ufkt <- n3fkt <- spfkt <- pgofkt <- K2fkt <- K3fkt <-
        sp1fkt <- sp2fkt <- sp3fkt  <- matrix(0, length(N2), length(RRGO))
  
  
  
  for(j in 1:length(RRGO)){
    
    RRgo <- RRGO[j]
    
    cl <-  parallel::makeCluster(getOption("cl.cores", num_cl)) #define cluster
    
    parallel::clusterExport(cl, c("pmvnorm", "dmvnorm", "prior_binary","Epgo_binary", "En3_binary_L",
                        "EPsProg_binary_L","Epgo_binary_L2", "En3_binary_L2",
                        "EPsProg_binary_L2","En3_binary_R", "EPsProg_binary_R", "Epgo_binary_R2", "En3_binary_R2",
                        "EPsProg_binary_R2", "t1", "t2", "t3", "alpha", "beta",
                        "steps1", "steps2", "stepm1", "stepm2", "stepl1", "stepl2",
                        "K", "N", "S", "fixed",
                        "c2", "c3", "c02", "c03",
                        "b1", "b2", "b3", "w", "RRgo", "Adj",
                        "p0", "p11", "p12", "in1", "in2"), envir=environment())
    
    if(strategy == 1){
      strat = "multipl."
      res <- parallel::parSapply(cl, N2, utility_binary_R, RRgo, Adj, w, p0, p11, p12, in1, in2,
                       alpha, beta, 
                       c2, c3, c02, c03, 
                       K, N, S,
                       steps1, stepm1, stepl1,
                       b1, b2, b3,
                       fixed)  
    }
    if(strategy == 2){
      strat = "add."
      res <- parallel::parSapply(cl, N2, utility_binary_L, RRgo, Adj, w, p0, p11, p12, in1, in2,
                       alpha, beta, 
                       c2, c3, c02, c03, 
                       K, N, S,
                       steps1, stepm1, stepl1,
                       b1, b2, b3,
                       fixed)  
    }
    if(strategy == 3){
      strat = "multipl2."
      res <- parallel::parSapply(cl, N2, utility_binary_R2, RRgo, Adj, w, p0, p11, p12, in1, in2,
                       alpha, beta, 
                       c2, c3, c02, c03, 
                       K, N, S,
                       steps1, stepm1, stepl1,
                       b1, b2, b3,
                       fixed)  
    }
    if(strategy == 4){
      strat = "add2."
      res <- parallel::parSapply(cl, N2, utility_binary_L2, RRgo, Adj, w, p0, p11, p12, in1, in2,
                       alpha, beta, 
                       c2, c3, c02, c03, 
                       K, N, S,
                       steps1, stepm1, stepl1,
                       b1, b2, b3,
                       fixed)  
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
  
  if(fixed){
    calresult <-  data.frame(Method= strat,
                          u = round(Eud,2), Adj = Adj, RRgo = RRGO[J], n2 = N2[I],
                          n3 = n3, n = N2[I] + n3,
                          pgo = round(pg,2), sProg = round(prob,2),
                          p0 = p0, p1 = p11, 
                          K = K, K2 = round(k2), K3 = round(k3),
                          sProg1 = round(prob1,2), sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                          steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                          alpha = alpha, beta = beta, c02 = c02,
                          c03 = c03, c2 = c2, c3 = c3, b1 = b1, b2 = b2, b3 = b3)  
  }else{
    calresult <-  data.frame(Method= strat,
                          u = round(Eud,2), Adj = Adj, RRgo = RRGO[J], n2 = N2[I],
                          n3 = n3, n = N2[I] + n3,
                          pgo = round(pg,2), sProg = round(prob,2),
                          w = w, p0 = p0, p11 = p11, p12 = p12, in1 = in1, in2 = in2,
                          K = K, K2 = round(k2), K3 = round(k3),
                          sProg1 = round(prob1,2), sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                          steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                          alpha = alpha, beta = beta, c02 = c02,
                          c03 = c03, c2 = c2, c3 = c3, b1 = b1, b2 = b2, b3 = b3)
  }
  
  calresults <- rbind(calresults, calresult)
  
    }
    
    index   <- which(calresults$u == max(calresults$u))
    
    result <- rbind(result, calresults[index,] ) 
  }
  
  
  comment(result) <-   c("\noptimization sequence RRgo:", RRGO,
                         "\noptimization sequence n2:", N2,
                         "\nset on date:", as.character(date),
                         "\nfinish date:", as.character(Sys.time()))
  close(pb)
  
  cat("", fill = TRUE)
  cat("", fill = TRUE)
  cat("Optimization result:", fill = TRUE)
  cat("", fill = TRUE)
  print(result)
  cat("", fill = TRUE)
  
  
  
  return(result)
}

