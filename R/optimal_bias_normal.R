#' Optimal phase II/III drug development planning when discounting phase II results with normally distributed endpoint
#'
#' The function \code{\link{optimal_bias_normal}} of the drugdevelopR package enables planning of phase II/III drug development programs with optimal sample size allocation and go/no-go decision rules including methods for discounting of phase II results for normally distributed endpoints (Preussler et. al, 2020). 
#' The discounting may be necessary as programs that proceed to phase III can be overoptimistic about the treatment effect (i.e. they are biased).
#' The assumed true treatment effects can be assumed fixed or modelled by a prior distribution.
#' The R Shiny application \href{https://web.imbi.uni-heidelberg.de/prior/}{prior} visualizes the prior distributions used in this package. 
#' Fast computing is enabled by parallel programming.
#' 
#' @name optimal_bias_normal
#' @inheritParams optimal_normal_generic
#' @inheritParams optimal_bias_generic
#' 
#' @return
#' `r optimal_return_doc(type = "normal", setting = "bias")`
#' 
#' @importFrom progressr progressor
#'
#' @examples
#' # Activate progress bar (optional)
#' \dontrun{progressr::handlers(global = TRUE)}
#' # Optimize
#' \donttest{
#' optimal_bias_normal(w=0.3,             # define parameters for prior
#'   Delta1 = 0.375, Delta2 = 0.625, in1=300, in2=600,    # (https://web.imbi.uni-heidelberg.de/prior/)
#'   a = 0.25, b = 0.75,
#'   n2min = 20, n2max = 100, stepn2 = 10,                # define optimization set for n2
#'   kappamin = 0.02, kappamax = 0.2, stepkappa = 0.02,   # define optimization set for kappa
#'   adj = "both",                                        # choose type of adjustment
#'   lambdamin = 0.2, lambdamax = 1, steplambda = 0.05,   # define optimization set for lambda
#'   alphaCImin = 0.025, alphaCImax = 0.5,
#'   stepalphaCI = 0.025,                                 # define optimization set for alphaCI
#'   alpha = 0.025, beta = 0.1,                           # drug development planning parameters
#'   c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,           # fixed and variable costs for phase II/III
#'   K = Inf, N = Inf, S = -Inf,                          # set constraints
#'   steps1 = 0,                                          # define lower boundary for "small"
#'   stepm1 = 0.5,                                        # "medium"
#'   stepl1 = 0.8,                                        # and "large" effect size categories
#'   b1 = 3000, b2 = 8000, b3 = 10000,                    # define expected benefits
#'   fixed = TRUE,                                        # true treatment effects are fixed/random
#'   num_cl = 1)                                          # number of coresfor parallelized computing 
#'   }                                        
#' @references
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences.
#'
#' @export

optimal_bias_normal <- function(w, Delta1, Delta2, in1, in2, a, b,
                           n2min, n2max, stepn2,
                           kappamin, kappamax, stepkappa,
                           adj = "both",
                           lambdamin = NULL, lambdamax = NULL, steplambda = NULL,
                           alphaCImin = NULL, alphaCImax = NULL, stepalphaCI = NULL,
                           alpha, beta, 
                           c2, c3, c02, c03, 
                           K = Inf, N = Inf, S = -Inf,
                           steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                           b1, b2, b3,
                           fixed = FALSE,  num_cl = 1){
  
  result <- NULL
  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- Inf
  
  
  date <- Sys.time()
  
  KAPPA    <- seq(kappamin, kappamax, stepkappa)
  N2      <- seq(n2min, n2max, stepn2)
  
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
    
    pb <- progressr::progressor(steps = length(ADJ)*length(KAPPA), label = "Optimization progress", message = "Optimization progress")
    pb(paste("Performing optimization for adjustment type", proz), class = "sticky", amount = 0)
    Adj <- NA_real_
    kappa <- NA_real_
    cl <-  parallel::makeCluster(getOption("cl.cores", num_cl)) #define cluster
    
    parallel::clusterExport(cl, c("pmvnorm", "dmvnorm", "dtnorm", "prior_normal","Epgo_normal", "En3_normal_L",
                                  "EPsProg_normal_L","Epgo_normal_L2", "En3_normal_L2",
                                  "EPsProg_normal_L2","En3_normal_R", "EPsProg_normal_R", "Epgo_normal_R2", "En3_normal_R2",
                                  "EPsProg_normal_R2", "alpha", "beta",
                                  "steps1", "steps2", "stepm1", "stepm2", "stepl1", "stepl2",
                                  "K", "N", "S", "fixed",
                                  "c2", "c3", "c02", "c03",
                                  "b1", "b2", "b3", "w", "kappa", "Adj",
                                  "Delta1", "Delta2", "in1", "in2", "a", "b"), envir=environment())
    
    for(l in 1:length(ADJ)){
      
      Adj <- ADJ[l]
      
      ufkt <- n3fkt <- spfkt <- pgofkt <- K2fkt <- K3fkt <-
        sp1fkt <- sp2fkt <- sp3fkt <- matrix(0, length(N2), length(KAPPA))
      
      for(j in 1:length(KAPPA)){
        
        kappa <- KAPPA[j]
        
        
        
        if(strategy == 1){
          strat = "multipl."
          res <- parallel::parSapply(cl, N2, utility_normal_R, kappa, Adj, w, Delta1, Delta2, in1, in2, a, b,
                           alpha, beta, 
                           c2, c3, c02, c03, 
                           K, N, S,
                           steps1, stepm1, stepl1,
                           b1, b2, b3,
                           fixed)  
        }
        if(strategy == 2){
          strat = "add."
          res <- parallel::parSapply(cl, N2, utility_normal_L, kappa, Adj, w, Delta1, Delta2, in1, in2, a, b,
                           alpha, beta, 
                           c2, c3, c02, c03, 
                           K, N, S,
                           steps1, stepm1, stepl1,
                           b1, b2, b3,
                           fixed)  
        }
        if(strategy == 3){
          strat = "multipl2."
          res <- parallel::parSapply(cl, N2, utility_normal_R2, kappa, Adj, w, Delta1, Delta2, in1, in2, a, b, 
                           alpha, beta,
                           c2, c3, c02, c03, 
                           K, N, S,
                           steps1, stepm1, stepl1,
                           b1, b2, b3,
                           fixed)  
        }
        if(strategy == 4){
          strat = "add2."
          res <- parallel::parSapply(cl, N2, utility_normal_L2, kappa, Adj, w, Delta1, Delta2, in1, in2, a, b,
                           alpha, beta, 
                           c2, c3, c02, c03, 
                           K, N, S,
                           steps1, stepm1, stepl1,
                           b1, b2, b3,
                           fixed)  
        }
        
        pb()
        
        
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
                              u = round(Eud,2), Adj = Adj, Kappa = KAPPA[J], n2 = N2[I],
                              n3 = n3, n = N2[I] + n3,
                              pgo = round(pg,2), sProg = round(prob,2),
                              Delta = Delta1,
                              K = K, K2 = round(k2), K3 = round(k3),
                              sProg1 = round(prob1,2), sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                              steps1 = steps1, stepm1 = stepm1, stepl1 = stepl1,
                              alpha = alpha, beta = beta, c02 = c02,
                              c03 = c03, c2 = c2, c3 = c3, b1 = b1, b2 = b2, b3 = b3)
      }else{
        calresult <-  data.frame(Method= strat,
                              u = round(Eud,2), Adj = Adj, Kappa = KAPPA[J], n2 = N2[I],
                              n3 = n3, n = N2[I] + n3,
                              pgo = round(pg,2), sProg = round(prob,2),
                              w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                              K = K, K2 = round(k2), K3 = round(k3),
                              sProg1 = round(prob1,2), sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                              steps1 = steps1, stepm1 = stepm1, stepl1 = stepl1,
                              alpha = alpha, beta = beta, c02 = c02,
                              c03 = c03, c2 = c2, c3 = c3, b1 = b1, b2 = b2, b3 = b3)
      }
      
      calresults <- rbind(calresults, calresult)
      
    }
    
    index   <- which(calresults$u == max(calresults$u))
    
    result <- rbind(result, calresults[index,] ) 
  }
  
  
  
  comment(result) <-   c("\noptimization sequence kappa:", KAPPA,
                         "\noptimization sequence n2:", N2,
                         "\nonset date:", as.character(date),
                         "\nfinish date:", as.character(Sys.time()))
  class(result) <- c("drugdevelopResult", class(result))
  
  parallel::stopCluster(cl)
  
  return(result)
  
}

