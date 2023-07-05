#' Optimal phase II/III drug development planning for multi-arm programs with 
#'  binary endpoint
#'
#' The \code{\link{optimal_multiarm_binary}} function enables planning of
#'  multi-arm phase II/III drug
#'  development programs with optimal sample size allocation and go/no-go 
#'  decision rules. For binary endpoints the treatment effect is measured by the 
#'  risk ratio (RR). So far, only three-arm trials with two treatments and one
#'  control are supported. The assumed true treatment effects can be assumed fixed
#'  or modelled by a prior distribution. The R Shiny application
#'  \href{https://web.imbi.uni-heidelberg.de/prior/}{prior} visualizes the prior 
#'  distributions used in this package. Fast computing is enabled by parallel 
#'  programming.
#' 
#' @name optimal_multiarm_binary
#' @inheritParams optimal_multiarm_generic
#' @inheritParams optimal_binary_generic
#' @param p0 assumed true rate of the control group
#' @param p11 assumed true rate of the treatment arm 1
#' @param p12 assumed true rate of treatment arm 2
#' 
#' @return
#' `r optimal_return_doc(type = "binary", setting = "multiarm")`
#' 
#' @examples
#' \donttest{optimal_multiarm_binary( p0 = 0.6, 
#'   p11 =  0.3, p12 = 0.5, 
#'   n2min = 20, n2max = 100, stepn2 = 4,               # define optimization set for n2
#'   rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.05,     # define optimization set for RRgo
#'   alpha = 0.05, beta = 0.1,                          # drug development planning parameters
#'   c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,           # fixed/variable costs for phase II/III
#'   K = Inf, N = Inf, S = -Inf,                        # set constraints
#'   steps1 = 1,                                        # define lower boundary for "small"
#'   stepm1 = 0.95,                                     # "medium"
#'   stepl1 = 0.85,                                     # and "large" effect size categories
#'   b1 = 1000, b2 = 2000, b3 = 3000,                   # define expected benefits 
#'   strategy = 1, num_cl = 1)                          # number of cores for parallelized computing 
#'   }
#' @references
#' IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at \href{https://www.iqwig.de/ueber-uns/methoden/methodenpapier/}{https://www.iqwig.de/ueber-uns/methoden/methodenpapier/}, assessed last 15.05.19.
#' 
#' @export
optimal_multiarm_binary <- function(p0, p11, p12, 
                             n2min, n2max, stepn2,
                             rrgomin, rrgomax, steprrgo,
                             alpha, beta,
                             c2, c3, c02, c03, 
                             K = Inf, N = Inf, S = -Inf,
                             steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                             b1, b2, b3,
                             strategy, num_cl = 1){
  
  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 0
  
  date <- Sys.time()
  
  RRGO <- seq(rrgomin, rrgomax, steprrgo)
  N2   <- seq(n2min, n2max, stepn2)
  
  if(strategy==1){STRATEGY = 1}
  if(strategy==2){STRATEGY = 2}
  if(strategy==3){STRATEGY = c(1, 2)}
  
  result <- NULL
  
  for(strategy in STRATEGY){
    
    ufkt <- spfkt <- pgofkt <- K2fkt <- K3fkt <-
      sp2fkt <- sp3fkt <- n3fkt <- matrix(0, length(N2), length(RRGO))
    
    cat("Optimization progress:", fill = TRUE)
    cat("", fill = TRUE)
    pb <- txtProgressBar(min = 0, max = length(RRGO), style = 3, label = "Optimization progess")
    
    for(j in 1:length(RRGO)){
      
      RRgo <- RRGO[j]
      
      cl <-  parallel::makeCluster(getOption("cl.cores", num_cl)) #define cluster
      
      parallel::clusterExport(cl, c("pmvnorm", "dmvnorm","qmvnorm","adaptIntegrate", "pgo_binary", "ss_binary", "Ess_binary",
                          "PsProg_binary", "alpha", "beta",
                          "steps1", "steps2", "stepm1", "stepm2", "stepl1", "stepl2",
                          "K", "N", "S", "strategy",
                          "c2", "c3", "c02", "c03",
                          "b1", "b2", "b3", "RRgo",
                          "p0", "p11", "p12"), envir = environment())
      
      
      res <- parallel::parSapply(cl, N2, utility_multiarm_binary, RRgo,
                       alpha,beta,p0,p11,p12,strategy,
                       c2,c02,c3,c03,K,N,S,
                       steps1, stepm1, stepl1,b1, b2, b3)
      
      setTxtProgressBar(title= "i", pb, j)
      parallel::stopCluster(cl)
      
      ufkt[, j]     <-  res[1, ]
      n3fkt[, j]    <-  res[2, ]
      spfkt[, j]    <-  res[3, ]
      pgofkt[, j]   <-  res[4, ]
      sp2fkt[, j]   <-  res[5, ]
      sp3fkt[, j]   <-  res[6, ]
      K2fkt[, j]    <-  res[7, ]
      K3fkt[, j]    <-  res[8, ]
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
    prob2 <- sp2fkt[I, J]
    prob3 <- sp3fkt[I, J]
    
    result <-  rbind(result, data.frame(Strategy = strategy,u = round(Eud,2), RRgo = RRGO[J], n2 = N2[I], 
                                        n3 = n3, n = N2[I] + n3,
                                        pgo = round(pg,2), sProg = round(prob,2),
                                        p0 = p0, p11 = p11, p12 = p12,
                                        K = K, N = N, S = S, K2 = round(k2), K3 = round(k3),
                                        sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                                        steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                                        alpha = alpha, beta = beta, 
                                        c02 = c02, c03 = c03, c2 = c2, c3 = c3, 
                                        b1 = b1, b2 = b2, b3 = b3))  
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

