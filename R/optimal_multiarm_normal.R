#' Optimal phase II/III drug development planning for multi-arm programs with 
#'  normally distributed endpoint
#'
#' The \code{\link{optimal_multiarm_normal}} function enables planning of
#'  multi-arm phase II/III drug development programs with optimal sample size 
#'  allocation and go/no-go decision rules. For normally distributed endpoints,
#'  the treatment effect is measured by the standardized difference in means
#'  (Delta). So far, only three-arm trials with two treatments and one control
#'  are supported. The assumed true treatment effects can be assumed fixed or
#'  modelled by a prior distribution. The R Shiny application
#'  \href{https://web.imbi.uni-heidelberg.de/prior/}{prior} visualizes the
#'  prior distributions used in this package. Fast computing is enabled by
#'  parallel programming.
#'  
#' @name optimal_multiarm_normal
#' @inheritParams optimal_multiarm_generic
#' @inheritParams optimal_normal_generic
#' @param Delta1 assumed true treatment effect as the standardized difference in
#'  means for treatment arm 1
#' @param Delta2 assumed true treatment effect as the standardized difference in
#'  means for treatment arm 2
#'  
#' @return
#' `r optimal_return_doc(type = "normal", setting = "multiarm")`
#' 
#' @importFrom progressr progressor
#' 
#' @examples
#' # Activate progress bar (optional)
#' \dontrun{progressr::handlers(global = TRUE)}
#' # Optimize
#' \donttest{
#' optimal_multiarm_normal(Delta1 = 0.375, Delta2 = 0.625,     
#'   n2min = 20, n2max = 100, stepn2 = 4,                 # define optimization set for n2
#'   kappamin = 0.02, kappamax = 0.2, stepkappa = 0.02,   # define optimization set for kappa
#'   alpha = 0.025, beta = 0.1,                           # drug development planning parameters
#'   c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,           # fixed/variable costs for phase II/III
#'   K = Inf, N = Inf, S = -Inf,                          # set constraints
#'   steps1 = 0,                                          # define lower boundary for "small"
#'   stepm1 = 0.5,                                        # "medium"
#'   stepl1 = 0.8,                                        # and "large" effect size categories
#'   b1 = 3000, b2 = 8000, b3 = 10000,                    # define expected benefits 
#'   strategy = 1,
#'   num_cl = 1)                                          # number of cores for parallelized computing 
#'   }
#' @references
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences.
#'
#' @export


optimal_multiarm_normal <- function(Delta1, Delta2,
                             n2min, n2max, stepn2,
                             kappamin, kappamax, stepkappa,
                             alpha, beta,
                             c2, c3, c02, c03, 
                             K = Inf, N = Inf, S = -Inf,
                             steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                             b1, b2, b3,
                             strategy, num_cl = 1){
  
  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- Inf
  
  date <- Sys.time()
  
  KAPPA <- seq(kappamin, kappamax, stepkappa)
  N2   <- seq(n2min, n2max, stepn2)
  
  if(strategy==1){STRATEGY = 1}
  if(strategy==2){STRATEGY = 2}
  if(strategy==3){STRATEGY = c(1, 2)}
  
  result <- NULL
  kappa <- NA_real_
  strategy <- NA_real_
  cl <-  parallel::makeCluster(getOption("cl.cores", num_cl)) #define cluster
  
  parallel::clusterExport(cl, c("pmvnorm", "dmvnorm","qmvnorm","adaptIntegrate", "pgo_normal", "ss_normal", "Ess_normal",
                                "PsProg_normal", "alpha", "beta",
                                "steps1", "steps2", "stepm1", "stepm2", "stepl1", "stepl2",
                                "K", "N", "S", "strategy",
                                "c2", "c3", "c02", "c03",
                                "b1", "b2", "b3", "KAPPA",
                                "Delta1", "Delta2"), envir = environment())
  
  for(strategy in STRATEGY){
    
    ufkt <- spfkt <- pgofkt <- K2fkt <- K3fkt <-
      sp2fkt <- sp3fkt <- n3fkt <- matrix(0, length(N2), length(KAPPA))
    
    pb <- progressr::progressor(steps = length(KAPPA),
                                label = "Optimization progress",
                                message = "Optimization progress")
    pb(paste("Performing optimization for strategy", strategy),
       class = "sticky", amount = 0)
    
    for(j in 1:length(KAPPA)){
      
      kappa <- KAPPA[j]
      
      
      
      
      res <- parallel::parSapply(cl, N2, utility_multiarm_normal, kappa,
                       alpha,beta, Delta1,Delta2,strategy,
                       c2,c02,c3,c03,K,N,S,
                       steps1, stepm1, stepl1,b1, b2, b3)
      
      pb()
      
      
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
    
    result <-  rbind(result, data.frame(Strategy = strategy,u = round(Eud,2), Kappa = KAPPA[J], n2 = N2[I], 
                                        n3 = n3, n = N2[I] + n3,
                                        pgo = round(pg,2), sProg = round(prob,2),
                                        Delta1 = Delta1, Delta2 = Delta2, 
                                        K = K, N = N, S = S, K2 = round(k2), K3 = round(k3),
                                        sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                                        steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                                        alpha = alpha, beta = beta, 
                                        c02 = c02, c03 = c03, c2 = c2, c3 = c3, 
                                        b1 = b1, b2 = b2, b3 = b3))  
  }
  
  
  comment(result) <-   c("\noptimization sequence Kappa:", KAPPA,
                         "\noptimization sequence n2:", N2,
                         "\nonset date:", as.character(date),
                         "\nfinish date:", as.character(Sys.time()))
  class(result) <- c("drugdevelopResult", class(result))
  
  parallel::stopCluster(cl)
  return(result)
  
}
