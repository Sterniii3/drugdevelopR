#' Optimal phase II/III drug development planning for multi-arm programs with 
#'  time-to-event endpoint
#'
#' The function \code{\link{optimal_multiarm}} of the drugdevelopR package 
#'  enables planning of multi-arm phase II/III drug development programs with 
#'  optimal sample size allocation and go/no-go decision rules 
#'  (Preussler et. al, 2019) for time-to-event endpoints. So far, only three-arm
#'   trials with two treatments and one control are supported. The assumed true
#'  treatment effects are assumed fixed (planning is also possible via
#'  user-friendly R Shiny App:
#'  \href{https://web.imbi.uni-heidelberg.de/multiarm/}{multiarm}). Fast
#'  computing is enabled by parallel programming.
#' 
#' @name optimal_multiarm
#' @inheritParams optimal_multiarm_generic
#' @inheritParams optimal_tte_generic
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param ec control arm event rate for phase II and III
#' @param beta type-II error rate for any pair, i.e. `1 - beta` is the
#'  (any-pair) power for calculation of the number of events for phase III
#' 
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom mvtnorm dmvnorm qmvnorm
#' @importFrom cubature adaptIntegrate
#' 
#' @return
#' `r optimal_return_doc(type = "tte", setting = "multiarm")`
#' 
#' @examples
#' \donttest{
#' # Activate progress bar (optional)
#' progressr::handlers(global = TRUE)
#' # Optimize
#' optimal_multiarm(hr1 = 0.75, hr2 = 0.80,    # define assumed true HRs 
#'   ec = 0.6,                                          # control arm event rate
#'   n2min = 30, n2max = 90, stepn2 = 6,                # define optimization set for n2
#'   hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05,     # define optimization set for HRgo
#'   alpha = 0.05, beta = 0.1,                          # drug development planning parameters
#'   c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,           # fixed/variable costs for phase II/III
#'   K = Inf, N = Inf, S = -Inf,                        # set constraints
#'   steps1 = 1,                                        # define lower boundary for "small"
#'   stepm1 = 0.95,                                     # "medium"
#'   stepl1 = 0.85,                                     # and "large" effect size categories
#'   b1 = 1000, b2 = 2000, b3 = 3000,                   # define expected benefit 
#'   strategy = 1,                                      # choose strategy: 1, 2 or 3
#'   num_cl = 1)                                        # number of cores for parallelized computing 
#'   }
#' 
#' @references
#' Preussler, S., Kirchner, M., Goette, H., Kieser, M. (2019). Optimal Designs for Multi-Arm Phase II/III Drug Development Programs. Submitted to peer-review journal.
#'
#' IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at \href{https://www.iqwig.de/ueber-uns/methoden/methodenpapier/}{https://www.iqwig.de/ueber-uns/methoden/methodenpapier/}, assessed last 15.05.19.
#' @export
optimal_multiarm <- function(hr1, hr2, ec,
                        n2min, n2max, stepn2,
                        hrgomin, hrgomax, stephrgo,
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

  HRGO <- seq(hrgomin, hrgomax, stephrgo)
  N2   <- seq(n2min, n2max, stepn2)

  if(strategy==1){STRATEGY = 1}
  if(strategy==2){STRATEGY = 2}
  if(strategy==3){STRATEGY = c(1, 2)}
  
  result <- NULL
  
  for(strategy in STRATEGY){
    
    ufkt <- spfkt <- pgofkt <- K2fkt <- K3fkt <-
      sp2fkt <- sp3fkt <- n3fkt <- matrix(0, length(N2), length(HRGO))
    
    pb <- progressr::progressor(steps = length(HRGO),
                                label = "Optimization progress",
                                message = "Optimization progress")
    pb(paste("Performing optimization for strategy", strategy),
       class = "sticky", amount = 0)
    
    for(j in 1:length(HRGO)){
      
      HRgo <- HRGO[j]
      
      cl <-  parallel::makeCluster(getOption("cl.cores", num_cl)) #define cluster
      
      parallel::clusterExport(cl, c("pmvnorm", "dmvnorm","qmvnorm","adaptIntegrate", "pgo_tte", "ss_tte", "Ess_tte",
                          "PsProg_tte", "alpha", "beta",
                          "steps1", "steps2", "stepm1", "stepm2", "stepl1", "stepl2",
                          "K", "N", "S", "strategy",
                          "c2", "c3", "c02", "c03",
                          "b1", "b2", "b3", "HRgo",
                          "hr1", "hr2", "ec"), envir = environment())
      
      
      res <- parallel::parSapply(cl, N2, utility_multiarm, HRgo,
                       alpha,beta,hr1,hr2,strategy,ec,
                       c2,c02,c3,c03,K,N,S,
                       steps1, stepm1, stepl1,b1, b2, b3)
      
      pb()
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
    
    result <-  rbind(result, data.frame(Strategy = strategy,u = round(Eud,2), HRgo = HRGO[J], n2 = N2[I], 
                          n3 = n3, n = N2[I] + n3,
                          pgo = round(pg,2), sProg = round(prob,2),
                          hr1 = hr1, hr2 = hr2, ec = ec,
                          K = K, N = N, S = S, K2 = round(k2), K3 = round(k3),
                          sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                          steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                          alpha = alpha, beta = beta, 
                          c02 = c02, c03 = c03, c2 = c2, c3 = c3, 
                          b1 = b1, b2 = b2, b3 = b3))  
  }


  comment(result) <-   c("\noptimization sequence HRgo:", HRGO,
                    "\noptimization sequence n2:", N2,
                    "\nonset date:", as.character(date),
                    "\nfinish date:", as.character(Sys.time()))
  
  return(result)
  
}
