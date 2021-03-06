#' Optimal phase II/III drug development planning for multi-arm programs
#'
#' The function \code{\link{optimal_multiarm}} of the drugdevelopR package enables planning of multi-arm phase II/III drug development programs with optimal sample size allocation and go/no-go decision rules (Preussler et. al, 2019). The assumed true treatment effects are assumed fixed (planning is also possible via user friendly R Shiny App: \href{https://web.imbi.uni-heidelberg.de/multiarm/}{multiarm}). Fast coputing is enabled by parallel programming.
#' 
#' @name optimal_multiarm
#' @param hr1 assumed true treatment effect on HR scale for treatment 1
#' @param hr2 assumed true treatment effect on HR scale for treatment 2
#' @param ec control arm event rate for phase II and III
#' @param n2min minimal total sample size in phase II, must be divisible by 3
#' @param n2max maximal total sample size in phase II, must be divisible by 3
#' @param stepn2 stepsize for the optimization over n2, must be divisible by 3
#' @param hrgomin minimal threshold value for the go/no-go decision rule
#' @param hrgomax maximal threshold value for the go/no-go decision rule
#' @param stephrgo stepsize for the optimization over HRgo
#' @param beta 1-beta (any-pair) power for calculation of the number of events for phase III
#' @param alpha one-sided significance level/ family wise error rate
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param steps1 lower boundary for effect size category "small" in HR scale, default: 1
#' @param stepm1 lower boundary for effect size category "medium" in HR scale = upper boundary for effect size category "small" in HR scale, default: 0.95
#' @param stepl1 lower boundary for effect size category "large" in HR scale = upper boundary for effect size category "medium" in HR scale, default: 0.85
#' @param b1 expected gain for effect size category "small"
#' @param b2 expected gain for effect size category "medium"
#' @param b3 expected gain for effect size category "large"
#' @param strategy choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)
#' @param num_cl number of clusters used for parallel computing, default: 1
#' @return
#' The output of the function \code{\link{optimal_multiarm}} is a data.frame containing the optimization results:
#' \describe{
#'   \item{Strategy}{Strategy, 1: "only best promising" or 2: "all promising"}
#'   \item{u}{maximal expected utility}
#'   \item{HRgo}{optimal threshold value for the decision rule to go to phase III}
#'   \item{n2}{optimal total sample size in phase II}
#'   \item{n3}{total expected sample size for phase III; rounded to next natural number}
#'   \item{n}{total sample size in the program; n = n2 + n3}
#'   \item{K}{maximal costs of the program}
#'   \item{pgo}{probability to go to phase III}
#'   \item{sProg}{probability of a successful program}
#'   \item{sProg1}{probability of a successful program with "small" treatment effect in Phase III}
#'   \item{sProg2}{probability of a successful program with "medium" treatment effect in Phase III}
#'   \item{sProg3}{probability of a successful program with "large" treatment effect in Phase III }
#'   \item{K2}{expected costs for phase II}
#'   \item{K3}{expected costs for phase III}
#'   }
#' and further input parameters.
#' res
#' Taking cat(comment()) of the data.frame object lists the used optimization sequences, start and finish date of the optimization procedure.
#' @examples
#' res <- optimal_multiarm(hr1 = 0.75, hr2 = 0.80, ec = 0.6,# define assumed true HRs and control arm event rate
#'   n2min = 30, n2max = 90, stepn2 = 6,                    # define optimization set for n2
#'   hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05,         # define optimization set for HRgo
#'   alpha = 0.05, beta = 0.1,                              # drug development planning parameters
#'   c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,               # define fixed and variable costs for phase II and III
#'   K = Inf, N = Inf, S = -Inf,                            # set maximal costs/ expected sample size for the program or minimal expected probability of a successful program
#'   steps1 = 1,                                            # define lower boundary for "small"
#'   stepm1 = 0.95,                                         # "medium"
#'   stepl1 = 0.85,                                         # and "large" treatment effect size categories as proposed by IQWiG (2016)
#'   b1 = 1000, b2 = 2000, b3 = 3000,                       # define expected benefit for a "small", "medium" and "large" treatment effect
#'   strategy = 1,                                          # choose strategy: 1, 2 or 3
#'   num_cl = 1)                                            # set number of cores used for parallelized computing (check maximum number possible with detectCores())
#' res
#' cat(comment(res))                                        # displays the optimization sequence, start and finish date of the optimization procedure.
#' @section drugdevelopR functions:
#' The drugdevelopR package provides the functions
#' \itemize{
#'   \item \code{\link{optimal_tte}},
#'   \item \code{\link{optimal_binary}} and
#'   \item \code{\link{optimal_normal}}
#' }
#' to plan optimal phase II/III drug development programs with
#' \itemize{
#'   \item time-to-event (treatment effect measured by hazard ratio (HR)),
#'   \item binary (treatment effect measured by risk ratio (RR)) or
#'   \item normally distributed (treatment effect measured by standardized difference in means (Delta))
#' }
#' endpoint, where the treatment effect is modelled by a \href{https://web.imbi.uni-heidelberg.de/prior/}{prior}. Optimal phase II/III drug development planning with fixed treatment effects can be done with the help of the R Shiny application \href{https://web.imbi.uni-heidelberg.de/basic/}{basic}. Extensions are 
#' \itemize{
#'   \item optimal planning of programs including methods for discounting of phase II results (function: \code{\link{optimal_bias}}, App: \href{https://web.imbi.uni-heidelberg.de/bias/}{bias}),
#'   \item optimal planning of programs with several phase III trials (function: \code{\link{optimal_multitrial}}, App: \href{https://web.imbi.uni-heidelberg.de/multitrial/}{multitrial}) and
#'   \item optimal planning of programs with multiple arms (function: \code{\link{optimal_multiarm}}, App: \href{https://web.imbi.uni-heidelberg.de/multiarm/}{multiarm}).
#' }
#' @references
#'Preussler, S., Kirchner, M., Goette, H., Kieser, M. (2019). Optimal Designs for Multi-Arm Phase II/III Drug Development Programs. Submitted to peer-review journal.
#'
#' IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at \href{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}, assessed last 15.05.19.
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
    
    cat("Optimization progress:", fill = TRUE)
    cat("", fill = TRUE)
    pb <- txtProgressBar(min = 0, max = length(HRGO), style = 3, label = "Optimization progess")
    
    for(j in 1:length(HRGO)){
      
      HRgo <- HRGO[j]
      
      cl <-  makeCluster(getOption("cl.cores", num_cl)) #define cluster
      
      clusterExport(cl, c("pmvnorm", "dmvnorm","qmvnorm","adaptIntegrate", "pgo", "ss", "Ess",
                          "PsProg", "alpha", "beta",
                          "steps1", "steps2", "stepm1", "stepm2", "stepl1", "stepl2",
                          "K", "N", "S", "strategy",
                          "c2", "c3", "c02", "c03",
                          "b1", "b2", "b3", "HRgo",
                          "hr1", "hr2", "ec"), envir = environment())
      
      
      res <- parSapply(cl, N2, utility_multiarm, HRgo,
                       alpha,beta,hr1,hr2,strategy,ec,
                       c2,c02,c3,c03,K,N,S,
                       steps1, stepm1, stepl1,b1, b2, b3)
      
      setTxtProgressBar(title= "i", pb, j)
      stopCluster(cl)
      
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
