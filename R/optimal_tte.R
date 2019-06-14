#' Optimal phase II/III drug development planning with time-to-event endpoint
#'
#' The drugdevelopR package enables planning of phase II/III drug development programs with optimal sample size allocation and go/no-go decision rules when assuming a prior distribution for the treatment effect. For time-to-event endpoints the treatment effect is measured by the hazard ratio (HR) and more information on the framework can be found in Kirchner et al. (2016). The R Shiny application \href{https://web.imbi.uni-heidelberg.de/prior/}{prior} visualizes the associated prior distributions used in this package.
#' @docType package
#' @name optimal_tte
#' @param w weight for mixture prior distribution
#' @param hr1 first assumed true treatment effect on HR scale for prior distribution
#' @param hr2 second assumed true treatment effect on HR scale for prior distribution
#' @param id1 amount of information for hr1 in terms of number of events
#' @param id2 amount of information for hr2 in terms of number of events
#' @param d2min minimal number of events for phase II
#' @param d2max maximal number of events for phase II
#' @param stepd2 stepsize for the optimization over d2
#' @param hrgomin minimal threshold value for the go/no-go decision rule
#' @param hrgomax maximal threshold value for the go/no-go decision rule
#' @param stephrgo stepsize for the optimization over HRgo
#' @param beta 1-beta power for calculation of the number of events for phase III by Schoenfeld (1981) formula
#' @param alpha significance level
#' @param xi2 event rate for phase II
#' @param xi3 event rate for phase III
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
#' @param gamma to model different populations in phase II and III choose gamma!=0, default: 0
#' @param skipII choose if skipping phase II is an option, default: FASLE
#' @param num_cl number of clusters used for parallel computing, default: 1
#' @return
#' The output of the function \code{\link{optimal_tte}} is a data.frame containing the optimization results:
#' \describe{
#'   \item{u}{maximal expected utility}
#'   \item{HRgo}{optimal threshold value for the decision rule to go to phase III}
#'   \item{d2}{optimal total number of events for phase II}
#'   \item{d3}{total expected number of events for phase III; rounded to next natural number}
#'   \item{d}{total expected number of events in the program; d = d2 + d3}
#'   \item{n2}{total sample size for phase II; rounded to the next even natural number}
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
#'   }
#' and further input parameters.
#'
#' Taking cat(comment()) of the data.frame object lists the used optimization sequences, start and finish date of the optimization procedure.
#' @examples
#' res <- optimal_tte(w = 0.3,                              # define parameters for prior
#'   hr1 = 0.69, hr2 = 0.88, id1 = 210, id2 = 420,          # (https://web.imbi.uni-heidelberg.de/prior/)
#'   d2min = 20, d2max = 100, stepd2 = 5,                   # define optimization set for n2
#'   hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05,         # define optimization set for HRgo
#'   alpha = 0.05, beta = 0.1, xi2 = 0.7, xi3 = 0.7,        # drug development planning parameters
#'   c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,               # define fixed and variable costs for phase II and III
#'   K = Inf, N = Inf, S = -Inf,                            # set maximal costs/ expected sample size for the program or minimal expected probability of a successful program
#'   steps1 = 1,                                            # define lower boundary for "small"
#'   stepm1 = 0.95,                                         # "medium"
#'   stepl1 = 0.85,                                         # and "large" treatment effect size categories as proposed by IQWiG (2016)
#'   b1 = 1000, b2 = 2000, b3 = 3000,                       # define expected benefit for a "small", "medium" and "large" treatment effect
#'   gamma = 0,                                             # assume different/same population structures in phase II and III
#'   skipII = FALSE,                                        # choose if skipping phase II would be an option
#'   num_cl = 1)                                            # set number of cores used for parallelized computing (check maximum number possible with detectCores())
#'
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
#' endpoint, where the treatment effect is modelled by a \href{https://web.imbi.uni-heidelberg.de/prior/}{prior}. Optimal phase II/III drug development planning with fixed treatment effects can be done with the help of the R Shiny application \href{https://web.imbi.uni-heidelberg.de/basic/}{basic}.
#' @references
#' Kirchner, M., Kieser, M., Goette, H., & Schueler, A. (2016). Utility-based optimization of phase II/III programs. Statistics in Medicine, 35(2), 305-316.
#'
#' IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at \href{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}, assessed last 15.05.19.
#'
#' Schoenfeld, D. (1981). The asymptotic properties of nonparametric tests for comparing survival distributions. Biometrika, 68(1), 316-319.
#'
#' @export
optimal_tte <- function(w,  hr1, hr2, id1, id2,
                        d2min, d2max, stepd2,
                        hrgomin, hrgomax, stephrgo,
                        alpha, beta, xi2, xi3,
                        c2, c3, c02, c03, 
                        K = Inf, N = Inf, S = -Inf,
                        steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                        b1, b2, b3,
                        gamma = 0,  fixed = FALSE,
                        skipII = FALSE,  num_cl = 1){

  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 0

  date <- Sys.time()

  if(skipII==TRUE){

    if(fixed==TRUE){
      median_prior = -log(hr1)
    }else{
      median_prior = round(quantile(box_tte(w, hr1, hr2, id1,id2),0.5),2)
      
      names(median_prior) = NULL  
    }

    res <- utility_skipII_tte(alpha = alpha, beta = beta, xi3 = xi3,
                              c03 = c03, c3 = c3,
                              b1 = b1, b2 = b2, b3 = b3,
                              median_prior = median_prior,
                              K = K, N = N, S = S,
                              steps1 = steps1, steps2 = steps2,
                              stepm1 = stepm1, stepm2 = stepm2,
                              stepl1 = stepl1, stepl2 = stepl2,
                              w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                              gamma = gamma)

    result_skipII <-  data.frame(u = round(res[1],2), median_prior_HR=round(exp(-median_prior),2),
                             HRgo = Inf, d2 = 0, d3 = res[2],
                             n2 = 0, n3 = res[3],
                             pgo = 1, sProg = round(res[4],2), K = K, N = N, S = S, K2 = 0, K3 = round(res[5]),
                             sProg1 = round(res[6],2), sProg2 = round(res[7],2), sProg3 = round(res[8],2),
                             steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                             alpha = alpha, beta = beta, xi3 = xi3, c02 = 0,
                             c03 = c03, c2 = 0, c3 = c3, b1 = b1, b2 = b2, b3 = b3,
                             w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2, gamma = gamma)

    cat("Result when skipping phase II:", fill = TRUE)
    cat("", fill = TRUE)
    print(result_skipII)
    cat("", fill = TRUE)
    cat("", fill = TRUE)
  }

  HRGO <- seq(hrgomin, hrgomax, stephrgo)
  D2   <- seq(d2min, d2max, stepd2)

  ufkt <- d3fkt <- spfkt <- pgofkt <- K2fkt <- K3fkt <-
   sp1fkt <- sp2fkt <- sp3fkt <- n2fkt <- n3fkt <- matrix(0, length(D2), length(HRGO))

  cat("Optimization progress:", fill = TRUE)
  cat("", fill = TRUE)
  pb <- txtProgressBar(min = 0, max = length(HRGO), style = 3, label = "Optimization progess")

  for(j in 1:length(HRGO)){

    HRgo <- HRGO[j]

    cl <-  makeCluster(getOption("cl.cores", num_cl)) #define cluster
    
    clusterExport(cl, c("pmvnorm", "dmvnorm", "prior_tte", "Epgo_tte", "Ed3_tte",
                        "EPsProg_tte", "alpha", "beta",
                        "steps1", "steps2", "stepm1", "stepm2", "stepl1", "stepl2",
                        "K", "N", "S", "gamma",
                        "xi2", "xi3", "c2", "c3", "c02", "c03",
                        "b1", "b2", "b3", "w", "HRgo",
                        "hr1", "hr2", "id1", "id2"), envir = environment())

    
    result <- parSapply(cl, D2, utility_tte, HRgo, w, hr1, hr2, id1, id2,
                        alpha, beta, xi2, xi3,
                        c2, c3, c02, c03, 
                        K, N, S,
                        steps1, stepm1, stepl1,
                        b1, b2, b3,
                        gamma, fixed)

    setTxtProgressBar(title= "i", pb, j)
    stopCluster(cl)

    ufkt[, j]      <-  result[1, ]
    d3fkt[, j]     <-  result[2, ]
    spfkt[, j]     <-  result[3, ]
    pgofkt[, j]    <-  result[4, ]
    K2fkt[, j]     <-  result[5, ]
    K3fkt[, j]     <-  result[6, ]
    sp1fkt[, j]    <-  result[7, ]
    sp2fkt[, j]    <-  result[8, ]
    sp3fkt[, j]    <-  result[9, ]
    n2fkt[, j]     <-  result[10, ]
    n3fkt[, j]     <-  result[11, ]

   }

  ind   <-  which(ufkt  ==  max(ufkt), arr.ind <-  TRUE)

  I <-  as.vector(ind[1, 1])
  J <-  as.vector(ind[1, 2])

  Eud   <- ufkt[I, J]
  d3    <- d3fkt[I, J]
  prob  <- spfkt[I, J]
  pg    <- pgofkt[I, J]
  k2    <- K2fkt[I, J]
  k3    <- K3fkt[I, J]
  prob1 <- sp1fkt[I, J]
  prob2 <- sp2fkt[I, J]
  prob3 <- sp3fkt[I, J]
  n2    <- n2fkt[I,J]
  n3    <- n3fkt[I,J]

  
  if(fixed){
    
    result <-  data.frame(u = round(Eud,2), HRgo = HRGO[J], d2 = D2[I], d3 = d3, d = D2[I] + d3,
                          n2 = n2, n3 = n3, n = n2 + n3,
                          pgo = round(pg,2), sProg = round(prob,2),
                          w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                          K = K, N = N, S = S, K2 = round(k2), K3 = round(k3),
                          sProg1 = round(prob1,2), sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                          steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                          alpha = alpha, beta = beta, xi2 = xi2, xi3 = xi3, c02 = c02,
                          c03 = c03, c2 = c2, c3 = c3, b1 = b1, b2 = b2, b3 = b3, gamma = gamma)
  }else{
    result <-  data.frame(u = round(Eud,2), HRgo = HRGO[J], d2 = D2[I], d3 = d3, d = D2[I] + d3,
                          n2 = n2, n3 = n3, n = n2 + n3,
                          pgo = round(pg,2), sProg = round(prob,2),
                          hr = hr1,
                          K = K, N = N, S = S, K2 = round(k2), K3 = round(k3),
                          sProg1 = round(prob1,2), sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                          steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                          alpha = alpha, beta = beta, xi2 = xi2, xi3 = xi3, c02 = c02,
                          c03 = c03, c2 = c2, c3 = c3, b1 = b1, b2 = b2, b3 = b3, gamma = gamma)
  }


  comment(result) <-   c("\noptimization sequence HRgo:", HRGO,
                    "\noptimization sequence d2:", D2,
                    "\nset on date:", as.character(date),
                    "\nfinish date:", as.character(Sys.time()))
  close(pb)

  cat("", fill = TRUE)
  cat("", fill = TRUE)
  cat("Optimization result:", fill = TRUE)
  cat("", fill = TRUE)

  print(result)

  if(skipII==TRUE){
   cat("", fill = TRUE)
   if(result_skipII[1]>result[1]){
     cat("Skipping phase II is the optimal option with respect to the maximal expected utility.", fill = TRUE)
   }else{
     cat("Skipping phase II is NOT the optimal option with respect to the maximal expected utility.", fill = TRUE)
   }
    
    return(list(result,result_skipII))
    
  }

  return(result)
  
}
