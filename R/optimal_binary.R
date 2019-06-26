#' Optimal phase II/III drug development planning with binary endpoint
#'
#' The drugdevelopR package enables planning of phase II/III drug development programs with optimal sample size allocation and go/no-go decision rules. For binary endpoints the treatment effect is measured by the risk ratio (RR).The assumed true treatment effects can be assumed fixed or modelled by a prior distribution. The R Shiny application \href{https://web.imbi.uni-heidelberg.de/prior/}{prior} visualizes the prior distributions used in this package. Fast coputing is enabled by parallel programming.
#' @docType package
#' @name optimal_binary
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
#' @param steps1 lower boundary for effect size category "small" in RR scale, default: 1
#' @param stepm1 lower boundary for effect size category "medium" in RR scale = upper boundary for effect size category "small" in RR scale, default: 0.95
#' @param stepl1 lower boundary for effect size category "large" in RR scale = upper boundary for effect size category "medium" in RR scale, default: 0.85
#' @param b1 expected gain for effect size category "small"
#' @param b2 expected gain for effect size category "medium"
#' @param b3 expected gain for effect size category "large"
#' @param gamma to model different populations in phase II and III choose gamma!=0, default: 0
#' @param fixed choose if true treatment effects are fixed or random, if TRUE p11 is used as fixed effect for p1
#' @param num_cl number of clusters used for parallel computing, default: 1
#' @return
#' The output of the function \code{\link{optimal_binary}} is a data.frame containing the optimization results:
#' \describe{
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
#' res <- optimal_binary(w = 0.3,                           # define parameters for prior
#'   p0 = 0.6, p11 =  0.3, p12 = 0.5, in1 = 30, in2 = 60,   # (https://web.imbi.uni-heidelberg.de/prior/)
#'   n2min = 20, n2max = 100, stepn2 = 4,                   # define optimization set for n2
#'   rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.05,         # define optimization set for RRgo
#'   alpha = 0.05, beta = 0.1,                              # drug development planning parameters
#'   c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,               # define fixed and variable costs for phase II and III,
#'   K = Inf, N = Inf, S = -Inf,                            # set maximal costs/ expected sample size for the program or minimal expected probability of a successful program
#'   steps1 = 1,                                            # define lower boundary for "small"
#'   stepm1 = 0.95,                                         # "medium"
#'   stepl1 = 0.85,                                         # and "large" treatment effect size categories as proposed by IQWiG (2016)
#'   b1 = 1000, b2 = 2000, b3 = 3000,                       # define expected benefit for a "small", "medium" and "large" treatment effect
#'   gamma = 0,                                             # assume different/same population structures in phase II and III
#'   fixed = FALSE,                                         # choose if true treatment effects are fixed or random
#'   skipII = FALSE,                                        # choose if skipping phase II would be an option
#'   num_cl = 1)                                            # set number of cores used for parallelized computing (check maximum number possible with detectCores())
#'
#' cat(comment(res))                                        # displays the optimization sequence, start and finish date of the optimization procedure.
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
#' endpoint, where the treatment effect is modelled by a \href{https://web.imbi.uni-heidelberg.de/prior/}{prior}. Optimal phase II/III drug development planning with fixed treatment effects can be done with the help of the R Shiny application \href{https://web.imbi.uni-heidelberg.de/basic/}{basic}.
#' @references
#' IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at \href{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}, assessed last 15.05.19.
#' @export
optimal_binary <- function(w, p0, p11, p12, in1, in2,
                        n2min, n2max, stepn2,
                        rrgomin, rrgomax, steprrgo,
                        alpha, beta, 
                        c2, c3, c02, c03, 
                        K = Inf, N = Inf, S = -Inf,
                        steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
                        b1, b2, b3,
                        gamma = 0, fixed = FALSE,
                        skipII = FALSE, num_cl = 1){

  date <- Sys.time()

  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 0

   if(skipII==TRUE){

     if(fixed){
       median_prior = p11
     }else{
       median_prior = round(quantile(box_binary(w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2),0.5),2)
       
       names(median_prior) = NULL 
     }
     

     res <- utility_skipII_binary(alpha = alpha, beta = beta,
                                  c03 = c03, c3 = c3,
                                  b1 = b1, b2 = b2, b3 = b3,
                                  p0 = p0, median_prior = median_prior,
                                  K = K, N = N, S = S,
                                  steps1 = steps1, steps2 = steps2,
                                  stepm1 = stepm1, stepm2 = stepm2,
                                  stepl1 = stepl1, stepl2 = stepl2,
                                  w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2,
                                  gamma = gamma, fixed = fixed)

     if(fixed){
       result_skipII <-  data.frame(u = round(res[1],2), RR=round(median_prior/p0,2),
                                    RRgo = Inf, n2 = 0, n3 = res[2],
                                    pgo = 1, sProg = round(res[3],2), K = K, K2 = 0, K3 = round(res[4]),
                                    sProg1 = round(res[5],2), sProg2 = round(res[6],2), sProg3 = round(res[7],2),
                                    steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                                    alpha = alpha, beta = beta, c02 = 0,
                                    c03 = c03, c2 = 0, c3 = c3, b1 = b1, b2 = b2, b3 = b3,
                                    p0 = p0, p1 = p11, gamma = gamma)  
     }else{
       result_skipII <-  data.frame(u = round(res[1],2), median_prior_RR=round(median_prior/p0,2),
                                    RRgo = Inf, n2 = 0, n3 = res[2],
                                    pgo = 1, sProg = round(res[3],2), K = K, K2 = 0, K3 = round(res[4]),
                                    sProg1 = round(res[5],2), sProg2 = round(res[6],2), sProg3 = round(res[7],2),
                                    steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                                    alpha = alpha, beta = beta, c02 = 0,
                                    c03 = c03, c2 = 0, c3 = c3, b1 = b1, b2 = b2, b3 = b3,
                                    w = w, p0 = p0, p11 = p11, p12 = p12, in1 = in1, in2 = in2, gamma = gamma)
     }
     


     cat("Result when skipping phase II:", fill = TRUE)
     cat("", fill = TRUE)
     print(result_skipII)
     cat("", fill = TRUE)
     cat("", fill = TRUE)

   }
  
   if(round(n2min/2) != n2min / 2) {
     n2min = n2min - 1
     cat(paste0("n2min must be equal number and is therefore set to: ", n2min), fill = TRUE)
     cat("", fill = TRUE)
   }
   if(round(n2max/2) != n2max / 2) {
     n2max = n2max + 1
     cat(paste0("n2max must be equal number and is therefore set to: ", n2max), fill = TRUE)
     cat("", fill = TRUE)
   }
   if(round(stepn2/2) != stepn2 / 2) {
     stepn2 = stepn2 + 1
     cat(paste0("stepn2 must be equal number and is therefore set to: ", stepn2), fill = TRUE)
     cat("", fill = TRUE)
   }

   HRGO <- seq(rrgomin, rrgomax, steprrgo)
   N2   <- seq(n2min, n2max, stepn2)

   ufkt <- spfkt <- pgofkt <- K2fkt <- K3fkt <-
     sp1fkt <- sp2fkt <- sp3fkt <- n2fkt <- n3fkt <- matrix(0, length(N2), length(HRGO))

   cat("Optimization progress:", fill = TRUE)
   cat("", fill = TRUE)
   pb <- txtProgressBar(min = 0, max = length(HRGO), style = 3, label = "Optimization progess")

   for(j in 1:length(HRGO)){

      RRgo <- HRGO[j]

      cl <-  makeCluster(getOption("cl.cores", num_cl)) #define cluster
      
      parallel::clusterExport(cl, c("pmvnorm", "dmvnorm", "prior_binary", "Epgo_binary", "En3_binary",
                          "EPsProg_binary","t1", "t2", "t3", "alpha", "beta",
                          "steps1", "steps2", "stepm1", "stepm2", "stepl1", "stepl2",
                          "K", "N", "S", "gamma", "fixed",
                          "c2", "c3", "c02", "c03",
                          "b1", "b2", "b3", "w", "RRgo",
                          "p0", "p11", "p12", "in1", "in2"), envir=environment())
      
      result <- parSapply(cl, N2, utility_binary, RRgo, w, p0, p11, p12, in1, in2,
                          alpha, beta, 
                          c2, c3, c02, c03, K, N, S,
                          steps1, stepm1, stepl1,
                          b1, b2, b3,
                          gamma, fixed)

      setTxtProgressBar(title= "i", pb, j)
      stopCluster(cl)

      ufkt[, j]      <-  result[1, ]
      n3fkt[, j]     <-  result[2, ]
      spfkt[, j]     <-  result[3, ]
      pgofkt[, j]    <-  result[4, ]
      K2fkt[, j]     <-  result[5, ]
      K3fkt[, j]     <-  result[6, ]
      sp1fkt[, j]    <-  result[7, ]
      sp2fkt[, j]    <-  result[8, ]
      sp3fkt[, j]    <-  result[9, ]

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
     result <-  data.frame(u = round(Eud,2), RRgo = HRGO[J], n2 = N2[I],
                           n3 = n3, n = N2[I] + n3,
                           pgo = round(pg,2), sProg = round(prob,2),
                           p0 = p0, p1 = p11, 
                           K = K, K2 = round(k2), K3 = round(k3),
                           sProg1 = round(prob1,2), sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                           steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                           alpha = alpha, beta = beta, c02 = c02,
                           c03 = c03, c2 = c2, c3 = c3, b1 = b1, b2 = b2, b3 = b3, gamma = gamma)  
   }else{
     result <-  data.frame(u = round(Eud,2), RRgo = HRGO[J], n2 = N2[I],
                           n3 = n3, n = N2[I] + n3,
                           pgo = round(pg,2), sProg = round(prob,2),
                           w = w, p0 = p0, p11 = p11, p12 = p12, in1 = in1, in2 = in2,
                           K = K, K2 = round(k2), K3 = round(k3),
                           sProg1 = round(prob1,2), sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                           steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                           alpha = alpha, beta = beta, c02 = c02,
                           c03 = c03, c2 = c2, c3 = c3, b1 = b1, b2 = b2, b3 = b3, gamma = gamma)
   }
   


   comment(result) <-   c("\noptimization sequence RRgo:", HRGO,
                      "\noptimization sequence n2:", N2,
                      "\nset on date:", as.character(date),
                      "\nfinish date:", as.character(Sys.time()))
   close(pb)

   cat("", fill = TRUE)
   cat("", fill = TRUE)
   cat("Optimization result:", fill = TRUE)
   cat("", fill = TRUE)

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

