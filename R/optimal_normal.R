#' Optimal phase II/III drug development planning with normally distributed endpoint
#'
#' The function \code{\link{optimal_normal}} of the \code{\link{drugdevelopR}}
#' package enables planning of phase II/III drug development programs with
#' optimal sample size allocation and go/no-go decision rules for normally
#' distributed endpoints. The treatment effect is measured by the standardized
#' difference in means. The assumed true treatment effects can be assumed to be 
#' fixed or modelled by a prior distribution. The R Shiny application
#' \href{https://web.imbi.uni-heidelberg.de/prior/}{prior} visualizes the prior
#' distributions used in this package. Fast computing is enabled by parallel
#'programming.
#' 
#' @name optimal_normal
#' @inheritParams optimal_normal_generic
#' @param skipII choose if skipping phase II is an option, default: FALSE; 
#' if TRUE, the program calculates the expected utility for the case when phase
#' II is skipped and compares it to the situation when phase II is not skipped.
#' The results are then returned as a list of two results lists, `res[[1]]`
#' being the results when including phase II and `res[[2]]` when skipping phase II.
#' 
#' @importFrom ggdmc dtnorm
#' @importFrom hbmem rtnorm
#' 
#' @return
#' The output of the function \code{\link{optimal_normal}} is a data.frame containing the optimization results:
#' \describe{
#'   \item{u}{maximal expected utility under the optimization constraints, i.e. the expected utility of the optimal sample size and threshold value}
#'   \item{Kappa}{optimal threshold value for the decision rule to go to phase III}
#'   \item{n2}{total sample size for phase II}
#'   \item{n3}{total sample size for phase III; rounded to the next even natural number}
#'   \item{n}{total sample size in the program; n = n2 + n3}
#'   \item{K}{maximal costs of the program}
#'   \item{pgo}{probability to go to phase III}
#'   \item{sProg}{probability of a successful program}
#'   \item{sProg1}{probability of a successful program with "small" treatment effect in phase III}
#'   \item{sProg2}{probability of a successful program with "medium" treatment effect in phase III}
#'   \item{sProg3}{probability of a successful program with "large" treatment effect in phase III }
#'   \item{K2}{expected costs for phase II}
#'   \item{K3}{expected costs for phase III}
#' }
#' and further input parameters.
#'
#' Taking `cat(comment())` of the data.frame object lists the used optimization sequences, start and finish date of the optimization procedure.
#' 
#' @examples
#' res <- optimal_normal(w=0.3,                             # define parameters for prior
#'   Delta1 = 0.375, Delta2 = 0.625, in1=300, in2=600,      # (https://web.imbi.uni-heidelberg.de/prior/)
#'   a = 0.25, b = 0.75,
#'   n2min = 20, n2max = 100, stepn2 = 4,                   # define optimization set for n2
#'   kappamin = 0.02, kappamax = 0.2, stepkappa = 0.02,     # define optimization set for kappa
#'   alpha = 0.05, beta = 0.1,                              # drug development planning parameters
#'   c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,             # define fixed and variable costs for phase II and III
#'   K = Inf, N = Inf, S = -Inf,                            # set maximal costs/ expected sample size for the program or minimal expected probability of a successful program
#'   steps1 = 0,                                            # define lower boundary for "small"
#'   stepm1 = 0.5,                                          # "medium"
#'   stepl1 = 0.8,                                          # and "large" treatment effect size categories as proposed by e.g. Cohen (1988)
#'   b1 = 3000, b2 = 8000, b3 = 10000,                      # define expected benefit for a "small", "medium" and "large" treatment effect
#'   gamma = 0,                                             # assume different/same population structures in phase II and III
#'   fixed = FALSE,                                         # choose if true treatment effects are fixed or random
#'   skipII = FALSE,                                        # choose if skipping phase II would be an option
#'   num_cl = 1)                                            # set number of cores used for parallelized computing (check maximum number possible with detectCores())
#' res
#' cat(comment(res))                                        # displays the optimization sequence, start and finish date of the optimization procedure
#' @references
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences.
#'
#' @export
optimal_normal <- function(w, Delta1, Delta2, in1, in2, a, b,
                        n2min, n2max, stepn2,
                        kappamin, kappamax, stepkappa,
                        alpha, beta, 
                        c2, c3, c02, c03, 
                        K = Inf, N = Inf, S = -Inf,
                        steps1 = 0, stepm1 = 0.5, stepl1 = 0.8,
                        b1, b2, b3,
                        gamma = 0,  fixed = FALSE,
                        skipII = FALSE,  num_cl = 1){

   date <- Sys.time()

   if(skipII){

     if(fixed){
       median_prior = Delta1
     }else{
       median_prior = round(quantile(box_normal(w = w, Delta1 = Delta1, Delta2 = Delta2,
                                                in1 = in1, in2 = in2, a = a, b = b),0.5),2)
       
       names(median_prior) = NULL  
     }
    
     res <- utility_skipII_normal(alpha = alpha, beta = beta,
                            c03 = c03, c3 = c3,
                            b1 = b1, b2 = b2, b3 = b3,
                            median_prior = median_prior,
                            K = K, N = N, S = S,
                            steps1 = steps1,
                            stepm1 = stepm1,
                            stepl1 = stepl1, 
                            w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                            gamma = gamma, fixed = fixed)
      # c(EU, n3, SP, K3, prob1, prob2, prob3)
     if(fixed){
       
       result_skipII <-  data.frame(u = round(res[1],2), Delta=round(median_prior,2),
                                    Kappa = -Inf, n2 = 0, n3 = res[2],
                                    pgo = 1, sProg = round(res[3],2), K = K, K2 = 0, K3 = round(res[4]),
                                    sProg1 = round(res[5],2), sProg2 = round(res[6],2), sProg3 = round(res[7],2),
                                    steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                                    alpha = alpha, beta = beta, c02 = 0,
                                    c03 = c03, c2 = 0, c3 = c3, b1 = b1, b2 = b2, b3 = b3,
                                    gamma = gamma)
       
     }else{
       result_skipII <-  data.frame(u = round(res[1],2), median_prior_Delta=round(median_prior,2),
                                    Kappa = -Inf, n2 = 0, n3 = res[2],
                                    pgo = 1, sProg = round(res[3],2), K = K, K2 = 0, K3 = round(res[4]),
                                    sProg1 = round(res[5],2), sProg2 = round(res[6],2), sProg3 = round(res[7],2),
                                    steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                                    alpha = alpha, beta = beta, c02 = 0,
                                    c03 = c03, c2 = 0, c3 = c3, b1 = b1, b2 = b2, b3 = b3,
                                    w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b=b,
                                    gamma = gamma)
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

   KAPPA <- seq(kappamin, kappamax, stepkappa)
   N2   <- seq(n2min, n2max, stepn2)

   ufkt <- spfkt <- pgofkt <- K2fkt <- K3fkt <-
     sp1fkt <- sp2fkt <- sp3fkt <- n2fkt <- n3fkt <- matrix(0, length(N2), length(KAPPA))

   cat("Optimization progress:", fill = TRUE)
   cat("", fill = TRUE)
   pb <- txtProgressBar(min = 0, max = length(KAPPA), style = 3, label = "Optimization progess")

   for(j in 1:length(KAPPA)){

      kappa <- KAPPA[j]

      cl <-  parallel::makeCluster(getOption("cl.cores", num_cl)) #define cluster

      parallel::clusterExport(cl, c("pmvnorm", "dmvnorm","dtnorm", "prior_normal", "Epgo_normal", "En3_normal",
                          "EPsProg_normal", "alpha", "beta",
                          "steps1", "stepm1", "stepl1", 
                          "K", "N", "S", "gamma", "fixed",
                          "c2", "c3", "c02", "c03",
                          "b1", "b2", "b3", "w", "kappa",
                          "Delta1", "Delta2", "in1", "in2", "a", "b"), envir=environment())

      result <- parallel::parSapply(cl, N2, utility_normal, kappa, w, Delta1, Delta2, in1, in2, a, b,
                          alpha, beta, 
                          c2, c3, c02, c03, 
                          K, N, S,
                          steps1, stepm1, stepl1,
                          b1, b2, b3,
                          gamma, fixed)
      
      setTxtProgressBar(title= "i", pb, j)
      parallel::stopCluster(cl)

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
     
     result <-  data.frame(u = round(Eud,2), Kappa = KAPPA[J], n2 = N2[I],
                           n3 = n3, n = N2[I] + n3,
                           pgo = round(pg,2), sProg = round(prob,2),
                           Delta = Delta1,
                           K = K, K2 = round(k2), K3 = round(k3),
                           sProg1 = round(prob1,2), sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                           steps1 = steps1, stepm1 = stepm1, stepl1 = stepl1,
                           alpha = alpha, beta = beta, c02 = c02,
                           c03 = c03, c2 = c2, c3 = c3, b1 = b1, b2 = b2, b3 = b3, gamma = gamma)
   }else{
     result <-  data.frame(u = round(Eud,2), Kappa = KAPPA[J], n2 = N2[I],
                           n3 = n3, n = N2[I] + n3,
                           pgo = round(pg,2), sProg = round(prob,2),
                           w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                           K = K, K2 = round(k2), K3 = round(k3),
                           sProg1 = round(prob1,2), sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                           steps1 = steps1, stepm1 = stepm1, stepl1 = stepl1,
                           alpha = alpha, beta = beta, c02 = c02,
                           c03 = c03, c2 = c2, c3 = c3, b1 = b1, b2 = b2, b3 = b3, gamma = gamma)
   }

   comment(result) <-   c("\noptimization sequence Kappa:", KAPPA,
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
   
   if(skipII==TRUE){
     cat("", fill = TRUE)
     if(result_skipII[1]>result[1]){
       cat("Skipping phase II is the optimal option with respect to the maximal expected utility.", fill = TRUE)
     }else{
       cat("Skipping phase II is NOT the optimal option with respect to the maximal expected utility.", fill = TRUE)
     }
     print(list(result,result_skipII))
     return(list(result,result_skipII))
     
   }

   return(result)
}
