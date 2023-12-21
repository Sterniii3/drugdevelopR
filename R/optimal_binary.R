#' Optimal phase II/III drug development planning with binary endpoint
#'
#' The \code{\link{optimal_binary}} function of the drugdevelopR package enables
#' planning of phase II/III drug development programs with optimal sample size 
#' allocation and go/no-go decision rules for binary endpoints. In this case,
#' the treatment effect is measured by the risk ratio (RR). The assumed true
#' treatment effects can be assumed to be fixed or modelled by a prior
#' distribution. The R Shiny application
#' \href{https://web.imbi.uni-heidelberg.de/prior/}{prior} visualizes the prior
#' distributions used in this package. Fast computing is enabled by parallel 
#' programming.
#' 
#' @name optimal_binary
#' 
#' @inheritParams optimal_binary_generic
#' @param skipII skipII choose if skipping phase II is an option, default: FALSE; 
#' if TRUE, the program calculates the expected utility for the case when phase
#' II is skipped and compares it to the situation when phase II is not skipped.
#' The results are then returned as a two-row data frame, `res[1, ]`
#' being the results when including phase II and `res[2, ]` when skipping phase II.
#' `res[2, ]` has an additional parameter, `res[2, ]$median_prior_RR`, which is
#' the assumed effect size used for planning the phase III study when the 
#' phase II is skipped.
#' 
#' @return
#' `r optimal_return_doc(type = "binary")` 
#' 
#' @importFrom progressr progressor
#'
#' @examples
#' # Activate progress bar (optional)
#' \dontrun{
#' progressr::handlers(global = TRUE)
#' }
#' # Optimize
#' \donttest{
#' optimal_binary(w = 0.3,                             # define parameters for prior
#'   p0 = 0.6, p11 =  0.3, p12 = 0.5,
#'   in1 = 30, in2 = 60,                               # (https://web.imbi.uni-heidelberg.de/prior/)
#'   n2min = 20, n2max = 100, stepn2 = 4,              # define optimization set for n2
#'   rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.05,    # define optimization set for RRgo
#'   alpha = 0.025, beta = 0.1,                        # drug development planning parameters
#'   c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,          # fixed and variable costs for phase II/III,
#'   K = Inf, N = Inf, S = -Inf,                       # set constraints
#'   steps1 = 1,                                       # define lower boundary for "small"
#'   stepm1 = 0.95,                                    # "medium"
#'   stepl1 = 0.85,                                    # and "large" treatment effect size categories
#'   b1 = 1000, b2 = 2000, b3 = 3000,                  # define expected benefits
#'   gamma = 0,                                        # population structures in phase II/III
#'   fixed = FALSE,                                    # true treatment effects are fixed/random
#'   skipII = FALSE,                                   # choose if skipping phase II is an option
#'   num_cl = 2)                                       # number of cores for parallelized computing
#'   }
#' @references
#' IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at \href{https://www.iqwig.de/ueber-uns/methoden/methodenpapier/}{https://www.iqwig.de/ueber-uns/methoden/methodenpapier/}, assessed last 15.05.19.
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
                                  steps1 = steps1, 
                                  stepm1 = stepm1, 
                                  stepl1 = stepl1, 
                                  w = w, p11 = p11, p12 = p12, in1 = in1, in2 = in2,
                                  gamma = gamma, fixed = fixed)

     if(fixed){
       result_skipII <-  data.frame(skipII = TRUE,
                                    u = round(res[1],2),
                                    RR=round(median_prior/p0,2),
                                    RRgo = Inf, n2 = 0, n3 = res[2],
                                    pgo = 1, sProg = round(res[3],2), K = K, K2 = 0, K3 = round(res[4]),
                                    sProg1 = round(res[5],2), sProg2 = round(res[6],2), sProg3 = round(res[7],2),
                                    steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                                    alpha = alpha, beta = beta, c02 = 0,
                                    c03 = c03, c2 = 0, c3 = c3, b1 = b1, b2 = b2, b3 = b3,
                                    p0 = p0, p1 = p11, gamma = gamma)  
     }else{
       result_skipII <-  data.frame(skipII = TRUE,
                                    u = round(res[1],2), median_prior_RR=round(median_prior/p0,2),
                                    RRgo = Inf, n2 = 0, n3 = res[2],
                                    pgo = 1, sProg = round(res[3],2), K = K, K2 = 0, K3 = round(res[4]),
                                    sProg1 = round(res[5],2), sProg2 = round(res[6],2), sProg3 = round(res[7],2),
                                    steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                                    alpha = alpha, beta = beta, c02 = 0,
                                    c03 = c03, c2 = 0, c3 = c3, b1 = b1, b2 = b2, b3 = b3,
                                    w = w, p0 = p0, p11 = p11, p12 = p12, in1 = in1, in2 = in2, gamma = gamma)
     }

   }
  
   if(round(n2min/2) != n2min / 2) {
     n2min = n2min - 1
     message(paste0("n2min must be an even number and is therefore set to: ", n2min))
   }
   if(round(n2max/2) != n2max / 2) {
     n2max = n2max + 1
     message(paste0("n2max must be an even number and is therefore set to: ", n2max))
   }
   if(round(stepn2/2) != stepn2 / 2) {
     stepn2 = stepn2 + 1
     message(paste0("stepn2 must be an even number and is therefore set to: ", stepn2))
   }

   HRGO <- seq(rrgomin, rrgomax, steprrgo)
   N2   <- seq(n2min, n2max, stepn2)

   ufkt <- spfkt <- pgofkt <- K2fkt <- K3fkt <-
     sp1fkt <- sp2fkt <- sp3fkt <- n2fkt <- n3fkt <- matrix(0, length(N2), length(HRGO))

   pb <- progressr::progressor(along = HRGO, label = "Optimization progress", message = "Optimization progress")
   pb("Performing optimization", class = "sticky", amount = 0)
   
   RRgo <- NA_real_
   cl <-  parallel::makeCluster(getOption("cl.cores", num_cl)) #define cluster
   
   parallel::clusterExport(cl, c("pmvnorm", "dmvnorm", "prior_binary", "Epgo_binary", "En3_binary",
                                 "EPsProg_binary","t1", "t2", "t3", "alpha", "beta",
                                 "steps1", "stepm1",  "stepl1", 
                                 "K", "N", "S", "gamma", "fixed",
                                 "c2", "c3", "c02", "c03",
                                 "b1", "b2", "b3", "w", "RRgo",
                                 "p0", "p11", "p12", "in1", "in2"), envir=environment())

   for(j in 1:length(HRGO)){

      RRgo <- HRGO[j]

      
      
      result <- parallel::parSapply(cl, N2, utility_binary, RRgo, w, p0, p11, p12, in1, in2,
                          alpha, beta, 
                          c2, c3, c02, c03, K, N, S,
                          steps1, stepm1, stepl1,
                          b1, b2, b3,
                          gamma, fixed)

      pb()
      

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
     result <-  data.frame(skipII = FALSE,
                           u = round(Eud,2), RRgo = HRGO[J], n2 = N2[I],
                           n3 = n3, n = N2[I] + n3,
                           pgo = round(pg,2), sProg = round(prob,2),
                           p0 = p0, p1 = p11, 
                           K = K, K2 = round(k2), K3 = round(k3),
                           sProg1 = round(prob1,2), sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                           steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                           alpha = alpha, beta = beta, c02 = c02,
                           c03 = c03, c2 = c2, c3 = c3, b1 = b1, b2 = b2, b3 = b3, gamma = gamma)  
   }else{
     result <-  data.frame(skipII = FALSE,
                           u = round(Eud,2), RRgo = HRGO[J], n2 = N2[I],
                           n3 = n3, n = N2[I] + n3,
                           pgo = round(pg,2), sProg = round(prob,2),
                           w = w, p0 = p0, p11 = p11, p12 = p12, in1 = in1, in2 = in2,
                           K = K, K2 = round(k2), K3 = round(k3),
                           sProg1 = round(prob1,2), sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                           steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                           alpha = alpha, beta = beta, c02 = c02,
                           c03 = c03, c2 = c2, c3 = c3, b1 = b1, b2 = b2, b3 = b3, gamma = gamma)
   }
   
   if(skipII){
     result <- merge(result,result_skipII, all = TRUE)
   }
   comment(result) <-   c("\noptimization sequence RRgo:", HRGO,
                          "\noptimization sequence n2:", N2,
                          "\nonset date:", as.character(date),
                          "\nfinish date:", as.character(Sys.time()))
   class(result) <- c("drugdevelopResult", class(result))
   
   parallel::stopCluster(cl)
   
   return(result)
}

