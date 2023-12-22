#' Optimal phase II/III drug development planning where several phase III trials are performed
#'
#' The `optimal_multitrial_binary` function enables planning of phase II/III
#' drug development programs with several phase III trials for the same
#' binary endpoint. The main output values are optimal sample size allocation
#' and go/no-go decision rules. For binary endpoints, the treatment effect is
#' measured by the risk ratio (RR).
#' 
#' The assumed true treatment effects can be assumed fixed or modelled by a
#' prior distribution. The R Shiny application
#' \href{https://web.imbi.uni-heidelberg.de/prior/}{prior} visualizes the prior 
#' distributions used in this package. 
#' 
#' Fast computing is enabled by parallel programming.
#' 
#' @name optimal_multitrial_binary
#' 
#' @inheritParams optimal_multitrial_generic
#' @inheritParams optimal_binary_generic
#' 
#' @inheritSection optimal_multitrial_generic Effect sizes
#' 
#' @return
#' `r optimal_return_doc(type = "binary", setting = "multitrial")`
#' 
#' @importFrom progressr progressor
#' 
#' @examples
#' # Activate progress bar (optional)
#' \dontrun{progressr::handlers(global = TRUE)}
#' # Optimize
#' \donttest{
#' optimal_multitrial_binary(w = 0.3,         # define parameters for prior
#'   p0 = 0.6, p11 =  0.3, p12 = 0.5,
#'   in1 = 30, in2 = 60,                             # (https://web.imbi.uni-heidelberg.de/prior/)
#'   n2min = 20, n2max = 100, stepn2 = 4,            # define optimization set for n2
#'   rrgomin = 0.7, rrgomax = 0.9, steprrgo = 0.05,  # define optimization set for RRgo
#'   alpha = 0.025, beta = 0.1,                      # drug development planning parameters
#'   c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,        # fixed and variable costs for phase II/III,
#'   K = Inf, N = Inf, S = -Inf,                     # set constraints
#'   b1 = 1000, b2 = 2000, b3 = 3000,                # expected benefit for a each effect size
#'   case = 1, strategy = TRUE,                      # chose Case and Strategy                                   
#'   fixed = TRUE,                                   # true treatment effects are fixed/random
#'   num_cl = 1)                                     # number of cores for parallelized computing
#'   }
#' @references
#' IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at \href{https://www.iqwig.de/ueber-uns/methoden/methodenpapier/}{https://www.iqwig.de/ueber-uns/methoden/methodenpapier/}, assessed last 15.05.19.
#' @export

optimal_multitrial_binary <- function(w, p0, p11, p12, in1, in2,
                               n2min, n2max, stepn2,
                               rrgomin, rrgomax, steprrgo,
                               alpha, beta, 
                               c2, c3, c02, c03, 
                               K = Inf, N = Inf, S = -Inf,
                               b1, b2, b3,
                               case, strategy = TRUE,
                               fixed = FALSE,  num_cl = 1){
  
  result <- result23 <- NULL
  
  # spezifications for one phase III trial
  steps1 = 1; stepm1 = 0.95; stepl1 = 0.85
  steps2  <- stepm1
  stepm2  <- stepl1
  stepl2  <- 0
  gamma   <- 0
  ymin <- -log(0.8)
  
  
  alpha_in <- alpha
  
  date <- Sys.time()
  
  RRGO <- seq(rrgomin, rrgomax, steprrgo)
  N2   <- seq(n2min, n2max, stepn2)
  
  if(!is.numeric(strategy)){
    if(case==1){
      # Strategy 1alpha vs. Strategy 1/2,
      STRATEGY = c(1, 2)
    }
    if(case==2){
      # Strategy 1alpha^2 vs. Strategy 2/2 vs. Strategy 2/3 vs. Strategy 2/2( + 1)
      STRATEGY = c(1, 2, 3, 23)
    }
    if(case==3){
      # Strategy 1alpha^3 vs. Strategy 3/3 vs. Strategy 3/4
      STRATEGY = c(1, 3, 4)
    }  
  }else{
    STRATEGY = strategy
  }
  Strategy <- NA_real_
  RRgo <- NA_real_
  cl <-  parallel::makeCluster(getOption("cl.cores", num_cl)) #define cluster
  parallel::clusterExport(cl, c("pmvnorm", "dmvnorm", "prior_binary", "Epgo_binary", "Epgo23_binary", "En3_binary",
                                "EPsProg_binary", "EPsProg2_binary", "EPsProg3_binary", "EPsProg4_binary", "EPsProg23_binary",
                                "alpha", "beta", "t1", "t2", "t3",
                                "steps1", "steps2", "stepm1", "stepm2", "stepl1", "stepl2",
                                "K", "N", "S", "fixed",
                                "c2", "c3", "c02", "c03",
                                "b1", "b2", "b3", "w", "RRgo", "ymin",
                                "p0", "p11", "p12", "in1", "in2"), envir = environment())
  
  for(Strategy in STRATEGY){
    
    ufkt <- spfkt <- pgofkt <- K2fkt <- K3fkt <-
      sp1fkt <- sp2fkt <- sp3fkt <- n2fkt <- n3fkt <- pgo3fkt <- 
      n33fkt<-  sp13fkt   <- sp23fkt   <- sp33fkt <- matrix(0, length(N2), length(RRGO))
    
    pb <- progressr::progressor(steps = length(STRATEGY)*length(RRGO),
                                label = "Optimization progress",
                                message = "Optimization progress")
    pb(paste("Performing optimization for strategy", Strategy),
       class = "sticky", amount = 0)
    
    
    
    for(j in 1:length(RRGO)){
      
      RRgo <- RRGO[j]
      
      
      
      ###################
      # Strategy 1alpha #
      ###################
      if(Strategy == 1){
        if(case==1){
          alpha <- alpha_in
        }
        if(case==2){
          alpha <- alpha_in^2
        }
        if(case==3){
          alpha <- alpha_in^3
        }  
      }else{
        alpha <- alpha_in
      }
      
      
      
      if(Strategy==1){
        res <- parallel::parSapply(cl, N2, utility_binary, RRgo, w, p0, p11, p12, in1, in2,
                         alpha, beta, 
                         c2, c3, c02, c03, K, N, S,
                         steps1, stepm1, stepl1,
                         b1, b2, b3,
                         gamma, fixed)  
      }
      if(Strategy==2){
        res <- parallel::parSapply(cl, N2, utility2_binary, RRgo, w, p0, p11, p12, in1, in2,
                         alpha, beta, 
                         c2, c3, c02, c03, 
                         K, N, S,
                         b1, b2, b3,
                         case, fixed)  
      }
      if(Strategy==3){
        res <- parallel::parSapply(cl, N2, utility3_binary, RRgo, w, p0, p11, p12, in1, in2,
                         alpha, beta, 
                         c2, c3, c02, c03,
                         K, N, S,
                         b1, b2, b3,
                         case, fixed)  
      }
      if(Strategy==23){
        res <- parallel::parSapply(cl, N2, utility23_binary, RRgo, w, p0, p11, p12, in1, in2,
                         alpha, beta, 
                         c2, c3, c02, c03, 
                         b1, b2, b3)  
      }
      if(Strategy==4){
        res <- parallel::parSapply(cl, N2, utility4_binary, RRgo, w, p0, p11, p12, in1, in2,
                         alpha, beta, 
                         c2, c3, c02, c03, 
                         K, N, S,
                         b1, b2, b3,
                         case, fixed)  
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
      
      if(Strategy==23){
        pgo3fkt[, j]    <-  res[10, ]
        n33fkt[, j]     <-  res[11, ]
        sp13fkt[, j]    <-  res[12, ]
        sp23fkt[, j]    <-  res[13, ]
        sp33fkt[, j]    <-  res[14, ]
      }
      
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
    
    if(Strategy==23){
      n33    <- n33fkt[I, J]
      pg3    <- pgo3fkt[I, J]
      prob13 <- sp13fkt[I, J]
      prob23 <- sp23fkt[I, J]
      prob33 <- sp33fkt[I, J]
      
    }else{
      n33    <- 0
      pg3    <- 0
      prob13 <- 0
      prob23 <- 0
      prob33 <- 0
    }
    
    
    
    if(!fixed){
      
      result <-  rbind(result, data.frame(Case = case, Strategy = Strategy, 
                                          u = round(Eud,2), RRgo = RRGO[J], n2 = N2[I],
                                          n3 = n3, n = N2[I] + n3,
                                          pgo = round(pg,2), sProg = round(prob,2),
                                          w = w, p0 = p0, p11 = p11, p12 = p12, in1 = in1, in2 = in2,
                                          K = K, N = N, S = S, K2 = round(k2), K3 = round(k3),
                                          sProg1 = round(prob1,2), sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                                          steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                                          pgo3 = round(pg3,2),  n33 = n33,
                                          sProg13 = round(prob13,2), sProg23 = round(prob23,2), sProg33 = round(prob33,2),
                                          alpha = alpha, beta = beta, c02 = c02,
                                          c03 = c03, c2 = c2, c3 = c3, b1 = b1, b2 = b2, b3 = b3))
    }else{
      
      result <-  rbind(result, data.frame(Case = case, Strategy = Strategy, 
                                          u = round(Eud,2), RRgo = RRGO[J], n2 = N2[I],
                                          n3 = n3, n = N2[I] + n3,
                                          pgo = round(pg,2), sProg = round(prob,2),
                                          p0 = p0, p1 = p11,
                                          K = K, N = N, S = S, K2 = round(k2), K3 = round(k3),
                                          sProg1 = round(prob1,2), sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                                          steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                                          alpha = alpha, beta = beta, c02 = c02,
                                          c03 = c03, c2 = c2, c3 = c3, b1 = b1, b2 = b2, b3 = b3))
    }
    
    
    
    comment(result) <-   c("\noptimization sequence RRgo:", RRGO,
                           "\noptimization sequence n2:", N2,
                           "\nonset date:", as.character(date),
                           "\nfinish date:", as.character(Sys.time()))
    
  }

  class(result) <- c("drugdevelopResult", class(result))
  parallel::stopCluster(cl)
  
  return(result)
  
}
