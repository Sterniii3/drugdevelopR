#' Optimal phase II/III drug development planning where several phase III trials are performed
#'
#' The `optimal_multitrial_normal` function enables planning of phase II/III
#' drug development programs with several phase III trials for
#' the same normally distributed endpoint. Its main output values are optimal
#' sample size allocation and go/no-go decision rules. For normally distributed 
#' endpoints, the treatment effect is measured by the standardized difference in
#' means (Delta). The assumed true treatment effects can be assumed fixed or
#' modelled by a prior distribution.
#' 
#' The R Shiny application \href{https://web.imbi.uni-heidelberg.de/prior/}{prior} 
#' visualizes the prior distributions used in this package. Fast computing is
#' enabled by parallel programming.
#' 
#' @name optimal_multitrial_normal
#' 
#' @inheritParams optimal_multitrial_generic
#' @inheritParams optimal_normal_generic
#' 
#' @inheritSection optimal_multitrial_generic Effect sizes
#' 
#' @return
#' `r optimal_return_doc(type = "normal", setting = "multitrial")`
#' 
#' @examples
#' # Activate progress bar (optional)
#' \dontrun{progressr::handlers(global = TRUE)}
#' # Optimize
#' \donttest{
#' optimal_multitrial_normal(w = 0.3,           # define parameters for prior
#'   Delta1 = 0.375, Delta2 = 0.625,
#'   in1 = 300, in2 = 600,                               # (https://web.imbi.uni-heidelberg.de/prior/)
#'   a = 0.25, b = 0.75,
#'   n2min = 20, n2max = 100, stepn2 = 4,                # define optimization set for n2
#'   kappamin = 0.02, kappamax = 0.2, stepkappa = 0.02,  # define optimization set for kappa
#'   alpha = 0.025, beta = 0.1,                          # drug development planning parameters
#'   c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,          # fixed and variable costs for phase II/III
#'   K = Inf, N = Inf, S = -Inf,                         # set constraints
#'   b1 = 3000, b2 = 8000, b3 = 10000,                   # expected benefit for each effect size                                         
#'   case = 1, strategy = TRUE,                          # chose Case and Strategy
#'   fixed = TRUE,                                       # true treatment effects are fixed/random
#'   num_cl = 1)                                         # number of cores for parallelized computing
#'   }
#' 
#' @references
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences.
#'
#' @export

optimal_multitrial_normal <- function(w, Delta1, Delta2, in1, in2, a, b,
                                      n2min, n2max, stepn2,
                                      kappamin, kappamax, stepkappa,
                                      alpha, beta, 
                                      c2, c3, c02, c03, 
                                      K = Inf, N = Inf, S = -Inf,
                                      b1, b2, b3,
                                      case, strategy = TRUE,
                                      fixed = FALSE,  num_cl = 1){
  
  result <- result23 <- NULL
  
  # spezifications for one phase III trial
  steps1 = 0; stepm1 = 0.5; stepl1 = 0.8
  steps2  <- stepm1
  stepm2  <- stepl1
  stepl2  <- Inf
  gamma   <- 0
  ymin <- 0.8
  
  
  alpha_in <- alpha
  
  date <- Sys.time()
  
  KAPPA <- seq(kappamin, kappamax, stepkappa)
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
  kappa <- NA_real_
  Strategy <- NA_real_
  cl <-  parallel::makeCluster(getOption("cl.cores", num_cl)) #define cluster
  parallel::clusterExport(cl, c("pmvnorm", "dmvnorm", "dtnorm", "prior_normal", "Epgo_normal", "Epgo23_normal", "En3_normal",
                                "EPsProg_normal", "EPsProg2_normal", "EPsProg3_normal", "EPsProg4_normal", "EPsProg23_normal",
                                "alpha", "beta",
                                "steps1", "steps2", "stepm1", "stepm2", "stepl1", "stepl2",
                                "K", "N", "S", "fixed", "case", "Strategy",
                                "b1", "b2", "b3", "w", "kappa",
                                "Delta1", "Delta2", "ymin", "in1", "in2", "a", "b" ), envir = environment())
  
  for(Strategy in STRATEGY){
    
    ufkt <- spfkt <- pgofkt <- K2fkt <- K3fkt <-
      sp1fkt <- sp2fkt <- sp3fkt <- n2fkt <- n3fkt <- pgo3fkt <- 
      n33fkt<-  sp13fkt   <- sp23fkt   <- sp33fkt <- matrix(0, length(N2), length(KAPPA))
    
    pb <- progressr::progressor(steps = length(STRATEGY)*length(KAPPA),
                                label = "Optimization progress",
                                message = "Optimization progress")
    pb(paste("Performing optimization for strategy", Strategy),
       class = "sticky", amount = 0)
    
    for(j in 1:length(KAPPA)){
      
      kappa <- KAPPA[j]
      
      
      
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
        res <- parallel::parSapply(cl, N2, utility_normal, kappa, w, Delta1, Delta2, in1, in2, a, b,
                         alpha, beta, 
                         c2, c3, c02, c03, 
                         K, N, S,
                         steps1, stepm1, stepl1,
                         b1, b2, b3,
                         gamma, fixed)  
      }
      if(Strategy==2){
        res <- parallel::parSapply(cl, N2, utility2_normal, kappa, w, Delta1, Delta2, in1, in2, a, b,
                         alpha, beta, 
                         c2, c3, c02, c03, 
                         K, N, S,
                         b1, b2, b3,
                         case, fixed)  
      }
      if(Strategy==3){
        res <- parallel::parSapply(cl, N2, utility3_normal, kappa, w, Delta1, Delta2, in1, in2, a, b,
                         alpha, beta, 
                         c2, c3, c02, c03, 
                         K, N, S,
                         b1, b2, b3,
                         case, fixed)  
      }
      if(Strategy==23){
        res <- parallel::parSapply(cl, N2, utility23_normal, kappa, w, Delta1, Delta2, in1, in2, a, b,
                         alpha, beta, 
                         c2, c3, c02, c03, 
                         b1, b2, b3)  
      }
      if(Strategy==4){
        res <- parallel::parSapply(cl, N2, utility4_normal, kappa, w, Delta1, Delta2, in1, in2, a, b,
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
      pg3    <- pgo3fkt[I, J]
      prob13 <- sp13fkt[I, J]
      prob23 <- sp23fkt[I, J]
      prob33 <- sp33fkt[I, J]
      n33    <- n33fkt[I,J]
    }else{
      pg3    <- 0
      prob13 <- 0
      prob23 <- 0
      prob33 <- 0
      n33    <- 0
    }
    
    
    if(!fixed){
      
      result <-  rbind(result, data.frame(Case = case, Strategy = Strategy, 
                                          u = round(Eud,2), Kappa = KAPPA[J], n2 = N2[I],
                                          n3 = n3, n = N2[I] + n3,
                                          pgo = round(pg,2), sProg = round(prob,2),
                                          w = w, Delta1 = Delta1, Delta2 = Delta2, in1 = in1, in2 = in2, a = a, b = b,
                                          K = K, N = N, S = S, K2 = round(k2), K3 = round(k3),
                                          sProg1 = round(prob1,2), sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                                          steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                                          pgo3 = round(pg3,2), n33 = n33,
                                          sProg13 = round(prob13,2), sProg23 = round(prob23,2), sProg33 = round(prob33,2),
                                          alpha = alpha, beta = beta, c02 = c02,
                                          c03 = c03, c2 = c2, c3 = c3, b1 = b1, b2 = b2, b3 = b3))
    }else{
      
      result <-  rbind(result, data.frame(Case = case, Strategy = Strategy, 
                                          u = round(Eud,2), Kappa = KAPPA[J], n2 = N2[I],
                                          n3 = n3, n = N2[I] + n3,
                                          pgo = round(pg,2), sProg = round(prob,2),
                                          Delta = Delta1,
                                          K = K, N = N, S = S, K2 = round(k2), K3 = round(k3),
                                          sProg1 = round(prob1,2), sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                                          steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                                          alpha = alpha, beta = beta, c02 = c02,
                                          c03 = c03, c2 = c2, c3 = c3, b1 = b1, b2 = b2, b3 = b3))
    }
    
    
    
    comment(result) <-   c("\noptimization sequence kappa:", KAPPA,
                           "\noptimization sequence n2:", N2,
                           "\nonset date:", as.character(date),
                           "\nfinish date:", as.character(Sys.time()))
    
  }
  class(result) <- c("drugdevelopResult", class(result))
  parallel::stopCluster(cl)
  
  return(result)
  
}
