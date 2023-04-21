#' Optimal phase II/III drug development planning where several phase III trials are performed for time-to-event endpoints
#'
#' The function \code{\link{optimal_multitrial}} of the drugdevelopR package enables planning of phase II/III drug development programs with time-to-event endpoints for programs with several phase III trials (Preussler et. al, 2019). 
#' Its main output values are the optimal sample size allocation and optimal go/no-go decision rules. 
#' The assumed true treatment effects can be assumed to be fixed (planning is then also possible via user friendly R Shiny App: \href{https://web.imbi.uni-heidelberg.de/multitrial/}{multitrial}) or can be modelled by a prior distribution. 
#' The R Shiny application \href{https://web.imbi.uni-heidelberg.de/prior/}{prior} visualizes the prior distributions used in this package. Fast computing is enabled by parallel programming.
#' 
#' @name optimal_multitrial
#' 
#' @inheritParams optimal_multitrial_generic
#' @inheritParams optimal_tte_generic
#' 
#' @inheritSection optimal_multitrial_generic Effect sizes
#' 
#' @return
#' `r optimal_return_doc(type = "tte", setting = "multitrial")`
#' 
#' 
#' @examples
#' res <- optimal_multitrial(w = 0.3,                       # define parameters for prior
#'   hr1 = 0.69, hr2 = 0.88, id1 = 210, id2 = 420,          # (https://web.imbi.uni-heidelberg.de/prior/)
#'   d2min = 20, d2max = 100, stepd2 = 5,                   # define optimization set for d2
#'   hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05,         # define optimization set for HRgo
#'   alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,       # drug development planning parameters
#'   c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,               # define fixed and variable costs for phase II and III
#'   K = Inf, N = Inf, S = -Inf,                            # set constraints
#'   b1 = 1000, b2 = 2000, b3 = 3000,                       # define expected benefit for a each effect size category
#'   case = 1, strategy = TRUE,                             # chose Case and Strategy
#'   fixed = TRUE,                                          # choose if true treatment effects are fixed or random
#'   num_cl = 1)                                            # set number of cores used for parallelized computing 
#' res
#' cat(comment(res))                                        # displays optimization sequence, start/ finish date of procedure.
#' 
#' @references
#' IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at \href{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}, assessed last 15.05.19.
#'
#' Preussler, S., Kieser, M., and Kirchner, M. (2019). Optimal sample size allocation and go/no-go decision rules for phase II/III programs where several phase III trials are performed. Biometrical Journal, 61(2), 357-378.
#'
#' Schoenfeld, D. (1981). The asymptotic properties of nonparametric tests for comparing survival distributions. Biometrika, 68(1), 316-319.
#'
#' @export
optimal_multitrial <- function(w,  hr1, hr2, id1, id2,
                        d2min, d2max, stepd2,
                        hrgomin, hrgomax, stephrgo,
                        alpha, beta, xi2, xi3,
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
  
  HRGO <- seq(hrgomin, hrgomax, stephrgo)
  D2   <- seq(d2min, d2max, stepd2)
  
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

  for(Strategy in STRATEGY){
    
    ufkt <- d3fkt <- spfkt <- pgofkt <- K2fkt <- K3fkt <-
      sp1fkt <- sp2fkt <- sp3fkt <- n2fkt <- n3fkt <- pgo3fkt <- 
      d33fkt <-  n33fkt<-  sp13fkt   <- sp23fkt   <- sp33fkt <- matrix(0, length(D2), length(HRGO))
    
    cat("", fill = TRUE)
    cat(paste("Case ", case,": Optimization progess for Strategy ", Strategy), fill = TRUE)
    cat("", fill = TRUE)
    pb <- txtProgressBar(min = 0, max = length(HRGO), style = 3, label = "Optimization progess")

    for(j in 1:length(HRGO)){
      
      HRgo <- HRGO[j]
      
      cl <-  parallel::makeCluster(getOption("cl.cores", num_cl)) #define cluster
    
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

      parallel::clusterExport(cl, c("pmvnorm", "dmvnorm", "prior_tte", "Epgo_tte", "Epgo23", "Ed3_tte",
                          "EPsProg_tte", "EPsProg2", "EPsProg3", "EPsProg4", "EPsProg23",
                          "alpha", "beta",
                          "steps1", "steps2", "stepm1", "stepm2", "stepl1", "stepl2",
                          "K", "N", "S","gamma", "fixed", "case", "Strategy",
                          "xi2", "xi3", "c2", "c3", "c02", "c03",
                          "b1", "b2", "b3", "w", "HRgo", "ymin",
                          "hr1", "hr2", "id1", "id2"), envir = environment())
      
      if(Strategy==1){
        res <- parallel::parSapply(cl, D2, utility_tte, HRgo, w, hr1, hr2, id1, id2,
                            alpha, beta, xi2, xi3,
                            c2, c3, c02, c03, 
                            K, N, S,
                            steps1, stepm1, stepl1,
                            b1, b2, b3,
                            gamma, fixed)  
      }
      if(Strategy==2){
        res <- parallel::parSapply(cl, D2, utility2, HRgo, w, hr1, hr2, id1, id2,
                            alpha, beta, xi2, xi3,
                            c2, c3, c02, c03, 
                            K, N, S,
                            b1, b2, b3,
                            case, fixed)  
      }
      if(Strategy==3){
        res <- parallel::parSapply(cl, D2, utility3, HRgo, w, hr1, hr2, id1, id2,
                            alpha, beta, xi2, xi3,
                            c2, c3, c02, c03, 
                            K, N, S,
                            b1, b2, b3,
                            case, fixed)  
      }
      if(Strategy==23){
        res <- parallel::parSapply(cl, D2, utility23, HRgo, w, hr1, hr2, id1, id2,
                         alpha, beta, xi2, xi3,
                         c2, c3, c02, c03, 
                         b1, b2, b3)  
      }
      if(Strategy==4){
        res <- parallel::parSapply(cl, D2, utility4, HRgo, w, hr1, hr2, id1, id2,
                            alpha, beta, xi2, xi3,
                            c2, c3, c02, c03, 
                            K, N, S,
                            b1, b2, b3,
                            case, fixed)  
      }
    
    setTxtProgressBar(title= "i", pb, j)
    parallel::stopCluster(cl)
  
    ufkt[, j]      <-  res[1, ]
    d3fkt[, j]     <-  res[2, ]
    spfkt[, j]     <-  res[3, ]
    pgofkt[, j]    <-  res[4, ]
    K2fkt[, j]     <-  res[5, ]
    K3fkt[, j]     <-  res[6, ]
    sp1fkt[, j]    <-  res[7, ]
    sp2fkt[, j]    <-  res[8, ]
    sp3fkt[, j]    <-  res[9, ]
    n2fkt[, j]     <-  res[10, ]
    n3fkt[, j]     <-  res[11, ]

    if(Strategy==23){
      pgo3fkt[, j]    <-  res[12, ]
      d33fkt[, j]     <-  res[13, ]
      n33fkt[, j]     <-  res[14, ]
      sp13fkt[, j]    <-  res[15, ]
      sp23fkt[, j]    <-  res[16, ]
      sp33fkt[, j]    <-  res[17, ]
    }
  
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

    if(Strategy==23){
      d33    <- d33fkt[I, J]
      pg3    <- pgo3fkt[I, J]
      prob13 <- sp13fkt[I, J]
      prob23 <- sp23fkt[I, J]
      prob33 <- sp33fkt[I, J]
      n33    <- n33fkt[I,J]
    }else{
      d33    <- 0
      pg3    <- 0
      prob13 <- 0
      prob23 <- 0
      prob33 <- 0
      n33    <- 0
    }
    
    
    close(pb)
    
    if(!fixed){

      result <-  rbind(result, data.frame(Case = case, Strategy = Strategy, 
                            u = round(Eud,2), HRgo = HRGO[J], d2 = D2[I], d3 = d3, d = D2[I] + d3,
                            n2 = n2, n3 = n3, n = n2 + n3,
                            pgo = round(pg,2), sProg = round(prob,2),
                            w = w, hr1 = hr1, hr2 = hr2, id1 = id1, id2 = id2,
                            K = K, N = N, S = S, K2 = round(k2), K3 = round(k3),
                            sProg1 = round(prob1,2), sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                            steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                            pgo3 = round(pg3,2), d33= d33, n33 = n33,
                            sProg13 = round(prob13,2), sProg23 = round(prob23,2), sProg33 = round(prob33,2),
                            alpha = alpha, beta = beta, xi2 = xi2, xi3 = xi3, c02 = c02,
                            c03 = c03, c2 = c2, c3 = c3, b1 = b1, b2 = b2, b3 = b3, gamma = gamma))
    }else{
    
      result <-  rbind(result, data.frame(Case = case, Strategy = Strategy, 
                            u = round(Eud,2), HRgo = HRGO[J], d2 = D2[I], d3 = d3, d = D2[I] + d3,
                            n2 = n2, n3 = n3, n = n2 + n3,
                            pgo = round(pg,2), sProg = round(prob,2),
                            hr = hr1,
                            K = K, N = N, S = S, K2 = round(k2), K3 = round(k3),
                            sProg1 = round(prob1,2), sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                            steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                            alpha = alpha, beta = beta, xi2 = xi2, xi3 = xi3, c02 = c02,
                            c03 = c03, c2 = c2, c3 = c3, b1 = b1, b2 = b2, b3 = b3, gamma = gamma))
    }
    

    
    comment(result) <-   c("\noptimization sequence HRgo:", HRGO,
                  "\noptimization sequence d2:", D2,
                  "\nset on date:", as.character(date),
                  "\nfinish date:", as.character(Sys.time()))

  }
  cat("", fill = TRUE)
  cat("", fill = TRUE)
  cat("Optimization result:", fill = TRUE)
  cat("", fill = TRUE)
  print(result)
  cat("", fill = TRUE)
  
  return(result)
  
}
