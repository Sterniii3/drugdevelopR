optimal_multiarm_normal <- function(Delta1, Delta2, ec,
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
  
  for(strategy in STRATEGY){
    
    ufkt <- spfkt <- pgofkt <- K2fkt <- K3fkt <-
      sp2fkt <- sp3fkt <- n3fkt <- matrix(0, length(N2), length(KAPPA))
    
    cat("Optimization progress:", fill = TRUE)
    cat("", fill = TRUE)
    pb <- txtProgressBar(min = 0, max = length(KAPPA), style = 3, label = "Optimization progess")
    
    for(j in 1:length(KAPPA)){
      
      kappa <- KAPPA[j]
      
      cl <-  makeCluster(getOption("cl.cores", num_cl)) #define cluster
      
      clusterExport(cl, c("pmvnorm", "dmvnorm","qmvnorm","adaptIntegrate", "pgo_normal", "ss_normal", "Ess_normal",
                          "PsProg_normal", "alpha", "beta",
                          "steps1", "steps2", "stepm1", "stepm2", "stepl1", "stepl2",
                          "K", "N", "S", "strategy",
                          "c2", "c3", "c02", "c03",
                          "b1", "b2", "b3", "KAPPA",
                          "Delta1", "Delta2", "ec"), envir = environment())
      
      
      res <- parSapply(cl, N2, utility_multiarm_normal, kappa,
                       alpha,beta, Delta1,Delta2,strategy,ec,
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
    
    result <-  rbind(result, data.frame(Strategy = strategy,u = round(Eud,2), kappa = KAPPA[J], n2 = N2[I], 
                                        n3 = n3, n = N2[I] + n3,
                                        pgo = round(pg,2), sProg = round(prob,2),
                                        Delta1 = Delta1, Delta2 = Delta2, ec = ec,
                                        K = K, N = N, S = S, K2 = round(k2), K3 = round(k3),
                                        sProg2 = round(prob2,2), sProg3 = round(prob3,2),
                                        steps1 = round(steps1,2), stepm1 = round(stepm1,2), stepl1 = round(stepl1,2),
                                        alpha = alpha, beta = beta, 
                                        c02 = c02, c03 = c03, c2 = c2, c3 = c3, 
                                        b1 = b1, b2 = b2, b3 = b3))  
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
  
  return(result)
  
}