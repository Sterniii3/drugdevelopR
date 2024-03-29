#' Optimal phase II/III drug development planning for programs with multiple
#' time-to-event endpoints
#'
#' The function \code{\link{optimal_multiple_tte}} of the drugdevelopR package
#' enables planning of phase II/III drug development programs with optimal
#' sample size allocation and go/no-go decision rules (Preussler et. al, 2019)
#' in a two-arm trial with two time-to-event endpoints.
#'
#' In this setting, the drug development program is defined to be successful if
#' it proceeds from phase II to phase III and at least one endpoint shows a
#' statistically significant treatment effect in phase III. For example,
#' this situation is found in oncology trials, where overall survival (OS)
#' and progression-free survival (PFS) are the two endpoints of interest.
#'
#' The gain of a successful program may differ according to the importance of
#' the endpoint that is significant. If endpoint 1 is significant (no matter
#' whether endpoint 2 is significant or not), then the gains `b11`, `b21`
#' and `b31` will be used for calculation of the utility. If only endpoint 2
#' is significant, then  `b12`, `b22` and `b32` will be used. This
#' also matches the oncology example, where OS (i.e. endpoint 1) implicates
#' larger expected gains than PFS alone (i.e. endpoint 2).
#'
#' Fast computing is enabled by parallel programming.
#' 
#' Monte Carlo simulations are applied for calculating utility, event count and
#' other operating characteristics in this setting. Hence, the results are affected
#' by random uncertainty. The extent of uncertainty is discussed in 
#' (Kieser et al. 2018).
#'
#' @name optimal_multiple_tte
#' @inheritParams optimal_multiple_generic
#' @inheritParams optimal_tte_generic
#' @param hr1 assumed true treatment effect on HR scale for endpoint 1 (e.g. OS)
#' @param hr2 assumed true treatment effect on HR scale for endpoint 2 (e.g. PFS)
#' @param id1 amount of information for hr1 in terms of number of events
#' @param id2 amount of information for hr2 in terms of number of events
#' @param beta type-II error rate for any pair, i.e. `1 - beta` is the (any-pair) power for calculation of the number of events for phase III
#' @param b11 expected gain for effect size category `"small"` if endpoint 1 is significant (and endpoint 2 may or may not be significant)
#' @param b21 expected gain for effect size category `"medium"` if endpoint 1 is significant (and endpoint 2 may or may not be significant)
#' @param b31 expected gain for effect size category `"large"` if endpoint 1 is significant (and endpoint 2 may or may not be significant)
#' @param b12 expected gain for effect size category `"small"` if endpoint 1 is not significant, but endpoint 2 is
#' @param b22 expected gain for effect size category `"medium"`if endpoint 1 is not significant, but endpoint 2 is
#' @param b32 expected gain for effect size category `"large"` if endpoint 1 is not significant, but endpoint 2 is
#'
#' @return
#' `r optimal_return_doc(type = "tte", setting = "multiple")`
#' 
#' @importFrom progressr progressor
#'
#' @examples
#' # Activate progress bar (optional)
#' \dontrun{progressr::handlers(global = TRUE)}
#' # Optimize
#' \donttest{
#' set.seed(123) # This function relies on Monte Carlo integration
#' optimal_multiple_tte(hr1 = 0.75,
#'   hr2 = 0.80, id1 = 210, id2 = 420,          # define assumed true HRs
#'   n2min = 30, n2max = 90, stepn2 = 6,        # define optimization set for n2
#'   hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05, # define optimization set for HRgo
#'   alpha = 0.025, beta = 0.1,                 # drug development planning parameters
#'   c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,   # fixed/variable costs for phase II/III
#'   K = Inf, N = Inf, S = -Inf,                # set constraints
#'   steps1 = 1,                                # define lower boundary for "small"
#'   stepm1 = 0.95,                             # "medium"
#'   stepl1 = 0.85,                             # and "large" effect size categories
#'   b11 = 1000, b21 = 2000, b31 = 3000,
#'   b12 = 1000, b22 = 1500, b32 = 2000,        # define expected benefits (both scenarios)
#'   rho = 0.6, fixed = TRUE,                   # correlation and treatment effect
#'   num_cl = 1)                                # number of cores for parallelized computing
#'  }
#'
#' @references
#' Kieser, M., Kirchner, M. Dölger, E., Götte, H. (2018).Optimal planning of phase II/III programs for clinical trials with multiple endpoints, Pharm Stat. 2018 Sep; 17(5):437-457.
#' 
#' Preussler, S., Kirchner, M., Goette, H., Kieser, M. (2019). Optimal Designs for Multi-Arm Phase II/III Drug Development Programs. Submitted to peer-review journal.
#'
#' IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at \href{https://www.iqwig.de/ueber-uns/methoden/methodenpapier/}{https://www.iqwig.de/ueber-uns/methoden/methodenpapier/}, assessed last 15.05.19.
#' @export

optimal_multiple_tte <- function(hr1,
                                 hr2,
                                 id1,
                                 id2,
                                 n2min,
                                 n2max,
                                 stepn2,
                                 hrgomin,
                                 hrgomax,
                                 stephrgo,
                                 alpha,
                                 beta,
                                 c2,
                                 c3,
                                 c02,
                                 c03,
                                 K = Inf,
                                 N = Inf,
                                 S = -Inf,
                                 b11,
                                 b21,
                                 b31,
                                 b12,
                                 b22,
                                 b32,
                                 steps1 = 1,
                                 stepm1 = 0.95,
                                 stepl1 = 0.85,
                                 rho,
                                 fixed = TRUE,
                                 num_cl = 1) {
  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 0
  
  rsamp <- NULL
  if(!fixed){
    rsamp <- get_sample_multiple_tte(hr1, hr2, id1, id2, rho)
  }
  
  date <- Sys.time()
  
  HRGO <- seq(hrgomin, hrgomax, stephrgo)
  N2   <- seq(n2min, n2max, stepn2)
  
  
  result <- NULL
  
  
  ufkt <- spfkt <- pgofkt <- K2fkt <- K3fkt <-
    sp2fkt <-
    sp3fkt <- n3fkt <- OSfkt <- matrix(0, length(N2), length(HRGO))
  
  pb <- progressr::progressor(steps = length(HRGO),
                              label = "Optimization progress",
                              message = "Optimization progress")
  pb("Performing optimization",
     class = "sticky",
     amount = 0)
  HRgo <- NA_real_
  cl <-
    parallel::makeCluster(getOption("cl.cores", num_cl)) #define cluster
  parallel::clusterExport(
    cl,
    c(
      "pnorm",
      "pmvnorm",
      "dnorm",
      "dmvnorm",
      "qnorm",
      "qmvnorm",
      "adaptIntegrate",
      "dbivanorm",
      "fmax",
      "pgo_multiple_tte",
      "pw",
      "Ess_multiple_tte",
      "EPsProg_multiple_tte",
      "os_tte",
      "alpha",
      "beta",
      "steps1",
      "steps2",
      "stepm1",
      "stepm2",
      "stepl1",
      "stepl2",
      "K",
      "N",
      "S",
      "c2",
      "c3",
      "c02",
      "c03",
      "b11",
      "b21",
      "b31",
      "b12",
      "b22",
      "b32",
      "HRgo",
      "hr1",
      "hr2",
      "id1",
      "id2",
      "rho",
      "fixed",
      "rsamp"
    ),
    envir = environment()
  )
  
  for (j in 1:length(HRGO)) {
    HRgo <- HRGO[j]
    

    res <-
      parallel::parSapply(
        cl,
        N2,
        utility_multiple_tte,
        HRgo = HRgo,
        alpha,
        beta,
        hr1,
        hr2,
        id1,
        id2,
        rho = rho,
        fixed = fixed,
        c2,
        c02,
        c3,
        c03,
        K,
        N,
        S,
        steps1,
        stepm1,
        stepl1,
        b11,
        b21,
        b31,
        b12,
        b22,
        b32,
        rsamp
      )
    
    pb()
    
    ufkt[, j]     <-  res[1,]
    n3fkt[, j]    <-  res[2,]
    spfkt[, j]    <-  res[3,]
    pgofkt[, j]   <-  res[4,]
    sp2fkt[, j]   <-  res[5,]
    sp3fkt[, j]   <-  res[6,]
    K2fkt[, j]    <-  res[7,]
    K3fkt[, j]    <-  res[8,]
    OSfkt[, j]    <-  res[9,]
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
  OS    <- OSfkt[I, J]
  
  if (!fixed) {
    result <-
      rbind(
        result,
        data.frame(
          u = round(Eud, 2),
          HRgo = HRGO[J],
          n2 = N2[I],
          n3 = n3,
          n = N2[I] + n3,
          pgo = round(pg, 2),
          sProg = round(prob, 2),
          hr1 = hr1,
          hr2 = hr2,
          id1 = id1,
          id2 = id2,
          rho = rho,
          K = K,
          N = N,
          S = S,
          K2 = round(k2),
          K3 = round(k3),
          sProg2 = round(prob2, 2),
          sProg3 = round(prob3, 2),
          OS = round(OS, 2),
          steps1 = round(steps1, 2),
          stepm1 = round(stepm1, 2),
          stepl1 = round(stepl1, 2),
          alpha = alpha,
          beta = beta,
          c02 = c02,
          c03 = c03,
          c2 = c2,
          c3 = c3,
          b11 = b11,
          b21 = b21,
          b31 = b31,
          b12 = b12,
          b22 = b22,
          b32 = b32
        )
      )
  } else{
    result <-
      rbind(
        result,
        data.frame(
          u = round(Eud, 2),
          HRgo = HRGO[J],
          n2 = N2[I],
          n3 = n3,
          n = N2[I] + n3,
          pgo = round(pg, 2),
          sProg = round(prob, 2),
          hr1 = hr1,
          hr2 = hr2,
          rho = rho,
          K = K,
          N = N,
          S = S,
          K2 = round(k2),
          K3 = round(k3),
          sProg2 = round(prob2, 2),
          sProg3 = round(prob3, 2),
          OS = round(OS, 2),
          steps1 = round(steps1, 2),
          stepm1 = round(stepm1, 2),
          stepl1 = round(stepl1, 2),
          alpha = alpha,
          beta = beta,
          c02 = c02,
          c03 = c03,
          c2 = c2,
          c3 = c3,
          b11 = b11,
          b21 = b21,
          b31 = b31,
          b12 = b12,
          b22 = b22,
          b32 = b32
        )
      )
  }
  
  
  
  comment(result) <-   c(
    "\noptimization sequence HRgo:",
    HRGO,
    "\noptimization sequence n2:",
    N2,
    "\nonset date:",
    as.character(date),
    "\nfinish date:",
    as.character(Sys.time())
  )
  parallel::stopCluster(cl)
  return(drugdevelopResult(result))
  
}
