#' Optimal phase II/III drug development planning for programs with multiple
#' normally distributed endpoints
#'
#' The function \code{\link{optimal_multiple_normal}} of the drugdevelopR
#' package enables planning of phase II/III drug development programs with
#'  optimal sample size allocation and go/no-go decision rules for two-arm
#'  trials with two normally  distributed endpoints and one control group
#'  (Preussler et. al, 2019).
#'
#'  For this setting, the drug development program is defined to be successful
#'  if it proceeds from phase II to phase III and all endpoints show a
#'  statistically significant treatment effect in phase III. For example, this
#'  situation is found in Alzheimer’s disease trials, where a drug should show
#'  significant results in improving cognition (cognitive endpoint) as well as
#'  in improving activities of daily living (functional endpoint).
#'
#'  The effect size categories small, medium and large are applied to both
#'  endpoints. In order to define an overall effect size from the two individual
#'  effect sizes, the function implements two different combination rules:
#'  * A strict rule (`relaxed = FALSE`) assigning a large overall effect in case
#'    both endpoints show an effect of large size, a small overall effect in
#'    case that at least one of the endpoints shows a small effect, and a medium
#'    overall effect otherwise, and
#'  * A relaxed rule (`relaxed = TRUE`) assigning a large overall effect if at
#'    least one of the endpoints shows a large effect, a small effect if both
#'    endpoints show a  small effect, and a medium overall effect otherwise.
#'
#'  Fast computing is enabled by parallel programming.
#'
#' @name optimal_multiple_normal
#' @inheritParams optimal_multiple_generic
#' @inheritParams optimal_normal_generic
#' @param Delta1 assumed true treatment effect for endpoint 1 measured as the
#'  difference in means
#' @param Delta2 assumed true treatment effect for endpoint 2 measured as the
#'   difference in means
#' @param in1 amount of information for Delta1 in terms of number of events
#' @param in2 amount of information for Delta2 in terms of number of events
#' @param sigma1 variance of endpoint 1
#' @param sigma2 variance of endpoint 2
#' @param relaxed relaxed or strict decision rule
#' @param beta type-II error rate for any pair, i.e. `1 - beta` is the (any-pair) power for calculation of the sample size for phase III
#'
#' @importFrom stats quantile rnorm
#' @importFrom MASS mvrnorm
#'
#' @return
#' `r optimal_return_doc(type = "normal", setting = "multiple")`
#'
#' @examples
#' # Activate progress bar (optional)
#' \dontrun{progressr::handlers(global = TRUE)}
#' # Optimize
#' \donttest{
#' set.seed(123) # This function relies on Monte Carlo integration
#' optimal_multiple_normal(Delta1 = 0.75,
#'   Delta2 = 0.80, in1=300, in2=600,                   # define assumed true HRs
#'   sigma1 = 8, sigma2= 12,                            # variances for both endpoints
#'   n2min = 30, n2max = 90, stepn2 = 10,               # define optimization set for n2
#'   kappamin = 0.05, kappamax = 0.2, stepkappa = 0.05, # define optimization set for HRgo
#'   alpha = 0.025, beta = 0.1,                         # planning parameters
#'   c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,           # fixed/variable costs: phase II/III
#'   K = Inf, N = Inf, S = -Inf,                        # set constraints
#'   steps1 = 0,                                        # define lower boundary for "small"
#'   stepm1 = 0.5,                                      # "medium"
#'   stepl1 = 0.8,                                      # and "large" effect size categories
#'   b1 = 1000, b2 = 2000, b3 = 3000,                   # define expected benefit
#'   rho = 0.5, relaxed = TRUE,                         # strict or relaxed rule
#'   fixed = TRUE,                                      # treatment effect
#'   num_cl = 1)                                        # parallelized computing
#'   }
#'
#' @references
#' Meinhard Kieser, Marietta Kirchner, Eva Dölger, Heiko Götte (2018). Optimal planning of phase II/III programs for clinical trials with multiple endpoints
#'
#' IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at \href{https://www.iqwig.de/ueber-uns/methoden/methodenpapier/}{https://www.iqwig.de/ueber-uns/methoden/methodenpapier/}, assessed last 15.05.19.
#' @export

optimal_multiple_normal <-
  function(Delta1,
           Delta2,
           in1,
           in2,
           sigma1,
           sigma2,
           n2min,
           n2max,
           stepn2,
           kappamin,
           kappamax,
           stepkappa,
           alpha,
           beta,
           c2,
           c3,
           c02,
           c03,
           K = Inf,
           N = Inf,
           S = -Inf,
           steps1 = 0,
           stepm1 = 0.5,
           stepl1 = 0.8,
           b1,
           b2,
           b3,
           rho,
           fixed,
           relaxed = FALSE,
           num_cl = 1) {
    date <- Sys.time()
    
    KAPPA <- seq(kappamin, kappamax, stepkappa)
    N2   <- seq(n2min, n2max, stepn2)
    
    
    result <- NULL
    
    
    
    ufkt <- spfkt <- pgofkt <- K2fkt <- K3fkt <-
      sp2fkt <-
      sp3fkt <- n3fkt <- matrix(0, length(N2), length(KAPPA))
    
    pb <- progressr::progressor(steps = length(KAPPA),
                                label = "Optimization progress",
                                message = "Optimization progress")
    pb("Performing optimization",
       class = "sticky",
       amount = 0)
    
    
    for (j in 1:length(KAPPA)) {
      kappa <- KAPPA[j]
      
      
      cl <-
        parallel::makeCluster(getOption("cl.cores", num_cl)) #define cluster
      
      parallel::clusterExport(
        cl,
        c(
          "pmvnorm",
          "pnorm",
          "dmvnorm",
          "dnorm",
          "qmvnorm",
          "qnorm",
          "dbivanorm",
          "max",
          "min",
          "pgo_multiple_normal",
          "Ess_multiple_normal",
          "EPsProg_multiple_normal",
          "posp_normal",
          "fmin",
          "alpha",
          "beta",
          "steps1",
          "stepm1",
          "stepl1",
          "K",
          "N",
          "S",
          "c2",
          "c3",
          "c02",
          "c03",
          "b1",
          "b2",
          "b3",
          "kappa",
          "integrate",
          "sapply",
          "Delta1",
          "Delta2",
          "in1",
          "in2",
          "sigma1",
          "sigma2",
          "rho",
          "fixed",
          "relaxed",
          "get_sample_multiple_normal"
        ),
        envir = environment()
      )
      
      
      
      res <-
        parallel::parSapply(
          cl,
          N2,
          utility_multiple_normal,
          kappa = kappa,
          alpha,
          beta,
          Delta1,
          Delta2,
          in1,
          in2,
          sigma1,
          sigma2,
          rho = rho,
          fixed = fixed,
          relaxed = relaxed,
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
          b1,
          b2,
          b3
        )
      
      pb()
      parallel::stopCluster(cl)
      
      ufkt[, j]     <-  res[1,]
      n3fkt[, j]    <-  res[2,]
      spfkt[, j]    <-  res[3,]
      pgofkt[, j]   <-  res[4,]
      sp2fkt[, j]   <-  res[5,]
      sp3fkt[, j]   <-  res[6,]
      K2fkt[, j]    <-  res[7,]
      K3fkt[, j]    <-  res[8,]
      
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
    
    if (!fixed) {
      result <-
        rbind(
          result,
          data.frame(
            u = round(Eud, 2),
            Kappa = KAPPA[J],
            n2 = N2[I],
            n3 = n3,
            n = N2[I] + n3,
            pgo = round(pg, 2),
            sProg = round(prob, 3),
            Delta1 = Delta1,
            Delta2 = Delta2,
            in1 = in1,
            in2 = in2,
            sigma1 = sigma1,
            sigma2 = sigma2,
            rho = rho,
            relaxed = relaxed,
            K = K,
            K2 = round(k2),
            K3 = round(k3),
            sProg1 = round(prob - prob2 -
                             prob3, 3),
            sProg2 = round(prob2, 3),
            sProg3 = round(prob3, 3),
            steps1 = round(steps1, 2),
            stepm1 = round(stepm1, 2),
            stepl1 = round(stepl1, 2),
            alpha = alpha,
            beta = beta,
            c02 = c02,
            c03 = c03,
            c2 = c2,
            c3 = c3,
            b1 = b1,
            b2 = b2,
            b3 = b3
          )
        )
    } else{
      result <-
        rbind(
          result,
          data.frame(
            u = round(Eud, 2),
            Kappa = KAPPA[J],
            n2 = N2[I],
            n3 = n3,
            n = N2[I] + n3,
            pgo = round(pg, 2),
            sProg = round(prob, 3),
            Delta1 = Delta1,
            Delta2 = Delta2,
            sigma1 = sigma1,
            sigma2 = sigma2,
            rho = rho,
            relaxed = relaxed,
            K = K,
            N = N,
            S = S,
            K2 = round(k2),
            K3 = round(k3),
            sProg1 = round(prob - prob2 -
                             prob3, 3),
            sProg2 = round(prob2, 3),
            sProg3 = round(prob3, 3),
            steps1 = round(steps1, 2),
            stepm1 = round(stepm1, 2),
            stepl1 = round(stepl1, 2),
            alpha = alpha,
            beta = beta,
            c02 = c02,
            c03 = c03,
            c2 = c2,
            c3 = c3,
            b1 = b1,
            b2 = b2,
            b3 = b3
          )
        )
    }
    
    
    
    comment(result) <-   c(
      "\noptimization sequence Kappa:",
      KAPPA,
      "\noptimization sequence n2:",
      N2,
      "\nonset date:",
      as.character(date),
      "\nfinish date:",
      as.character(Sys.time())
    )
    class(result) <- c("drugdevelopResult", class(result))
    return(result)
    
  }
