# prior distribution for theta

#' Prior distribution for time-to-event outcomes
#'
#' If we do not assume the treatment effects to be fixed, i.e. `fixed = FALSE`,
#' the function `prior_tte` allows us to model the treatment effect following a prior distribution.
#' For more details concerning the definition of a prior distribution, see the \href{https://sterniii3.github.io/drugdevelopR/articles/Introduction-to-drugdevelopR.html}{vignette on priors}
#' as well as the \href{https://web.imbi.uni-heidelberg.de/prior/}{Shiny app}.
#'
#' @param x integration variable
#' @param w weight for mixture prior distribution
#' @param hr1 first assumed true treatment effect on HR scale for prior distribution
#' @param hr2 second assumed true treatment effect on HR scale for prior distribution
#' @param id1 amount of information for `hr1` in terms of number of events
#' @param id2 amount of information for `hr2` in terms of number of events
#' @return The output of the functions `Epgo_tte()` is the expected number of participants in phase III with conservative decision rule and sample size calculation.
#' @examples res <- prior_tte(x = 0.5, w = 0.5, hr1 = 0.69, hr2 = 0.88, id1 = 240, id2 = 420)
#' @export
#' @keywords internal
prior_tte<-function(x, w, hr1, hr2, id1, id2){
    w * dnorm(x, -log(hr1), sqrt(4/id1)) + 
    (1 - w) * dnorm(x, -log(hr2), sqrt(4/id2))
}

# 10000 realizations of the prior distribution
box_tte<-function(w, hr1, hr2, id1, id2){
  w * rnorm(1000000, -log(hr1),sqrt(4/id1)) + 
    (1 - w) * rnorm(1000000, -log(hr2), sqrt(4/id2))
}

# expected probability to go to phase III

#' Expected probability to go to phase III for time-to-event outcomes
#' 
#' @param HRgo threshold value for the go/no-go decision rule
#' @param d2 total number of events for phase II; must be even number
#' @param w weight for mixture prior distribution
#' @param hr1 first assumed true treatment effect on HR scale for prior distribution
#' @param hr2 second assumed true treatment effect on HR scale for prior distribution
#' @param id1 amount of information for `hr1` in terms of number of events
#' @param id2 amount of information for `hr2` in terms of number of events
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `hr1` is used as fixed effect
#' @return The output of the functions `Epgo_tte()` is the expected probability to go to phase III.
#' @examples res <- Epgo_tte(HRgo = 0.8, d2 = 50,  
#'                                 w = 0.3, hr1 = 0.69, hr2 = 0.81, 
#'                                 id1 = 280, id2 = 420, fixed = FALSE)
#' @export
#' @keywords internal
Epgo_tte <-  function(HRgo, d2, w, hr1, hr2, id1, id2, fixed){
  if(!fixed){
    return(  
      integrate(function(x){
        sapply(x, function(x){
          pnorm((log(HRgo) + x)/sqrt(4/d2))*
            prior_tte(x, w, hr1, hr2, id1, id2)
        })
      },  - Inf, Inf)$value
    ) 
  }else{
    return(
      pnorm((log(HRgo) - log(hr1))/sqrt(4/d2))
    )
  }
}

# expected number of events for phase III 
# in before phase II perspective
#' Expected sample size for phase III for time-to-event outcomes
#' 
#' @param HRgo threshold value for the go/no-go decision rule
#' @param d2 total events for phase II; must be even number
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param w weight for mixture prior distribution
#' @param hr1 first assumed true treatment effect on HR scale for prior distribution
#' @param hr2 second assumed true treatment effect on HR scale for prior distribution
#' @param id1 amount of information for `hr1` in terms of number of events
#' @param id2 amount of information for `hr2` in terms of number of events
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `hr1` is used as fixed effect
#' @return The output of the the functions `Ed3_tte` is the expected number of events in phase III. 
#' @examples res <-  Ed3_tte(HRgo = 0.8, d2 = 50,
#'                         alpha = 0.025, beta = 0.1, w = 0.3, 
#'                         hr1 =  0.69, hr2 = 0.81, 
#'                         id1 = 280, id2 = 420, fixed = FALSE)
#' @export
#' @keywords internal
Ed3_tte <-  function(HRgo, d2, alpha, beta, 
                     w, hr1, hr2, id1, id2, fixed){
  if(!fixed){
    return(  
      integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            ((4*(qnorm(1-alpha)+qnorm(1-beta))^2)/(y^2))*
              dnorm(y,
                  mean = x,
                  sd = sqrt(4/d2))*
            prior_tte(x, w, hr1, hr2, id1, id2)
          }, -log(HRgo), Inf)$value
        })
      },  - Inf, Inf)$value
    )
  }else{
    return(  
      integrate(function(y){
        ((4*(qnorm(1-alpha)+qnorm(1-beta))^2)/(y^2))*
            dnorm(y,
                  mean = -log(hr1),
                  sd = sqrt(4/d2)) 
      }, -log(HRgo), Inf)$value
    )
  }
}

# expected probability of a successful program

#' Expected probability of a successful program for time-to-event outcomes
#' 
#' @param HRgo threshold value for the go/no-go decision rule
#' @param d2 total events for phase II; must be even number
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param step1 lower boundary for effect size
#' @param step2 upper boundary for effect size
#' @param w weight for mixture prior distribution
#' @param hr1 first assumed true treatment effect on HR scale for prior distribution
#' @param hr2 second assumed true treatment effect on HR scale for prior distribution
#' @param id1 amount of information for `hr1` in terms of number of events
#' @param id2 amount of information for `hr2` in terms of number of events
#' @param gamma difference in treatment effect due to different population structures in phase II and III
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `hr1` is used as fixed effect
#' @return  The output of the functions `EPsProg_tte()` is the expected probability of a successful program.
#' @examples res <- EPsProg_tte(HRgo = 0.8, d2 = 50, 
#'                            alpha = 0.025, beta = 0.1, 
#'                            step1 = 1, step2 = 0.95, 
#'                            w = 0.3, hr1 = 0.69, hr2 = 0.81,
#'                            id1 = 280, id2 = 420,
#'                            gamma = 0, fixed = FALSE)
#' @export
#' @keywords internal

EPsProg_tte <-  function(HRgo, d2, alpha, beta, 
                         step1, step2, 
                         w, hr1, hr2, id1, id2, 
                         gamma, fixed){

  c = (qnorm(1 - alpha) + qnorm(1 - beta))^2

  if(!fixed){
    return(
      integrate(function(x){
        sapply(x, function(x){
          integrate(function(y){
            (pnorm(qnorm(1-alpha)-log(step2)/(sqrt(y^2/c)),
                    mean = (x+gamma)/(sqrt(y^2/c)),
                    sd = 1) -
                pnorm(qnorm(1-alpha)-log(step1)/(sqrt(y^2/c)),
                      mean = (x+gamma)/(sqrt(y^2/c)),
                      sd = 1) )*
                  dnorm(y,
                      mean = x,
                      sd = sqrt(4/d2))*
                    prior_tte(x, w, hr1, hr2, id1, id2)
          }, -log(HRgo), Inf)$value
        })
      },  - Inf, Inf)$value
    )
  }else{
    return(
      integrate(function(y){
        (pnorm(qnorm(1-alpha)-log(step2)/(sqrt(y^2/c)),
                mean = (-log(hr1)+gamma)/(sqrt(y^2/c)),
                sd = 1) -
            pnorm(qnorm(1-alpha)-log(step1)/(sqrt(y^2/c)),
                  mean = (-log(hr1)+gamma)/(sqrt(y^2/c)),
                  sd = 1))*
              dnorm(y,
                mean = -log(hr1),
                sd = sqrt(4/d2)) 
      },  - log(HRgo), Inf)$value
    )
  }
}

# utility function

#' Utility function for time-to-event outcomes.
#' 
#' The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters and the expected probability of a successful program. 
#' The utility is in a further step maximized by the `optimal_tte()` function.
#' @param HRgo threshold value for the go/no-go decision rule
#' @param d2 total events for phase II; must be even number
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param xi2 event rate for phase II
#' @param xi3 event rate for phase III
#' @param w weight for mixture prior distribution
#' @param hr1 first assumed true treatment effect on HR scale for prior distribution
#' @param hr2 second assumed true treatment effect on HR scale for prior distribution
#' @param id1 amount of information for `hr1` in terms of number of events
#' @param id2 amount of information for `hr2` in terms of number of events
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param steps1 lower boundary for effect size category `"small"` in RR scale, default: 1
#' @param stepm1 lower boundary for effect size category `"medium"` in RR scale = upper boundary for effect size category "small" in RR scale, default: 0.95
#' @param stepl1 lower boundary for effect size category `"large"` in RR scale = upper boundary for effect size category "medium" in RR scale, default: 0.85
#' @param b1 expected gain for effect size category `"small"`
#' @param b2 expected gain for effect size category `"medium"`
#' @param b3 expected gain for effect size category `"large"`
#' @param gamma difference in treatment effect due to different population structures in phase II and III
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `hr1` is used as fixed effect
#' @return The output of the functions `utility_tte()` is the expected utility of the program.
#' @examples res <- utility_tte(d2 = 50, HRgo = 0.8, w = 0.3, 
#'                                  hr1 =  0.69, hr2 = 0.81, 
#'                                  id1 = 280, id2 = 420, xi2 = 0.7, xi3 = 0.7,
#'                                  alpha = 0.025, beta = 0.1,
#'                                  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,
#'                                  K = Inf, N = Inf, S = -Inf,
#'                                  steps1 = 1, stepm1 = 0.95, stepl1 = 0.85,
#'                                  b1 = 1000, b2 = 2000, b3 = 3000, 
#'                                  gamma = 0, fixed = TRUE)
#' @export
#' @keywords internal
utility_tte <-  function(d2, HRgo, w, hr1, hr2, id1, id2,
                         alpha, beta, xi2, xi3,
                         c2, c3, c02, c03, 
                         K, N, S,
                         steps1, stepm1, stepl1,
                         b1, b2, b3,
                         gamma, fixed){

  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 0
  
  d3  <- Ed3_tte(HRgo = HRgo, d2 = d2, alpha = alpha, 
                 beta = beta, w = w, hr1 = hr1, hr2 = hr2,
                 id1 = id1, id2 = id2, fixed = fixed)

  # sample size is rounded up to next even natural number
  n2  <- ceiling(d2*(1/xi2))
  if(round(n2/2) != n2 / 2) {n2 <- n2 + 1}

  n3  <- ceiling(d3 * (1/xi3))
  if(round(n3/2) != n3 / 2) {n3 <- n3 + 1}

  # expected number of events is rounded to natural number
  d3  <- ceiling(d3)
  
  if(n2+n3>N){
    return(c(-9999, -9999, -9999, -9999, -9999, -9999, 
             -9999, -9999, -9999, -9999, -9999))
  }else{
    pg    <-  Epgo_tte(HRgo = HRgo, d2 = d2, 
                       w = w, hr1 = hr1, hr2 = hr2, 
                       id1 = id1, id2 = id2,
                       fixed = fixed)
    
    K2    <-  c02 + c2 * n2         # cost phase II
    K3    <-  c03 * pg + c3 * n3    # cost phase III
    
    if(K2+K3>K){
      return(c(-9999, -9999, -9999, -9999, -9999, -9999, 
               -9999, -9999, -9999, -9999, -9999))
    }else{
      # probability of a successful program:
      # small, medium and large effect size
      prob1 <- EPsProg_tte(HRgo = HRgo, d2 = d2,
                           alpha = alpha, beta = beta,
                           step1 = steps1, step2 =  steps2,
                           w = w, hr1 = hr1, hr2 = hr2, 
                           id1 = id1, id2 = id2, 
                           gamma = gamma, fixed = fixed)
      prob2 <- EPsProg_tte(HRgo = HRgo, d2 = d2, 
                           alpha = alpha, beta = beta,
                           step1 =  stepm1, step2 =  stepm2,
                           w = w, hr1 = hr1, hr2 = hr2, 
                           id1 = id1, id2 = id2, 
                           gamma = gamma, fixed = fixed)
      prob3 <- EPsProg_tte(HRgo = HRgo, d2 = d2, 
                           alpha = alpha, beta = beta,
                           step1 =  stepl1, step2 = stepl2,
                           w = w, hr1 = hr1, hr2 = hr2, 
                           id1 = id1, id2 = id2, 
                           gamma = gamma, fixed = fixed)
      
      SP    <-  prob1 + prob2 + prob3
      
      if(SP<S){
        return(c(-9999, -9999, -9999, -9999, -9999, 
                 -9999, -9999, -9999, -9999, -9999, -9999))
      }else{
        G     <-  b1 * prob1 + b2 * prob2 + b3 * prob3 
        EU    <-  - K2 - K3 + G

        return(
          c(EU, d3, SP, pg, K2, K3, prob1, prob2, prob3, n2, n3)
        )
      } 
    } 
  }
}

#################
# skip phase II #
#################

# number of events for phase III based on median_prior
#' Expected probability to go to phase III for time-to-event outcomes
#' 
#' If choosing `skipII = TRUE`, the program calculates the expected utility for the case when phase
#' II is skipped and compares it to the situation when phase II is not skipped.
#'  This function calculates the expected sample size for phase III for time-to-event outcomes using a median prior. 
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III 
#' @param median_prior the median_prior is given as -log(hr1), the assumed true treatment effect
#' @return The output of the functions `d3_skipII_tte()` is the expected number of events in phase III when skipping phase II.
#' @examples res <- d3_skipII_tte(alpha = 0.05, beta = 0.1, median_prior = 0.35)
#' @export
#' @keywords internal
d3_skipII_tte <-function(alpha, beta, median_prior){
  return(
    (4*(qnorm(1-alpha)+qnorm(1-beta))^2)/(median_prior^2)
  )
}

# expected probability of a successful program 
# based on median_prior
EPsProg_skipII_tte <-function(alpha, beta, step1, step2, 
                              median_prior, w, hr1, hr2, 
                              id1, id2, gamma, fixed){

  c=(qnorm(1-alpha)+qnorm(1-beta))^2

  if(!fixed){
    return(
      integrate(function(x){
        sapply(x,function(x){
          (pnorm(qnorm(1-alpha)-
                   log(step2)/(sqrt(median_prior^2/c)),
                  mean=(x+gamma)/(sqrt(median_prior^2/c)),
                  sd=1)-
              pnorm(qnorm(1-alpha)-
                      log(step1)/(sqrt(median_prior^2/c)),
                    mean=(x+gamma)/(sqrt(median_prior^2/c)),
                    sd=1))*
                prior_tte(x, w, hr1, hr2, id1, id2)
        })
      }, -Inf, Inf)$value
    )  
  }else{
    return(
      pnorm(qnorm(1-alpha)- 
              log(step2)/(sqrt(median_prior^2/c)),
            mean=(-log(hr1)+gamma)/(sqrt(median_prior^2/c)),
            sd=1)-
        pnorm(qnorm(1-alpha)- 
                log(step1)/(sqrt(median_prior^2/c)),
              mean=(-log(hr1)+gamma)/(sqrt(median_prior^2/c)),
              sd=1) 
    )
  }
}

#utility function
utility_skipII_tte <-function(alpha, beta, xi3, c03, c3, 
                              b1, b2, b3, median_prior,
                              K, N, S, 
                              steps1, stepm1, stepl1,
                              w, hr1, hr2, id1, id2, 
                              gamma, fixed){

  steps2 <- stepm1
  stepm2 <- stepl1
  stepl2 <- 0
  
  d3  <- d3_skipII_tte(alpha = alpha, beta = beta, 
                       median_prior = median_prior)

  n3  <- ceiling(d3*(1/xi3))
  if(round(n3/2) != n3 / 2) {n3 = n3 + 1}

  d3  <- ceiling(d3)
  
  if(n3>N){
    return(c(-9999, -9999, -9999, -9999,
             -9999, -9999, -9999, -9999))
  }else{
    K2  <- 0 
    K3  <- c03 + c3*n3
    if(K2+K3>K){
      return(c(-9999, -9999, -9999, -9999,
               -9999, -9999, -9999, -9999))
    }else{

      # probability of a successful program:
      # small, medium, large effect size
      prob1 <- EPsProg_skipII_tte(alpha = alpha, beta = beta, 
                                  step1 = steps1, 
                                  step2 = steps2,
                                  median_prior = median_prior, 
                                  w = w, hr1 = hr1, hr2 = hr2, 
                                  id1 = id1, id2 = id2,
                                  gamma = gamma, fixed = fixed)
      prob2 <- EPsProg_skipII_tte(alpha = alpha, beta = beta, 
                                  step1 = stepm1, 
                                  step2 = stepm2,
                                  median_prior = median_prior, 
                                  w = w, hr1 = hr1, hr2 = hr2, 
                                  id1 = id1, id2 = id2,
                                  gamma = gamma, fixed = fixed)
      prob3 <- EPsProg_skipII_tte(alpha = alpha, beta = beta, 
                                  step1 = stepl1, 
                                  step2 = stepl2,
                                  median_prior = median_prior, 
                                  w = w, hr1 = hr1, hr2 = hr2, 
                                  id1 = id1, id2 = id2,
                                  gamma = gamma, fixed = fixed)
      
      SP    <-  prob1 + prob2 + prob3
      
      if(SP<S){
        return(c(-9999, -9999, -9999, -9999,
                 -9999, -9999, -9999, -9999))
        
      }else{
        G     <- b1 * prob1 + b2 * prob2 + b3 * prob3 
        EU    <-  - K2 - K3 + G
        
        return(
          c(EU, d3, n3, SP, K3, prob1, prob2, prob3)
        )
      }
    }
  }
}

