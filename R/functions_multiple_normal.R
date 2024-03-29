############################################################################################
## R-Code for maximizing expected utility
# Utility-based optimization of phase II/III programs with two endpoints (continuous, normally distributed)
# Success: both endpoints significant in Phase III
# Two normally distributed endpoints
# endpoints: ADAS-cog, ADCS-ADL (Alzheimer's Disease)
#
# Author: Marietta Kirchner
# Date: 13.12.2016
############################################################################################


#'Density for the minimum of two normally distributed random variables
#'
#' The function `fmin()` will return the value of f(z), which is the value of the density function of the
#' minimum of two normally distributed random variables.
#'
#' Z= min(X,Y) with X ~ N(mu1,sigma1^2), Y ~ N(mu2,sigma2^2)
#'
#' f(z)=f1(z)+f2(z)
#'@param y integral variable
#'@param mu1 mean of second endpoint
#'@param mu2 mean of first endpoint
#'@param sigma1 standard deviation of first endpoint
#'@param sigma2 standard deviation of second endpoint
#'@param rho correlation between the two endpoints
#'@return  The function `fmin()` will return the value of f(z), which is the value of the density function of the
#'minimum of two normally distributed random variables.
#'@keywords internal
#' @export
fmin <- function (y, mu1, mu2, sigma1, sigma2, rho)
{
  t1 <- dnorm(y, mean = mu1, sd = sigma1)
  tt <- rho * (y - mu1) / (sigma1 * sqrt(1 - rho * rho))
  tt <- tt - (y - mu2) / (sigma2 * sqrt(1 - rho * rho))
  t1 <- t1 * pnorm(tt)
  t2 <- dnorm(y, mean = mu2, sd = sigma2)
  tt <- rho * (y - mu2) / (sigma2 * sqrt(1 - rho * rho))
  tt <- tt - (y - mu1) / (sigma1 * sqrt(1 - rho * rho))
  t2 <- t2 * pnorm(tt)
  return(t1 + t2)
}

#'Density of the bivariate normal distribution
#'@param x integral variable
#'@param y integral variable
#'@param mu1 mean of second endpoint
#'@param mu2 mean of first endpoint
#'@param sigma1 standard deviation of first endpoint
#'@param sigma2 standard deviation of second endpoint
#'@param rho correlation between the two endpoints
#'@return The Function `dbivanorm()` will return the density of a bivariate normal distribution.
#'
#'@name dbivanorm
#'@keywords internal
#' @export
dbivanorm <- function(x, y, mu1, mu2, sigma1, sigma2, rho) {
  covariancemat <-
    matrix(c(
      sigma1,
      rho * sqrt(sigma1) * sqrt(sigma2),
      rho * sqrt(sigma1) * sqrt(sigma2),
      sigma2
    ), ncol = 2)
  ff <-
    mvtnorm::dmvnorm(cbind(x, y), mean = c(mu1, mu2), sigma = covariancemat)
  return(ff)
}


#' Probability to go to phase III for multiple endpoints with normally distributed outcomes
#'
#' This function calculated the probability that we go to phase III, i.e. that results of phase II are promising enough to
#' get a successful drug development program. Successful means that both endpoints show a statistically significant positive treatment effect in phase III.
#' @param kappa threshold value for the go/no-go decision rule; vector for both endpoints
#' @param n2 total sample size for phase II; must be even number
#' @param Delta1 assumed true treatment effect given as difference in means for endpoint 1
#' @param Delta2 assumed true treatment effect given as difference in means for endpoint 2
#' @param in1 amount of information for Delta1 in terms of sample size
#' @param in2 amount of information for Delta2 in terms of sample size
#' @param sigma1 standard deviation of first endpoint
#' @param sigma2 standard deviation of second endpoint
#' @param fixed choose if true treatment effects are fixed or random, if TRUE Delta1 is used as fixed effect
#' @param rho correlation between the two endpoints
#' @param rsamp sample data set for Monte Carlo integration
#' @return The output of the function `pgo_multiple_normal()` is the probability to go to phase III.
#' 
#' @keywords internal
#' @export
pgo_multiple_normal <- function(kappa,
                                n2,
                                Delta1,
                                Delta2,
                                in1,
                                in2,
                                sigma1,
                                sigma2,
                                fixed,
                                rho,
                                rsamp) {
  Kappa <- c(kappa * sigma1, kappa * sigma2)
  var1 <- (4 * sigma1 ^ 2) / n2  #variance of effect for endpoint 1
  var2 <- (4 * sigma2 ^ 2) / n2  #variance of effect for endpoint 2
  covmat <-
    matrix(c(
      var1,
      rho * sqrt(var1) * sqrt(var2),
      rho * sqrt(var1) * sqrt(var2),
      var2
    ), ncol = 2) #covariance-Matrix of c(Delta2,Delta2)
  
  
  if (fixed) {
    return(mvtnorm::pmvnorm(
      lower = Kappa,
      upper = c(Inf, Inf),
      mean = c(Delta1, Delta2),
      sigma = covmat
    )[1])
  }
  
  else {
    nsim <- nrow(rsamp)
    pgo_vector <- vector(length = nsim)
    
    for (m in 1:nsim) {
      Delta1_prior <- rsamp[m, 1]
      Delta2_prior <- rsamp[m, 2]
      
      pgo_vector[m] <-
        mvtnorm::pmvnorm(
          lower = Kappa,
          upper = c(Inf, Inf),
          mean = c(Delta1_prior, Delta2_prior),
          sigma = covmat
        )[1]
    }
    
    return(1 / nsim * sum(pgo_vector))
    
  }
  
  #    return(integrate(function(u){
  #      sapply(u,function(u){
  #        integrate(function(v){
  #          sapply(v,function(v){
  #            (mvtnorm::pmvnorm(lower=Kappa, upper=c(Inf,Inf), mean=c(u,v),sigma=covmat)[1])*dbivanorm(u,v,Delta1,Delta2, 4/in1, 4/in2, rho)
  #          })
  #        },0,Inf )$value
  #      })
  #   },0,Inf )$value)
  
  
  
}


#' Expected sample size for phase III for multiple endpoints with normally distributed outcomes
#'
#' Given phase II results are promising enough to get the "go"-decision to go to phase III this function now calculates the expected sample size for phase III.
#' The results of this function are necessary for calculating the utility of the program, which is then in a further step maximized by the `optimal_multiple_normal()` function
#' @param kappa threshold value for the go/no-go decision rule; vector for both endpoints
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta 1-beta power for calculation of sample size for phase III
#' @param Delta1 assumed true treatment effect given as difference in means for endpoint 1
#' @param Delta2 assumed true treatment effect given as difference in means for endpoint 2
#' @param in1 amount of information for Delta1 in terms of sample size
#' @param in2 amount of information for Delta2 in terms of sample size
#' @param sigma1 standard deviation of first endpoint
#' @param sigma2 standard deviation of second endpoint
#' @param fixed choose if true treatment effects are fixed or random, if TRUE Delta1 is used as fixed effect
#' @param rho correlation between the two endpoints
#' @param rsamp sample data set for Monte Carlo integration
#' @return the output of the function Ess_multiple_normal is the expected number of participants in phase III
#' 
#' @keywords internal
#' @export
Ess_multiple_normal <-
  function(kappa,
           n2,
           alpha,
           beta,
           Delta1,
           Delta2,
           in1,
           in2,
           sigma1,
           sigma2,
           fixed,
           rho,
           rsamp) {
    Kappa <- c(kappa * sigma1, kappa * sigma2)
    Sigma <- c(sigma1, sigma2)
    r <-
      c(4 * Sigma[1] ^ 2, 4 * Sigma[2] ^ 2) #(r1,r2) known constant for endpoint i
    var1 <- r[1] / n2 #variance of effect for endpoint 1
    var2 <- r[2] / n2 #variance of effect for endpoint 2
    covmat <-
      matrix(c(
        var1,
        rho * sqrt(var1) * sqrt(var2),
        rho * sqrt(var1) * sqrt(var2),
        var2
      ), ncol = 2) #covariance-Matrix of c(true1,true2)
    
    c <- qnorm(1 - alpha) + qnorm(1 - beta)
    
    
    if (fixed)  {
      return(integrate(function(x) {
        sapply(x, function(x) {
          (4 * (qnorm(1 - alpha) + qnorm(1 - beta)) ^ 2 / x ^ 2) * fmin(x, Delta1, Delta2, var1, var2, rho)
        })
      }, kappa, Inf)$value)
    }
    
    else   {
      nsim <- nrow(rsamp)
      Ess_vector <- vector(length = nsim)
      
      for (m in 1:nsim) {
        Delta1_prior <- rsamp[m, 1]
        Delta2_prior <- rsamp[m, 2]
        
        Ess_vector[m] <- integrate(function(x) {
          sapply(x, function(x) {
            (4 * (qnorm(1 - alpha) + qnorm(1 - beta)) ^ 2 / x ^ 2) * fmin(x, Delta1_prior, Delta2_prior, var1, var2, rho)
          })
        }, kappa, Inf)$value
        
      }
      
      return(1 / nsim * sum(Ess_vector))
      
      
      #    return(integrate(function(u){ sapply(u,function(u){
      #      integrate(function(v){sapply(v,function(v){
      #        integrate(function(y){ sapply(y,function(y){
      #          integrate(function(x){ sapply(x,function(x){
      #            max((r[1]*(qnorm(1-alpha)+qnorm(1-beta))^2/x^2),(r[2]*(qnorm(1-alpha)+qnorm(1-beta))^2/y^2))*dbivanorm(x,y,u,v,var1,var2,rho)*dbivanorm(u,v,Delta1,Delta2,4/in1, 4/in2,rho)
      #          })
      #          },Kappa[1],Inf)$value
      #        })
      #        },Kappa[2],Inf)$value
      #      })
      #      }, -Inf,Inf )$value
      #    })
      #   },-Inf,Inf )$value)
      
    }
  }

#' Probability of a successful program, when going to phase III for multiple endpoint with normally distributed outcomes
#'
#' After getting the "go"-decision to go to phase III, i.e. our results of phase II are over the predefined threshold `kappa`, this function
#' calculates the probability, that our program is successfull, i.e. that both endpoints show a statistically significant positive treatment effect in phase III.
#' @param kappa threshold value for the go/no-go decision rule;
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param Delta1 assumed true treatment effect given as difference in means for endpoint 1
#' @param Delta2 assumed true treatment effect given as difference in means for endpoint 2
#' @param in1 amount of information for `Delta1` in terms of sample size
#' @param in2 amount of information for `Delta2` in terms of sample size
#' @param sigma1 standard deviation of first endpoint
#' @param sigma2 standard deviation of second endpoint
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `Delta1` is used as fixed effect
#' @param rho correlation between the two endpoints
#' @param rsamp sample data set for Monte Carlo integration
#' @return The output of the function `posp_normal()` is the probability of a successful program, when going to phase III.
#' 
#' @keywords internal
#' @export
posp_normal <-
  function(kappa,
           n2,
           alpha,
           beta,
           Delta1,
           Delta2,
           sigma1,
           sigma2,
           in1,
           in2,
           fixed,
           rho,
           rsamp) {
    Kappa <- c(kappa * sigma1, kappa * sigma2)
    Sigma <- c(sigma1, sigma2)
    r <-
      c(4 * Sigma[1] ^ 2, 4 * Sigma[2] ^ 2) #(r1,r2) known constant for endpoint i
    var1 <- r[1] / n2 #variance of effect for endpoint 1
    var2 <- r[2] / n2 #variance of effect for endpoint 2
    covmat <-
      matrix(c(
        var1,
        rho * sqrt(var1) * sqrt(var2),
        rho * sqrt(var1) * sqrt(var2),
        var2
      ), ncol = 2) #covariance-Matrix of c(true1,true2)
    covmatt3 <- matrix(c(1, rho, rho, 1), ncol = 2)
    
    c <- (qnorm(1 - alpha) + qnorm(1 - beta))^2
    
    if (fixed) {
      return(integrate(function(y) {
        sapply(y, function(y) {
          integrate(function(x) {
            sapply(x, function(x) {
              mvtnorm::pmvnorm(
                lower = c(qnorm(1 - alpha), qnorm(1 - alpha)),
                upper = c(Inf, Inf),
                mean = c(Delta1 / sqrt(r[1] / max((r[1] * c / x ^ 2), (r[2] *
                                                                         c / y ^ 2)
                )), Delta2 / sqrt(r[2] / max((r[1] * c / x ^ 2), (r[2] * c / y ^ 2)
                ))),
                sigma = covmatt3
              )[1] * dbivanorm(x, y, Delta1, Delta2, var1, var2, rho)
            })
          }, Kappa[1], Inf)$value #x
        })
      }, Kappa[2], Inf)$value)
    }
    
    else {
      nsim <- nrow(rsamp)
      
      for (m in 1:nsim) {
        Delta1_prior <- rsamp[m, 1]
        Delta2_prior <- rsamp[m, 2]
        posp_vector <- vector(length = nsim)
        
        posp_vector[m] <- integrate(function(y) {
          sapply(y, function(y) {
            integrate(function(x) {
              sapply(x, function(x) {
                mvtnorm::pmvnorm(
                  lower = c(qnorm(1 - alpha), qnorm(1 - alpha)),
                  upper = c(Inf, Inf),
                  mean = c(
                    Delta1_prior / sqrt(r[1] / max((
                      r[1] * c / x ^ 2
                    ), (
                      r[2] * c / y ^ 2
                    ))),
                    Delta2_prior / sqrt(r[2] / max((
                      r[1] * c / x ^ 2
                    ), (
                      r[2] * c / y ^ 2
                    )))
                  ),
                  sigma = covmatt3
                )[1] * dbivanorm(x, y, Delta1_prior, Delta2_prior, var1, var2, rho)
              })
            }, Kappa[1], Inf)$value #x
          })
        }, Kappa[2], Inf)$value
        
      }
      
      return(sum(posp_vector))
      
      
      #   return (integrate(function(u){ sapply(u,function(u){
      #      integrate(function(v){sapply(v,function(v){
      #        integrate(function(y){ sapply(y,function(y){
      #          integrate(function(x){ sapply(x,function(x){
      #            mvtnorm::pmvnorm(lower=c(qnorm(1-alpha),qnorm(1-alpha)),
      #                    upper=c(Inf,Inf),
      #                    mean=c(u/sqrt(r[1]/max((r[1]*c/x^2),(r[2]*c/y^2))),v/sqrt(r[2]/max((r[1]*c/x^2),(r[2]*c/y^2)))),
      #                    sigma=covmatt3)[1]*dbivanorm(x,y,u,v,var1,var2,rho)*dbivanorm(u,v,Delta1,Delta2,4/in1,4/in2,rho)
      #          })
      #          }, Kappa[1],Inf)$value #x
      #        })
      #        }, Kappa[2],Inf)$value #y
      #      })
      #      },-Inf,Inf )$value
      #    })
      #    },-Inf,Inf )$value)
    }
    
    
  }

#E(n3|GO)
#  expn3go_normal<-Ess_multiple_normal(kappa, n2, alpha, beta, Delta1, Delta2, in1, in2, sigma1, sigma2, fixed, rho)/pgo_normal(kappa, n2, Delta1, Delta2, in1, in2, sigma1, sigma2, fixed, rho)




#' Expected probability of a successful program for multiple endpoints and normally distributed outcomes
#'
#' This function calculates the probability that our drug development program is successful.
#' Successful is defined as both endpoints showing a statistically significant positive treatment effect in phase III.
#' @param kappa threshold value for the go/no-go decision rule;
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param Delta1 assumed true treatment effect given as difference in means for endpoint 1
#' @param Delta2 assumed true treatment effect given as difference in means for endpoint 2
#' @param in1 amount of information for `Delta1` in terms of sample size
#' @param in2 amount of information for `Delta2` in terms of sample size
#' @param sigma1 standard deviation of first endpoint
#' @param sigma2 standard deviation of second endpoint
#' @param step11 lower boundary for effect size for first endpoint
#' @param step12 lower boundary for effect size for second endpoint
#' @param step21 upper boundary for effect size for first endpoint
#' @param step22 upper boundary for effect size for second endpoint
#' @param fixed choose if true treatment effects are fixed or random, if TRUE then `Delta1` is used as fixed effect
#' @param rho correlation between the two endpoints
#' @param rsamp sample data set for Monte Carlo integration
#' @return The output of the function `EPsProg_multiple_normal()` is the expected probability of a successfull program, when going to phase III.
#' 
#' @keywords internal
#' @export
EPsProg_multiple_normal <-
  function(kappa,
           n2,
           alpha,
           beta,
           Delta1,
           Delta2,
           sigma1,
           sigma2,
           step11,
           step12,
           step21,
           step22,
           in1,
           in2,
           fixed,
           rho,
           rsamp) {
    Kappa <- c(kappa * sigma1, kappa * sigma2)
    Sigma <- c(sigma1, sigma2)
    r <-
      c(4 * Sigma[1] ^ 2, 4 * Sigma[2] ^ 2) #(r1,r2) known constant for endpoint i
    var1 <- r[1] / n2 #variance of effect for endpoint 1
    var2 <- r[2] / n2 #variance of effect for endpoint 2
    covmat <-
      matrix(c(
        var1,
        rho * sqrt(var1) * sqrt(var2),
        rho * sqrt(var1) * sqrt(var2),
        var2
      ), ncol = 2) #covariance-Matrix of c(true1,true2)
    covmatt3 <- matrix(c(1, rho, rho, 1), ncol = 2)
    
    
    c <- (qnorm(1 - alpha) + qnorm(1 - beta)) ^ 2
    
    if (fixed) {
      return(integrate(function(y) {
        sapply(y, function(y) {
          integrate(function(x) {
            sapply(x, function(x) {
              mvtnorm::pmvnorm(
                lower = c(
                  qnorm(1 - alpha) + step11 * sqrt(max((
                    r[1] * c / x ^ 2
                  ), (
                    r[2] * c / y ^ 2
                  ))) / 2,
                  qnorm(1 - alpha) + step12 * sqrt(max((
                    r[1] * c / x ^ 2
                  ), (
                    r[2] * c / y ^ 2
                  ))) / 2
                ),
                upper = c(
                  qnorm(1 - alpha) + step21 * sqrt(max((
                    r[1] * c / x ^ 2
                  ), (
                    r[2] * c / y ^ 2
                  ))) / 2,
                  qnorm(1 - alpha) + step22 * sqrt(max((
                    r[1] * c / x ^ 2
                  ), (
                    r[2] * c / y ^ 2
                  ))) / 2
                ),
                mean = c(Delta1 / sqrt(r[1] / max((r[1] * c / x ^ 2), (r[2] *
                                                                         c / y ^ 2)
                )), Delta2 / sqrt(r[2] / max((r[1] * c / x ^ 2), (r[2] * c / y ^ 2)
                ))),
                sigma = covmatt3
              )[1] * dbivanorm(x, y, Delta1, Delta2, var1, var2, rho)
            })
          }, Kappa[1], Inf)$value
        })
      }, Kappa[2], Inf)$value)
    }
    
    
    else {
      nsim <- nrow(rsamp)
      
      for (m in 1:nsim) {
        Delta1_prior <- rsamp[m, 1]
        Delta2_prior <- rsamp[m, 2]
        EPsProg_vector <- vector(length = nsim)
        
        EPsProg_vector[m] <- integrate(function(y) {
          sapply(y, function(y) {
            integrate(function(x) {
              sapply(x, function(x) {
                mvtnorm::pmvnorm(
                  lower = c(
                    qnorm(1 - alpha) + step11 * sqrt(max((
                      r[1] * c / x ^ 2
                    ), (
                      r[2] * c / y ^ 2
                    ))) / 2,
                    qnorm(1 - alpha) + step12 * sqrt(max((
                      r[1] * c / x ^ 2
                    ), (
                      r[2] * c / y ^ 2
                    ))) / 2
                  ),
                  upper = c(
                    qnorm(1 - alpha) + step21 * sqrt(max((
                      r[1] * c / x ^ 2
                    ), (
                      r[2] * c / y ^ 2
                    ))) / 2,
                    qnorm(1 - alpha) + step22 * sqrt(max((
                      r[1] * c / x ^ 2
                    ), (
                      r[2] * c / y ^ 2
                    ))) / 2
                  ),
                  mean = c(
                    Delta1_prior / sqrt(r[1] / max((
                      r[1] * c / x ^ 2
                    ), (
                      r[2] * c / y ^ 2
                    ))),
                    Delta2_prior / sqrt(r[2] / max((
                      r[1] * c / x ^ 2
                    ), (
                      r[2] * c / y ^ 2
                    )))
                  ),
                  sigma = covmatt3
                )[1] * dbivanorm(x, y, Delta1_prior, Delta2_prior, var1, var2, rho)
              })
            }, Kappa[1], Inf)$value
          })
        }, Kappa[2], Inf)$value
        
      }
      
      return(sum(EPsProg_vector))
      
      #    return(integrate(function(u){ sapply(u,function(u){
      #      integrate(function(v){sapply(v,function(v){
      #        integrate(function(y){ sapply(y,function(y){
      #          integrate(function(x){ sapply(x,function(x){
      #            mvtnorm::pmvnorm(lower=c(qnorm(1-alpha)+step11*sqrt(max((r[1]*c/x^2),(r[2]*c/y^2)))/2,
      #                            qnorm(1-alpha)+step12*sqrt(max((r[1]*c/x^2),(r[2]*c/y^2)))/2),
      #                    upper=c(qnorm(1-alpha)+step21*sqrt(max((r[1]*c/x^2),(r[2]*c/y^2)))/2,
      #                            qnorm(1-alpha)+step22*sqrt(max((r[1]*c/x^2),(r[2]*c/y^2)))/2),
      #                    mean=c(u/sqrt(r[1]/max((r[1]*c/x^2),(r[2]*c/y^2))),v/sqrt(r[2]/max((r[1]*c/x^2),(r[2]*c/y^2)))),
      #                    sigma=covmatt3)[1]*dbivanorm(x,y,u,v,var1,var2,rho)*dbivanorm(u,v,Delta1,Delta2,4/in1,4/in2,rho)
      #          })
      #          },Kappa[1],Inf)$value
      #        })
      #        },Kappa[2],Inf)$value
      #      })
      #      },-Inf,Inf )$value
      #    })
      #    },-Inf,Inf )$value)
      
    }
  }

#' Utility function for multiple endpoints with normally distributed outcomes.
#'
#' The utility function calculates the expected utility of our drug development program and is given as gains minus costs and depends on the parameters and the expected probability of a successful program.
#' The utility is in a further step maximized by the `optimal_multiple_normal()` function.
#' @param kappa threshold value for the go/no-go decision rule; vector for both endpoints
#' @param n2 total sample size for phase II; must be even number
#' @param alpha significance level
#' @param beta `1-beta` power for calculation of sample size for phase III
#' @param Delta1 assumed true treatment effect given as difference in means for endpoint 1
#' @param Delta2 assumed true treatment effect given as difference in means for endpoint 2
#' @param in1 amount of information for `Delta1` in terms of sample size
#' @param in2 amount of information for `Delta2` in terms of sample size
#' @param sigma1 standard deviation of first endpoint
#' @param sigma2 standard deviation of second endpoint
#' @param c2 variable per-patient cost for phase II
#' @param c3 variable per-patient cost for phase III
#' @param c02 fixed cost for phase II
#' @param c03 fixed cost for phase III
#' @param K constraint on the costs of the program, default: Inf, e.g. no constraint
#' @param N constraint on the total expected sample size of the program, default: Inf, e.g. no constraint
#' @param S constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint
#' @param steps1 lower boundary for effect size category `"small"` in HR scale, default: 1
#' @param stepm1 lower boundary for effect size category `"medium"` in HR scale = upper boundary for effect size category `"small"` in HR scale, default: 0.95
#' @param stepl1 lower boundary for effect size category `"large"` in HR scale = upper boundary for effect size category `"medium"` in HR scale, default: 0.85
#' @param b1 expected gain for effect size category `"small"`
#' @param b2 expected gain for effect size category `"medium"`
#' @param b3 expected gain for effect size category `"large"`
#' @param fixed choose if true treatment effects are fixed or random, if TRUE `Delta1` is used as fixed effect
#' @param rho correlation between the two endpoints
#' @param relaxed relaxed or strict decision rule
#' @return The output of the function `utility_multiple_normal()` is the expected utility of the program.
#' 
#' @keywords internal
#' @export
utility_multiple_normal <- function(kappa,
                                    n2,
                                    alpha,
                                    beta,
                                    Delta1,
                                    Delta2,
                                    in1,
                                    in2,
                                    sigma1,
                                    sigma2,
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
                                    b3,
                                    fixed,
                                    rho,
                                    relaxed,
                                    rsamp) {

  n3 <-
    Ess_multiple_normal(
      kappa = kappa,
      n2 = n2,
      alpha = alpha ,
      beta = beta,
      Delta1 = Delta1,
      Delta2 = Delta2,
      in1 = in1,
      in2 = in2,
      sigma1 = sigma1,
      sigma2 = sigma2,
      fixed = fixed,
      rho = rho,
      rsamp = rsamp
    )
  n3 <- ceiling(n3)
  
  
  POSP = posp_normal(
    kappa = kappa,
    n2 = n2,
    alpha = alpha,
    beta = beta,
    Delta1 = Delta1,
    Delta2 = Delta2,
    in1,
    in2,
    sigma1 = sigma1,
    sigma2 = sigma2,
    fixed = fixed,
    rho = rho,
    rsamp = rsamp
  )
  
  
  if (n2 + n3 > N) {
    return(c(-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999))
    
  } else{
    pgo   = pgo_multiple_normal(
      kappa = kappa,
      n2 = n2,
      Delta1 = Delta1,
      Delta2 = Delta2,
      in1 = in1,
      in2 = in2,
      sigma1 = sigma1,
      sigma2 = sigma2,
      fixed = fixed,
      rho = rho,
      rsamp = rsamp
    )
    
    K2    <-  c02 + c2 * n2  #cost phase II
    K3    <-  c03 * pgo + c3 * n3  #cost phase III
    
    if (K2 + K3 > K) {
      return(c(-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999))
      #output: expected utility Eud, En3, EsP, Epgo
      
    } else{
      # probability of a successful program; small, medium, large effect size
      prob11 = EPsProg_multiple_normal(
        kappa = kappa ,
        n2 = n2,
        alpha = alpha,
        beta = beta,
        Delta1 = Delta1,
        Delta2 = Delta2,
        sigma1 = sigma1,
        sigma2 = sigma2,
        step11 = steps1 ,
        step12 = stepm1  ,
        step21 = stepm1,
        step22 = Inf,
        in1 = in1,
        in2 = in2,
        fixed = fixed,
        rho = rho,
        rsamp = rsamp
      )
      prob21 = EPsProg_multiple_normal(
        kappa = kappa ,
        n2 = n2,
        alpha = alpha,
        beta = beta,
        Delta1 = Delta1,
        Delta2 = Delta2,
        sigma1 = sigma1,
        sigma2 = sigma2,
        step11 = stepm1,
        step12 = steps1 ,
        step21 = Inf ,
        step22 = stepm1,
        in1 = in1,
        in2 = in2,
        fixed = fixed,
        rho = rho,
        rsamp = rsamp
      )
      prob31 = EPsProg_multiple_normal(
        kappa = kappa ,
        n2 = n2,
        alpha = alpha,
        beta = beta,
        Delta1 = Delta1,
        Delta2 = Delta2,
        sigma1 = sigma1,
        sigma2 = sigma2,
        step11 = steps1,
        step12 = steps1,
        step21 = stepm1,
        step22 = stepm1,
        in1 = in1,
        in2 = in2,
        fixed = fixed,
        rho = rho,
        rsamp = rsamp
      )
      
      proba1 <- prob11 + prob21 + prob31
      proba3 = EPsProg_multiple_normal(
        kappa = kappa ,
        n2 = n2,
        alpha = alpha,
        beta = beta,
        Delta1 = Delta1,
        Delta2 = Delta2,
        sigma1 = sigma1,
        sigma2 = sigma2,
        step11 = stepl1,
        step12 = stepl1,
        step21 = Inf,
        step22 = Inf,
        in1 = in1,
        in2 = in2,
        fixed = fixed,
        rho = rho,
        rsamp = rsamp
      )
      proba2 <- POSP - proba1 - proba3
      
      prob13 = EPsProg_multiple_normal(
        kappa = kappa ,
        n2 = n2,
        alpha = alpha,
        beta = beta,
        Delta1 = Delta1,
        Delta2 = Delta2,
        sigma1 = sigma1,
        sigma2 = sigma2,
        step11 = stepl1,
        step12 = steps1,
        step21 = Inf,
        step22 = stepl1,
        in1 = in1,
        in2 = in2,
        fixed = fixed,
        rho = rho,
        rsamp = rsamp
      )
      prob23 = EPsProg_multiple_normal(
        kappa = kappa ,
        n2 = n2,
        alpha = alpha,
        beta = beta,
        Delta1 = Delta1,
        Delta2 = Delta2,
        sigma1 = sigma1,
        sigma2 = sigma2,
        step11 = steps1,
        step12 = stepl1,
        step21 = stepl1,
        step22 = Inf,
        in1 = in1,
        in2 = in2,
        fixed = fixed ,
        rho = rho,
        rsamp = rsamp
      )
      prob33 = EPsProg_multiple_normal(
        kappa = kappa ,
        n2 = n2,
        alpha = alpha,
        beta = beta,
        Delta1 = Delta1,
        Delta2 = Delta2,
        sigma1 = sigma1,
        sigma2 = sigma2,
        step11 = stepl1,
        step12 = stepl1,
        step21 = Inf,
        step22 = Inf,
        in1 = in1,
        in2 = in2,
        fixed = fixed,
        rho = rho,
        rsamp = rsamp
      )
      
      probb3 <- prob13 + prob23 + prob33
      probb1 = EPsProg_multiple_normal(
        kappa = kappa ,
        n2 = n2,
        alpha = alpha,
        beta = beta,
        Delta1 = Delta1,
        Delta2 = Delta2,
        sigma1 = sigma1,
        sigma2 = sigma2,
        step11 = steps1,
        step12 = steps1,
        step21 = stepm1,
        step22 = stepm1,
        in1 = in1,
        in2 = in2,
        fixed = fixed,
        rho = rho,
        rsamp = rsamp
      )
      probb2 <- POSP - probb1 - probb3
      
      if (relaxed == "TRUE") {
        prob1 <- probb1
        prob2 <- probb2
        prob3 <- probb3
      } else {
        prob1 <- proba1
        prob2 <- proba2
        prob3 <- proba3
      }
      
      
      
      SP    = prob1 + prob2 + prob3                         # probability of a successful program
      
      
      if (SP < S) {
        return(c(-9999,-9999,-9999,-9999,-9999,-9999,-9999,-9999))
        
      } else{
        G     = b1 * prob1 + b2 * prob2 + b3 * prob3                      # gain
        
        EU    = -K2 - K3 + G                                            # total expected utility
        
        SP2 = prob2
        SP3 = prob3
        
        return(c(EU, n3, SP, pgo, SP2, SP3, K2, K3))
        
      }
    }
  }
  
}

#' Generate sample for Monte Carlo integration in the multiple setting
#' 
#' @param Delta1 assumed true treatment effect given as difference in means for endpoint 1
#' @param Delta2 assumed true treatment effect given as difference in means for endpoint 2
#' @param in1 amount of information for `Delta1` in terms of sample size
#' @param in2 amount of information for `Delta2` in terms of sample size
#' @param rho correlation between the two endpoints
#' 
#' @return a randomly generated data frame
#' @keywords internal
get_sample_multiple_normal <-
  function(Delta1, Delta2, in1, in2, rho) {
    mu_prior <- c(Delta1, Delta2) # true treatment effect theta
    Sigma_prior <- matrix(c(
      4 / in1,
      rho * sqrt(4 / in1) * sqrt(4 / in2),
      rho * sqrt(4 / in1) * sqrt(4 / in2),
      4 / in2
    ), ncol = 2)
    nsim <- 100
    rsamp <-
      MASS::mvrnorm(
        n = nsim,
        mu_prior,
        Sigma_prior,
        tol = 1e-6,
        empirical = FALSE,
        EISPACK = FALSE
      )
    return(rsamp)
  }
