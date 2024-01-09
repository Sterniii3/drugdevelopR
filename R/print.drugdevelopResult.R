#' Printing a drugdevelopResult Object
#' 
#' Displays details about the optimal results from a \code{drugdevelopResult} object.
#'
#' @param x Data frame of class \code{drugdevelopResult}.
#' @param sequence logical, print optimization sequence (default = `FALSE`)?
#' @param ... Further arguments.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Activate progress bar (optional)
#' \dontrun{progressr::handlers(global = TRUE)}
#' # Optimize
#' res <- optimal_normal(w=0.3,                       
#'   Delta1 = 0.375, Delta2 = 0.625, in1=300, in2=600,  
#'   a = 0.25, b = 0.75,
#'   n2min = 20, n2max = 100, stepn2 = 4,               
#'   kappamin = 0.02, kappamax = 0.2, stepkappa = 0.02, 
#'   alpha = 0.025, beta = 0.1,                       
#'   c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,   
#'   K = Inf, N = Inf, S = -Inf,                        
#'   steps1 = 0,                                      
#'   stepm1 = 0.5,                                      
#'   stepl1 = 0.8,                                      
#'   b1 = 3000, b2 = 8000, b3 = 10000,                  
#'   gamma = 0,                                         
#'   fixed = FALSE,                                    
#'   skipII = FALSE,                                    
#'   num_cl = 1)
#' # Print results
#' print(res)                                  
#' }
#' @keywords internal
print.drugdevelopResult <- function(x, sequence = FALSE, ...) {
  if(nrow(x) > 1){
    # Skip phase II option
    if("skipII" %in% names(x)){
      cat("Optimization result including phase II:\n")
      print_drugdevelopResult_helper(x[which(!x$skipII),], ...)
      cat("\nResult when skipping phase II:\n")
      print_drugdevelopResult_helper(x[which(x$skipII),], ...)
      
      cat("\n")
      if(x[which(x$skipII), "u"]>x[which(!x$skipII), "u"]){
        cat("Skipping phase II is the optimal option with respect to the maximal expected utility.\n")
      }else{
        cat("Skipping phase II is NOT the optimal option with respect to the maximal expected utility.\n")
      }
    }
    # Bias adjustment
    if("Method" %in% names(x)){
      if("multipl." %in% x$Method){
        cat("Optimization result with multiplicative adjustment of the treatment effect:\n")
        print_drugdevelopResult_helper(x[which(x$Method == "multipl."),], ...)
        cat("\n")
      }
      if("add." %in% x$Method){
        cat("Optimization result with additive adjustment of the treatment effect:\n")
        print_drugdevelopResult_helper(x[which(x$Method == "add."),], ...)
        cat("\n")
      }
      if("multipl2." %in% x$Method){
        cat("Optimization result with multiplicative adjustment of treatment effect and decision rule:\n")
        print_drugdevelopResult_helper(x[which(x$Method == "multipl2."),], ...)
        cat("\n")
      }
      if("add2." %in% x$Method){
        cat("Optimization result with additive adjustment of treatment effect and decision rule:\n")
        print_drugdevelopResult_helper(x[which(x$Method == "add2."),], ...)
        cat("\n")
      }
    }
    # Multi-trial
    if("Case" %in% names(x) & "Strategy" %in% names(x)){
      for(i in (1:nrow(x))){
        cat("Optimization result with ", x[i,"Case"], " significant trial(s) needed,",
            " strategy ", x[i,"Strategy"], ":", "\n", sep = "")
        print_drugdevelopResult_helper(x[i, ], ...)
        cat("\n")
      }
    }
    # Multi-arm
    if(!("Case" %in% names(x)) & "Strategy" %in% names(x)){
      for(i in (1:nrow(x))){
        if(x[i,"Strategy"] == 1){
          cat("Optimization result where only the most promising candidate continues:\n", sep = "")
        }
        if(x[i,"Strategy"] == 2){
          cat("Optimization result where all promising candidates continue:\n", sep = "")
        }
        print_drugdevelopResult_helper(x[i, ], ...)
        cat("\n")
      }
    }
  }
  # General option
  else{
    cat("Optimization result:\n")
    print_drugdevelopResult_helper(x, ...)
  }
  if(sequence){
    cat(attr(x, "comment"))
  }
}

#' Helper function for printing a drugdevelopResult Object
#'
#' @param x Data frame
#' @param ... Further arguments.
#'
#' @return No return value, called for printing to the console using `cat()`
#' @keywords internal
print_drugdevelopResult_helper <- function(x, ...){
  cat(" Utility: ", x$u, "\n", sep = "")
  if("Adj" %in% names(x)){
    cat(" Bias adjustment parameter: ", x$Adj, "\n", sep = "")
  }
  cat(" Sample size:\n")
  cat("   phase II: ", x$n2, ", phase III: ", x$n3, ", total: ", x$n, 
      "\n", sep = "")
  if("d2" %in% names(x)){
    cat(" Expected number of events:\n")
    cat("   phase II: ", x$d2, ", phase III: ", x$d3, ", total: ", x$d, 
        "\n", sep = "")
  }
  if("xi3" %in% names(x)){
    
    if("xi2" %in% names(x)) {
      cat(" Assumed event rate:\n")
      cat("   phase II: ", x$xi2, ", phase III: ", x$xi3, 
          "\n", sep = "")
    } else {
      cat(" Assumed event rate in phase III:", x$xi3, 
          "\n", sep = "")
    }

  }
  cat(" Probability to go to phase III: ", x$pgo, "\n", sep = "")
  cat(" Total cost:\n")
  cat("   phase II: ", x$K2, ", phase III: ", x$K3, 
      ", cost constraint: ", x$K, "\n", sep = "")
  cat(" Fixed cost:\n")
  cat("   phase II: ", x$c02, ", phase III: ", x$c03, "\n", sep = "")
  cat(" Variable cost per patient:\n")
  cat("   phase II: ", x$c2, ", phase III: ", x$c3, "\n", sep = "")
  if(!("b11" %in% names(x))){
    cat(" Effect size categories (expected gains):\n")
    cat("  small: ", x$steps1, " (", x$b1, "),", " medium: ", x$stepm1, 
        " (", x$b2, "),", " large: ", x$stepl1, " (", x$b3, ")\n", sep = "")
  }
  if("b11" %in% names(x)){
    cat(" Effect size categories:\n")
    cat("  small: ", x$steps1, " medium: ", x$stepm1, " large: ", x$stepl1,
        "\n", sep = "")
    cat(" Expected gains if endpoint 1 is significant:\n")
    cat("  small: ", x$b11, " medium: ", x$b21, " large: ", x$b31,
        "\n", sep = "")
    cat(" Expected gains if only endpoint 2 is significant:\n")
    cat("  small: ", x$b12, " medium: ", x$b22, " large: ", x$b32,
        "\n", sep = "")
  }
  cat(" Success probability: ", x$sProg, "\n", sep = "")
  if("sProg1" %in% names(x)){
    cat(" Success probability by effect size:\n")
    cat("   small: ", x$sProg1, ", medium: ", x$sProg2, 
        ", large: ", x$sProg3, "\n", sep = "")
  } else if("sProg2" %in% names(x)){
    cat(" Success probability for a trial with:\n")
    cat("   two arms in phase III: ", x$sProg2, 
        ", three arms in phase III: ", x$sProg3, "\n", sep = "")
  }
  if("OS" %in% names(x)){
    cat(" Probability of endpoint 1 being significant in phase III: ", x$OS, "\n", sep = "")
  }
  if("OP" %in% names(x)){
    cat(" Probability for at least one endpoint being significant: ", x$OP, "\n", sep = "")
  }
  cat(" Significance level: ", x$alpha, "\n", sep = "")
  cat(" Targeted power: ", 1-x$beta, "\n", sep = "")
  if("Kappa" %in% names(x)){
    cat(" Decision rule threshold: ", x$Kappa, " [Kappa] \n", sep = "")
  }
  if("RRgo" %in% names(x)){
    cat(" Decision rule threshold: ", x$RRgo, " [RRgo] \n", sep = "")
  }
  if("HRgo" %in% names(x)){
    cat(" Decision rule threshold: ", x$HRgo, " [HRgo] \n", sep = "")
  }
  if("Delta" %in% names(x)){
    cat(" Assumed true effect: ", x$Delta, " [Delta] \n", sep = "")
  }
  if("p1" %in% names(x)){
    cat(" Assumed true effect:\n")
    cat("  rate in the control: ", x$p0, ", rate in the treatment group: ", x$p1,
        "\n", sep = "")
  }
  if("hr" %in% names(x)){
    cat(" Assumed true effect: ", x$hr, " [hr] \n", sep = "")
  }
  if("median_prior_Delta" %in% names(x)){
    if(!is.na(x$median_prior_Delta)){
      cat(" Assumed true effect for planning the phase III trial (median of prior): ", 
          x$median_prior_Delta, "\n", sep = "")
    }
  }
  if("median_prior_RR" %in% names(x)){
    if(!is.na(x$median_prior_RR)){
      cat(" Assumed true effect for planning the phase III trial (median of prior): ", 
          x$median_prior_RR, "\n", sep = "")
    }
  }
  if("median_prior_HR" %in% names(x)){
    if(!is.na(x$median_prior_HR)){
      cat(" Assumed true effect for planning the phase III trial (median of prior): ", 
          x$median_prior_HR, "\n", sep = "")
    }
  }
  if("Delta1" %in% names(x) &
     !(!("Case" %in% names(x)) & "Strategy" %in% names(x)) &
     !("rho" %in% names(x))){
    cat(" Parameters of the prior distribution: \n", sep = "")
    cat("   Delta1: ", x$Delta1, ", Delta2: ", x$Delta2, ", in1: ", x$in1,
        ", in2: ", x$in2, ",\n   a: ", x$a, ", b: ", x$b, ", w: ", x$w, "\n", sep = "")
  }
  if("Delta1" %in% names(x) & 
     (!("Case" %in% names(x)) & "Strategy" %in% names(x))){
    cat(" Assumed true effects [Delta]: \n", sep = "")
    cat("   treatment 1: ", x$Delta1, ", treatment 2: ",
        x$Delta2, "\n", sep = "")
  }
  if("Delta1" %in% names(x) & 
     ("rho" %in% names(x))){
    cat(" Assumed true effects [Delta]: \n", sep = "")
    cat("   endpoint 1: ", x$Delta1, ", endpoint 2: ",
        x$Delta2, "\n", sep = "")
    cat(" Correlation between endpoints: ", x$rho, "\n", sep = "")
  }
  if("p11" %in% names(x) & !(!("Case" %in% names(x)) & "Strategy" %in% names(x))){
    cat(" Parameters of the prior distribution: \n", sep = "")
    cat("   p0: ", x$p0,",  p11: ", x$p11, ", p12: ", x$p12, ", in1: ", x$in1,
        ", in2: ", x$in2, ", w: ", x$w, "\n", sep = "")
  }
  if("p11" %in% names(x) & (!("Case" %in% names(x)) & "Strategy" %in% names(x))){
    cat(" Assumed true rates: \n", sep = "")
    cat("   control group: ", x$p0,",  treatment 1: ", x$p11, ", treatment 2: ",
        x$p12,"\n", sep = "")
  }
  if("hr1" %in% names(x) & !("ec" %in% names(x)) & !("rho" %in% names(x))){
    cat(" Parameters of the prior distribution: \n", sep = "")
    cat("   hr1: ", x$hr1,",  hr2: ", x$hr2, ", id1: ", x$id1,
        ", id2: ", x$id2, ", w: ", x$w, "\n", sep = "")
  }
  if("hr1" %in% names(x) & "ec" %in% names(x)){
    cat(" Assumed true effects [HR]: \n", sep = "")
    cat("   treatment 1: ", x$hr1,",  treatment 2: ", x$hr2, "\n", sep = "")
    cat(" Control arm event rate: ", x$ec, "\n", sep = "")
  }
  if("hr1" %in% names(x) & ("rho" %in% names(x))){
    cat(" Assumed true effects [HR]: \n", sep = "")
    cat("   endpoint 1: ", x$hr1,",  endpoint2 2: ", x$hr2, "\n", sep = "")
    cat(" Correlation between endpoints: ", x$rho, "\n", sep = "")
  }
  if("gamma" %in% names(x)){
    if(!is.null(x$gamma)){
      cat(" Treatment effect offset between phase II and III: ", x$gamma, 
          " [gamma] \n", sep = "")
    }
  }
}
