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
#' res <- optimal_normal(w=0.3,                       
#'   Delta1 = 0.375, Delta2 = 0.625, in1=300, in2=600,  
#'   a = 0.25, b = 0.75,
#'   n2min = 20, n2max = 100, stepn2 = 4,               
#'   kappamin = 0.02, kappamax = 0.2, stepkappa = 0.02, 
#'   alpha = 0.05, beta = 0.1,                       
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
#' print(res)                                  
#' }
#' 
#' 
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
      if(sequence){
        cat(attr(res, "comment"))
      }
    }
    
  }
  # General option
  else{
    cat("Optimization result:\n")
    print_drugdevelopResult_helper(x, ...)
    if(sequence){
      cat(attr(res, "comment"))
    }
  }
  

}

#' Helper function for printing a drugdevelopResult Object
#'
#' @param x Data frame
#' @param ... Further arguments.
#'
#' @keywords internal
print_drugdevelopResult_helper <- function(x, ...){
  cat(" Utility: ", x$u, "\n", sep = "")
  cat(" Sample size:\n")
  cat("   Phase II: ", x$n2, ", phase III: ", x$n3, ", total: ", x$n, 
      "\n", sep = "")
  if("d2" %in% names(x)){
    cat(" Expected number of events:\n")
    cat("   Phase II: ", x$d2, ", phase III: ", x$d3, ", total: ", x$d, 
        "\n", sep = "")
  }
  if("xi3" %in% names(x)){
    
    if("xi2" %in% names(x)) {
      cat(" Assumed event rate:\n")
      cat("   Phase II: ", x$xi2, ", phase III: ", x$xi3, ", total: ", x$d, 
          "\n", sep = "")
    } else {
      cat(" Assumed event rate in phase III:", x$xi3, 
          "\n", sep = "")
    }

  }
  cat(" Probability to go to phase III: ", x$pgo, "\n", sep = "")
  cat(" Total cost:\n")
  cat("   Phase II: ", x$K2, ", phase III: ", x$K3, 
      ", cost constraint: ", x$K, "\n", sep = "")
  cat(" Fixed cost:\n")
  cat("   Phase II: ", x$c02, ", phase III: ", x$c03, "\n", sep = "")
  cat(" Variable cost per patient:\n")
  cat("   Phase II: ", x$c2, ", phase III: ", x$c3, "\n", sep = "")
  cat(" Effect size categories (expected gains):\n")
  cat("  small: ", x$steps1, " (", x$b1, "),", " medium: ", x$stepm1, 
      " (", x$b2, "),", " large: ", x$stepl1, " (", x$b3, ")\n", sep = "")
  cat(" Success probability: ", x$sProg, "\n", sep = "")
  cat(" Success probability for effect size:\n")
  cat("   small: ", x$sProg1, ", medium: ", x$sProg2, 
      ", large: ", x$sProg3, "\n", sep = "")
  if("Kappa" %in% names(x)){
    cat(" Decision rule threshold: ", x$Kappa, " (Kappa) \n", sep = "")
  }
  if("RRgo" %in% names(x)){
    cat(" Decision rule threshold: ", x$RRgo, " (RRgo) \n", sep = "")
  }
  if("HRgo" %in% names(x)){
    cat(" Decision rule threshold: ", x$RRgo, " (HRgo) \n", sep = "")
  }
  if("Delta" %in% names(x)){
    cat(" Assumed fixed effect: ", x$Delta, " (Delta) \n", sep = "")
  }
  if("p1" %in% names(x)){
    cat(" Assumed fixed effect:\n")
    cat("  Rate in the control:", x$p0, ", rate in the treatment group:", x$p1,
        "\n", sep = "")
  }
  if("hr" %in% names(x)){
    cat(" Assumed fixed effect: ", x$Delta, " (Delta) \n", sep = "")
  }
  if("median_prior_Delta" %in% names(x)){
    if(!is.na(x$median_prior_Delta)){
      cat(" Assumed fixed effect for planning the phase III trial (median of prior): ", 
          x$median_prior_Delta, "\n", sep = "")
    }
  }
  if("median_prior_RR" %in% names(x)){
    if(!is.na(x$median_prior_RR)){
      cat(" Assumed fixed effect for planning the phase III trial (median of prior): ", 
          x$median_prior_RR, "\n", sep = "")
    }
  }
  if("median_prior_HR" %in% names(x)){
    if(!is.na(x$median_prior_HR)){
      cat(" Assumed fixed effect for planning the phase III trial (median of prior): ", 
          x$median_prior_HR, "\n", sep = "")
    }
  }
  if("Delta1" %in% names(x)){
    cat(" Parameters of the prior distribution: \n", sep = "")
    cat("   Delta1: ", x$Delta1, ", Delta2: ", x$Delta2, ", in1: ", x$in1,
        ", in2: ", x$in2, ", a: ", x$a, ", b: ", x$b, ", w: ", x$w, "\n", sep = "")
  }
  if("p11" %in% names(x)){
    cat(" Parameters of the prior distribution: \n", sep = "")
    cat("   p0: ", x$p0,",  p11: ", x$p11, ", p12: ", x$p12, ", in1: ", x$in1,
        ", in2: ", x$in2, ", w: ", x$w, "\n", sep = "")
  }
  if("hr1" %in% names(x)){
    cat(" Parameters of the prior distribution: \n", sep = "")
    cat("   hr1: ", x$hr1,",  hr2: ", x$hr2, ", id1: ", x$id1,
        ", id2: ", x$id2, ", w: ", x$w, "\n", sep = "")
  }
  if("gamma" %in% names(x)){
    if(!is.null(x$gamma)){
      cat(" Treatment effect offset between phase II and III: ", x$gamma, 
          " (gamma) \n", sep = "")
    }
  }
}