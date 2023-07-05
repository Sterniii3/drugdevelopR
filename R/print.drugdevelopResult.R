#' Printing a drugdevelopResult Object
#' 
#' Displays details about the optimal results from a \code{drugdevelopResult} object.
#'
#' @param x Data frame of class \code{drugdevelopResult}.
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
#'  print(res)                                  
#' }
#' 
#' 
print.drugdevelopResult <- function(x, ...) {
  details <- match.arg(details)
  detail.df <- get_details(x, details, max_n)
  
  if (details == "low") {
    x_list <- list(
      Point_Estimate_n = x$sample_size,
      Minimum_Sufficient_n = detail.df,
      Message = x$exit.mes
    )
  } else {
    x_list <- list(
      Point_Estimate_n = x$sample_size,
      Details = format(detail.df, digits = digits),
      Message = x$exit.mes
    )
  }
  if (invisible) {
    invisible(x_list)
  } else {
    print(x_list, ...)
  }
}