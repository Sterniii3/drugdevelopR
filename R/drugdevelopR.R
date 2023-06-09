#' Utility based optimal phase II/III drug development planning 
#'
#' The drugdevelopR package enables utility based planning of phase II/III drug development programs with optimal sample size allocation and go/no-go decision rules. The assumed true treatment effects can be assumed fixed (planning is then also possible via user friendly R Shiny App: \href{https://web.imbi.uni-heidelberg.de/drugdevelopR/}{drugdevelopR}) or modelled by a prior distribution. The R Shiny application \href{https://web.imbi.uni-heidelberg.de/prior/}{prior} visualizes the prior distributions used in this package. Fast computing is enabled by parallel programming.
#' @docType package
#' @name drugdevelopR
#' @section drugdevelopR package and R Shiny App:
#' The drugdevelopR package provides the functions to plan optimal phase II/III drug development programs with
#' \itemize{
#'   \item time-to-event endpoint (\code{\link{optimal_tte}}),
#'   \item binary endpoint (\code{\link{optimal_binary}}) and
#'   \item normally distributed endpoint (\code{\link{optimal_normal}}),
#' }
#' where the treatment effect is assumed fixed or modelled by a \href{https://web.imbi.uni-heidelberg.de/prior/}{prior}. In these settings, optimal phase II/III drug development planning with fixed assumed treatment effects can also be done with the help of the R Shiny application \href{https://web.imbi.uni-heidelberg.de/basic/}{basic}. Extensions to the basic setting are 
#' \itemize{
#'   \item optimal planning of programs including methods for discounting of phase II results (function: \code{\link{optimal_bias}}, App: \href{https://web.imbi.uni-heidelberg.de/bias/}{bias}),
#'   \item optimal planning of programs with several phase III trials (function: \code{\link{optimal_multitrial}}, App: \href{https://web.imbi.uni-heidelberg.de/multitrial/}{multitrial}) and
#'   \item optimal planning of programs with multiple arms (function: \code{\link{optimal_multiarm}}, App: \href{https://web.imbi.uni-heidelberg.de/multiarm/}{multiarm}). 
#' }
#' The R Shiny App \href{https://web.imbi.uni-heidelberg.de/drugdevelopR/}{drugdevelopR} serves as homepage, navigating the different parts of drugdevelopR via links.
#' @references
#' Kirchner, M., Kieser, M., Goette, H., & Schueler, A. (2016). Utility-based optimization of phase II/III programs. Statistics in Medicine, 35(2), 305-316.
#' 
#' Preussler, S., Kieser, M., and Kirchner, M. (2019). Optimal sample size allocation and go/no-go decision rules for phase II/III programs where several phase III trials are performed. Biometrical Journal, 61(2), 357-378.
#'
#' Preussler, S., Kirchner, M., Goette, H., Kieser, M. (2019). Optimal designs for phase II/III drug development programs including methods for discounting of phase II results. Submitted to peer-review journal.
#' 
#' Preussler, S., Kirchner, M., Goette, H., Kieser, M. (2019). Optimal designs for multi-arm Phase II/III drug development programs. Submitted to peer-review journal.
#'
#' 
#' @keywords internal
#' @export
drugdevelopR <- function(){
  cat("The package drugdevelopR enables utility based optimal phase II/III drug development planning")
}  
