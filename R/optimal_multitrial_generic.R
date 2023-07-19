#' Generic function for optimizing multi-trial programs
#' 
#' @name optimal_multitrial_generic
#' 
#' @section Effect sizes:
#' In other settings, the definition of "small", "medium" and "large" effect
#' sizes can be user-specified using the input parameters `steps1`, `stepm1` and
#' `stepl1`. Due to the complexity of the multitrial setting, this feature is
#' not included for this setting. Instead, the effect sizes were set to
#' to predefined values as explained under sProg1, sProg2 and sProg3 in the 
#' *Value* section.
#' 
#' @param case choose case: "at least 1, 2 or 3 significant trials needed for approval"
#' @param strategy choose strategy: "conduct 1, 2, 3 or 4 trials in order to achieve the case's goal"; TRUE calculates all strategies of the selected `case` 
#' 
#' @keywords internal
optimal_multitrial_generic <- function(case, strategy = TRUE){
  # This function is only used for documentation.
  # Hence, it contains no code.
}