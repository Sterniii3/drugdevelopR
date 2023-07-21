#' Construct a drugdevelopResult object from a data frame
#'
#' This is a short wrapper for adding the "drugdevelopR" string to the list
#' of S3 classes of a data frame.
#'
#' @param x data frame
#'
#' @export
#' @keywords internal
drugdevelopResult <- function(x, ...) {
  class(x) <- c("drugdevelopResult", class(x))
  return(x)
}
