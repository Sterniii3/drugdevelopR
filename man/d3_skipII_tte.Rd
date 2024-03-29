% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_tte.R
\name{d3_skipII_tte}
\alias{d3_skipII_tte}
\title{Expected probability to go to phase III for time-to-event outcomes}
\usage{
d3_skipII_tte(alpha, beta, median_prior)
}
\arguments{
\item{alpha}{significance level}

\item{beta}{\code{1-beta} power for calculation of sample size for phase III}

\item{median_prior}{the median_prior is given as -log(hr1), the assumed true treatment effect}
}
\value{
The output of the functions \code{d3_skipII_tte()} is the expected number of events in phase III when skipping phase II.
}
\description{
If choosing \code{skipII = TRUE}, the program calculates the expected utility for the case when phase
II is skipped and compares it to the situation when phase II is not skipped.
This function calculates the expected sample size for phase III for time-to-event outcomes using a median prior.
}
\examples{
res <- d3_skipII_tte(alpha = 0.05, beta = 0.1, median_prior = 0.35)
}
\keyword{internal}
