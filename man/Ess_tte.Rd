% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_multiarm.R
\name{Ess_tte}
\alias{Ess_tte}
\title{Expected sample size for phase III for multiarm programs with time-to-event outcomes}
\usage{
Ess_tte(HRgo, n2, alpha, beta, ec, hr1, hr2, strategy, case)
}
\arguments{
\item{HRgo}{threshold value for the go/no-go decision rule}

\item{n2}{total sample size for phase II; must be even number}

\item{alpha}{significance level}

\item{beta}{1-beta power for calculation of sample size for phase III}

\item{ec}{control arm event rate for phase II and III}

\item{hr1}{assumed true treatment effect on HR scale for treatment 1}

\item{hr2}{assumed true treatment effect on HR scale for treatment 2}

\item{strategy}{choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)}

\item{case}{different cases: 1 ("nogo"), 21 (treatment 1 is promising, treatment 2 is not), 22 (treatment 2 is promising, treatment 1 is not), 31 (both treatments are promising, treatment 1 is better), 32 (both treatments are promising, treatment 2 is better)}
}
\value{
the function Ess_tte() returns the expected sample size for phase III when going to phase III
}
\description{
Given phase II results are promising enough to get the "go"-decision to go to phase III this function now calculates the expected sample size for phase III given the cases and strategies listed below.
The results of this function are necessary for calculating the utility of the program, which is then in a further step maximized by the \code{optimal_multiarm()} function
}
\examples{
res <- Ess_tte(HRgo = 0.8 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
                            ec = 0.6, hr1 = 0.7, hr2 = 0.8, strategy = 2, case = 21)
}
\keyword{internal}
