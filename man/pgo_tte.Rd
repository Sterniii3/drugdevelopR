% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_multiarm.R
\name{pgo_tte}
\alias{pgo_tte}
\title{Probability to go to phase III for multiarm programs with time-to-event outcomes}
\usage{
pgo_tte(HRgo, n2, ec, hr1, hr2, strategy, case)
}
\arguments{
\item{HRgo}{threshold value for the go/no-go decision rule}

\item{n2}{total sample size in phase II, must be divisible by 3}

\item{ec}{control arm event rate for phase II and III}

\item{hr1}{assumed true treatment effect on HR scale for treatment 1}

\item{hr2}{assumed true treatment effect on HR scale for treatment 2}

\item{strategy}{choose Strategy: 1 ("only best promising"), 2 ("all promising")}

\item{case}{different cases: 1 ("nogo"), 21 (treatment 1 is promising, treatment 2 is not), 22 (treatment 2 is promising, treatment 1 is not), 31 (both treatments are promising, treatment 1 is better), 32 (both treatments are promising, treatment 2 is better)}
}
\value{
The function pgo_tte() returns the probability to go to phase III.
}
\description{
Given our parameters this function calculates the probability to go to phase III after the second phase was conducted. The considered strategies are as follows:
\itemize{
\item 
\enumerate{
\item Strategy: Only best promising treatment goes to phase III
}
}
}
\examples{
res <- pgo_tte(HRgo = 0.8, n2 = 48 , ec = 0.6, hr1 = 0.7, hr2 = 0.8, strategy = 2, case = 31)
}
\keyword{internal}
