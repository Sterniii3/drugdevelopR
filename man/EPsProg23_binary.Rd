% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_multitrial_binary.R
\name{EPsProg23_binary}
\alias{EPsProg23_binary}
\title{Expected probability of a successful program deciding between two or three phase III trials for a binary distributed outcome}
\usage{
EPsProg23_binary(
  RRgo,
  n2,
  alpha,
  beta,
  w,
  p0,
  p11,
  p12,
  in1,
  in2,
  case,
  size,
  ymin
)
}
\arguments{
\item{RRgo}{threshold value for the go/no-go decision rule}

\item{n2}{total sample size for phase II; must be even number}

\item{alpha}{significance level}

\item{beta}{\code{1-beta} power for calculation of sample size for phase III}

\item{w}{weight for mixture prior distribution}

\item{p0}{assumed true rate of control group}

\item{p11}{assumed true rate of treatment group}

\item{p12}{assumed true rate of treatment group}

\item{in1}{amount of information for \code{p11} in terms of sample size}

\item{in2}{amount of information for \code{p12} in terms of sample size}

\item{case}{choose case: "at least 1, 2 or 3 significant trials needed for approval"}

\item{size}{size category \code{"small"}, \code{"medium"} or \code{"large"}}

\item{ymin}{assumed minimal clinical relevant effect}
}
\value{
The output of the function \code{EPsProg23_binary()} is the expected probability of a successful program.
}
\description{
The function \code{EPsProg23_binary()} calculates the expected probability of a successful program
with a normally distributed outcome. This function follows a special decision rule in order to determine
whether two or three phase III trials should be conducted. First, two phase III trials are performed. Depending
on their success, the decision for a third phase III trial is made:
\itemize{
\item If both trials are successful, no third phase III trial will be conducted.
\item If only one of the two trials is successful and the other trial has a treatment effect that points in the same direction,
a third phase III trial will be conducted with a sample size of N3 = N3(ymin), which depends on an assumed minimal clinical relevant effect (\code{ymin}).
The third trial then has to be significant at level \code{alpha}
\item If only one of the two trials is successful and the treatment effect of the other points in opposite direction or
if none of the two trials are successful, then no third trial is performed and the drug development development program is not successful.
In the utility function, this will lead to a utility of -9999.
}
}
\examples{
res <- EPsProg23_binary(RRgo = 0.8, n2 = 50,  alpha = 0.025, beta = 0.1, 
                                 w = 0.6,  p0 = 0.3, p11 =  0.3, p12 = 0.5, 
                                 in1 = 300, in2 = 600, case = 2, size = "small",
                                 ymin = 0.5)
}
\keyword{internal}
