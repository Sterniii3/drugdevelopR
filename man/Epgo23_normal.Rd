% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_multitrial_normal.R
\name{Epgo23_normal}
\alias{Epgo23_normal}
\title{Expected probability to do third phase III trial}
\usage{
Epgo23_normal(kappa, n2, alpha, beta, a, b, w, Delta1, Delta2, in1, in2)
}
\arguments{
\item{kappa}{threshold value for the go/no-go decision rule}

\item{n2}{total sample size for phase II; must be even number}

\item{alpha}{significance level}

\item{beta}{\code{1-beta} power for calculation of sample size for phase III}

\item{a}{lower boundary for the truncation}

\item{b}{upper boundary for the truncation}

\item{w}{weight for mixture prior distribution}

\item{Delta1}{assumed true treatment effect for standardized difference in means}

\item{Delta2}{assumed true treatment effect for standardized difference in means}

\item{in1}{amount of information for \code{Delta1} in terms of sample size}

\item{in2}{amount of information for \code{Delta2} in terms of sample size}
}
\value{
The output of the function \code{Epgo23_normal()} is the probability to a third phase III trial.
}
\description{
In the setting of Case 2: Strategy 2/2( + 1); at least two trials significant (and the
treatment effect of the other one at least showing in the same direction) this function calculates the probability that a third phase III trial is necessary.
}
\examples{
\donttest{Epgo23_normal(kappa = 0.1, n2 = 50, w = 0.3, alpha = 0.025, beta = 0.1, a = 0.25, b=0.75,
                               Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600)}
}
\keyword{internal}
