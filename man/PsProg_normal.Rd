% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_multiarm_normal.R
\name{PsProg_normal}
\alias{PsProg_normal}
\title{Probability of a successful program for multiarm programs with normally distributed outcomes}
\usage{
PsProg_normal(
  kappa,
  n2,
  alpha,
  beta,
  Delta1,
  Delta2,
  step1,
  step2,
  strategy,
  case
)
}
\arguments{
\item{kappa}{threshold value for the go/no-go decision rule}

\item{n2}{total sample size for phase II; must be even number}

\item{alpha}{significance level}

\item{beta}{1-beta power for calculation of sample size for phase III}

\item{Delta1}{assumed true treatment effect for standardized difference in means}

\item{Delta2}{assumed true treatment effect for standardized difference in means}

\item{step1}{lower boundary for effect size}

\item{step2}{upper boundary for effect size}

\item{strategy}{choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)}

\item{case}{different cases: 1 ("nogo"), 21 (treatment 1 is promising, treatment 2 is not), 22 (treatment 2 is promising, treatment 1 is not), 31 (both treatments are promising, treatment 1 is better), 32 (both treatments are promising, treatment 2 is better)}
}
\value{
The function PsProg_normal() returns the probability of a successful program.
}
\description{
Given we get the "go"-decision in phase II, this functions now calculates the probability that the results of the confirmatory trial (phase III) are significant, i.e. we have a statistically relevant positive effect of the treatment.
}
\examples{
res <- PsProg_normal(kappa = 0.1 ,n2 = 50 ,alpha = 0.05, beta = 0.1,
                            Delta1 = 0.375, Delta2 = 0.625,  step1 = 0, step2 = 0.5,
                            strategy = 3, case = 31)
}
\keyword{internal}
