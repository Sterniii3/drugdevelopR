% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_multitrial_normal.R
\name{EPsProg_multitrial_normal}
\alias{EPsProg_multitrial_normal}
\alias{EPsProg2_normal}
\alias{EPsProg3_normal}
\alias{EPsProg4_normal}
\title{Expected probability of a successful program for multitrial programs with normally distributed outcomes}
\usage{
EPsProg2_normal(
  kappa,
  n2,
  alpha,
  beta,
  w,
  Delta1,
  Delta2,
  in1,
  in2,
  a,
  b,
  case,
  size,
  fixed
)

EPsProg3_normal(
  kappa,
  n2,
  alpha,
  beta,
  w,
  Delta1,
  Delta2,
  in1,
  in2,
  a,
  b,
  case,
  size,
  fixed
)

EPsProg4_normal(
  kappa,
  n2,
  alpha,
  beta,
  w,
  Delta1,
  Delta2,
  in1,
  in2,
  a,
  b,
  case,
  size,
  fixed
)
}
\arguments{
\item{kappa}{threshold value for the go/no-go decision rule}

\item{n2}{total sample size for phase II; must be even number}

\item{alpha}{significance level}

\item{beta}{\code{1-beta} power for calculation of sample size for phase III}

\item{w}{weight for mixture prior distribution}

\item{Delta1}{assumed true treatment effect for standardized difference in means}

\item{Delta2}{assumed true treatment effect for standardized difference in means}

\item{in1}{amount of information for \code{Delta1} in terms of sample size}

\item{in2}{amount of information for \code{Delta2} in terms of sample size}

\item{a}{lower boundary for the truncation}

\item{b}{upper boundary for the truncation}

\item{case}{choose case: "at least 1, 2 or 3 significant trials needed for approval"}

\item{size}{size category \code{"small"}, \code{"medium"} or \code{"large"}}

\item{fixed}{choose if true treatment effects are fixed or random}
}
\value{
The output of the function \code{EPsProg2_normal()}, \code{EPsProg3_normal()} and \code{EPsProg4_normal()} is the expected probability of a successful program when performing several phase III trials (2, 3 or 4 respectively).
}
\description{
These functions calculate the expected probability of a successful program given the parameters.
Each function represents a specific strategy, e.g. the function \code{EpsProg3_normal()} calculates the expected probability if three phase III trials are performed.
The parameter case specifies how many of the trials have to be successful, i.e. how many trials show a significantly relevant positive treatment effect.
}
\details{
The following cases can be investigated by the software:
\itemize{
\item Two phase III trials
\itemize{
\item Case 1: Strategy 1/2; at least one trial significant, the treatment effect of the other one at least showing in the same direction
\item Case 2: Strategy 2/2; both trials significant
}
\item Three phase III trials
\itemize{
\item Case 2: Strategy 2/3; at least two trials significant, the treatment effect of the other one at least showing in the same direction
\item Case 3: Strategy 3/3; all trials significant
}
\item Four phase III trials
\itemize{
\item Case 3: Strategy 3/4; at least three trials significant, the treatment effect of the other one at least showing in the same direction
}
}
}
\examples{
\donttest{EPsProg2_normal(kappa = 0.1, n2 = 50, alpha = 0.025, beta = 0.1, w = 0.3,
                                 Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
                                 a = 0.25, b = 0.75, 
                                 case = 2, size = "small", fixed = FALSE)}
          \donttest{EPsProg3_normal(kappa = 0.1, n2 = 50, alpha = 0.025, beta = 0.1, w = 0.3,
                                 Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
                                 a = 0.25, b = 0.75, 
                                 case = 2, size = "small", fixed = TRUE)}
          \donttest{EPsProg4_normal(kappa = 0.1, n2 = 50, alpha = 0.025, beta = 0.1, w = 0.3,
                                 Delta1 = 0.375, Delta2 = 0.625, in1 = 300, in2 = 600, 
                                 a = 0.25, b = 0.75, 
                                 case = 3, size = "small", fixed = TRUE)}                      
}
\keyword{internal}
