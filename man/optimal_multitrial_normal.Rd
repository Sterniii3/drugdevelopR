% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimal_multitrial_normal.R
\name{optimal_multitrial_normal}
\alias{optimal_multitrial_normal}
\title{Optimal phase II/III drug development planning where several phase III trials are performed}
\usage{
optimal_multitrial_normal(
  w,
  Delta1,
  Delta2,
  in1,
  in2,
  a,
  b,
  n2min,
  n2max,
  stepn2,
  kappamin,
  kappamax,
  stepkappa,
  alpha,
  beta,
  c2,
  c3,
  c02,
  c03,
  K = Inf,
  N = Inf,
  S = -Inf,
  b1,
  b2,
  b3,
  case,
  strategy = TRUE,
  fixed = FALSE,
  num_cl = 1
)
}
\arguments{
\item{w}{weight for
\href{https://web.imbi.uni-heidelberg.de/prior/}{mixture prior distribution}}

\item{Delta1}{assumed true prior treatment effect measured as the
standardized difference in means, see
\href{https://web.imbi.uni-heidelberg.de/prior/}{here} for details}

\item{Delta2}{assumed true prior treatment effect measured as the
standardized difference in means, see
\href{https://web.imbi.uni-heidelberg.de/prior/}{here} for details}

\item{in1}{amount of information for \code{Delta1} in terms of sample size, see
\href{https://web.imbi.uni-heidelberg.de/prior/}{here}
for details}

\item{in2}{amount of information for \code{Delta2} in terms of sample size, see
\href{https://web.imbi.uni-heidelberg.de/prior/}{here}
for details}

\item{a}{lower boundary for the truncation of the \href{https://web.imbi.uni-heidelberg.de/prior/}{prior distribution}}

\item{b}{upper boundary for the truncation of the \href{https://web.imbi.uni-heidelberg.de/prior/}{prior distribution}}

\item{n2min}{minimal total sample size for phase II; must be an even number}

\item{n2max}{maximal total sample size for phase II, must be an even number}

\item{stepn2}{step size for the optimization over n2; must be an even number}

\item{kappamin}{minimal threshold value kappa for the go/no-go decision rule}

\item{kappamax}{maximal threshold value  kappa for the go/no-go decision rule}

\item{stepkappa}{step size for the optimization over the threshold value kappa}

\item{alpha}{one-sided significance level}

\item{beta}{type II error rate; i.e. \code{1 - beta} is the power for calculation of the sample size for phase III}

\item{c2}{variable per-patient cost for phase II in 10^5 $}

\item{c3}{variable per-patient cost for phase III in 10^5 $}

\item{c02}{fixed cost for phase II in 10^5 $}

\item{c03}{fixed cost for phase III in 10^5 $}

\item{K}{constraint on the costs of the program, default: Inf, e.g. no constraint}

\item{N}{constraint on the total expected sample size of the program, default: Inf, e.g. no constraint}

\item{S}{constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint}

\item{b1}{expected gain for effect size category "small" in 10^5 $}

\item{b2}{expected gain for effect size category "medium" in 10^5 $}

\item{b3}{expected gain for effect size category "large" in 10^5 $}

\item{case}{choose case: "at least 1, 2 or 3 significant trials needed for approval"}

\item{strategy}{choose strategy: "conduct 1, 2, 3 or 4 trials in order to achieve the case's goal"; TRUE calculates all strategies of the selected \code{case}}

\item{fixed}{choose if true treatment effects are fixed or following a prior distribution, if TRUE \code{Delta1} is used as fixed effect}

\item{num_cl}{number of clusters used for parallel computing, default: 1}
}
\value{
The output of the function is a \code{data.frame} object containing the optimization results:
\describe{\item{Case}{Case: "number of significant trials needed"}\item{Strategy}{Strategy: "number of trials to be conducted in order to achieve the goal of the case"}\item{u}{maximal expected utility under the optimization constraints, i.e. the expected utility of the optimal sample size and threshold value}\item{Kappa}{optimal threshold value for the decision rule to go to phase III}\item{n2}{total sample size for phase II; rounded to the next even natural number}
\item{n3}{total sample size for phase III; rounded to the next even natural number}
\item{n}{total sample size in the program; n = n2 + n3}
\item{K}{maximal costs of the program (i.e. the cost constraint, if it is set or the sum K2+K3 if no cost constraint is set)}
\item{pgo}{probability to go to phase III}
\item{sProg}{probability of a successful program}\item{sProg1}{probability of a successful program with "small" treatment effect in phase III (lower boundary in HR scale is set to 0, as proposed by Cohen (1988))}
\item{sProg2}{probability of a successful program with "medium" treatment effect in phase III (lower boundary in HR scale is set to 0.5, as proposed Cohen (1988))}
\item{sProg3}{probability of a successful program with "large" treatment effect in phase III (lower boundary in HR scale is set to 0.8, as proposed Cohen (1988))}\item{K2}{expected costs for phase II}
\item{K3}{expected costs for phase III}}
and further input parameters. Taking \code{cat(comment())} of the
data frame lists the used optimization sequences, start and
finish time of the optimization procedure. The attribute
\code{attr(,"trace")} returns the utility values of all parameter
combinations visited during optimization.
}
\description{
The \code{optimal_multitrial_normal} function enables planning of phase II/III
drug development programs with several phase III trials for
the same normally distributed endpoint. Its main output values are optimal
sample size allocation and go/no-go decision rules. For normally distributed
endpoints, the treatment effect is measured by the standardized difference in
means (Delta). The assumed true treatment effects can be assumed fixed or
modelled by a prior distribution.
}
\details{
The R Shiny application \href{https://web.imbi.uni-heidelberg.de/prior/}{prior}
visualizes the prior distributions used in this package. Fast computing is
enabled by parallel programming.
}
\section{Effect sizes}{

In other settings, the definition of "small", "medium" and "large" effect
sizes can be user-specified using the input parameters \code{steps1}, \code{stepm1} and
\code{stepl1}. Due to the complexity of the multitrial setting, this feature is
not included for this setting. Instead, the effect sizes were set to
to predefined values as explained under sProg1, sProg2 and sProg3 in the
\emph{Value} section.
}

\examples{
# Activate progress bar (optional)
\dontrun{progressr::handlers(global = TRUE)}
# Optimize
\donttest{
optimal_multitrial_normal(w = 0.3,           # define parameters for prior
  Delta1 = 0.375, Delta2 = 0.625,
  in1 = 300, in2 = 600,                               # (https://web.imbi.uni-heidelberg.de/prior/)
  a = 0.25, b = 0.75,
  n2min = 20, n2max = 100, stepn2 = 4,                # define optimization set for n2
  kappamin = 0.02, kappamax = 0.2, stepkappa = 0.02,  # define optimization set for kappa
  alpha = 0.025, beta = 0.1,                          # drug development planning parameters
  c2 = 0.675, c3 = 0.72, c02 = 15, c03 = 20,          # fixed and variable costs for phase II/III
  K = Inf, N = Inf, S = -Inf,                         # set constraints
  b1 = 3000, b2 = 8000, b3 = 10000,                   # expected benefit for each effect size                                         
  case = 1, strategy = TRUE,                          # chose Case and Strategy
  fixed = TRUE,                                       # true treatment effects are fixed/random
  num_cl = 1)                                         # number of cores for parallelized computing
  }

}
\references{
Cohen, J. (1988). Statistical power analysis for the behavioral sciences.
}
