% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimal_multitrial.R
\name{optimal_multitrial}
\alias{optimal_multitrial}
\title{Optimal phase II/III drug development planning where several phase III trials are performed for time-to-event endpoints}
\usage{
optimal_multitrial(
  w,
  hr1,
  hr2,
  id1,
  id2,
  d2min,
  d2max,
  stepd2,
  hrgomin,
  hrgomax,
  stephrgo,
  alpha,
  beta,
  xi2,
  xi3,
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
\item{w}{weight for mixture prior distribution, see
\href{https://web.imbi.uni-heidelberg.de/prior/}{this Shiny application}
for the choice of weights}

\item{hr1}{first assumed true treatment effect on HR scale for \href{https://web.imbi.uni-heidelberg.de/prior/}{prior distribution}}

\item{hr2}{second assumed true treatment effect on HR scale for \href{https://web.imbi.uni-heidelberg.de/prior/}{prior distribution}}

\item{id1}{amount of information for \code{hr1} in terms of number of events}

\item{id2}{amount of information for \code{hr2} in terms of number of events}

\item{d2min}{minimal number of events for phase II}

\item{d2max}{maximal number of events for phase II}

\item{stepd2}{step size for the optimization over d2}

\item{hrgomin}{minimal threshold value for the go/no-go decision rule}

\item{hrgomax}{maximal threshold value for the go/no-go decision rule}

\item{stephrgo}{step size for the optimization over HRgo}

\item{alpha}{one-sided significance level}

\item{beta}{type II error rate; i.e. \code{1 - beta} is the power for calculation of the number of events for phase III by Schoenfeld's formula (Schoenfeld 1981)}

\item{xi2}{assumed event rate for phase II, used for calculating the sample size of phase II via \code{n2 = d2/xi2}}

\item{xi3}{event rate for phase III, used for calculating the sample size of phase III in analogy to \code{xi2}}

\item{c2}{variable per-patient cost for phase II in 10^5 $.}

\item{c3}{variable per-patient cost for phase III in 10^5 $.}

\item{c02}{fixed cost for phase II in 10^5 $.}

\item{c03}{fixed cost for phase III in 10^5 $.}

\item{K}{constraint on the costs of the program, default: Inf, e.g. no constraint}

\item{N}{constraint on the total expected sample size of the program, default: Inf, e.g. no constraint}

\item{S}{constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint}

\item{b1}{expected gain for effect size category "small"}

\item{b2}{expected gain for effect size category "medium"}

\item{b3}{expected gain for effect size category "large"}

\item{case}{choose case: "at least 1, 2 or 3 significant trials needed for approval"}

\item{strategy}{choose strategy: "conduct 1, 2, 3 or 4 trials in order to achieve the case's goal"; TRUE calculates all strategies of the selected \code{case}}

\item{fixed}{choose if true treatment effects are fixed or random, if TRUE hr1 is used as a fixed effect and hr2 is ignored}

\item{num_cl}{number of clusters used for parallel computing, default: 1}
}
\value{
The output of the function is a \code{data.frame} object containing the optimization results:
\describe{\item{Case}{Case: "number of significant trials needed"}\item{Strategy}{Strategy: "number of trials to be conducted in order to achieve the goal of the case"}\item{u}{maximal expected utility under the optimization constraints, i.e. the expected utility of the optimal sample size and threshold value}\item{HRgo}{optimal threshold value for the decision rule to go to phase III}
\item{d2}{optimal total number of events for phase II}
\item{d3}{total expected number of events for phase III; rounded to next natural number}
\item{d}{total expected number of events in the program; d = d2 + d3}\item{n2}{total sample size for phase II; rounded to the next even natural number}
\item{n3}{total sample size for phase III; rounded to the next even natural number}
\item{n}{total sample size in the program; n = n2 + n3}
\item{K}{maximal costs of the program (i.e. the cost constraint, if it is set or the sum K2+K3 if no cost constraint is set)}
\item{pgo}{probability to go to phase III}
\item{sProg}{probability of a successful program}\item{sProg1}{probability of a successful program with "small" treatment effect in phase III (lower boundary in HR scale is set to 1, as proposed by IQWiG (2016))}
\item{sProg2}{probability of a successful program with "medium" treatment effect in phase III (lower boundary in HR scale is set to 0.95, as proposed by IQWiG (2016))}
\item{sProg3}{probability of a successful program with "large" treatment effect in phase III (lower boundary in HR scale is set to 0.85, as proposed by IQWiG (2016))}\item{K2}{expected costs for phase II}
\item{K3}{expected costs for phase III}}
and further input parameters. Taking \code{cat(comment())} of the
data frame lists the used optimization sequences, start and
finish time of the optimization procedure. The attribute
\code{attr(,"trace")} returns the utility values of all parameter
combinations visited during optimization.
}
\description{
The function \code{\link{optimal_multitrial}} of the drugdevelopR package enables planning of phase II/III drug development programs with time-to-event endpoints for programs with several phase III trials (Preussler et. al, 2019).
Its main output values are the optimal sample size allocation and optimal go/no-go decision rules.
The assumed true treatment effects can be assumed to be fixed (planning is then also possible via user friendly R Shiny App: \href{https://web.imbi.uni-heidelberg.de/multitrial/}{multitrial}) or can be modelled by a prior distribution.
The R Shiny application \href{https://web.imbi.uni-heidelberg.de/prior/}{prior} visualizes the prior distributions used in this package. Fast computing is enabled by parallel programming.
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
optimal_multitrial(w = 0.3,                # define parameters for prior
  hr1 = 0.69, hr2 = 0.88, id1 = 210, id2 = 420,     # (https://web.imbi.uni-heidelberg.de/prior/)
  d2min = 20, d2max = 100, stepd2 = 5,              # define optimization set for d2
  hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05,    # define optimization set for HRgo
  alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7,  # drug development planning parameters
  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,          # fixed and variable costs for phase II/III
  K = Inf, N = Inf, S = -Inf,                       # set constraints
  b1 = 1000, b2 = 2000, b3 = 3000,                  # expected benefit for each effect size
  case = 1, strategy = TRUE,                        # chose Case and Strategy
  fixed = TRUE,                                     # true treatment effects are fixed/random
  num_cl = 1)                                       # number of cores for parallelized computing
}
}
\references{
IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at \href{https://www.iqwig.de/ueber-uns/methoden/methodenpapier/}{https://www.iqwig.de/ueber-uns/methoden/methodenpapier/}, assessed last 15.05.19.

Preussler, S., Kieser, M., and Kirchner, M. (2019). Optimal sample size allocation and go/no-go decision rules for phase II/III programs where several phase III trials are performed. Biometrical Journal, 61(2), 357-378.

Schoenfeld, D. (1981). The asymptotic properties of nonparametric tests for comparing survival distributions. Biometrika, 68(1), 316-319.
}
