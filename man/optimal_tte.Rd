% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimal_tte.R
\name{optimal_tte}
\alias{optimal_tte}
\title{Optimal phase II/III drug development planning with time-to-event endpoint}
\format{
data.frame containing the optimization results (see Value)
}
\usage{
optimal_tte(
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
  steps1 = 1,
  stepm1 = 0.95,
  stepl1 = 0.85,
  b1,
  b2,
  b3,
  gamma = 0,
  fixed = FALSE,
  skipII = FALSE,
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

\item{steps1}{lower boundary for effect size category "small" in HR scale, default: 1}

\item{stepm1}{lower boundary for effect size category "medium" in HR scale = upper boundary for effect size category "small" in HR scale, default: 0.95}

\item{stepl1}{lower boundary for effect size category "large" in HR scale = upper boundary for effect size category "medium" in HR scale, default: 0.85}

\item{b1}{expected gain for effect size category "small"}

\item{b2}{expected gain for effect size category "medium"}

\item{b3}{expected gain for effect size category "large"}

\item{gamma}{to model different populations in phase II and III choose \code{gamma != 0}, default: 0}

\item{fixed}{choose if true treatment effects are fixed or random, if TRUE hr1 is used as a fixed effect and hr2 is ignored}

\item{skipII}{choose if skipping phase II is an option, default: FALSE;
if TRUE, the program calculates the expected utility for the case when phase
II is skipped and compares it to the situation when phase II is not skipped.
The results are then returned as a two-row data frame, \code{res[1, ]}
being the results when including phase II and \code{res[2, ]} when skipping phase II.
\code{res[2, ]} has an additional parameter, \code{res[2, ]$median_prior_HR}, which is
the assumed hazards ratio used for planning the phase III study when the
phase II is skipped. It is calculated as the exponential function of the
median of the prior function.}

\item{num_cl}{number of clusters used for parallel computing, default: 1}
}
\value{
The output of the function is a \code{data.frame} object containing the optimization results:
\describe{\item{u}{maximal expected utility under the optimization constraints, i.e. the expected utility of the optimal sample size and threshold value}\item{HRgo}{optimal threshold value for the decision rule to go to phase III}
\item{d2}{optimal total number of events for phase II}
\item{d3}{total expected number of events for phase III; rounded to next natural number}
\item{d}{total expected number of events in the program; d = d2 + d3}\item{n2}{total sample size for phase II; rounded to the next even natural number}
\item{n3}{total sample size for phase III; rounded to the next even natural number}
\item{n}{total sample size in the program; n = n2 + n3}
\item{K}{maximal costs of the program (i.e. the cost constraint, if it is set or the sum K2+K3 if no cost constraint is set)}
\item{pgo}{probability to go to phase III}
\item{sProg}{probability of a successful program}\item{sProg1}{probability of a successful program with "small" treatment effect in phase III}
\item{sProg2}{probability of a successful program with "medium" treatment effect in phase III}
\item{sProg3}{probability of a successful program with "large" treatment effect in phase III}\item{K2}{expected costs for phase II}
\item{K3}{expected costs for phase III}}
and further input parameters. Taking \code{cat(comment())} of the
data frame lists the used optimization sequences, start and
finish date of the optimization procedure.
}
\description{
The function \code{\link{optimal_tte}} of the \code{\link{drugdevelopR}}
package enables planning of phase II/III drug development programs with optimal
sample size allocation and go/no-go decision rules for time-to-event endpoints
(Kirchner et al., 2016). The assumed true treatment effects can be assumed to
be fixed or modelled by
a prior distribution. When assuming fixed true treatment effects, planning can
also be done with the user-friendly R Shiny app
\href{https://web.imbi.uni-heidelberg.de/basic/}{basic}.
The app \href{https://web.imbi.uni-heidelberg.de/prior/}{prior} visualizes
the prior distributions used in this package. Fast computing is enabled by
parallel programming.
}
\examples{
# Activate progress bar (optional)
\dontrun{
progressr::handlers(global = TRUE)
}
# Optimize
\donttest{
optimal_tte(w = 0.3,                    # define parameters for prior
  hr1 = 0.69, hr2 = 0.88, id1 = 210, id2 = 420,   # (https://web.imbi.uni-heidelberg.de/prior/)
  d2min = 20, d2max = 100, stepd2 = 5,            # define optimization set for d2
  hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05,  # define optimization set for HRgo
  alpha = 0.025, beta = 0.1, xi2 = 0.7, xi3 = 0.7, # drug development planning parameters
  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,        # fixed/variable costs for phase II/III
  K = Inf, N = Inf, S = -Inf,                     # set constraints
  steps1 = 1,                                     # define lower boundary for "small"
  stepm1 = 0.95,                                  # "medium"
  stepl1 = 0.85,                                  # and "large" treatment effect size categories
  b1 = 1000, b2 = 2000, b3 = 3000,                # expected benefit for each effect size category
  gamma = 0,                                      # population structures in phase II/III
  fixed = FALSE,                                  # true treatment effects are fixed/random
  skipII = FALSE,                                 # skipping phase II 
  num_cl = 1)                                     # number of cores for parallelized computing 
}
}
\references{
Kirchner, M., Kieser, M., Goette, H., & Schueler, A. (2016). Utility-based optimization of phase II/III programs. Statistics in Medicine, 35(2), 305-316.

IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at \href{https://www.iqwig.de/ueber-uns/methoden/methodenpapier/}{https://www.iqwig.de/ueber-uns/methoden/methodenpapier/}, assessed last 15.05.19.

Schoenfeld, D. (1981). The asymptotic properties of nonparametric tests for comparing survival distributions. Biometrika, 68(1), 316-319.
}
\seealso{
\code{\link{optimal_binary}}, \code{\link{optimal_normal}}, \code{\link{optimal_bias}}, \code{\link{optimal_multitrial}} and \code{\link{optimal_multiarm}}
}
