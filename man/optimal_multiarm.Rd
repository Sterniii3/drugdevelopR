% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimal_multiarm.R
\name{optimal_multiarm}
\alias{optimal_multiarm}
\title{Optimal phase II/III drug development planning for multi-arm programs}
\usage{
optimal_multiarm(hr1, hr2, ec, n2min, n2max, stepn2, hrgomin, hrgomax,
  stephrgo, alpha, beta, c2, c3, c02, c03, K = Inf, N = Inf,
  S = -Inf, steps1 = 1, stepm1 = 0.95, stepl1 = 0.85, b1, b2, b3,
  strategy, num_cl = 1)
}
\arguments{
\item{hr1}{assumed true treatment effect on HR scale for treatment 1}

\item{hr2}{assumed true treatment effect on HR scale for treatment 2}

\item{ec}{control arm event rate for phase II and III}

\item{n2min}{minimal total sample size in phase II, must be divisible by 3}

\item{n2max}{maximal total sample size in phase II, must be divisible by 3}

\item{stepn2}{stepsize for the optimization over n2, must be divisible by 3}

\item{hrgomin}{minimal threshold value for the go/no-go decision rule}

\item{hrgomax}{maximal threshold value for the go/no-go decision rule}

\item{stephrgo}{stepsize for the optimization over HRgo}

\item{alpha}{one-sided significance level/ family wise error rate}

\item{beta}{1-beta (any-pair) power for calculation of the number of events for phase III}

\item{c2}{variable per-patient cost for phase II}

\item{c3}{variable per-patient cost for phase III}

\item{c02}{fixed cost for phase II}

\item{c03}{fixed cost for phase III}

\item{K}{constraint on the costs of the program, default: Inf, e.g. no constraint}

\item{N}{constraint on the total expected sample size of the program, default: Inf, e.g. no constraint}

\item{S}{constraint on the expected probability of a successful program, default: -Inf, e.g. no constraint}

\item{steps1}{lower boundary for effect size category "small" in HR scale, default: 1}

\item{stepm1}{lower boundary for effect size category "medium" in HR scale = upper boundary for effect size category "small" in HR scale, default: 0.95}

\item{stepl1}{lower boundary for effect size category "large" in HR scale = upper boundary for effect size category "medium" in HR scale, default: 0.85}

\item{b1}{expected gain for effect size category "small"}

\item{b2}{expected gain for effect size category "medium"}

\item{b3}{expected gain for effect size category "large"}

\item{strategy}{choose Strategy: 1 ("only best promising"), 2 ("all promising") or 3 (both)}

\item{num_cl}{number of clusters used for parallel computing, default: 1}
}
\value{
The output of the function \code{\link{optimal_multiarm}} is a data.frame containing the optimization results:
\describe{
  \item{Strategy}{Strategy, 1: "only best promising" or 2: "all promising"}
  \item{u}{maximal expected utility}
  \item{HRgo}{optimal threshold value for the decision rule to go to phase III}
  \item{n2}{optimal total sample size in phase II}
  \item{n3}{total expected sample size for phase III; rounded to next natural number}
  \item{n}{total sample size in the program; n = n2 + n3}
  \item{K}{maximal costs of the program}
  \item{pgo}{probability to go to phase III}
  \item{sProg}{probability of a successful program}
  \item{sProg1}{probability of a successful program with "small" treatment effect in Phase III}
  \item{sProg2}{probability of a successful program with "medium" treatment effect in Phase III}
  \item{sProg3}{probability of a successful program with "large" treatment effect in Phase III }
  \item{K2}{expected costs for phase II}
  \item{K3}{expected costs for phase III}
  }
and further input parameters.
res
Taking cat(comment()) of the data.frame object lists the used optimization sequences, start and finish date of the optimization procedure.
}
\description{
The function \code{\link{optimal_multiarm}} of the drugdevelopR package enables planning of multi-arm phase II/III drug development programs with optimal sample size allocation and go/no-go decision rules (Preussler et. al, 2019). The assumed true treatment effects are assumed fixed (planning is also possible via user friendly R Shiny App: \href{https://web.imbi.uni-heidelberg.de/multiarm/}{multiarm}). Fast coputing is enabled by parallel programming.
}
\section{drugdevelopR functions}{

The drugdevelopR package provides the functions
\itemize{
  \item \code{\link{optimal_tte}},
  \item \code{\link{optimal_binary}} and
  \item \code{\link{optimal_normal}}
}
to plan optimal phase II/III drug development programs with
\itemize{
  \item time-to-event (treatment effect measured by hazard ratio (HR)),
  \item binary (treatment effect measured by risk ratio (RR)) or
  \item normally distributed (treatment effect measured by standardized difference in means (Delta))
}
endpoint, where the treatment effect is modelled by a \href{https://web.imbi.uni-heidelberg.de/prior/}{prior}. Optimal phase II/III drug development planning with fixed treatment effects can be done with the help of the R Shiny application \href{https://web.imbi.uni-heidelberg.de/basic/}{basic}. Extensions are 
\itemize{
  \item optimal planning of programs including methods for discounting of phase II results (function: \code{\link{optimal_bias}}, App: \href{https://web.imbi.uni-heidelberg.de/bias/}{bias}),
  \item optimal planning of programs with several phase III trials (function: \code{\link{optimal_multitrial}}, App: \href{https://web.imbi.uni-heidelberg.de/multitrial/}{multitrial}) and
  \item optimal planning of programs with multiple arms (function: \code{\link{optimal_multiarm}}, App: \href{https://web.imbi.uni-heidelberg.de/multiarm/}{multiarm}).
}
}

\examples{
res <- optimal_multiarm(hr1 = 0.75, hr2 = 0.80, ec = 0.6,# define assumed true HRs and control arm event rate
  n2min = 30, n2max = 90, stepn2 = 6,                    # define optimization set for n2
  hrgomin = 0.7, hrgomax = 0.9, stephrgo = 0.05,         # define optimization set for HRgo
  alpha = 0.05, beta = 0.1,                              # drug development planning parameters
  c2 = 0.75, c3 = 1, c02 = 100, c03 = 150,               # define fixed and variable costs for phase II and III
  K = Inf, N = Inf, S = -Inf,                            # set maximal costs/ expected sample size for the program or minimal expected probability of a successful program
  steps1 = 1,                                            # define lower boundary for "small"
  stepm1 = 0.95,                                         # "medium"
  stepl1 = 0.85,                                         # and "large" treatment effect size categories as proposed by IQWiG (2016)
  b1 = 1000, b2 = 2000, b3 = 3000,                       # define expected benefit for a "small", "medium" and "large" treatment effect
  strategy = 1,                                          # choose strategy: 1, 2 or 3
  num_cl = 1)                                            # set number of cores used for parallelized computing (check maximum number possible with detectCores())
res
cat(comment(res))                                        # displays the optimization sequence, start and finish date of the optimization procedure.
}
\references{
Preussler, S., Kirchner, M., Goette, H., Kieser, M. (2019). Optimal Designs for Multi-Arm Phase II/III Drug Development Programs. Submitted to peer-review journal.

IQWiG (2016). Allgemeine Methoden. Version 5.0, 10.07.2016, Technical Report. Available at \href{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}{https://www.iqwig.de/de/methoden/methodenpapier.3020.html}, assessed last 15.05.19.
}
