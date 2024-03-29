% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_multitrial.R
\name{EPsProg_multitrial}
\alias{EPsProg_multitrial}
\alias{EPsProg2}
\alias{EPsProg3}
\alias{EPsProg4}
\title{Expected probability of a successful program for multitrial programs in a time-to-event setting}
\usage{
EPsProg2(HRgo, d2, alpha, beta, w, hr1, hr2, id1, id2, case, size, fixed)

EPsProg3(HRgo, d2, alpha, beta, w, hr1, hr2, id1, id2, case, size, fixed)

EPsProg4(HRgo, d2, alpha, beta, w, hr1, hr2, id1, id2, case, size, fixed)
}
\arguments{
\item{HRgo}{threshold value for the go/no-go decision rule}

\item{d2}{total number of events in phase II}

\item{alpha}{significance level}

\item{beta}{1-beta power for calculation of sample size for phase III}

\item{w}{weight for mixture prior distribution}

\item{hr1}{first assumed true treatment effect on HR scale for prior distribution}

\item{hr2}{second assumed true treatment effect on HR scale for prior distribution}

\item{id1}{amount of information for hr1 in terms of number of events}

\item{id2}{amount of information for hr2 in terms of number of events}

\item{case}{choose case: "at least 1, 2 or 3 significant trials needed for approval"}

\item{size}{size category "small", "medium" or "large"}

\item{fixed}{choose if true treatment effects are fixed or random}
}
\value{
The output of the function EPsProg2(), EPsProg3() and EPsProg4() is the expected probability of a successful program when performing several phase III trials (2, 3 or 4 respectively)
}
\description{
These functions calculate the expected probability of a successful program given the parameters.
Each function represents a specific strategy, e.g. the function EpsProg3() calculates the expected probability if three phase III trials are performed.
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
\donttest{EPsProg2(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, 
                                 w = 0.3, hr1 =  0.69, hr2 = 0.81, 
                                 id1 = 210, id2 = 420, case = 2, size = "small",
                                 fixed = FALSE)}
          \donttest{EPsProg3(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, 
                                 w = 0.3, hr1 =  0.69, hr2 = 0.81, 
                                 id1 = 210, id2 = 420, case = 2, size = "small",
                                 fixed = TRUE)}
          \donttest{EPsProg4(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, 
                                 w = 0.3, hr1 =  0.69, hr2 = 0.81, 
                                 id1 = 210, id2 = 420, case = 3, size = "small",
                                 fixed = TRUE)}
}
\keyword{internal}
