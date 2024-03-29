% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_multitrial.R
\name{EPsProg23}
\alias{EPsProg23}
\title{Expected probability of a successful program deciding between two or three phase III trials in a time-to-event setting}
\usage{
EPsProg23(HRgo, d2, alpha, beta, w, hr1, hr2, id1, id2, case, size, ymin)
}
\arguments{
\item{HRgo}{threshold value for the go/no-go decision rule}

\item{d2}{total number of events in phase II}

\item{alpha}{significance level}

\item{beta}{\code{1-beta} power for calculation of sample size for phase III}

\item{w}{weight for mixture prior distribution}

\item{hr1}{first assumed true treatment effect on HR scale for prior distribution}

\item{hr2}{second assumed true treatment effect on HR scale for prior distribution}

\item{id1}{amount of information for \code{hr1} in terms of number of events}

\item{id2}{amount of information for \code{hr2} in terms of number of events}

\item{case}{choose case: "at least 1, 2 or 3 significant trials needed for approval"}

\item{size}{size category \code{"small"}, \code{"medium"} or \code{"large"}}

\item{ymin}{assumed minimal clinical relevant effect}
}
\value{
The output of the function \code{EPsProg23()} is the expected probability of a successful program.
}
\description{
The function \code{EPsProg23()} calculates the expected probability of a successful program
in a time-to-event setting. This function follows a special decision rule in order to determine
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
res <- EPsProg23(HRgo = 0.8, d2 = 50,  alpha = 0.025, beta = 0.1, 
                                  w = 0.3, hr1 =  0.69, hr2 = 0.81, 
                                  id1 = 280, id2 = 420, case = 2, size = "small",
                                  ymin = 0.5)
}
\keyword{internal}
