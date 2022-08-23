\name{soccer}
\alias{soccer}
\docType{data}
\title{German Bundesliga data for the seasons 2008-2010}
\description{
The \code{soccer} data contains different covariables for the teams that played in the first Germna soccer division, the Bundesliga, in the seasons 2007/2008 until 2009/2010.
}
\usage{data(soccer)}
\format{
  A data frame with 54 observations on the following 16 variables.
  \describe{
    \item{\code{pos}}{the final league rank of a soccer team at the end of the season}
    \item{\code{team}}{soccer teams}
    \item{\code{points}}{number of the points a team has earned during the season}
    \item{\code{transfer.spendings}}{the amount (in Euro) that a team has spent for new players at the start of the season}
    \item{\code{transfer.receits}}{the amount (in Euro) that a team has earned for the selling of players at the start of the season}
    \item{\code{yellow.card}}{number of the yellow cards a team has received during the season}
    \item{\code{yellow.red.card}}{number of the yellow-red cards a team has received during the season}
    \item{\code{red.card}}{number of the red cards a team has received during the season}
    \item{\code{unfair.score}}{unfairness score which is derived by the number of yellow cards (1 unfairness point), yellow-red cards (2 unfairness points) and red cards (3 unfairness points) a team has received during the season}
    \item{\code{ave.unfair.score}}{average unfairness score per match}
    \item{\code{ball.possession}}{average percentage of ball possession per match}
    \item{\code{tackles}}{average percentage of head-to-head duels won per match}
    \item{\code{capacity}}{capacity of the team's soccer stadium}
    \item{\code{total.attend}}{total attendance of a soccer team for the whole season}
    \item{\code{ave.attend}}{average attendance of a soccer team  per match}
    \item{\code{sold.out}}{number of stadium sold outs during a season}
  }
  
}



\references{
Groll, A. and G. Tutz (2012). Regularization for generalized additive mixed models by likelihood-based boosting. \emph{Methods of Information in Medicine.} 51(2), 168-177.

Groll, A. and G. Tutz (2014). 
Variable selection for generalized linear mixed models by
L1-penalized estimation. \emph{Statistics and Computing} 24(2), 137--154.

We are grateful to Jasmin Abedieh for providing the German Bundesliga data, which were part of her bachelor thesis.
}
\seealso{
  \code{\link{glmmLasso}}, \code{\link[GMMBoost]{bGLMM}}, \code{\link[GMMBoost]{bGAMM}}.
}
\keyword{datasets}
