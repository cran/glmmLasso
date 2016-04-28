\name{knee}
\alias{knee}
\docType{data}
\title{Clinical pain study on knee data}
\description{
The \code{knee} data set illustrates the effect of a medical spray on the pressure pain in the knee due to sports injuries.
}
\usage{data(knee)}
\format{
  A data frame with 381 patients, each with three replicates, and the following 7 variables:
  \describe{
    \item{\code{pain}}{the magnitude of pressure pain in the knee given in 5 categories (1: lowest pain; 5: strongest pain).}
    \item{\code{time}}{the number of replication}
    \item{\code{id}}{number of patient}
    \item{\code{th}}{the therapy (1: spray; 0: placebo)}
    \item{\code{age}}{age of the patient in years}
    \item{\code{sex}}{sex of the patient (1: male; 0: female)}
    \item{\code{pain.start}}{the magnitude of pressure pain in the knee at the beginning of the study}
  }
  
}



\references{
Tutz, G. (2000). \emph{Die Analyse kategorialer Daten - eine anwendungsorientierte Einfuehrung
in Logit-Modellierung und kategoriale Regression}. Muenchen: Oldenbourg Verlag.

Tutz, G. and A. Groll (2012). Likelihood-based boosting in binary and ordinal random effects models. 
Journal of Computational and Graphical Statistics, \bold{22(2)}: 356-378}
\seealso{
  \code{\link[GMMBoost]{OrdinalBoost}}, \code{\link{glmmLasso}}.
}
\keyword{datasets}
