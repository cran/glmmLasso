\name{cumulative}
\alias{cumulative}

\title{Family Object for Ordinal Regression with Cumulative Probabilities }
\description{
  Provides necessary family components to fit a proportional odds regression model to an ordered response
  based on the corresponding (multivariate) binary design representation.
}
\usage{
cumulative()
}

\value{

  \item{linkinv}{function: the inverse of the link function as a function of eta. Solely the logit link is implemented, hence, the response function \eqn{h(\eta)=exp(\eta)/(1+exp(\eta))} is used.}
  
  \item{deriv.mat}{function: derivative function as a function of the mean (not of eta as normally).}

  \item{SigmaInv}{function: the inverse of the variance as a function of the mean.}
  
  \item{family}{character: the family name.}
  
  \item{multivariate}{Logical. Is always set to TRUE if the family is used.}

}
\details{For a response variable \eqn{Y} 
  with ordered values \eqn{1,2,\ldots,M+1} the design of the corresponding (multivariate) binary response 
  representation is automatically created by the \link{glmmLasso} function. The result is  
  a linear predictor matrix \eqn{\eta} with \eqn{n} rows and \eqn{M} columns. 
  
  Based on this \eqn{(n x M)} predictor matrix \eqn{\eta} or on the 
  corresponding \eqn{(n x M)} matrix \eqn{\mu} the below mentioned family components
  can be calculated.
  
  Solely the logit link is implemented, hence, a proportional odds model is fitted.
}

\references{
Agresti, A. (2013)
\emph{Categorical Data Analysis},
3rd ed. Hoboken, NJ, USA: Wiley.

Dobson, A. J. and Barnett, A. (2008)
\emph{An Introduction to Generalized Linear Models},
3rd ed. Boca Raton: Chapman & Hall/CRC Press.


McCullagh, P. and Nelder, J. A. (1989)
\emph{Generalized Linear Models}, 2nd ed. London: Chapman & Hall.

Simonoff, J. S. (2003)
\emph{Analyzing Categorical Data},
New York: Springer-Verlag.

Tutz, G. (2012)
\emph{Regression for Categorical Data},
Cambridge University Press. 


Yee, T. W. and Wild, C. J. (1996)
Vector generalized additive models.
\emph{Journal of the Royal Statistical Society, Series B, Methodological},
\bold{58}, 481--493.

}
\author{
Andreas Groll  \email{groll@math.lmu.de}
}

\seealso{
    \code{\link{acat}},
    \code{\link{glmmLasso}},
    \link{knee}
}

\examples{
\dontrun{
data(knee)

knee[,c(2,4:6)]<-scale(knee[,c(2,4:6)],center=TRUE,scale=TRUE)
knee<-data.frame(knee)

## fit adjacent category model
glm.obj <- glmmLasso(pain ~ time + th + age + sex, rnd = NULL,  
        family = cumulative(), data = knee, lambda=10,
        switch.NR=TRUE, control=list(print.iter=TRUE)) 

summary(glm.obj)
}}
\keyword{ordinal regression}

