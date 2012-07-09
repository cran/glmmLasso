\name{glmmLassoControl}
\alias{glmmLassoControl}
\concept{glmmLassoControl}
\title{Control Values for \code{glmmLasso} fit}
\description{
  The values supplied in the function call replace the defaults and a list with all possible arguments is returned. The returned list is used as the \code{control} argument to the \code{bGLMM} function.
}

\usage{
glmmLassoControl(nue=1,lin="(Intercept)",start=NULL,q_start=NULL, phi_start = 1, steps=2000,method="EM", 
                 overdispersion=FALSE, epsilon=1e-7, maxIter=1000, print.iter=FALSE, method.final="REML")
} 
    
\arguments{
  \item{nue}{weakness of the learner. Choose 0 < nue =< 1. Default is 1.}  
  \item{lin}{a vector specifying fixed effects, which are excluded from shrinkage.}
  \item{start}{a vector containing starting values for fixed and random effects of suitable length. Default is a vector full of zeros.}
  \item{q_start}{a scalar or matrix of suitable dimension, specifying starting values for the random-effects variance-covariance matrix. Default is a scalar 0.1 or diagonal matrix with 0.1 in the diagonal.}
  \item{phi_start}{a scalar specifying the the starting value of the scale parameter. Default is 1.}
  \item{steps}{the number of interations. Default is 2000.}
  \item{method}{two methods for the computation of the random-effects variance-covariance parameter estimates can be chosen, an EM-type estimate and an REML-type estimate. The REML-type estimate uses the \code{bobyqa} function for optimization.
                Default is \code{EM}.}
 \item{overdispersion}{logical scalar. If \code{FALSE}, no scale parameter is derived, if \code{TRUE}, in each boosting iteration a scale parameter is estimated by use of Pearson residuals. 
       This can be used to fit overdispersed Poisson models. Default is \code{FALSE}.}
  \item{epsilon}{controls the speed of convergence. Default is 1e-7.}
  \item{maxIter}{the number of interations for the final Fisher scoring reestimation procedure. Default is 1000.}
    \item{print.iter}{logical. Should the number of interations be printed?. Default is FALSE.}
      \item{method.final}{two methods for the computation of the random-effects variance-covariance parameter estimates  
                    for the final Fisher scoring reestimation procedure  can be chosen, an EM-type estimate and an REML-type estimate. The REML-type estimate uses the \code{bobyqa} function for optimization.
                Default is \code{REML}.}


}

\value{
  a list with components for each of the possible arguments.
}

\author{
Andreas Groll \email{andreas.groll@stat.uni-muenchen.de}
}

\seealso{
  \code{\link{glmmLasso}}, \code{\link[minqa]{bobyqa}}
}

\examples{
# and REML estimates for random effects covariance parameters
# and lighten the convergence criterion 
glmmLassoControl(method="REML", epsilon=1e-5)
}
