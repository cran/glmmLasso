\name{glmmLassoControl}
\alias{glmmLassoControl}
\concept{glmmLassoControl}
\title{Control Values for \code{glmmLasso} fit}
\description{
  The values supplied in the function call replace the defaults and a list with all possible arguments is returned. The returned list is used as the \code{control} argument to the \code{bGLMM} function.
}

\usage{

glmmLassoControl(nue=1,index=NULL,smooth=NULL, start=NULL, q_start=NULL, 
                 phi_start=1, steps=1000, method="EM", overdispersion=FALSE,     
                 epsilon=1e-5, maxIter=200, print.iter=FALSE, 
                 print.iter.final=FALSE, method.final="EM", 
                 eps.final=1e-5, Q.fac=5)
} 
    
\arguments{
  \item{nue}{weakness of the learner. Choose 0 < nue =< 1. Default is 1.}  
  \item{index}{vector which defines the grouping of the variables. Components sharing the same number build a group and factor variables get a single number (and are automatically treated as a group). Non-penalized coefficients are marked with NA.}
  \item{smooth}{a list specifying the formula of the smooth terms, together with the number of basis functions \code{nbasis}, the degree of the B-splines \code{spline.degree}, the order of differences that is used for penalization \code{diff.ord} and finally a correspodning penalty parameter \code{penal}.}
  \item{start}{a vector containing starting values for fixed and random effects of suitable length. Default is a vector full of zeros.}
  \item{q_start}{a scalar or matrix of suitable dimension, specifying starting values for the random-effects variance-covariance matrix. Default is a scalar 0.1 or diagonal matrix with 0.1 in the diagonal, depending on the dimension of the random effects.}
  \item{phi_start}{a scalar specifying the starting value of the scale parameter. Default is 1.}
  \item{steps}{the number of iterations. Default is 1000.}
  \item{method}{two methods for the computation of the random-effects variance-covariance parameter estimates can be chosen, an EM-type estimate and an REML-type estimate. The REML-type estimate uses the \code{nlminb} or the \code{bobyqa} function for optimization, depending on the dimension of the random effects. Default is \code{EM}.}
 \item{overdispersion}{logical scalar. If \code{FALSE}, no scale parameter is derived, if \code{TRUE}, in each boosting iteration a scale parameter is estimated by use of Pearson residuals. 
       This can be used e.g. to fit overdispersed Poisson models. Default is \code{FALSE}.
       If the Gaussian family is used, overdispersion is automatically set \code{TRUE}.}
  \item{epsilon}{controls the speed of convergence. Default is 1e-5.}
  \item{maxIter}{the number of iterations for the final Fisher scoring reestimation procedure. Default is 200.}
    \item{print.iter}{logical. Should the number of interations be printed? Default is FALSE.}
    \item{print.iter.final}{logical. Should the number of interations in the final re-estimation step be printed? Default is FALSE.}
      \item{method.final}{two methods for the computation of the random-effects variance-covariance parameter estimates  
                    for the final Fisher scoring reestimation procedure  can be chosen, an EM-type estimate and an REML-type estimate. The REML-type estimate uses the \code{bobyqa} function for optimization.
                Default is \code{EM}.}
\item{eps.final}{controls the speed of convergence in the final re-estimation. Default is 1e-5.}
\item{Q.fac}{Factor which controls the interval on which is searched for the optimal parameters of the random-effects variance-covariance matrix, if method.final="REML". Default is 5.}
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
# Use REML estimates for random effects covariance parameters
# and lighten the convergence criterion 
glmmLassoControl(method="REML", epsilon=1e-4)
}
