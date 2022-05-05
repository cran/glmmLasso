\name{glmmLassoControl}
\alias{glmmLassoControl}
\concept{glmmLassoControl}
\title{Control Values for \code{glmmLasso} fit}
\description{
  The values supplied in the function call replace the defaults and a list with all possible arguments is returned. The returned list is used as the \code{control} argument to the \code{glmmLasso} function.
}

\usage{

glmmLassoControl(nue=1,index=NULL,smooth=NULL, start=NULL, q_start=NULL, 
                 center = TRUE, standardize = TRUE, steps=1000, 
                 method="EM", overdispersion=FALSE,     
                 epsilon=1e-4, maxIter=200, print.iter=FALSE, 
                 print.iter.final=FALSE, method.final="EM", 
                 eps.final=1e-4, Q.fac=5, complexity="hat.matrix",...)
} 
    
\arguments{
  \item{nue}{weakness of the learner. Choose 0 < nue =< 1. Default is 1.}  
  \item{index}{vector which defines the grouping of the variables. Components sharing the same number build a group and factor variables get a single number (and are automatically treated as a group). Non-penalized coefficients are marked with NA.}
  \item{smooth}{a list specifying the formula of the smooth terms, together with the number of basis functions \code{nbasis}, the degree of the B-splines \code{spline.degree}, the order of differences that is used for penalization \code{diff.ord} and finally a correspodning penalty parameter \code{penal}.}
  \item{start}{a vector containing starting values for fixed and random effects of suitable length. Default is a vector full of zeros.}
  \item{q_start}{a scalar or matrix of suitable dimension, specifying starting values for the random-effects variance-covariance matrix. Default is a scalar 0.1 or diagonal matrix with 0.1 in the diagonal, depending on the dimension of the random effects.}
  \item{center}{logical. If true, the columns of the design matrix will be
    centered (except a possible intercept column).}
  \item{standardize}{logical. If true, the design matrix will be
    blockwise orthonormalized such that for each block \eqn{X^TX = n 1}
         (*after* possible centering).}
  \item{steps}{the number of iterations. Default is 1000.}
  \item{method}{two methods for the computation of the random-effects variance-covariance parameter estimates can be chosen, an EM-type estimate and an REML-type estimate. The REML-type estimate uses the \code{nlminb} or the \code{bobyqa} function for optimization, depending on the dimension of the random effects. Default is \code{EM}.}
 \item{overdispersion}{logical scalar. If \code{FALSE}, no scale parameter is derived, if \code{TRUE}, in each             iteration a scale parameter is estimated by use of Pearson residuals. 
       This can be used e.g. to fit overdispersed Poisson models. Default is \code{FALSE}.
       If the Gaussian family is used, overdispersion is automatically set \code{TRUE}.}
  \item{epsilon}{controls the speed of convergence. Default is 1e-4.}
  \item{maxIter}{the number of iterations for the final Fisher scoring re-estimation procedure. Default is 200.}
    \item{print.iter}{logical. Should the number of iterations be printed? Default is FALSE.}
    \item{print.iter.final}{logical. Should the number of iterations in the final re-estimation step be printed? Default is FALSE.}
      \item{method.final}{two methods for the computation of the random-effects variance-covariance parameter estimates  
                    for the final Fisher scoring re-estimation procedure  can be chosen, an EM-type estimate and an REML-type estimate. The REML-type estimate uses the \code{bobyqa} function for optimization.
                Default is \code{EM}.}
\item{eps.final}{controls the speed of convergence in the final re-estimation. Default is 1e-4.}
\item{Q.fac}{Factor which controls the interval on which is searched for the optimal parameters of the random-effects variance-covariance matrix, if method.final="REML". Default is 5.}
\item{complexity}{Character which determines how the model complexity is computed. Default is "hat.matrix", which sums up the trace of the corresponding hat matrix. Alternatively, simply the number of estimated (non-zero) parameters can be used by setting complexity="non-zero".}
\item{...}{Futher arguments to be passed.}
}

\value{
  a list with components for each of the possible arguments.
}

\author{
Andreas Groll \email{groll@math.lmu.de}
}

\seealso{
  \code{\link{glmmLasso}}, \code{\link[minqa]{bobyqa}}
}

\examples{
# Use REML estimates for random effects covariance parameters
# and lighten the convergence criterion 
glmmLassoControl(method="REML", epsilon=1e-4)
}
