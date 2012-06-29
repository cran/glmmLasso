\name{glmmLasso}
\alias{glmmLasso}
\alias{glmmLasso}
\docType{package}
\title{
Variable selection for generalized linear mixed models by
L1-penalized estimation.
}
\description{
Variable selection for generalized linear mixed models by
L1-penalized estimation.
}
\details{
\tabular{ll}{
Package: \tab glmmLasso\cr
Type: \tab Package\cr
Version: \tab 1.0.3\cr
Date: \tab 2012-06-29\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
for loading a dataset type data(nameofdataset)
}

\usage{
glmmLasso(fix=formula, rnd=formula, data, lambda, family = NULL, control = list())
}     
\arguments{
  \item{fix}{a two-sided linear formula object describing the
    fixed-effects part of the model, with the response on the left of a
    \code{~} operator and the terms, separated by \code{+} operators, on
    the right. For categorical covariables use \code{as.factor(.)} in the formula. 
    Note, that the corresponding dummies are treated as a group and are updated blockwise}  
  \item{rnd}{a two-sided linear formula object describing the
    random-effects part of the model, with the grouping factor on the left of a
    \code{~} operator and the random terms, separated by \code{+} operators, on
    the right.}
  \item{data}{the data frame containing the variables named in
    \code{formula}.}
  \item{lambda}{the penalty parameter that controls the shrinkage of fixed terms and controls the variable selection.
  The optimal penalty parameter is a tuning parameter of the procedure that has to be determined, 
 e.g. by use of information criteria or cross validation. (See details or the quick demo for an example.)} 
  \item{family}{
    a GLM family, see \code{\link[stats]{glm}} and
    \code{\link[stats]{family}}. If \code{family} is missing then a
    linear mixed model is fit; otherwise a generalized linear mixed
    model is fit.}
  \item{control}{a list of control values for the estimation algorithm to replace the default values returned by the function \code{bGLMMControl}. Defaults to an empty list.}
}
\value{Generic functions such as \code{print}, \code{predict} and \code{summary} have methods to show the results of the fit.
The \code{predict} function uses also estimates of random effects for prediction, if possible (i.e. for known subjects of the grouping factor).

   \item{call}{a list containing an image of the \code{bGLMM} call that produced the object.}  
  \item{coefficients}{a vector containing the estimated fixed effects}
  \item{ranef}{a vector containing the estimated random effects.}
  \item{StdDev}{a scalar or matrix containing the estimates of the random effects standard deviation or variance-covariance parameters, respectively.}
  \item{fitted.values}{a vector of fitted values.}
  \item{phi}{estimated scale parameter, if \code{overdispersion=TRUE} is used. Otherwise, it is equal to one.}  
  \item{Deltamatrix}{a matrix containing the estimates of fixed and random effects (columns) for each boosting iteration (rows).}
  \item{Q_long}{a list containing the estimates of the random effects standard deviation or variance-covariance parameters, respectively, for each boosting iteration.}
  \item{fixerror}{a vector with standrad errors for the fixed effects.}
  \item{ranerror}{a vector with standrad errors for the random effects.}
  \item{aic}{AIC: The negative of twice the log-likelihood plus twice the corresponding degrees of freedom. The corresponding degrees of freedom are determined by the sum of nonzero coefficients corresponding to fixed
effects plus the number of random effects covariance parameters that have to be  estimated.}
  \item{bic}{BIC: The negative of twice the log-likelihood plus the product of the logarithm of the overall number of observations
  with the corresponding degrees of freedom. The corresponding degrees of freedom are determined by the sum of nonzero coefficients corresponding to fixed
effects plus the number of random effects covariance parameters that have to be  estimated.}
}



\author{
Andreas Groll
}

\references{
Groll, A. and G. Tutz (2011). Variable selection for generalized linear mixed models by
L1-penalized estimation. Technical Report \bold{108}, Ludwig-Maximilians-University.

Goeman, J. J. (2010). L1 Penalized Estimation in the Cox Proportional Hazards Model.
\emph{Biometrical Journal} 52, 70--84.
}


\seealso{
\code{\link{glmmLassoControl},\link{soccer}}
}
\examples{

data("soccer")

## center all metric variables
soccer[,c(4,5,9:16)]<-scale(soccer[,c(4,5,9:16)],center=TRUE,scale=TRUE)
soccer<-data.frame(soccer)

## linear mixed models
lm1 <- glmmLasso(points ~ transfer.spendings + I(transfer.spendings^2)
       + ave.unfair.score + transfer.receits + ball.possession
       + tackles + ave.attend + sold.out, rnd = list(team=~1), data = soccer, lambda=400)
      
lm2 <- glmmLasso(points~transfer.spendings + I(transfer.spendings^2)
       + ave.unfair.score + transfer.receits + ball.possession
       + tackles + ave.attend + sold.out, rnd = list(team=~1 + ave.attend), lambda=500, 
       data = soccer, control = list(steps=100, lin="ave.attend", method="REML"))

## linear mixed models with categorical covariates
lm3 <- glmmLasso(points ~ transfer.spendings + I(transfer.spendings^2)
       + as.factor(red.card) + as.factor(yellow.red.card) 
       + transfer.receits + ball.possession + tackles + ave.attend
       + sold.out, rnd = list(team=~1), data = soccer, lambda=500)


## generalized linear mixed model
glm1 <- glmmLasso(points~transfer.spendings  + I(transfer.spendings^2)
        + ave.unfair.score + transfer.receits + ball.possession
        + tackles + ave.attend + sold.out, rnd = list(team=~1),  
        family = poisson(link = log), data = soccer, lambda=500,
        control = list(start=c(5,rep(0,31)),overdispersion=TRUE))

# see also demo("glmmLasso-soccer")
}
\keyword{
Lasso, Shrinkage, Variable selection, Generalized linear mixed model
}


