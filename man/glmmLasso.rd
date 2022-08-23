\name{glmmLasso}
\alias{glmmLasso}
\docType{package}
\title{
Variable Selection for Generalized Linear Mixed Models by L1-Penalized Estimation.}
\description{
A variable selection approach for generalized linear mixed models by L1-penalized estimation is provided.
}
\details{
The \code{glmmLasso} algorithm is a gradient ascent algorithm designed for generalized linear mixed models, which incorporates variable selection by L1-penalized estimation. In a final re-estimation step a model the includes only the variables corresponding to the non-zero fixed effects is fitted by simple Fisher scoring. For both the main algorithm as well as for the final re-estimation Fisher scoring 
two methods for the computation of the random-effects variance-covariance parameter estimates can be chosen, an EM-type estimate and an REML-type estimate.

\tabular{ll}{
Package: \tab glmmLasso\cr
Type: \tab Package\cr
Version: \tab 1.6.1\cr
Date: \tab 2022-10-05\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
for loading a dataset type data(nameofdataset)
}

\usage{
glmmLasso(fix=formula, rnd=formula, data, lambda, family = gaussian(link="identity"), 
          switch.NR=FALSE, final.re=FALSE, control = list())
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
    the right; aternatively, the random effects design matrix can be given directly (with suitable column names). If set to NULL, no random effects are included.}
  \item{data}{the data frame containing the variables named in
    \code{formula}.}
  \item{lambda}{the penalty parameter that controls the shrinkage of fixed terms and controls the variable selection.
  The optimal penalty parameter is a tuning parameter of the procedure that has to be determined, 
 e.g. by use of information criteria or cross validation. (See details or the quick demo for an example.)} 
  \item{family}{
    a GLM family, see \code{\link[stats]{glm}} and
    \code{\link[stats]{family}}. Also ordinal response models can be fitted: use \code{family=\link{acat}()} and \code{family=\link{cumulative}()} for the fitting of an adjacent category or cumulative model, respectively. If \code{family} is missing then a
    linear mixed model is fit; otherwise a generalized linear mixed
    model is fit.}
  \item{switch.NR}{logical. Should the algorithm swith to a Newton-Raphson update step, when reasonable? Default is FALSE.}
  \item{final.re}{logical. Should the final Fisher scoring re-estimation be performed? Default is FALSE.}
  \item{control}{a list of control values for the estimation algorithm to replace the default values returned by the function \code{\link{glmmLassoControl}}. Defaults to an empty list.}
}
\value{Generic functions such as \code{print}, \code{predict}, \code{plot} and \code{summary} have methods to show the results of the fit.
The \code{predict} function returns predictions on the scale of the response variable and uses also estimates of random effects for prediction, if possible (i.e. for known subjects of the grouping factor). The \code{plot} function 
plots the smooth terms, if any have been specified. 

   \item{call}{a list containing an image of the \code{glmmLasso} call that produced the object.}  
  \item{coefficients}{a vector containing the estimated fixed effects. By default the covariates are standardized/centered within the procedure (see \code{\link{glmmLassoControl}}), so the coefficients are transformed back to the original scale.}
  \item{smooth}{a vector containing the estimated spline coefficients, if smooth terms have been specified.}  
  \item{ranef}{a vector containing the estimated random effects.}
  \item{StdDev}{a scalar or matrix containing the estimates of the random effects standard deviation or variance-covariance parameters, respectively.}
  \item{fitted.values}{a vector of fitted values.}
  \item{phi}{estimated scale parameter, if \code{overdispersion=TRUE} is used. Otherwise, it is equal to one.}  
  \item{Deltamatrix}{a matrix containing the estimates of fixed and random effects (columns) for each iteration (rows) of the main algorithm (i.e. before the final re-estimation step is performed, see details).}
  \item{Q_long}{a list containing the estimates of the random effects variance-covariance parameters for each iteration of the main algorithm.}
  \item{fixerror}{a vector with standrad errors for the fixed effects.}
  \item{ranerror}{a vector with standrad errors for the random effects.}
  \item{aic}{AIC: The negative of twice the log-likelihood plus twice the corresponding degrees of freedom. The corresponding degrees of freedom are determined by the sum of nonzero coefficients corresponding to fixed
effects plus the number of random effects covariance parameters that have to be  estimated.}
  \item{bic}{BIC: The negative of twice the log-likelihood plus the product of the logarithm of the overall number of observations
  with the corresponding degrees of freedom. The corresponding degrees of freedom are determined by the sum of nonzero coefficients corresponding to fixed
effects plus the number of random effects covariance parameters that have to be  estimated.}
  \item{conv.step}{number of iterations until the main algorithm has converged.}
}



\author{
Andreas Groll  \email{groll@statistik.tu-dortmund.de}
}

\references{
Groll, A. and G. Tutz (2014). 
Variable selection for generalized linear mixed models by
L1-penalized estimation. \emph{Statistics and Computing} 24(2), 137--154.

Goeman, J. J. (2010). L1 Penalized Estimation in the Cox Proportional Hazards Model.
\emph{Biometrical Journal} 52, 70--84.
}


\seealso{
\code{\link{glmmLassoControl},\link{soccer},\link{knee}}
}
\examples{
\dontrun{
data("soccer")

soccer[,c(4,5,9:16)]<-scale(soccer[,c(4,5,9:16)],center=TRUE,scale=TRUE)
soccer<-data.frame(soccer)

## linear mixed model
lm1 <- glmmLasso(points ~ transfer.spendings + ave.unfair.score 
       + ball.possession + tackles 
       + ave.attend + sold.out, rnd = list(team=~1), 
       lambda=50, data = soccer)
      
summary(lm1)

## similar linear model without random effects
lm1b <- glmmLasso(points ~ transfer.spendings + ave.unfair.score 
       + ball.possession + tackles 
       + ave.attend + sold.out, rnd = NULL, 
       lambda=50, data = soccer)
      
summary(lm1b)


## linear mixed model with slope on ave.attend;  
## the coefficient of ave.attend is not penalized;
lm2 <- glmmLasso(points~transfer.spendings + ave.unfair.score 
      + ball.possession + tackles + ave.attend 
       + sold.out, rnd = list(team=~1 + ave.attend), lambda=10, 
       data = soccer, control = list(index=c(1,2,3,4,NA,5), 
       method="REML",print.iter=TRUE))

summary(lm2)

## linear mixed model with categorical covariates
## and final Fisher scoring re-estimation step
lm3 <- glmmLasso(points ~ transfer.spendings + as.factor(red.card)  
       + as.factor(yellow.red.card) + ball.possession 
       + tackles + ave.attend + sold.out, rnd = list(team=~1), 
       data = soccer, lambda=100, final.re=TRUE,
       control = list(print.iter=TRUE,print.iter.final=TRUE))

summary(lm3)

## generalized linear mixed model
## with starting values
glm1 <- glmmLasso(points~transfer.spendings  
        + ave.unfair.score + sold.out 
        + tackles + ave.attend + ball.possession, rnd = list(team=~1),  
        family = poisson(link = log), data = soccer, lambda=100, 
        control = list(print.iter=TRUE,start=c(1,rep(0,29)),q_start=0.7)) 

summary(glm1)

## generalized linear mixed model with a smooth term
glm2 <- glmmLasso(points~ + ave.unfair.score + ave.attend 
        + ball.possession + tackles  + sold.out, 
        rnd = list(team=~1),  family = poisson(link = log), 
        data = soccer, lambda=100, control = list(smooth=
        list(formula=~-1 + transfer.spendings, nbasis=7, 
        spline.degree=3, diff.ord=2, penal=5), 
        print.iter=TRUE)) 
 
summary(glm2)
 
plot(glm2,plot.data=FALSE)        

#####################
#####################
#####################

data(knee)

knee[,c(2,4:6)]<-scale(knee[,c(2,4:6)],center=TRUE,scale=TRUE)
knee<-data.frame(knee)

## fit cumulative model
glm3 <- glmmLasso(pain ~ time + th + age + sex, rnd = NULL,  
        family = cumulative(), data = knee, lambda=10,
        control=list(print.iter=TRUE)) 

summary(glm3)

## fit adjacent category model
glm4 <- glmmLasso(pain ~ time + th + age + sex, rnd = NULL,  
        family = acat(), data = knee, lambda=10,
        control=list(print.iter=TRUE)) 

summary(glm4)


# see also demo("glmmLasso-soccer")
}}
\concept{Lasso}
\concept{Shrinkage}
\concept{Variable selection}
\concept{Generalized linear mixed model}


