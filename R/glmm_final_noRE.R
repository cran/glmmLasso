glmm_final_noRE<-function(y,X,K,Delta_start,steps=1000,family,overdispersion,phi,
                     nue=1,print.iter.final=FALSE,flushit,eps.final=1e-5)
{
  #browser()
  ## Print stuff.
  ia <- if(flushit) interactive() else FALSE
  
N<-length(y)
lin<-ncol(as.matrix(X))
Eta<-X%*%Delta_start

if(is.null(family$multivariate)){
  D<-family$mu.eta(Eta)
  Mu<-family$linkinv(Eta)
  SigmaInv <- 1/family$variance(Mu)
}else{
  Mu <- family$linkinv(Eta, K)
  D <- family$deriv.mat(Eta, K)
  SigmaInv <- family$SigmaInv(Mu, K)
}


if(print.iter.final)
  #     message()
{
  cat(if(ia) "\r" else NULL)
  cat("\nFinal Re-estimation Iteration  1")
  if(.Platform$OS.type != "unix" & ia) flush.console()
}


Z_alles<-X

Delta<-matrix(0,steps,lin)
Eta.ma<-matrix(0,steps+1,N)
Eta.ma[1,]<-Eta

l=1
opt<-steps

if(is.null(family$multivariate)){
  D <- drop(D);SigmaInv <- drop(SigmaInv)
  score_vec <- t(Z_alles)%*%((y-Mu)*D*SigmaInv)
  F_gross<-t(Z_alles)%*%(Z_alles*D*SigmaInv*D)
}else{
  score_vec<-RcppEigenProd1(Z_alles, D, SigmaInv, y, Mu)
  W_opt <- RcppEigenProd2(D, SigmaInv)
  F_gross <- t(Z_alles)%*%(W_opt%*%Z_alles)
}

InvFisher<-try(chol2inv(chol(F_gross)),silent=TRUE)
if(inherits(InvFisher, "try-error"))
InvFisher<-solve(F_gross)  

half.index<-0
solve.test<-FALSE
Delta_r<-InvFisher%*%score_vec

######### big while loop for testing if the update leads to Fisher matrix which can be inverted
while(!solve.test)
{  
  if(half.index>50)
    half.index<-Inf
  
Delta[1,]<-Delta_start+nue*(0.5^half.index)*Delta_r

Eta<-Z_alles%*%Delta[1,]

if(is.null(family$multivariate)){
  D<-family$mu.eta(Eta)
  Mu<-family$linkinv(Eta)
  SigmaInv <- 1/family$variance(Mu)
}else{
  Mu <- family$linkinv(Eta, K)
  D <- family$deriv.mat(Eta, K)
  SigmaInv <- family$SigmaInv(Mu, K)
}


if(is.null(family$multivariate)){
  D <- drop(D);SigmaInv <- drop(SigmaInv)
  score_vec <- t(Z_alles)%*%((y-Mu)*D*SigmaInv)
  F_gross<-t(Z_alles)%*%(Z_alles*D*SigmaInv*D)
}else{
  score_vec<-RcppEigenProd1(Z_alles, D, SigmaInv, y, Mu)
  W_opt <- RcppEigenProd2(D, SigmaInv)
  F_gross <- t(Z_alles)%*%(W_opt%*%Z_alles)
}


InvFisher<-try(chol2inv(chol(F_gross)),silent=TRUE)

  if(inherits(InvFisher, "try-error"))
  InvFisher<-try(solve(F_gross),silent=TRUE)  

  if(inherits(InvFisher, "try-error"))
  {
    half.index<-half.index+1  
  }else{
    solve.test<-TRUE 
  }
}

Eta.ma[2,]<-Eta

###############################################################################################################################################
################################################################### Main Iterations ###################################################################
eps<-eps.final*sqrt(length(Delta_r))

for (l in 2:steps)
{
  
  if(print.iter.final)
    #  message("Iteration ",l)
  {
    cat(if(ia) "\r" else if(l > 1) "\n" else NULL)
    cat(paste("Final Re-estimation Iteration ",l))
    if(.Platform$OS.type != "unix" & ia) flush.console()
  }
  
  
half.index<-0
solve.test<-FALSE

Delta_r<-InvFisher%*%score_vec
######### big while loop for testing if the update leads to Fisher matrix which can be inverted
while(!solve.test)
{  

if(half.index>50)
 half.index<-Inf

Delta[l,]<-Delta[l-1,]+nue*(0.5^half.index)*Delta_r
Eta<-Z_alles%*%Delta[l,]

if(is.null(family$multivariate)){
  D<-family$mu.eta(Eta)
  Mu<-family$linkinv(Eta)
  SigmaInv <- 1/family$variance(Mu)
}else{
  Mu <- family$linkinv(Eta, K)
  D <- family$deriv.mat(Eta, K)
  SigmaInv <- family$SigmaInv(Mu, K)
}


if(is.null(family$multivariate)){
  D <- drop(D);SigmaInv <- drop(SigmaInv)
  score_vec <- t(Z_alles)%*%((y-Mu)*D*SigmaInv)
  F_gross<-t(Z_alles)%*%(Z_alles*D*SigmaInv*D)
}else{
  score_vec<-RcppEigenProd1(Z_alles, D, SigmaInv, y, Mu)
  W_opt <- RcppEigenProd2(D, SigmaInv)
  F_gross <- t(Z_alles)%*%(W_opt%*%Z_alles)
}


InvFisher<-try(chol2inv(chol(F_gross)),silent=TRUE)
  if(inherits(InvFisher, "try-error"))
  InvFisher<-try(solve(F_gross),silent=TRUE)  
  if(inherits(InvFisher, "try-error"))
  {
    half.index<-half.index+1  
  }else{
    solve.test<-TRUE 
  }
}

Eta.ma[l+1,]<-Eta

finish<-(sqrt(sum((Eta.ma[l,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l,])^2))<eps)
finish2<-(sqrt(sum((Eta.ma[l-1,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l-1,])^2))<eps)

if(finish ||  finish2) 
  break
  
}

######## Final calculation

opt<-l

  if(is.null(family$multivariate)){
    W_opt <- D*SigmaInv*D
    FinalHat<-(Z_alles*sqrt(W_opt))%*%InvFisher%*%t(Z_alles*sqrt(W_opt))
  }else{
    W_inv_t <- RcppEigenSpChol(W_opt)
    FinalHat <- RcppEigenProd3(W_inv_t, Z_alles, InvFisher)
  }

  complexity<-sum(diag(FinalHat))

if(overdispersion)
  phi<-(sum((y-Mu)^2/family$variance(Mu)))/(N-complexity)

Deltafinal<-Delta[l,]
Standard_errors<-InvFisher


ret.obj<-list()
ret.obj$opt<-opt
ret.obj$Delta<-Deltafinal
ret.obj$Standard_errors<-Standard_errors
ret.obj$phi<-phi
ret.obj$complexity<-complexity
return(ret.obj)
}
