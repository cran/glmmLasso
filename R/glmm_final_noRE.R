glmm_final_noRE<-function(y,X,Delta_start,steps=1000,family,overdispersion,phi,
                     nue=1,print.iter.final=FALSE,eps.final=1e-5)
{
N<-length(y)
lin<-ncol(as.matrix(X))
Eta<-X%*%Delta_start
Mu<-as.vector(family$linkinv(Eta))
Sigma<-as.vector(family$variance(Mu))
if(overdispersion)
Sigma<-Sigma*phi
D<-as.vector(family$mu.eta(Eta))

if(print.iter.final)
  message("Final Re-estimation Iteration ", 1)
#print(paste("Final Re-estimation Iteration ", 1,sep=""))

Z_alles<-X

Delta<-matrix(0,steps,lin)
Eta.ma<-matrix(0,steps+1,N)
Eta.ma[1,]<-Eta

l=1
opt<-steps

score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)
F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)

InvFisher<-try(chol2inv(chol(F_gross)),silent=TRUE)
if(class(InvFisher)=="try-error")
InvFisher<-solve(F_gross)  

half.index<-0
solve.test<-FALSE
Delta_r<-InvFisher%*%score_vec

######### big while loop for testing if the update leads to Fisher matrix which can be inverted
while(!solve.test)
{  

solve.test2<-FALSE  
while(!solve.test2)
{  
Delta[1,]<-Delta_start+nue*(0.5^half.index)*Delta_r

Eta<-Z_alles%*%Delta[1,]

Mu<-as.vector(family$linkinv(Eta))
Sigma<-as.vector(family$variance(Mu))
D<-as.vector(family$mu.eta(Eta))


if (overdispersion)
{  
F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)
InvFisher<-try(chol2inv(chol(F_gross)),silent=TRUE)
  if(class(InvFisher)=="try-error")
  InvFisher<-try(solve(F_gross),silent=TRUE)  
  if(class(InvFisher)=="try-error")
  {
    half.index<-half.index+1  
  }else{
    solve.test2<-TRUE 
}}else{
  solve.test2<-TRUE
}}


if(overdispersion)
{
FinalHat<-(Z_alles*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher%*%t(Z_alles*sqrt(D*1/Sigma*D*1/Sigma))#E-Uu
phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(FinalHat)))
Sigma<-Sigma*phi
}


score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)
F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)

InvFisher<-try(chol2inv(chol(F_gross)),silent=TRUE)

  if(class(InvFisher)=="try-error")
  InvFisher<-try(solve(F_gross),silent=TRUE)  

  if(class(InvFisher)=="try-error")
  {
    half.index<-half.index+1  
  }else{
    solve.test<-TRUE 
  }
}

Eta.ma[2,]<-Eta

y_dach<-as.vector(family$linkinv(Eta))
Dev_neu<-sum(family$dev.resids(y,y_dach,wt=rep(1,N))^2)

###############################################################################################################################################
################################################################### Main Iterations ###################################################################
eps<-eps.final*sqrt(length(Delta_r))

for (l in 2:steps)
{
  
if(print.iter.final)
  message("Final Re-estimation Iteration ", l)
#print(paste("Final Re-estimation Iteration ", l,sep=""))

half.index<-0
solve.test<-FALSE

Delta_r<-InvFisher%*%score_vec
######### big while loop for testing if the update leads to Fisher matrix which can be inverted
first.time<-FALSE
while(!solve.test)
{  

solve.test2<-FALSE  
while(!solve.test2)
{  

if(half.index>50)
 half.index<-Inf

Delta[l,]<-Delta[l-1,]+nue*(0.5^half.index)*Delta_r

Eta<-Z_alles%*%Delta[l,]
Mu<-as.vector(family$linkinv(Eta))
Sigma<-as.vector(family$variance(Mu))
D<-as.vector(family$mu.eta(Eta))

if (overdispersion)
{  
  F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)
  InvFisher<-try(chol2inv(chol(F_gross)),silent=TRUE)
  if(class(InvFisher)=="try-error")
    InvFisher<-try(solve(F_gross),silent=TRUE)  
  if(class(InvFisher)=="try-error")
  {
    half.index<-half.index+1  
  }else{
    if(!first.time)
    half.index.final<-half.index
    solve.test2<-TRUE 
    first.time<-TRUE
  }}else{
    if(!first.time)
    half.index.final<-half.index
    solve.test2<-TRUE 
    first.time<-TRUE
  }}


if(overdispersion)
{
FinalHat<-(Z_alles*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher%*%t(Z_alles*sqrt(D*1/Sigma*D*1/Sigma))#E-Uu
phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(FinalHat)))
Sigma<-Sigma*phi
}



score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)
F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)

InvFisher<-try(chol2inv(chol(F_gross)),silent=TRUE)
  if(class(InvFisher)=="try-error")
  InvFisher<-try(solve(F_gross),silent=TRUE)  
  if(class(InvFisher)=="try-error")
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

#browser()

opt<-l


if(!solve.test)
{
Delta[l,]<-Delta[l-1,]+nue*(0.5^half.index.final)*Delta_r

Eta<-Z_alles%*%Delta[l,]
Mu<-as.vector(family$linkinv(Eta))
Sigma<-as.vector(family$variance(Mu))
D<-as.vector(family$mu.eta(Eta))

if(overdispersion)
{  
  F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)
  InvFisher<-try(chol2inv(chol(F_gross)),silent=TRUE)
  if(class(InvFisher)=="try-error")
    InvFisher<-try(solve(F_gross),silent=TRUE)  
}

  FinalHat<-(Z_alles*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher%*%t(Z_alles*sqrt(D*1/Sigma*D*1/Sigma))
  complexity<-sum(diag(FinalHat))

if(overdispersion)
{
  phi<-(sum((y-Mu)^2/Mu))/(N-complexity)
  Sigma<-Sigma*phi
}

}

FinalHat<-(Z_alles*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher%*%t(Z_alles*sqrt(D*1/Sigma*D*1/Sigma))
complexity<-sum(diag(FinalHat))

Deltafinal<-Delta[l,]
Standard_errors<-sqrt(diag(InvFisher))

#FinalHat<-(Z_alles*sqrt(D*1/Sigma*D))%*%Inv_F_opt%*%t(Z_alles*sqrt(D*1/Sigma*D))

ret.obj<-list()
ret.obj$opt<-opt
ret.obj$Delta<-Deltafinal
ret.obj$Standard_errors<-Standard_errors
ret.obj$phi<-phi
ret.obj$complexity<-complexity
return(ret.obj)
}
