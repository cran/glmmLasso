glmm_final<-function(y,X,W,k,q_start,Delta_start,s,steps=1000,family,method,overdispersion,phi)
{
N<-length(y)
lin<-ncol(as.matrix(X))
n<-length(k)
Eta<-cbind(X,W)%*%Delta_start
Mu<-as.vector(family$linkinv(Eta))
Sigma<-as.vector(family$variance(Mu))
if(overdispersion)
Sigma<-Sigma*phi
D<-as.vector(family$mu.eta(Eta))
W0_inv<-D*1/Sigma*D


Z_alles<-cbind(X,W)

if(s==1)
{
P1<-c(rep(0,lin),rep((q_start^(-1)),n*s))
P1<-diag(P1)
}else{
P1<-matrix(0,lin+n*s,lin+n*s)
for(jf in 1:n)
P1[(lin+(jf-1)*s+1):(lin+jf*s),(lin+(jf-1)*s+1):(lin+jf*s)]<-chol2inv(chol(q_start))
}

if(overdispersion)
{
E<-diag(N)
M0<-Z_alles%*%chol2inv(chol((t(Z_alles)%*%(Z_alles*W0_inv))+P1))%*%t(Z_alles*W0_inv)
phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(M0)))
}

Delta<-matrix(0,steps,(lin+s*n))
Delta[1,]<-Delta_start

Q<-list()
Q[[1]]<-q_start

score_vec<-rep(0,(lin+s*n))
D<-as.vector(family$mu.eta(Eta))
Delta_r<-rep(0,lin+s*n)

l=1
opt<-steps

score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-P1%*%Delta[1,]
F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1

InvFisher<-chol2inv(chol(F_gross))
Delta_r<-InvFisher%*%score_vec
Delta[1,]<-Delta[1,]+Delta_r

Eta<-Z_alles%*%Delta[1,]

Mu<-as.vector(family$linkinv(Eta))
Sigma<-as.vector(family$variance(Mu))
D<-as.vector(family$mu.eta(Eta))

if (method=="EM")
{
F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
InvFisher<-chol2inv(chol(F_gross))
############################# Q updaten ################
Q1<-InvFisher[(lin+1):(lin+s),(lin+1):(lin+s)]+Delta[1,(lin+1):(lin+s)]%*%t(Delta[1,(lin+1):(lin+s)])
for (i in 2:n)
Q1<-Q1+InvFisher[(lin+(i-1)*s+1):(lin+i*s),(lin+(i-1)*s+1):(lin+i*s)]+Delta[1,(lin+(i-1)*s+1):(lin+i*s)]%*%t(Delta[1,(lin+(i-1)*s+1):(lin+i*s)])
Q1<-1/n*Q1
}else{
Eta_tilde<-Eta+(y-Mu)*1/D

Betadach<-Delta[1,1:lin]

if(s==1)
{
optim.obj<-bobyqa(sqrt(q_start),likelihood_bobyqa,D=D,Sigma=Sigma,X=X,X_aktuell=X,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W,lower = 1e-14, upper=20)
Q1<-as.matrix(optim.obj$par)^2
}else{
q_start_vec<-c(diag(q_start),q_start[lower.tri(q_start)])
up1<-min(20,50*max(q_start_vec))#[(s+1):(s*(s+1)*0.5)]))
upp<-rep(up1,length(q_start_vec))
low<-c(rep(0,s),rep(-up1,0.5*(s^2-s)))
kkk_vec<-c(rep(-1,s),rep(0.5,0.5*(s^2-s)))
optim.obj<-bobyqa(q_start_vec,likelihood,D=D,Sigma=Sigma,X=X,X_aktuell=X,Eta_tilde=Eta_tilde,Betadach=Betadach,W=W,n=n,s=s,k=k,lower=low,upper=upp)
Q1<-matrix(0,s,s)
Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
Q1<-Q1+t(Q1)
diag(Q1)<-(optim.obj$par[1:s])

#### Check for positive definitness ########
      for (ttt in 0:100)
      {
      Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
      Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
      Q_solvetest<-try(solve(Q1))
         if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
         break
      }
}}

Q[[2]]<-Q1

if(overdispersion)
{
Uu<-(E-(Z_alles*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher%*%t(Z_alles*sqrt(D*1/Sigma*D*1/Sigma)))%*%(E-M0)
FinalHat<-E-Uu
phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(FinalHat)))
Sigma<-Sigma*phi
}


###############################################################################################################################################
################################################################### Boost ###################################################################
for (l in 2:steps)
{

  if(s==1)
  {
  P1<-c(rep(0,lin),rep((Q1^(-1)),n*s))
  P1<-diag(P1)
  }else{
  P1<-matrix(0,lin+n*s,lin+n*s)
  for(jf in 1:n)
  P1[(lin+(jf-1)*s+1):(lin+jf*s),(lin+(jf-1)*s+1):(lin+jf*s)]<-chol2inv(chol(Q1))
  }

  score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-P1%*%Delta[l-1,]
  F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1

InvFisher<-chol2inv(chol(F_gross))
Delta_r<-InvFisher%*%score_vec
Delta[l,]<-Delta[l-1,]+Delta_r
Eta<-Z_alles%*%Delta[l,]
Mu<-as.vector(family$linkinv(Eta))
Sigma<-as.vector(family$variance(Mu))
D<-as.vector(family$mu.eta(Eta))

if (method=="EM")
{
F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
InvFisher<-chol2inv(chol(F_gross))
############################# Q update ################
Q1<-InvFisher[(lin+1):(lin+s),(lin+1):(lin+s)]+Delta[l,(lin+1):(lin++s)]%*%t(Delta[l,(lin+1):(lin+s)])
for (i in 2:n)
Q1<-Q1+InvFisher[(lin+(i-1)*s+1):(lin+i*s),(lin+(i-1)*s+1):(lin+i*s)]+Delta[l,(lin+(i-1)*s+1):(lin+i*s)]%*%t(Delta[l,(lin+(i-1)*s+1):(lin+i*s)])

Q1<-1/n*Q1
}else{
Eta_tilde<-Eta+(y-Mu)*1/D

Betadach<-Delta[l,1:lin]

if(s==1)
{
optim.obj<-bobyqa(sqrt(Q1),likelihood_bobyqa,D=D,Sigma=Sigma,X=X,X_aktuell=X,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = 1e-12, upper = 20)
Q1<-as.matrix(optim.obj$par)^2
}else{
Q1_vec<-c(diag(Q1),Q1[lower.tri(Q1)])
optim.obj<-bobyqa(Q1_vec,likelihood,D=D,Sigma=Sigma,X=X,X_aktuell=X,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W,lower=low,upper=upp)

Q1<-matrix(0,s,s)
Q1[lower.tri(Q1)]<-optim.obj$par[(s+1):(s*(s+1)*0.5)]
Q1<-Q1+t(Q1)
diag(Q1)<-(optim.obj$par[1:s])

#### Check for positiv definitness ########
for (ttt in 0:100)
      {
      Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
      Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
       Q_solvetest<-try(solve(Q1))
         if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
         break
      }
}}

Q[[l+1]]<-Q1

if(overdispersion)
{
Uu<-(E-(Z_alles*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher%*%t(Z_alles*sqrt(D*1/Sigma*D*1/Sigma)))%*%Uu
FinalHat<-E-Uu
phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(FinalHat)))
}


kritval<-sqrt(sum((Delta[l-1,]-Delta[l,])^2))/sqrt(sum(Delta[l-1,]^2))
if(kritval<1e-6)
break

if(l>2)
{
kritval2<-sqrt(sum((Delta[l-2,]-Delta[l,])^2))/sqrt(sum(Delta[l-2,]^2))
if(kritval2<1e-6)
break
}}

opt<-l
Deltafinal<-Delta[l,]
Q_final<-Q[[l+1]]

  if(s==1)
  {
  P1<-c(rep(0,lin),rep((Q_final^(-1)),n*s))
  P1<-diag(P1)
  }else{
  P1<-matrix(0,lin+n*s,lin+n*s)
  for(jf in 1:n)
  P1[(lin+(jf-1)*s+1):(lin+jf*s),(lin+(jf-1)*s+1):(lin+jf*s)]<-chol2inv(chol(Q_final))
  }
F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
Inv_F_opt<-chol2inv(chol(F_gross))

Standard_errors<-sqrt(diag(Inv_F_opt))

ret.obj=list()
ret.obj$opt<-opt
ret.obj$Delta<-Deltafinal
ret.obj$Q<-Q_final
ret.obj$Standard_errors<-Standard_errors
ret.obj$phi<-phi
return(ret.obj)
}
