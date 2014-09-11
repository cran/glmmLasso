glmm_final_multi_random_smooth<-function(y,X,Phi,W,k,penal.vec,q_start,Delta_start,s,n,steps=1000,
                                         family,method,overdispersion,phi,nue=1,rnd.len,
                                         print.iter.final=FALSE,eps.final=1e-5,
                                         Q.min=1e-13,Q.max=20,Q.fac=5)
{
dim.smooth<-dim(Phi)[2]  
N<-length(y)
lin<-ncol(as.matrix(X))
Eta<-cbind(X,Phi,W)%*%Delta_start
Mu<-as.vector(family$linkinv(Eta))
Sigma<-as.vector(family$variance(Mu))
if(overdispersion)
  Sigma<-Sigma*phi
D<-as.vector(family$mu.eta(Eta))
W0_inv<-D*1/Sigma*D

if(print.iter.final)
  message("Final Re-estimation Iteration ", 1)
#print(paste("Final Re-estimation Iteration ", 1,sep=""))


Z_alles<-cbind(X,Phi,W)

if(all(s==1))
{
P1<-c(rep(0,lin),penal.vec,rep(diag(q_start)^(-1),n))
P1<-diag(P1)
}else{
P1<-matrix(0,lin+dim.smooth+n*s,lin+dim.smooth+n%*%s)
diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
inv.act<-chol2inv(chol(q_start[1:s[1],1:s[1]]))
for(jf in 1:n[1])
P1[(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1]),(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1])]<-inv.act

     for (zu in 2:rnd.len)
     {
     inv.act<-chol2inv(chol(q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
     for(jf in 1:n[zu])
     P1[(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
     (lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-inv.act
     }
}


Delta<-matrix(0,steps,(lin+s%*%n))
Eta.ma<-matrix(0,steps+1,N)
Eta.ma[1,]<-Eta

Q<-list()
Q[[1]]<-q_start

l=1
opt<-steps

score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-P1%*%Delta[1,]
F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1

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

if (method=="EM" || overdispersion)
{  
  F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
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

if (method=="EM")
{
############################# Q updaten ################
   Q1<-matrix(0,sum(s),sum(s))
   Q1[1:s[1],1:s[1]]<-InvFisher[(lin+dim.smooth+1):(lin+dim.smooth+s[1]),(lin+dim.smooth+1):(lin+dim.smooth+s[1])]+Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+s[1])]%*%t(Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+s[1])])
   for (i in 2:n[1])
   Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher[(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1]),(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])]+Delta[1,(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])]%*%t(Delta[1,(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])])
   Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]

     for (zu in 2:rnd.len)
     {
     Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher[(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
     for (i in 2:n[zu])
     Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher[(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
     Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-1/n[zu]*Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
     }

}else{
Eta_tilde<-Eta+(y-Mu)*1/D

Betadach<-Delta[1,1:(lin+dim.smooth)]

   if(all(s==1))
   {
   q_start_vec<-diag(q_start)
   upp<-rep(Q.fac*Q.max,sum(s))
   low<-rep((1/Q.fac)*Q.min,sum(s))
   optim.obj<-try(bobyqa(sqrt(q_start_vec),likelihood_diag,D=D,Sigma=Sigma,X=cbind(X,Phi),X_aktuell=cbind(X,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
   Q1<-diag(optim.obj$par)^2
   }else{
   q_start_vec<-c(diag(q_start)[1:s[1]],q_start[1:s[1],1:s[1]][lower.tri(q_start[1:s[1],1:s[1]])])
   up1<-Q.fac*Q.max
   low<-c(rep(0,s[1]),rep(-up1,0.5*(s[1]^2-s[1])))

     for (zu in 2:rnd.len)
     {
     q_start_vec<-c(q_start_vec,c(diag(q_start)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])],q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]))
     up1<-Q.fac*Q.max
     low<-c(low,c(rep(0,s[zu]),rep(-up1,0.5*(s[zu]^2-s[zu]))))
     }
     upp<-rep(up1,length(q_start_vec))
     optim.obj<-try(bobyqa(q_start_vec,likelihood_block,D=D,Sigma=Sigma,X=cbind(X,Phi),X_aktuell=cbind(X,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
     optim.vec<-optim.obj$par
     
          Q1<-matrix(0,sum(s),sum(s))
     diag(Q1)[1:s[1]]<-optim.vec[1:s[1]]
     if(s[1]>1)
     Q1[1:s[1],1:s[1]][lower.tri(Q1[1:s[1],1:s[1]])]<-optim.vec[(s[1]+1):(s[1]*(s[1]+1)*0.5)]
     optim.vec<-optim.vec[-c(1:(s[1]*(s[1]+1)*0.5))]
     
     for (zu in 2:rnd.len)
     {
     diag(Q1)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-optim.vec[1:s[zu]]
     if(s[zu]>1)
     Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]<-optim.vec[(s[zu]+1):(s[zu]*(s[zu]+1)*0.5)]
     optim.vec<-optim.vec[-c(1:(s[zu]*(s[zu]+1)*0.5))]
     }

     #### Check for positive definitness ########
      for (ttt in 0:100)
      {
      Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
      Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
      Q_solvetest<-try(solve(Q1))
         if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
         break
      }
   }
}

Q[[2]]<-Q1

if(overdispersion)
{  
FinalHat<-(Z_alles*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher%*%t(Z_alles*sqrt(D*1/Sigma*D*1/Sigma))
phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(FinalHat)))
Sigma<-Sigma*phi
}

P1.old.temp<-P1

if(all(s==1))
{
  P1<-c(rep(0,lin),penal.vec,rep(diag(Q1)^(-1),n))
  P1<-diag(P1)
}else{
  P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
  diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
  inv.act<-chol2inv(chol(Q1[1:s[1],1:s[1]]))
  for(jf in 1:n[1])
    P1[(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1]),(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1])]<-inv.act
  
  for (zu in 2:rnd.len)
  {
    inv.act<-chol2inv(chol(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
    for(jf in 1:n[zu])
      P1[(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
         (lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-inv.act
  }
}

score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-P1%*%Delta[1,]
F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1

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


y_dach<-as.vector(family$linkinv(Eta))
Dev_neu<-sum(family$dev.resids(y,y_dach,wt=rep(1,N))^2)

Eta.ma[2,]<-Eta

###############################################################################################################################################
################################################################### Boost ###################################################################
eps<-eps.final*sqrt(length(Delta_r))

for (l in 2:steps)
{

if(print.iter.final)
  message("Final Re-estimation Iteration ", l)
#print(paste("Final Re-estimation Iteration ", l,sep=""))


half.index<-0
solve.test<-FALSE

P1.old<-P1
P1.old2<-P1.old.temp

Delta_r<-InvFisher%*%score_vec
######### big while loop for testing if the update leads to Fisher matrix which can be inverted
while(!solve.test)
{  
  
solve.test2<-FALSE  
while(!solve.test2)
{ 
if(half.index>500)
{
half.index<-Inf;P1.old<-P1.old2
}
Delta[l,]<-Delta[l-1,]+nue*(0.5^half.index)*Delta_r
Eta<-Z_alles%*%Delta[l,]
Mu<-as.vector(family$linkinv(Eta))
Sigma<-as.vector(family$variance(Mu))
D<-as.vector(family$mu.eta(Eta))

if (method=="EM" || overdispersion)
{  
  F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1.old
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

if (method=="EM")
{
############################# Q update ################
   Q1<-matrix(0,sum(s),sum(s))
   Q1[1:s[1],1:s[1]]<-InvFisher[(lin+dim.smooth+1):(lin+dim.smooth+s[1]),(lin+dim.smooth+1):(lin+dim.smooth+s[1])]+Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+s[1])]%*%t(Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+s[1])])
   for (i in 2:n[1])
   Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher[(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1]),(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])]+Delta[l,(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])]%*%t(Delta[1,(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])])
   Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]

     for (zu in 2:rnd.len)
     {
     Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher[(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[l,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
     for (i in 2:n[zu])
     Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher[(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[l,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
     Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-1/n[zu]*Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
     }

}else{
Eta_tilde<-Eta+(y-Mu)*1/D

Betadach<-Delta[l,1:(lin+dim.smooth)]

   if(all(s==1))
   {
   Q1_vec<-diag(Q1)
   up1<-max(max(upp),Q.fac*max(Q1))
   min1<-min(min(upp),(1/Q.fac)*min(Q1))
   upp<-rep(up1,sum(s))
   low<-rep(min1,sum(s))
   optim.obj<-try(bobyqa(sqrt(Q1_vec),likelihood_diag,D=D,Sigma=Sigma,X=cbind(X,Phi),X_aktuell=cbind(X,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
   Q1<-diag(optim.obj$par)^2
   }else{
   Q1_vec<-c(diag(Q1)[1:s[1]],Q1[1:s[1],1:s[1]][lower.tri(Q1[1:s[1],1:s[1]])])
   up1<-max(up1,Q.fac*max(Q1))  
   low<-c(rep(0,s[1]),rep(-up1,0.5*(s[1]^2-s[1])))
   
     for (zu in 2:rnd.len)
     {
     Q1_vec<-c(Q1_vec,c(diag(Q1)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])],Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]))
     up1<-max(up1,Q.fac*max(Q1))  
     low<-c(low,c(rep(0,s[zu]),rep(-up1,0.5*(s[zu]^2-s[zu]))))
     }   
     upp<-rep(up1,length(Q1_vec))
     optim.obj<-try(bobyqa(Q1_vec,likelihood_block,D=D,Sigma=Sigma,X=cbind(X,Phi),X_aktuell=cbind(X,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
     optim.vec<-optim.obj$par
     
     Q1<-matrix(0,sum(s),sum(s))
     diag(Q1)[1:s[1]]<-optim.vec[1:s[1]]
     if(s[1]>1)
     Q1[1:s[1],1:s[1]][lower.tri(Q1[1:s[1],1:s[1]])]<-optim.vec[(s[1]+1):(s[1]*(s[1]+1)*0.5)]
     optim.vec<-optim.vec[-c(1:(s[1]*(s[1]+1)*0.5))]
     
     for (zu in 2:rnd.len)
     {
     diag(Q1)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-optim.vec[1:s[zu]]
     if(s[zu]>1)
     Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]<-optim.vec[(s[zu]+1):(s[zu]*(s[zu]+1)*0.5)]
     optim.vec<-optim.vec[-c(1:(s[zu]*(s[zu]+1)*0.5))]
     }

     #### Check for positive definitness ########
      for (ttt in 0:100)
      {
      Q1[lower.tri(Q1)]<-((0.5)^ttt)*Q1[lower.tri(Q1)]
      Q1[upper.tri(Q1)]<-((0.5)^ttt)*Q1[upper.tri(Q1)]
      Q_solvetest<-try(solve(Q1))
         if(all (eigen(Q1)$values>0) & class(Q_solvetest)!="try-error")
         break
      }
   }

}

Q[[l+1]]<-Q1

y_dach<-as.vector(family$linkinv(Eta))

FinalHat<-(Z_alles*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher%*%t(Z_alles*sqrt(D*1/Sigma*D*1/Sigma))
complexity<-sum(diag(FinalHat))
if(overdispersion)
{
  phi<-(sum((y-Mu)^2/Mu))/(N-complexity)
Sigma<-Sigma*phi
}

P1.old.temp<-P1.old

if(all(s==1))
{
  P1<-c(rep(0,lin),penal.vec,rep(diag(Q1)^(-1),n))
  P1<-diag(P1)
}else{
  P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
  diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
  inv.act<-chol2inv(chol(Q1[1:s[1],1:s[1]]))
  for(jf in 1:n[1])
    P1[(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1]),(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1])]<-inv.act
  
  for (zu in 2:rnd.len)
  {
    inv.act<-chol2inv(chol(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
    for(jf in 1:n[zu])
      P1[(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
         (lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-inv.act
  }
}

score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-P1%*%Delta[l,]
F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1

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

FinalHat<-(Z_alles*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher%*%t(Z_alles*sqrt(D*1/Sigma*D*1/Sigma))
complexity<-sum(diag(FinalHat))


#print(paste("Final Iteration =", l,sep=""))


opt<-l
Deltafinal<-Delta[l,]
Q_final<-Q[[l+1]]

Standard_errors<-sqrt(diag(InvFisher))

## compute ranef part of loglik
if(all(s==1))
{
  P1.ran<-rep(diag(Q_final)^(-1),n)
  P1.ran<-diag(P1.ran)
}else{
  P1.ran<-matrix(0,n%*%s,n%*%s)
  inv.act<-chol2inv(chol(Q_final[1:s[1],1:s[1]]))
  for(jf in 1:n[1])
    P1.ran[( (jf-1)*s[1]+1):( jf*s[1]),( (jf-1)*s[1]+1):( jf*s[1])]<-inv.act
  
  for (zu in 2:rnd.len)
  {
    inv.act<-chol2inv(chol(Q_final[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
    for(jf in 1:n[zu])
      P1.ran[( n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):( n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
             ( n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):( n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-inv.act
  }
}

ranef.logLik<--0.5*t(Deltafinal[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)])%*%P1.ran%*%Deltafinal[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]

#FinalHat<-(Z_alles*sqrt(D*1/Sigma*D))%*%Inv_F_opt%*%t(Z_alles*sqrt(D*1/Sigma*D))

ret.obj=list()
ret.obj$ranef.logLik<-ranef.logLik
ret.obj$opt<-opt
ret.obj$Delta<-Deltafinal
ret.obj$Q<-Q_final
ret.obj$Standard_errors<-Standard_errors
ret.obj$phi<-phi
ret.obj$complexity<-complexity
return(ret.obj)
}
