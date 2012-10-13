glmm_final_multi_random<-function(y,X,W,k,q_start,Delta_start,s,n,steps=1000,family,method,overdispersion,phi,nue=1,rnd.len,print.iter.final=FALSE)
{
N<-length(y)
lin<-ncol(as.matrix(X))
Eta<-cbind(X,W)%*%Delta_start
Mu<-as.vector(family$linkinv(Eta))
Sigma<-as.vector(family$variance(Mu))
if(overdispersion)
Sigma<-Sigma*phi
D<-as.vector(family$mu.eta(Eta))
W0_inv<-D*1/Sigma*D

if(print.iter.final)
print(paste("Final Re-estimation Iteration ", 1,sep=""))


Z_alles<-cbind(X,W)

if(all(s==1))
{
P1<-c(rep(0,lin),rep(diag(q_start)^(-1),n))
P1<-diag(P1)
}else{
P1<-matrix(0,lin+n%*%s,lin+n%*%s)
inv.act<-chol2inv(chol(q_start[1:s[1],1:s[1]]))
for(jf in 1:n[1])
P1[(lin+(jf-1)*s[1]+1):(lin+jf*s[1]),(lin+(jf-1)*s[1]+1):(lin+jf*s[1])]<-inv.act

     for (zu in 2:rnd.len)
     {
     inv.act<-chol2inv(chol(q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
     for(jf in 1:n[zu])
     P1[(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
     (lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-inv.act
     }
}

if(overdispersion)
{
E<-diag(N)
M0<-Z_alles%*%chol2inv(chol((t(Z_alles)%*%(Z_alles*W0_inv))+P1))%*%t(Z_alles*W0_inv)
phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(M0)))
}

Delta<-matrix(0,steps,(lin+s%*%n))
Delta[1,]<-Delta_start

Q<-list()
Q[[1]]<-q_start

score_vec<-rep(0,(lin+s%*%n))
D<-as.vector(family$mu.eta(Eta))
Delta_r<-rep(0,lin+s%*%n)

l=1
opt<-steps

score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-P1%*%Delta[1,]
F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1

InvFisher<-chol2inv(chol(F_gross))
Delta_r<-InvFisher%*%score_vec
Delta[1,]<-Delta[1,]+nue*Delta_r

Eta<-Z_alles%*%Delta[1,]

Mu<-as.vector(family$linkinv(Eta))
Sigma<-as.vector(family$variance(Mu))
D<-as.vector(family$mu.eta(Eta))

if (method=="EM")
{
F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
InvFisher<-chol2inv(chol(F_gross))
############################# Q updaten ################
   Q1<-matrix(0,sum(s),sum(s))
   Q1[1:s[1],1:s[1]]<-InvFisher[(lin+1):(lin+s[1]),(lin+1):(lin+s[1])]+Delta[1,(lin+1):(lin+s[1])]%*%t(Delta[1,(lin+1):(lin+s[1])])
   for (i in 2:n[1])
   Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher[(lin+(i-1)*s[1]+1):(lin+i*s[1]),(lin+(i-1)*s[1]+1):(lin+i*s[1])]+Delta[1,(lin+(i-1)*s[1]+1):(lin+i*s[1])]%*%t(Delta[1,(lin+(i-1)*s[1]+1):(lin+i*s[1])])
   Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]

     for (zu in 2:rnd.len)
     {
     Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher[(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
     for (i in 2:n[zu])
     Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher[(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
     Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-1/n[zu]*Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
     }

}else{
Eta_tilde<-Eta+(y-Mu)*1/D

Betadach<-Delta[1,1:lin]

   if(all(s==1))
   {
   q_start_vec<-diag(q_start)
   upp<-rep(20,sum(s))
   low<-rep(1e-13,sum(s))
   optim.obj<-try(bobyqa(sqrt(q_start_vec),likelihood_diag,D=D,Sigma=Sigma,X=X,X_aktuell=X,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
   Q1<-diag(optim.obj$par)^2
   }else{
   q_start_vec<-c(diag(q_start)[1:s[1]],q_start[1:s[1],1:s[1]][lower.tri(q_start[1:s[1],1:s[1]])])
   up1<-min(20,50*max(q_start_vec))
   low<-c(rep(0,s[1]),rep(-up1,0.5*(s[1]^2-s[1])))

     for (zu in 2:rnd.len)
     {
     q_start_vec<-c(q_start_vec,c(diag(q_start)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])],q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]))
     up1<-min(20,50*max(q_start_vec))
     low<-c(low,c(rep(0,s[zu]),rep(-up1,0.5*(s[zu]^2-s[zu]))))
     }
     upp<-rep(up1,length(q_start_vec))
     optim.obj<-try(bobyqa(q_start_vec,likelihood_block,D=D,Sigma=Sigma,X=X,X_aktuell=X,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
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
Uu<-(E-(Z_alles*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher%*%t(Z_alles*sqrt(D*1/Sigma*D*1/Sigma)))%*%(E-M0)
FinalHat<-E-Uu
phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(FinalHat)))
Sigma<-Sigma*phi
}


###############################################################################################################################################
################################################################### Boost ###################################################################
for (l in 2:steps)
{

if(print.iter.final)
print(paste("Final Re-estimation Iteration ", l,sep=""))

if(all(s==1))
{
P1<-c(rep(0,lin),rep(diag(Q1)^(-1),n))
P1<-diag(P1)
}else{
P1<-matrix(0,lin+n%*%s,lin+n%*%s)
inv.act<-chol2inv(chol(Q1[1:s[1],1:s[1]]))
for(jf in 1:n[1])
P1[(lin+(jf-1)*s[1]+1):(lin+jf*s[1]),(lin+(jf-1)*s[1]+1):(lin+jf*s[1])]<-inv.act

     for (zu in 2:rnd.len)
     {
     inv.act<-chol2inv(chol(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
     for(jf in 1:n[zu])
     P1[(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
     (lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-inv.act
     }
}

  score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-P1%*%Delta[l-1,]
  F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1

InvFisher<-chol2inv(chol(F_gross))
Delta_r<-InvFisher%*%score_vec
Delta[l,]<-Delta[l-1,]+nue*Delta_r
Eta<-Z_alles%*%Delta[l,]
Mu<-as.vector(family$linkinv(Eta))
Sigma<-as.vector(family$variance(Mu))
D<-as.vector(family$mu.eta(Eta))

if (method=="EM")
{
F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
InvFisher<-chol2inv(chol(F_gross))
############################# Q update ################
   Q1<-matrix(0,sum(s),sum(s))
   Q1[1:s[1],1:s[1]]<-InvFisher[(lin+1):(lin+s[1]),(lin+1):(lin+s[1])]+Delta[1,(lin+1):(lin+s[1])]%*%t(Delta[l,(lin+1):(lin+s[1])])
   for (i in 2:n[1])
   Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher[(lin+(i-1)*s[1]+1):(lin+i*s[1]),(lin+(i-1)*s[1]+1):(lin+i*s[1])]+Delta[l,(lin+(i-1)*s[1]+1):(lin+i*s[1])]%*%t(Delta[1,(lin+(i-1)*s[1]+1):(lin+i*s[1])])
   Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]

     for (zu in 2:rnd.len)
     {
     Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher[(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
     for (i in 2:n[zu])
     Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher[(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
     Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-1/n[zu]*Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
     }

}else{
Eta_tilde<-Eta+(y-Mu)*1/D

Betadach<-Delta[l,1:lin]

   if(all(s==1))
   {
   Q1_vec<-diag(Q1)
   optim.obj<-try(bobyqa(sqrt(Q1_vec),likelihood_diag,D=D,Sigma=Sigma,X=X,X_aktuell=X,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
   Q1<-diag(optim.obj$par)^2
   }else{
   Q1_vec<-c(diag(Q1)[1:s[1]],Q1[1:s[1],1:s[1]][lower.tri(Q1[1:s[1],1:s[1]])])

     for (zu in 2:rnd.len)
     Q1_vec<-c(Q1_vec,c(diag(Q1)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])],Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]))

     optim.obj<-try(bobyqa(Q1_vec,likelihood_block,D=D,Sigma=Sigma,X=X,X_aktuell=X,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
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

#print(paste("Final Iteration =", l,sep=""))


opt<-l
Deltafinal<-Delta[l,]
Q_final<-Q[[l+1]]

  
if(all(s==1))
{
P1<-c(rep(0,lin),rep(diag(Q_final)^(-1),n))
P1<-diag(P1)
}else{
P1<-matrix(0,lin+n%*%s,lin+n%*%s)
inv.act<-chol2inv(chol(Q_final[1:s[1],1:s[1]]))
for(jf in 1:n[1])
P1[(lin+(jf-1)*s[1]+1):(lin+jf*s[1]),(lin+(jf-1)*s[1]+1):(lin+jf*s[1])]<-inv.act

     for (zu in 2:rnd.len)
     {
     inv.act<-chol2inv(chol(Q_final[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
     for(jf in 1:n[zu])
     P1[(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
     (lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-inv.act
     }
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
