gradient.lasso.block<-function(score.beta,b,lambda.b,block)       
{
p<-length(b)
grad.lasso<-rep(0,p)

lambda_vec<-rep(0,length(b))
lambda_vec[1:block[1]]<-lambda.b*sqrt(block[1])
for (i in 2:length(block))
lambda_vec[(block[i-1]+1):sum(block[1:i])]<-lambda.b*sqrt(block[i])

group.sum<-rep(0,length(block))
group.sum[1]<-sqrt(sum(b[1:block[1]]^2))
for (i in 2:length(block))
group.sum[i]<-sqrt(sum(b[(sum(block[1:(i-1)])+1):sum(block[1:i])]^2))

if(group.sum[1]!=0)
{
grad.lasso[1:block[1]]<-score.beta[1:block[1]]-lambda_vec[1]*(b[1:block[1]]/group.sum[1])
}else{

  if(group.sum[1]==0 & sqrt(sum(score.beta[1:block[1]]^2))>lambda_vec[1] )
  {
  grad.lasso[1:block[1]]<-score.beta[1:block[1]]-lambda_vec[1]*(score.beta[1:block[1]]/sqrt(sum(score.beta[1:block[1]]^2)))
  }else{
  grad.lasso[1:block[1]]<-0
}}

for (i in 2:length(block))
{
if(group.sum[i]!=0)
{
grad.lasso[(sum(block[1:(i-1)])+1):sum(block[1:i])]<-score.beta[(sum(block[1:(i-1)])+1):sum(block[1:i])]-lambda_vec[i]*(b[(sum(block[1:(i-1)])+1):sum(block[1:i])]/group.sum[i])
}else{

  if(group.sum[i]==0 & sqrt(sum(score.beta[(sum(block[1:(i-1)])+1):sum(block[1:i])]^2))>lambda_vec[i] )
  {
  grad.lasso[(sum(block[1:(i-1)])+1):sum(block[1:i])]<-score.beta[(sum(block[1:(i-1)])+1):sum(block[1:i])]-lambda_vec[i]*(score.beta[(sum(block[1:(i-1)])+1):sum(block[1:i])]/sqrt(sum(score.beta[(sum(block[1:(i-1)])+1):sum(block[1:i])]^2)))
  }else{
  grad.lasso[(sum(block[1:(i-1)])+1):sum(block[1:i])]<-0
}}
}
return(grad.lasso)
}

gradient.lasso<-function(score.beta,b,lambda.b)
{
p<-length(b)
grad.lasso<-rep(0,p)

b.isnt.0<-b!=0
grad.lasso[b.isnt.0]<-score.beta[b.isnt.0]-lambda.b*sign(b[b.isnt.0])
b.is.0<-(b==0 & abs(score.beta)>=lambda.b)
grad.lasso[b.is.0]<-score.beta[b.is.0]-lambda.b*sign(score.beta[b.is.0])
b.is.0<-(b==0 & abs(score.beta)<lambda.b)
grad.lasso[b.is.0]<-0
return(grad.lasso)
}


t.change<-function(grad,b)
{
a<-(sign(b)==-sign(grad)&sign(grad)!=0)
rate<--b[a]/grad[a]
if(all(a==FALSE))rate<-Inf
return(min(rate))
}

l2norm<-function(vec)
{
l2length<-sqrt(sum(vec^2))
if(l2length>0){vec.norm<-vec/l2length}
if(l2length==0){vec.norm<-vec}
return(list("length"=l2length,"normed"=vec.norm))
}                


#############################################################################################################################################               
###################################################### main Lasso function #################################################################
#  fix=formula(points ~ transfer.spendings + I(transfer.spendings^2)+ ave.unfair.score + transfer.receits + ball.possession+ ave.unfair.score + transfer.receits + ball.possession+ tackles + ave.attend + sold.out);rnd=list(team=~1 + ave.attend);data=soccer;lambda=100;family=gaussian(link = "identity");control = list(steps=500, lin="ave.attend", method="REML", sel.method="bic")
# fix=formula(y~x1+x2+x3+x4+x5);rnd=list(Person=~1);data=Hirst;lambda=4;family=poisson(link=log);control=list()
est.glmmLasso<-function(fix,rnd,data,lambda,family=gaussian(link = "identity"),control=list())
{                      
y <- model.response(model.frame(fix, data))

X <- model.matrix(fix, data)

very.old.names<-attr(terms(fix),"term.labels")

old.names<-attr(X,"dimnames")[[2]]

rndformula <- as.character(rnd)

trmsrnd <- terms(rnd[[1]])
newrndfrml <- "~ -1"
newrndfrml <- paste(newrndfrml,  if(attr(trmsrnd, "intercept")) names(rnd)[1] else "", sep=" + ")

if(length(attr(trmsrnd, "variables"))>1)
{
newrndfrml <- paste(newrndfrml,  
         paste(sapply(attr(trmsrnd,"term.labels"), function(lbl){
                     paste(lbl, names(rnd)[1], sep=":")
                 }), collapse=" + "), sep="+") }


W_start <- model.matrix(formula(newrndfrml), data)

rnlabels<-terms(formula(newrndfrml))
random.labels<-attr(rnlabels,"term.labels")
k<-table(data[,colnames(data)==(names(rnd)[1])])   
n<-length(k)
s<-dim(W_start)[2]/n

if(s>1)
{
W<-W_start[,seq(from=1,to=1+(s-1)*n,by=n)]
for (i in 2:n)
W<-cbind(W,W_start[,seq(from=i,to=i+(s-1)*n,by=n)])
}else{
W<-W_start
}

control<-do.call(glmmLassoControl, control)

if(control$print.iter)
print(paste("Iteration ", 1,sep=""))


if(sum(substr(control$lin,1,9)=="as.factor")>0)
{
group0<-substr(control$lin,1,9)=="as.factor"
spl0<-strsplit(control$lin[group0],"\\)")

control$lin<-control$lin[!group0]
for(ur in 1:length(spl0))
control$lin<-c(control$lin,old.names[is.element(substr(old.names,1,nchar(spl0[[ur]])),spl0[[ur]])])
}

if(length(control$lin)>1)
{
lin.out<-!is.element(very.old.names,control$lin[2:length(control$lin)])
very.old.names<-very.old.names[lin.out]
}

group<-substr(very.old.names,1,9)=="as.factor"

if(sum(group)>0)
{
spl<-strsplit(very.old.names[group],"\\(")

categ.names<-character()
for(uz in 1:length(spl))
categ.names<-c(categ.names,spl[[uz]][2])

spl2<-strsplit(categ.names,"\\)")
categ.names2<-character()
for(uz in 1:length(spl2))
categ.names2<-c(categ.names2,spl2[[uz]])

block<-numeric()
posi<-1
for(ip in 1:length(group))
{
  if(!group[ip])
  {
  block<-c(block,1)
  }else{
  block<-c(block,length(levels(as.factor(data[,categ.names2[posi]])))-1)
  posi<-posi+1
  }
}}else{
block<-rep(1,length(group))  
}

BLOCK<-FALSE
if(!all(block==1))
BLOCK<-TRUE

if(length(control$start)==0)
control$start<-c(rep(0,(length(control$lin)+n*s)))

if(length(control$q_start)==0)
{
control$q_start<-rep(0.1,s)
if(s>1)
control$q_start<-diag(control$q_start)
}

beta_null<-control$start[1:length(control$lin)]
ranef_null<-control$start[(length(control$lin)+1):(length(control$lin)+n*s)]
q_start<-control$q_start

N<-length(y)
if (all (X==0))
{
X<-rep(1,N)
X<-as.matrix(X)
}
# add Intercept 
if (!all (X[,1]==1))
X<-cbind(1,X)

lin<-dim(X)[2]

Z_fastalles<-X

no.sel<-is.element(attr(X,"dimnames")[[2]],control$lin)

U<-X[,!no.sel]
X<-as.matrix(X[,no.sel])

q<-dim(X)[2]
phi<-1

if(lin>1)
{
Eta_start<-X%*%beta_null+W%*%ranef_null
}else{
Eta_start<-rep(beta_null,N)+W%*%ranef_null
}

D<-as.vector(family$mu.eta(Eta_start))
Mu<-as.vector(family$linkinv(Eta_start))
Sigma<-as.vector(family$variance(Mu))


lin0<-sum(beta_null!=0)
if(s==1)
{
Q_start<-diag(q_start^2,s)
p_start<-c(rep(0,lin0),rep((q_start^2)^(-1),n*s))
P1<-diag(p_start)
}else{
Q_start<-q_start
P1<-matrix(0,lin0+n*s,lin0+n*s)
for(jf in 1:n)
P1[(lin0+(jf-1)*s+1):(lin0+jf*s),(lin0+(jf-1)*s+1):(lin0+jf*s)]<-chol2inv(chol(q_start))
}

if(control$overdispersion)
{  
W0_inv<-D*1/Sigma*D
Z_start<-cbind(X[,beta_null!=0],W)
M0<-Z_start%*%chol2inv(chol((t(Z_start)%*%(Z_start*W0_inv))+P1))%*%t(Z_start*W0_inv)
phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(M0)))
Sigma<-Sigma*phi
E<-diag(N)
}

Z_alles<-cbind(X,U,W)
########################################################## some definitions ################################################
Delta<-matrix(0,control$steps,(lin+s*n))
Delta[1,1:q]<-beta_null
Delta[1,(lin+1):(lin+s*n)]<-t(ranef_null)

Delta_start<-Delta[1,]

active_old<-!is.element(Delta[1,],0)

Q<-list()
Q[[1]]<-Q_start

l=1

if(s==1)
{
P1<-c(rep(0,lin),rep((Q_start^(-1)),n*s))
score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-t(t(P1)*Delta[1,])
}else{
Q_inv<-chol2inv(chol(Q_start))
P1<-matrix(0,lin+n*s,lin+n*s)
for(j in 1:n)
P1[(lin+(j-1)*s+1):(lin+j*s),(lin+(j-1)*s+1):(lin+j*s)]<-Q_inv
score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-P1%*%Delta[1,]
}

if (BLOCK)
{
grad.1<-gradient.lasso.block(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
}else{
grad.1<-gradient.lasso(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
}

score_vec<-c(score_vec[1:q],grad.1,score_vec[(lin+1):(lin+s*n)])

if(s==1)
{
F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+diag(P1)
}else{
F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
}

t_edge<-t.change(grad=score_vec,b=Delta[1,])

grad.2<-t(score_vec)%*%F_gross%*%score_vec

t_opt<-l2norm(score_vec)$length/grad.2
nue<-control$nue
Delta[1,]<-Delta[1,]+min(t_opt,t_edge)*nue*score_vec
Eta<-Z_alles%*%Delta[1,]
Mu<-as.vector(family$linkinv(Eta))
Sigma<-as.vector(family$variance(Mu))
D<-as.vector(family$mu.eta(Eta))

active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,n*s))
Z_aktuell<-Z_alles[,active]
lin_akt<-sum(!is.element(Delta[1,1:lin],0))

if(s==1)
{
P_akt<-c(rep(0,lin_akt),rep((Q_start^(-1)),n*s))
F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
}else{
P_akt<-matrix(0,lin_akt+n*s,lin_akt+n*s)
for(jf in 1:n)
P_akt[(lin_akt+(jf-1)*s+1):(lin_akt+jf*s),(lin_akt+(jf-1)*s+1):(lin_akt+jf*s)]<-Q_inv
F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
}

betaact<-Delta[1,active]
InvFisher<-chol2inv(chol(F_gross))

if (control$method=="EM")
{   
############################# Q update ################
Q1<-InvFisher[(lin_akt+1):(lin_akt+s),(lin_akt+1):(lin_akt+s)]+Delta[1,(lin+1):(lin+s)]%*%t(Delta[1,(lin+1):(lin+s)])
for (i in 2:n)
Q1<-Q1+InvFisher[(lin_akt+(i-1)*s+1):(lin_akt+i*s),(lin_akt+(i-1)*s+1):(lin_akt+i*s)]+Delta[1,(lin+(i-1)*s+1):(lin+i*s)]%*%t(Delta[1,(lin+(i-1)*s+1):(lin+i*s)])
Q1<-1/n*Q1
}else{
Eta_tilde<-Eta+(y-Mu)*1/D
Betadach<-Delta[1,1:(lin)]     
aktuell_vec<-!is.element(Delta[1,1:(lin)],0)
X_aktuell<-Z_fastalles[,aktuell_vec]

if(s==1)
{
optim.obj<-bobyqa(sqrt(Q_start),likelihood_bobyqa,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = 1e-14, upper = 20)
Q1<-as.matrix(optim.obj$par)^2
}else{
q_start_vec<-c(diag(q_start),q_start[lower.tri(q_start)])
up1<-min(20,50*max(q_start_vec))#[(s+1):(s*(s+1)*0.5)]))
upp<-rep(up1,length(q_start_vec))
low<-c(rep(0,s),rep(-up1,0.5*(s^2-s)))
kkk_vec<-c(rep(-1,s),rep(0.5,0.5*(s^2-s)))
optim.obj<-try(bobyqa(q_start_vec,likelihood,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp))
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

if(control$overdispersion)
{
Uu<-(E-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma)))%*%(E-M0)
FinalHat<-E-Uu
phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(FinalHat)))
Sigma<-Sigma*phi
}


y_dach<-as.vector(family$linkinv(Eta))
penal_neu<-lambda*sum(abs(Delta[l,(q+1):(lin)]))
Dev_neu<-sum(family$dev.resids(y,y_dach,wt=rep(1,N))^2)


vorz<-F
NRstep<-F
###############################################################################################################################################
################################################################### Boost ###################################################################
if(control$steps!=1)
{
for (l in 2:control$steps)
{
if(control$print.iter)
print(paste("Iteration ", l,sep=""))


if(s==1)
{
P1<-c(rep(0,lin),rep((Q1^(-1)),n*s))

score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-t(t(P1)*Delta[l-1,])
}else{
Q_inv<-solve(Q1)

P1<-matrix(0,lin+n*s,lin+n*s)

for(j in 1:n)
P1[(lin+(j-1)*s+1):(lin+j*s),(lin+(j-1)*s+1):(lin+j*s)]<-Q_inv
score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-P1%*%Delta[1,]
}

if (BLOCK)
{
grad.1<-gradient.lasso.block(score.beta=score_vec[(q+1):lin],b=Delta[l-1,(q+1):lin],lambda.b=lambda,block=block)
}else{
grad.1<-gradient.lasso(score.beta=score_vec[(q+1):lin],b=Delta[l-1,(q+1):lin],lambda.b=lambda)
}
score_vec<-c(score_vec[1:q],grad.1,score_vec[(lin+1):(lin+s*n)])

if(s==1)
{
F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+diag(P1)
}else{
F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
}

t_edge<-t.change(grad=score_vec,b=Delta[l-1,])

grad.2<-t(score_vec)%*%F_gross%*%score_vec

t_opt<-l2norm(score_vec)$length/grad.2
tryNR<- (t_opt<t_edge) && !(all(active_old==active)  && !NRstep)  && (lin_akt>q+1) 

if(tryNR) 
{
NRDelta<-Delta[l-1,active]+nue*InvFisher%*%score_vec[active]
vorz<-all(sign(NRDelta[(q+1):lin_akt])==sign(betaact[(q+1):lin_akt]))  
}   

if(!vorz)
tryNR<-F

if(tryNR)
{
Delta[l,active]<- NRDelta
NRstep<-T
}else{
Delta[l,]<-Delta[l-1,]+min(t_opt,t_edge)*nue*score_vec
}

Eta<-Z_alles%*%Delta[l,]
Mu<-as.vector(family$linkinv(Eta))
Sigma<-as.vector(family$variance(Mu))
D<-as.vector(family$mu.eta(Eta))

active_old<-active
active<-c(rep(T,q),!is.element(Delta[l,(q+1):lin],0),rep(T,n*s))
Z_aktuell<-Z_alles[,active]
lin_akt<-sum(!is.element(Delta[l,1:lin],0))

if(s==1)
{
P_akt<-c(rep(0,lin_akt),rep((Q1^(-1)),n*s))
F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
}else{
P_akt<-matrix(0,lin_akt+n*s,lin_akt+n*s)
for(jf in 1:n)
P_akt[(lin_akt+(jf-1)*s+1):(lin_akt+jf*s),(lin_akt+(jf-1)*s+1):(lin_akt+jf*s)]<-Q_inv
F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
}

betaact<-Delta[l,active]
InvFisher<-chol2inv(chol(F_gross))
      
if (control$method=="EM")
{        
############################# Q updaten ################
Q1<-InvFisher[(lin_akt+1):(lin_akt+s),(lin_akt+1):(lin_akt+s)]+Delta[l,(lin+1):(lin+s)]%*%t(Delta[l,(lin+1):(lin+s)])

for (i in 2:n)
Q1<-Q1+InvFisher[(lin_akt+(i-1)*s+1):(lin_akt+i*s),(lin_akt+(i-1)*s+1):(lin_akt+i*s)]+Delta[l,(lin+(i-1)*s+1):(lin+i*s)]%*%t(Delta[l,(lin+(i-1)*s+1):(lin+i*s)])
Q1<-1/n*Q1
}else{
Eta_tilde<-Eta+(y-Mu)*1/D

Betadach<-Delta[l,1:(lin)]

aktuell_vec<-!is.element(Delta[l,1:(lin)],0)
X_aktuell<-Z_fastalles[,aktuell_vec]

if(s==1)
{        
optim.obj<-try(bobyqa(sqrt(Q1),likelihood_bobyqa,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = 1e-12, upper = 20))
Q1<-as.matrix(optim.obj$par)^2
}else{
Q1_vec<-c(diag(Q1),Q1[lower.tri(Q1)])
optim.obj<-try(bobyqa(Q1_vec,likelihood,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp))

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

if(control$overdispersion)
{
Uu<-(E-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma)))%*%Uu
FinalHat<-E-Uu
phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(FinalHat)))
Sigma<-Sigma*phi
}


Q[[l+1]]<-Q1

y_dach<-as.vector(family$linkinv(Eta))

penal_old<-penal_neu

penal_neu<-lambda*sum(abs(Delta[l,(q+1):(lin)]))
Dev_old<-Dev_neu
Dev_neu<-sum(family$dev.resids(y,y_dach,wt=rep(1,N))^2)
###
finish1<-(2*abs(Dev_neu-Dev_old)/(2*abs(Dev_neu-penal_neu)+0.1)<control$epsilon)
finish2<-(2*abs(penal_neu-penal_old)/(2*abs(Dev_neu-penal_neu)+0.1)<control$epsilon)

if((finish1 && finish2) || (all(grad.1 == 0) ))
break
}}


Delta_neu<-Delta[l,]
Eta_opt<-Z_alles%*%Delta_neu
Mu_opt<-as.vector(family$linkinv(Eta_opt))
Sigma_opt<-as.vector(family$variance(Mu_opt))    
D_opt<-as.vector(family$mu.eta(Eta_opt))
Qfinal<-Q[[l+1]]

aaa<-!is.element(Delta_neu[1:(lin)],0)

glmm_fin<-try(glmm_final(y,Z_fastalles[,aaa],W,k,q_start=Qfinal,Delta_start=Delta_neu[c(aaa,rep(T,n*s))],s,steps=control$maxIter,family=family,method=control$method,overdispersion=control$overdispersion,phi=phi))

if(class(glmm_fin)=="try-error" || glmm_fin$opt>(control$maxIter-5))
{
glmm_fin<-try(glmm_final(y,Z_fastalles[,aaa],W,k,q_start=q_start,Delta_start=Delta_start[c(aaa,rep(T,n*s))],s,steps=2000,family=family,method=control$method,overdispersion=control$overdispersion,phi=phi))

if(class(glmm_fin)=="try-error" || glmm_fin$opt>1990)
{
cat("Warning:\n")
cat("Final Fisher scoring reestimation did not converge!")
}}

if(class(glmm_final)=="try-error")
break

Delta_neu2<-Delta_neu
Delta_neu2[c(aaa,rep(T,n*s))]<-glmm_fin$Delta

Delta_neu<-Delta_neu2

Standard_errors<-rep(0,length(Delta_neu))
Standard_errors[c(aaa,rep(T,n*s))]<-glmm_fin$Standard_errors


Qfinal<-glmm_fin$Q

phi<-glmm_fin$phi

Eta_opt<-Z_alles%*%Delta_neu
Mu_opt<-as.vector(family$linkinv(Eta_opt))

if(s==1)
Qfinal<-sqrt(Qfinal)

if(!is.matrix(Qfinal))
Qfinal<-as.matrix(Qfinal)
colnames(Qfinal)<-random.labels
rownames(Qfinal)<-random.labels

if(dim(X)[2]==1)
{
names(Delta_neu)[1]<-"(Intercept)"
names(Standard_errors)[1]<-"(Intercept)"
}else{
names(Delta_neu)[1:dim(X)[2]]<-colnames(X)
names(Standard_errors)[1:dim(X)[2]]<-colnames(X)
}

if(lin>1)
{
names(Delta_neu)[(dim(X)[2]+1):lin]<-colnames(U)
names(Standard_errors)[(dim(X)[2]+1):lin]<-colnames(U)
}

permut<-numeric()
for(ys in 1:length(old.names))
permut<-c(permut,match(old.names[ys],names(Delta_neu)))

Delta_neu[1:length(old.names)]<-Delta_neu[permut]
Standard_errors[1:length(old.names)]<-Standard_errors[permut]

names(Delta_neu)[1:length(old.names)]<-old.names
names(Standard_errors)[1:length(old.names)]<-old.names


Delta[,1:length(old.names)]<-Delta[,permut]

names(Delta_neu)[(lin+1):(lin+n*s)]<-colnames(W)
names(Standard_errors)[(lin+1):(lin+n*s)]<-colnames(W)
colnames(Delta)<-c(old.names,colnames(W))

ret.obj=list()
ret.obj$Deltamatrix<-Delta
ret.obj$ranef<-Delta_neu[(lin+1):(lin+n*s)]
ret.obj$coefficients<-Delta_neu[1:(lin)]
ret.obj$fixerror<-Standard_errors[1:(lin)]
ret.obj$ranerror<-Standard_errors[(lin+1):(lin+n*s)]
ret.obj$Q_long<-Q
ret.obj$Q<-Qfinal
ret.obj$y_hat<-Mu_opt
ret.obj$phi<-phi
ret.obj$family<-family
ret.obj$fix<-fix
ret.obj$newrndfrml<-newrndfrml
ret.obj$subject<-names(rnd)[1]
ret.obj$data<-data
return(ret.obj)
}



glmmLasso <- function(fix=formula, rnd=formula, data, lambda, family=NULL, control=list()) UseMethod("glmmLasso")

glmmLasso.formula <- function(fix,rnd,...)
{
est <- est.glmmLasso(fix,rnd,...)
est$fitted.values <- est$y_hat
est$StdDev <- est$Q
est$call <- match.call()
class(est) <- "glmmLasso"
est
}


print.glmmLasso <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\nFixed Effects:\n")
cat("\nCoefficients:\n")
print(x$coefficients)
cat("\nRandom Effects:\n")
cat("\nStdDev:\n")
print(x$StdDev)
}




summary.glmmLasso <- function(object, ...)
{
se <- object$fixerror
zval <- coefficients(object) / se
TAB <- cbind(Estimate = coefficients(object),
StdErr = se,
z.value = zval,
p.value = 2*pnorm(-abs(zval)))
res <- list(call=object$call,
coefficients=TAB,StdDev=object$StdDev)
class(res) <- "summary.glmmLasso"
res
}


print.summary.glmmLasso <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\n")
cat("\nFixed Effects:\n")
cat("\nCoefficients:\n")
printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
cat("\nRandom Effects:\n")
cat("\nStdDev:\n")
print(x$StdDev)
}


predict.glmmLasso <- function(object,newdata=NULL,...)
{
if(is.null(newdata))
{
y<-fitted(object)
}else{                   
family<-object$family

if(length(object$factor.names>0))
{
for (i in 1:length(object$factor.names))
{
newdata[,object$factor.names[i]]<-factor(newdata[,object$factor.names[i]],levels=object$factor.list[[i]])
}}

X <- model.matrix(object$fix, newdata)

subj.new<-levels(as.factor(newdata[,object$subject]))
subj.old<-levels(as.factor(object$data[,object$subject]))
subj.test<-is.element(subj.new,subj.old)
subj.ok<-subj.new[subj.test]

krit.random<-!all(!is.element(subj.new,subj.old))
if(krit.random)
{
W_start <- model.matrix(formula(object$newrndfrml), newdata)
}else{
W_start <- NULL
}

rnlabels<-terms(formula(object$newrndfrml))
random.labels<-attr(rnlabels,"term.labels")
s<-length(random.labels)
k<-table(newdata[,colnames(newdata)==(object$subject)])   
n<-length(k)

if(s>1)
{for (i in 2:s)
subj.test<-cbind(subj.test,subj.test)
subj.test<-as.vector(t(subj.test))
}


if(krit.random)
{
if(s>1)
{
W<-W_start[,seq(from=1,to=1+(s-1)*n,by=n)]
for (i in 2:n)
W<-cbind(W,W_start[,seq(from=i,to=i+(s-1)*n,by=n)])
}else{
W<-W_start
}
y<- as.vector(family$linkinv(X%*%object$coef))
rand.ok<-is.element(newdata[,object$subject],subj.ok)
W.neu<-W[,subj.test]
y[rand.ok]<- family$linkinv(cbind(X,W.neu)[rand.ok,]%*%c(object$coef,object$ranef[match(colnames(W.neu),names(object$ranef))]))
}else{
W<-NULL
y<- as.vector(family$linkinv(X%*%object$coef))
}
}
y
}



