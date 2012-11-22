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

est.glmmLasso<-function(fix,rnd,data,lambda,family=gaussian(link = "identity"),control=list())
{  
very.old.names<-attr(terms(fix),"term.labels")
                    
y <- model.response(model.frame(fix, data))

X <- model.matrix(fix, data)

if(dim(X)[2]==1)
{
if(colnames(X)=="(Intercept)")
stop("No terms to select! Use glmer, glmmPQL or glmmML!")
}

if(all(substr(very.old.names,1,9)!="as.factor"))
very.old.names<-colnames(X)[colnames(X)!="(Intercept)"]


fix.part<-very.old.names[1]
if(length(very.old.names)>1)
{
for(iu in 2:length(very.old.names))
fix.part<-paste(fix.part,very.old.names[iu],sep=" + ")
}
fix<-formula(paste(attr(terms(fix),"variables")[[2]],fix.part,sep=" ~ "))

old.names<-attr(X,"dimnames")[[2]]

rnd.len<-length(rnd)


if(rnd.len==1)
{
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

}else{     
rndformula <- list()
newrndfrml <- list()
n <- numeric()
s <- numeric()
k<-NULL
random.labels<-list()
W<-NULL

     for (zu in 1:rnd.len)
     {
     rndformula[[zu]] <- as.character(rnd[[zu]])

     trmsrnd <- terms(rnd[[zu]])
     newrndfrml[[zu]] <- "~ -1"
     newrndfrml[[zu]] <- paste(newrndfrml[[zu]],  if(attr(trmsrnd, "intercept")) names(rnd)[zu] else "", sep=" + ")

     if(length(attr(trmsrnd, "variables"))>1)
     {
     newrndfrml[[zu]] <- paste(newrndfrml[[zu]],  
     paste(sapply(attr(trmsrnd,"term.labels"), function(lbl){
                         paste(lbl, names(rnd)[zu], sep=":")
                      }), collapse=" + "), sep="+") }
     W_start <- model.matrix(formula(newrndfrml[[zu]]), data)
     
     
     rnlabels<-terms(formula(newrndfrml[[zu]]))
     random.labels[[zu]]<-attr(rnlabels,"term.labels")
     k1<-table(data[,colnames(data)==(names(rnd)[zu])])   
     n[zu]<-length(k1)
     s[zu]<-dim(W_start)[2]/n[zu]

     if(s[zu]>1)
     {
     W2<-W_start[,seq(from=1,to=1+(s[zu]-1)*n[zu],by=n[zu])]
     for (i in 2:n[zu])
     W2<-cbind(W2,W_start[,seq(from=i,to=i+(s[zu]-1)*n[zu],by=n[zu])])
     }else{
     W2<-W_start
     }
     W<-cbind(W,W2)
     k<-c(k,k1)
     }
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
control$start<-c(rep(0,(length(control$lin)+n%*%s)))

if(length(control$q_start)==0)
{
control$q_start<-rep(0.1,sum(s))
if(sum(s)>1)
control$q_start<-diag(control$q_start)
}

beta_null<-control$start[1:length(control$lin)]
ranef_null<-control$start[(length(control$lin)+1):(length(control$lin)+n%*%s)]
q_start<-control$q_start

N<-length(y)
if (all (X==0))
{
X<-rep(1,N)
X<-as.matrix(X)
}
# add Intercept 
if (!all (X[,1]==1))
{
cat("Warning:\n")
cat("Intercept has to be incorporated and is added!\n")
X<-cbind(1,X)
colnames(X)[1]<-"(Intercept)"
old.names<-c("(Intercept)",old.names)
}

lin<-dim(X)[2]

Z_fastalles<-X

no.sel<-is.element(attr(X,"dimnames")[[2]],control$lin)


U<-as.matrix(X[,!no.sel])
colnames(U)<-old.names[!no.sel]
X<-as.matrix(X[,no.sel])
colnames(X)<-old.names[no.sel]


q<-dim(X)[2]
phi<-control$phi_start

if(lin>1)
{
Eta_start<-X%*%beta_null+W%*%ranef_null
}else{
Eta_start<-rep(beta_null,N)+W%*%ranef_null
}

D<-as.vector(family$mu.eta(Eta_start))
Mu<-as.vector(family$linkinv(Eta_start))
Sigma<-as.vector(family$variance(Mu))


if(rnd.len==1)
{
lin0<-sum(beta_null!=0)
if(s==1)
{
Q_start<-diag(q_start,s)
p_start<-c(rep(0,lin0),rep((q_start)^(-1),n*s))
P1<-diag(p_start)
}else{
Q_start<-q_start
P1<-matrix(0,lin0+n*s,lin0+n*s)
inv.act<-chol2inv(chol(q_start))
for(jf in 1:n)
P1[(lin0+(jf-1)*s+1):(lin0+jf*s),(lin0+(jf-1)*s+1):(lin0+jf*s)]<-inv.act
}
}else{
lin0<-sum(beta_null!=0)
if(all(s==1))
{
Q_start<-diag(diag(q_start),sum(s))
p_start<-c(rep(0,lin0),rep(diag(q_start)^(-1),n))
P1<-diag(p_start)
}else{
Q_start<-matrix(0,sum(s),sum(s))
Q_start[1:s[1],1:s[1]]<-q_start[1:s[1],1:s[1]]
P1<-matrix(0,lin0+n%*%s,lin0+n%*%s)
inv.act<-chol2inv(chol(q_start[1:s[1],1:s[1]]))
for(jf in 1:n[1])
P1[(lin0+(jf-1)*s[1]+1):(lin0+jf*s[1]),(lin0+(jf-1)*s[1]+1):(lin0+jf*s[1])]<-inv.act

     for (zu in 2:rnd.len)
     {
     Q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
     inv.act<-chol2inv(chol(q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
     for(jf in 1:n[zu])
     P1[(lin0+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin0+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
     (lin0+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin0+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-inv.act
     }
}
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
Delta<-matrix(0,control$steps,(lin+n%*%s))
Delta[1,1:q]<-beta_null
Delta[1,(lin+1):(lin+n%*%s)]<-t(ranef_null)

Delta_start<-Delta[1,]

active_old<-!is.element(Delta[1,],0)

Q<-list()
Q[[1]]<-Q_start

l=1

if(rnd.len==1)
{
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
}else{
   if(all(s==1))
   {
   P1<-c(rep(0,lin),rep(diag(Q_start)^(-1),n))
   score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-t(t(P1)*Delta[1,])
   }else{
   Q_inv<-list()
   Q_inv[[1]]<-chol2inv(chol(Q_start[1:s[1],1:s[1]]))
   P1<-matrix(0,lin+n%*%s,lin+n%*%s)
   for(jf in 1:n[1])
   P1[(lin+(jf-1)*s[1]+1):(lin+jf*s[1]),(lin+(jf-1)*s[1]+1):(lin+jf*s[1])]<-Q_inv[[1]]

     for (zu in 2:rnd.len)
     {
     Q_inv[[zu]]<-chol2inv(chol(Q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
     for(jf in 1:n[zu])
     P1[(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
     (lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
     }

   score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-P1%*%Delta[1,]
   }
}




if (BLOCK)
{
grad.1<-gradient.lasso.block(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
}else{
grad.1<-gradient.lasso(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
}

score_vec<-c(score_vec[1:q],grad.1,score_vec[(lin+1):(lin+n%*%s)])

if(all(s==1))
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

active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,n%*%s))
Z_aktuell<-Z_alles[,active]
lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))


if(rnd.len==1)
{
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
}else{
   if(all(s==1))
   {
   P_akt<-c(rep(0,lin_akt),rep(diag(Q_start)^(-1),n))
   F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
   }else{
   P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
   for(jf in 1:n[1])
   P_akt[(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1]),(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1])]<-Q_inv[[1]]

     for (zu in 2:rnd.len)
     {
     for(jf in 1:n[zu])
     P_akt[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
     (lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
     }
     F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
   }
}



betaact<-Delta[1,active]
InvFisher<-chol2inv(chol(F_gross))

if (control$method=="EM")
{   
############################# Q update ################
   if(rnd.len==1)
   {
   Q1<-InvFisher[(lin_akt+1):(lin_akt+s),(lin_akt+1):(lin_akt+s)]+Delta[1,(lin+1):(lin+s)]%*%t(Delta[1,(lin+1):(lin+s)])
   for (i in 2:n)
   Q1<-Q1+InvFisher[(lin_akt+(i-1)*s+1):(lin_akt+i*s),(lin_akt+(i-1)*s+1):(lin_akt+i*s)]+Delta[1,(lin+(i-1)*s+1):(lin+i*s)]%*%t(Delta[1,(lin+(i-1)*s+1):(lin+i*s)])
   Q1<-1/n*Q1
   }else{
   Q1<-matrix(0,sum(s),sum(s))
   Q1[1:s[1],1:s[1]]<-InvFisher[(lin_akt+1):(lin_akt+s[1]),(lin_akt+1):(lin_akt+s[1])]+Delta[1,(lin+1):(lin+s[1])]%*%t(Delta[1,(lin+1):(lin+s[1])])
   for (i in 2:n[1])
   Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher[(lin_akt+(i-1)*s[1]+1):(lin_akt+i*s[1]),(lin_akt+(i-1)*s[1]+1):(lin_akt+i*s[1])]+Delta[1,(lin+(i-1)*s[1]+1):(lin+i*s[1])]%*%t(Delta[1,(lin+(i-1)*s[1]+1):(lin+i*s[1])])
   Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]

     for (zu in 2:rnd.len)
     {
     Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
     for (i in 2:n[zu])
     Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
     Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-1/n[zu]*Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
     }
   }
}else{
Eta_tilde<-Eta+(y-Mu)*1/D
Betadach<-Delta[1,1:(lin)]     
aktuell_vec<-!is.element(Delta[1,1:(lin)],0)
X_aktuell<-Z_fastalles[,aktuell_vec]


if(rnd.len==1)
{

   if(s==1)
   {
   upp<-min(20,50*Q_start)
   low<-1e-14
   optim.obj<-bobyqa(sqrt(Q_start),likelihood_bobyqa,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp)
   Q1<-as.matrix(optim.obj$par)^2
   }else{
   q_start_vec<-c(diag(q_start),q_start[lower.tri(q_start)])
   up1<-min(20,50*max(q_start_vec))#[(s+1):(s*(s+1)*0.5)]))
   upp<-rep(up1,length(q_start_vec))
   low<-c(rep(0,s),rep(-up1,0.5*(s^2-s)))
#   kkk_vec<-c(rep(-1,s),rep(0.5,0.5*(s^2-s)))
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
      }   
}else{
   if(all(s==1))
   {
   q_start_vec<-diag(q_start)
   upp<-rep(min(20,50*diag(q_start)),sum(s))
   low<-rep(1e-14,sum(s))
   optim.obj<-try(bobyqa(sqrt(q_start_vec),likelihood_diag,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
   Q1<-diag(optim.obj$par)^2
   }else{
   q_start_vec<-c(diag(q_start)[1:s[1]],q_start[1:s[1],1:s[1]][lower.tri(q_start[1:s[1],1:s[1]])])
   up1<-min(20,50*max(q_start_vec))
   low<-c(rep(0,s[1]),rep(-up1,0.5*(s[1]^2-s[1])))
#   kkk_vec<-c(rep(-1,s[1]),rep(0.5,0.5*(s[1]^2-s[1])))

     for (zu in 2:rnd.len)
     {
     q_start_vec<-c(q_start_vec,c(diag(q_start)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])],q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]))
     up1<-min(20,50*max(q_start_vec))
     low<-c(low,c(rep(0,s[zu]),rep(-up1,0.5*(s[zu]^2-s[zu]))))
#     kkk_vec<-c(kkk_vec,c(rep(-1,s[1]),rep(0.5,0.5*(s[1]^2-s[1]))))
     }
     upp<-rep(up1,length(q_start_vec))
     optim.obj<-try(bobyqa(q_start_vec,likelihood_block,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
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



if(rnd.len==1)
{
if(s==1)
{
P1<-c(rep(0,lin),rep((Q1^(-1)),n*s))
score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-t(t(P1)*Delta[1,])
}else{
Q_inv<-solve(Q1)
P1<-matrix(0,lin+n*s,lin+n*s)
for(j in 1:n)
P1[(lin+(j-1)*s+1):(lin+j*s),(lin+(j-1)*s+1):(lin+j*s)]<-Q_inv
score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-P1%*%Delta[l-1,]
}
}else{
   if(all(s==1))
   {
   P1<-c(rep(0,lin),rep(diag(Q1)^(-1),n))
   score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-t(t(P1)*Delta[l-1,])
   }else{
   Q_inv<-list()
   Q_inv[[1]]<-chol2inv(chol(Q1[1:s[1],1:s[1]]))
   P1<-matrix(0,lin+n%*%s,lin+n%*%s)
   for(jf in 1:n[1])
   P1[(lin+(jf-1)*s[1]+1):(lin+jf*s[1]),(lin+(jf-1)*s[1]+1):(lin+jf*s[1])]<-Q_inv[[1]]

     for (zu in 2:rnd.len)
     {
     Q_inv[[zu]]<-chol2inv(chol(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
     for(jf in 1:n[zu])
     P1[(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
     (lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
     }

   score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-P1%*%Delta[l-1,]
   }
}


if (BLOCK)
{
grad.1<-gradient.lasso.block(score.beta=score_vec[(q+1):lin],b=Delta[l-1,(q+1):lin],lambda.b=lambda,block=block)
}else{
grad.1<-gradient.lasso(score.beta=score_vec[(q+1):lin],b=Delta[l-1,(q+1):lin],lambda.b=lambda)
}
score_vec<-c(score_vec[1:q],grad.1,score_vec[(lin+1):(lin+n%*%s)])

if(all(s==1))
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
active<-c(rep(T,q),!is.element(Delta[l,(q+1):lin],0),rep(T,n%*%s))
Z_aktuell<-Z_alles[,active]
lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))



if(rnd.len==1)
{
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
}else{
   if(all(s==1))
   {
   P_akt<-c(rep(0,lin_akt),rep(diag(Q1)^(-1),n))
   F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
   }else{
   P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
   for(jf in 1:n[1])
   P_akt[(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1]),(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1])]<-Q_inv[[1]]

     for (zu in 2:rnd.len)
     {
     for(jf in 1:n[zu])
     P_akt[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
     (lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
     }
     F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
   }
}


betaact<-Delta[l,active]
InvFisher<-chol2inv(chol(F_gross))
      
if (control$method=="EM")
{        
############################# Q update ################
   if(rnd.len==1)
   {
   Q1<-InvFisher[(lin_akt+1):(lin_akt+s),(lin_akt+1):(lin_akt+s)]+Delta[l,(lin+1):(lin+s)]%*%t(Delta[l,(lin+1):(lin+s)])
   for (i in 2:n)
   Q1<-Q1+InvFisher[(lin_akt+(i-1)*s+1):(lin_akt+i*s),(lin_akt+(i-1)*s+1):(lin_akt+i*s)]+Delta[l,(lin+(i-1)*s+1):(lin+i*s)]%*%t(Delta[l,(lin+(i-1)*s+1):(lin+i*s)])
   Q1<-1/n*Q1
   }else{
   Q1<-matrix(0,sum(s),sum(s))
   Q1[1:s[1],1:s[1]]<-InvFisher[(lin_akt+1):(lin_akt+s[1]),(lin_akt+1):(lin_akt+s[1])]+Delta[l,(lin+1):(lin+s[1])]%*%t(Delta[1,(lin+1):(lin+s[1])])
   for (i in 2:n[1])
   Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher[(lin_akt+(i-1)*s[1]+1):(lin_akt+i*s[1]),(lin_akt+(i-1)*s[1]+1):(lin_akt+i*s[1])]+Delta[l,(lin+(i-1)*s[1]+1):(lin+i*s[1])]%*%t(Delta[l,(lin+(i-1)*s[1]+1):(lin+i*s[1])])
   Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]

     for (zu in 2:rnd.len)
     {
     Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
     for (i in 2:n[zu])
     Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
     Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-1/n[zu]*Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
     }
   }  
}else{
Eta_tilde<-Eta+(y-Mu)*1/D

Betadach<-Delta[l,1:(lin)]

aktuell_vec<-!is.element(Delta[l,1:(lin)],0)
X_aktuell<-Z_fastalles[,aktuell_vec]


if(rnd.len==1)
{

   if(s==1)
   {
   if(Q1<1e-14)
   low<-0

   optim.obj<-bobyqa(sqrt(Q1),likelihood_bobyqa,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp)
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
   }   
}else{
   if(all(s==1))
   {
   Q1_vec<-diag(Q1)
   optim.obj<-try(bobyqa(sqrt(Q1_vec),likelihood_diag,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
   Q1<-diag(optim.obj$par)^2
   }else{
     Q1_vec<-c(diag(Q1)[1:s[1]],Q1[1:s[1],1:s[1]][lower.tri(Q1[1:s[1],1:s[1]])])

     for (zu in 2:rnd.len)
     Q1_vec<-c(Q1_vec,c(diag(Q1)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])],Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]))
     
     optim.obj<-try(bobyqa(Q1_vec,likelihood_block,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
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

if(rnd.len==1)
{
glmm_fin<-try(glmm_final(y,Z_fastalles[,aaa],W,k,q_start=q_start,Delta_start=Delta_start[c(aaa,rep(T,n%*%s))],s,steps=control$maxIter,family=family,method=control$method.final,overdispersion=control$overdispersion,phi=control$phi,print.iter.final=control$print.iter.final))
}else{
glmm_fin<-try(glmm_final_multi_random(y,Z_fastalles[,aaa],W,k,q_start=q_start,Delta_start=Delta_start[c(aaa,rep(T,n%*%s))],s,n,steps=control$maxIter,family=family,method=control$method.final,overdispersion=control$overdispersion,phi=control$phi,rnd.len=rnd.len,print.iter.final=control$print.iter.final))
}


if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
{
cat("Warning:\n")
cat("Final Fisher scoring reestimation did not converge!\n")
}


Delta_neu2<-Delta_neu
Standard_errors<-rep(NA,length(Delta_neu))

if(class(glmm_fin)!="try-error")
{
Delta_neu2[c(aaa,rep(T,n%*%s))]<-glmm_fin$Delta
Standard_errors[c(aaa,rep(T,n%*%s))]<-glmm_fin$Standard_errors
Qfinal<-glmm_fin$Q
phi<-glmm_fin$phi
}

Delta_neu<-Delta_neu2

Eta_opt<-Z_alles%*%Delta_neu
Mu_opt<-as.vector(family$linkinv(Eta_opt))


if(rnd.len==1)
{
if(s==1)
Qfinal<-sqrt(Qfinal)


if(!is.matrix(Qfinal))
Qfinal<-as.matrix(Qfinal)
colnames(Qfinal)<-random.labels
rownames(Qfinal)<-random.labels
}else{
Qfinal_old<-Qfinal
Qfinal<-list()
Qfinal[[1]]<-as.matrix(Qfinal_old[1:s[1],1:s[1]])
colnames(Qfinal[[1]])<-random.labels[[1]]
rownames(Qfinal[[1]])<-random.labels[[1]]

     for (zu in 2:rnd.len)
     {
     Qfinal[[zu]]<-as.matrix(Qfinal_old[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])
     colnames(Qfinal[[zu]])<-random.labels[[zu]]
     rownames(Qfinal[[zu]])<-random.labels[[zu]]
     }
}

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

names(Delta_neu)[(lin+1):(lin+n%*%s)]<-colnames(W)
names(Standard_errors)[(lin+1):(lin+n%*%s)]<-colnames(W)
colnames(Delta)<-c(old.names,colnames(W))

aic<-NaN
bic<-NaN


if(rnd.len==1)
{
complexity<-0.5*(s*(s+1))
}else{
complexity<-0.5*(s[1]*(s[1]+1))
for(zu in 2:rnd.len)
complexity<-complexity+0.5*(s[zu]*(s[zu]+1))
}

if (family$family=="poisson")
{
aic<--2*sum(y*log(Mu_opt)-Mu_opt)+2*(sum(Delta_neu[1:(lin)]!=0)+complexity)
bic<--2*sum(y*log(Mu_opt)-Mu_opt)+log(N)*(sum(Delta_neu[1:(lin)]!=0)+complexity)
}

if (family$family=="binomial")
{
aic<--2*sum(y*log(Mu_opt)+(1-y)*log(1-Mu_opt))+2*(sum(Delta_neu[1:(lin)]!=0)+complexity)
bic<--2*sum(y*log(Mu_opt)+(1-y)*log(1-Mu_opt))+log(N)*(sum(Delta_neu[1:(lin)]!=0)+complexity)
}


if (family$family=="gaussian")
{
aic<--2*sum(y*Mu_opt-0.5*(Mu_opt^2))+2*(sum(Delta_neu[1:(lin)]!=0)+complexity)
bic<--2*sum(y*Mu_opt-0.5*(Mu_opt^2))+log(N)*(sum(Delta_neu[1:(lin)]!=0)+complexity)
}



ret.obj=list()
ret.obj$aic<-aic
ret.obj$bic<-bic
ret.obj$Deltamatrix<-Delta
ret.obj$ranef<-Delta_neu[(lin+1):(lin+n%*%s)]
ret.obj$coefficients<-Delta_neu[1:(lin)]
ret.obj$fixerror<-Standard_errors[1:(lin)]
ret.obj$ranerror<-Standard_errors[(lin+1):(lin+n%*%s)]
ret.obj$Q_long<-Q
ret.obj$Q<-Qfinal
ret.obj$y_hat<-Mu_opt
ret.obj$phi<-phi
ret.obj$family<-family
ret.obj$fix<-fix
ret.obj$newrndfrml<-newrndfrml
ret.obj$subject<-names(rnd)
ret.obj$data<-data
ret.obj$rnd.len<-rnd.len
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
printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
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
rnd.len<-object$rnd.len                 
family<-object$family

if(length(object$factor.names>0))
{
for (i in 1:length(object$factor.names))
newdata[,object$factor.names[i]]<-factor(newdata[,object$factor.names[i]],levels=object$factor.list[[i]])
}

X <- model.matrix(object$fix, newdata)

if(rnd.len==1)
{
subj.new<-levels(as.factor(newdata[,object$subject]))
subj.old<-levels(as.factor(object$data[,object$subject]))
subj.test<-is.element(subj.new,subj.old)
subj.ok<-subj.new[subj.test]

krit.random<-!all(!is.element(subj.new,subj.old))

if(krit.random)
{
W_start <- model.matrix(formula(object$newrndfrml), newdata)

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

subj.test<-subj.test.long

if(s>1)
{
W<-W_start[,seq(from=1,to=1+(s-1)*n,by=n)]
for (i in 2:n)
W<-cbind(W,W_start[,seq(from=i,to=i+(s-1)*n,by=n)])
}else{
W<-W_start
}

y<- as.vector(family$linkinv(X[,is.element(colnames(X),names(object$coef))]%*%object$coef[is.element(names(object$coef),colnames(X))]))
rand.ok<-is.element(newdata[,object$subject],subj.ok)
W.neu<-W[,subj.test]
if(dim(X)[1]!=1)
{
y[rand.ok]<- family$linkinv(cbind(X[,is.element(colnames(X),names(object$coef))],W.neu)[rand.ok,]%*%c(object$coef[is.element(names(object$coef),colnames(X))],object$ranef[match(colnames(W.neu),names(object$ranef))]))}else{
y[rand.ok]<- family$linkinv(c(X[,is.element(colnames(X),names(object$coef))],W.neu)[rand.ok]%*%c(object$coef[is.element(names(object$coef),colnames(X))],object$ranef[match(names(W.neu),names(object$ranef))]))}
}else{
W<-NULL
y<- as.vector(family$linkinv(X[,is.element(colnames(X),names(object$coef))]%*%object$coef[is.element(names(object$coef),colnames(X))]))
}

}else{

rnlabels<-list()
random.labels<-list()
s<-numeric()
k<-NULL   
n<-numeric()
W<- NULL
subj.test.long<-numeric()
subj.ok<-character()
krit.random<-logical()
subj.ok<-list()
W.single <- list()

for(zu in 1:rnd.len)
{
subj.new<-levels(as.factor(newdata[,object$subject[zu]]))
subj.old<-levels(as.factor(object$data[,object$subject[zu]]))
subj.test<-is.element(subj.new,subj.old)
subj.ok[[zu]]<-subj.new[subj.test]

krit.random[zu]<-!all(!is.element(subj.new,subj.old))

if(krit.random[zu])
{
rnlabels[[zu]]<-terms(formula(object$newrndfrml[[zu]]))
random.labels[[zu]]<-attr(rnlabels[[zu]],"term.labels")
s[zu]<-length(random.labels[[zu]])
k1<-table(newdata[,colnames(newdata)==(object$subject[zu])])
k<-c(k,k1)   
n[zu]<-length(k1)

W_start <- model.matrix(formula(object$newrndfrml[[zu]]), newdata)

     if(s[zu]>1)
     {
     W2<-W_start[,seq(from=1,to=1+(s[zu]-1)*n[zu],by=n[zu])]
     for (i in 2:n[zu])
     W2<-cbind(W2,W_start[,seq(from=i,to=i+(s[zu]-1)*n[zu],by=n[zu])])
     }else{
     W2<-W_start
     }
     W<-cbind(W,W2)
     W.single[[zu]]<-W2

if(s[zu]>1)
{for (i in 2:s[zu])
subj.test<-cbind(subj.test,subj.test)
subj.test<-as.vector(t(subj.test))
}
subj.test.long<-c(subj.test.long,subj.test)
}}


dim.W.single<-rep(0,rnd.len+1)
for(zu in 1:rnd.len)
dim.W.single[zu+1]<-dim(W.single[[zu]])[2]

if(!all(!krit.random))
{
rand.ok<-matrix(0,dim(newdata)[1],rnd.len)
for(zu in 1:rnd.len)
rand.ok[,zu]<-is.element(newdata[,object$subject[zu]],subj.ok[[zu]])

W.rnd<-matrix(0,dim(W)[1],dim(W)[2])
for(ur in 1:dim(newdata)[1])
{
  for (zu in 1:rnd.len)
  {
  if(rand.ok[ur,zu]==1)
  W.rnd[ur,sum(dim.W.single[1:zu])+1:sum(dim.W.single[zu+1])]<-W.single[[zu]][ur,]
  }
}

W.neu<-W.rnd[,as.logical(subj.test.long)]
if(!is.matrix(W.neu))
W.neu<-t(as.matrix(W.neu))
colnames(W.neu)<-colnames(W)[as.logical(subj.test.long)]
  
  if(dim(X)[1]!=1)
  {
  y<- family$linkinv(cbind(X[,is.element(colnames(X),names(object$coef))],W.neu)%*%c(object$coef[is.element(names(object$coef),colnames(X))],object$ranef[match(colnames(W.neu),names(object$ranef))]))
  }else{
  y<- family$linkinv(c(X[,is.element(colnames(X),names(object$coef))],W.neu)%*%c(object$coef[is.element(names(object$coef),colnames(X))],object$ranef[match(colnames(W.neu),names(object$ranef))]))
  }

}else{
y<- as.vector(family$linkinv(X[,is.element(colnames(X),names(object$coef))]%*%object$coef[is.element(names(object$coef),colnames(X))]))
}}
}
y
}



