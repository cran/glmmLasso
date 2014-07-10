gradient.lasso.block<-function(score.beta,b,lambda.b,block)       
{
  p<-length(b)
  grad.lasso<-rep(0,p)

  block2 <- rep(block,block)
  lambda_vec<-rep(lambda.b,length(b))
  lambda_vec <- lambda_vec*sqrt(block2)
    
  group.sum<-rep(0,length(block))
  group.sum[1]<-sqrt(sum(b[1:block[1]]^2))
  if(length(block)>1)
  {  
  for (i in 2:length(block))
    group.sum[i]<-sqrt(sum(b[(sum(block[1:(i-1)])+1):sum(block[1:i])]^2))
  }
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
  
  if(length(block)>1)
  {  
  for (i in 2:length(block))
  {
    if(group.sum[i]!=0)
    {
      grad.lasso[(sum(block[1:(i-1)])+1):sum(block[1:i])]<-score.beta[(sum(block[1:(i-1)])+1):sum(block[1:i])]-lambda_vec[sum(block[1:(i-1)])+1]*(b[(sum(block[1:(i-1)])+1):sum(block[1:i])]/group.sum[i])
    }else{
      
      if(group.sum[i]==0 & sqrt(sum(score.beta[(sum(block[1:(i-1)])+1):sum(block[1:i])]^2))>lambda_vec[sum(block[1:(i-1)])+1])
      {
        grad.lasso[(sum(block[1:(i-1)])+1):sum(block[1:i])]<-score.beta[(sum(block[1:(i-1)])+1):sum(block[1:i])]-lambda_vec[sum(block[1:(i-1)])+1]*(score.beta[(sum(block[1:(i-1)])+1):sum(block[1:i])]/sqrt(sum(score.beta[(sum(block[1:(i-1)])+1):sum(block[1:i])]^2)))
      }else{
        grad.lasso[(sum(block[1:(i-1)])+1):sum(block[1:i])]<-0
      }}
  }}
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
  rate <- rep(Inf, length(a))
  rate[a] <--b[a]/grad[a]
  ret.obj<-list()
  ret.obj$min.rate<-min(rate)
  ret.obj$whichmin <- which.min(rate)
  return(ret.obj)
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

est.glmmLasso<-function(fix,rnd,data,lambda,family=gaussian(link = "identity"),
                        final.re=FALSE,switch.NR=T,control=list())
{  
  
  
 if(grepl("\\*", fix[3]))
 stop("Usage of '*' not allowed in formula! Please specify the corresponding variables separately.")  
   
  fix.old<-fix
  ic.dummy<-attr(terms(fix),"intercept")  
  y <- model.response(model.frame(fix, data))
  very.old.names<-attr(terms(fix),"term.labels")
    
      
  
  if(ic.dummy!=1 && sum(substr(very.old.names,1,9)=="as.factor")>0){
    fix.help <- update(fix,~ .+1)
    orig.names <- colnames(model.matrix(fix.help, data))[-1]
  }else{
    orig.names <- colnames(model.matrix(fix, data))
  }

  control<-do.call(glmmLassoControl, control)
      
  if(!is.null(control$index))
  {
    order.vec<-order(control$index)
    very.old.names<-very.old.names[order.vec]
    control$index<-control$index[order.vec]
  }else{
    control$index<-1:length(very.old.names)
  }
  
  if(length(control$index)!=length(very.old.names))
    stop("Length of vector defining the grouping of the variables doesn't match with 
         the formula!")

  attr(control$index,"names")<-very.old.names
 
  fix<-formula(paste("y~-1+",paste(very.old.names,collapse="+"))) 
  
  if(ic.dummy==1)
  {
    fix<-update(fix,~ .+1) 
    control$index<-c(NA,control$index)
    names(control$index)[1]<-"(Intercept)"
  }
  
 
  index.new<-c()
  fac.variab<-logical()
  for(i in 1:length(control$index))
  {
    if(!grepl("as.factor",names(control$index)[i]))
    {
      index.new<-c(index.new,control$index[i]) 
      fac.variab<-c(fac.variab,F)
    }else{
      if(!grepl("\\:", names(control$index)[i]))
      {  
      fac.name<-strsplit(strsplit(names(control$index)[i],"\\(")[[1]][2],"\\)")[[1]][1]
      }else{
      fac.name<-unlist(strsplit(unlist(strsplit(unlist(strsplit(names(control$index)[i],"\\(")),"\\)")),"\\:"))
      fac.name<-paste(fac.name[2],":",fac.name[length(fac.name)],sep="")
      }
      if(!grepl("\\:", fac.name))
      {
        length.fac<-length(levels(as.factor(data[,fac.name])))-1
      }else{
        length.fac<-(length(levels(data[,strsplit(fac.name,":")[[1]][1]]))-1)*(length(levels(data[,strsplit(fac.name,":")[[1]][2]]))-1)
      }
      index.new<-c(index.new,rep(control$index[i],length.fac))
      fac.variab<-c(fac.variab,rep(T,length.fac))
    }
  }
      
if(ic.dummy!=1 && sum(substr(very.old.names,1,9)=="as.factor")>0){
  fix.help <- update(fix,~ .+1)
  X <- model.matrix(fix.help, data)[,-1]
}else{
  X <- model.matrix(fix, data)  
}

if(any(fac.variab))
  X[,fac.variab]<-scale(X[,fac.variab])

if(dim(X)[2]==1)
{
  if(colnames(X)=="(Intercept)")
    stop("No terms to select! Use glmer, glmmPQL or glmmML!")
}


#browser()

old.names<-attr(X,"dimnames")[[2]]


if(control$print.iter)
  print(paste("Iteration ", 1,sep=""))


if(is.list(rnd))
{
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
                        }), collapse=" + "), sep="+") 
  }
  
  
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
subject.names<-names(rnd)  
}else{
  W<-rnd
  if(!is.null(control$W.index))
  {
    s<-length(unique(control$W.index))
    n<-ncol(W)/s
  }else{
    n<-ncol(W)
    s<-1
  }
  
  
  newrndfrml<-NULL  
  random.labels<-attr(W,"W.name")   
  k<-NULL
  rnd.len <- 1
  attr(rnd,"names")<-colnames(rnd)
  subject.names<-colnames(rnd)  
}



very.old.names<-very.old.names[!is.na(control$index)]

block<-as.numeric(table(index.new[!is.na(index.new)]))

BLOCK<-FALSE
if(!all(block==1))
  BLOCK<-TRUE

lin<-ncol(X)

if(is.null(control$start))
  control$start<-c(rep(0,(lin+n%*%s)))

if(is.null(control$q_start))
{
  control$q_start<-rep(0.1,sum(s))
  if(sum(s)>1)
    control$q_start<-diag(control$q_start)
}

q_start<-control$q_start

N<-length(y)


beta_null<-control$start[1:lin]
attr(beta_null,"names")<-orig.names
beta_null<-beta_null[colnames(X)]

ranef_null<-control$start[(lin+1):(lin+n%*%s)]

Z_fastalles<-X
  
if(!control$overdispersion && family$family=="gaussian")
control$overdispersion<-T
  
#######################################################################  
######################## allow switch to Newton Raphson ###############
#######################################################################  
if(switch.NR)
{
#######################################################################  
###########################  1. No Smooth #############################  
#######################################################################  
if(is.null(control$smooth))
{  
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
    }else{
      Q_start<-q_start
    }
  }else{
    lin0<-sum(beta_null!=0)
    if(all(s==1))
    {
      Q_start<-diag(diag(q_start),sum(s))
    }else{
      Q_start<-matrix(0,sum(s),sum(s))
      Q_start[1:s[1],1:s[1]]<-q_start[1:s[1],1:s[1]]
      
      for (zu in 2:rnd.len)
      {
        Q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
      }
    }
  }
    
  U<-X[,!is.na(index.new),drop=FALSE]
  X<-X[,is.na(index.new),drop=FALSE]
  
  final.names<-c(colnames(X),colnames(U))
  
  q<-dim(X)[2]

  
  Z_alles<-cbind(X,U,W)
  ########################################################## some definitions ################################################
  Delta<-matrix(0,control$steps,(lin+n%*%s))
  Delta[1,1:lin]<-beta_null[final.names]
  Delta[1,(lin+1):(lin+n%*%s)]<-t(ranef_null)
  
  Eta.ma<-matrix(0,control$steps+1,N)
  Eta.ma[1,]<-Eta_start
  
  
  control$epsilon<-control$epsilon*sqrt(dim(Delta)[2])
  
  Delta_start<-Delta[1,]
  
  active_old<-!is.element(Delta[1,],0)
  
  Q_inv<-NULL
  Q_inv.old.temp<-NULL
  Q_inv.start<-NULL

  
  ## start value for overdispersion
  if(control$overdispersion)
  {
    if(!is.null(control$phi_start))
    {  
      phi<-control$phi_start
    }else{
      if(is.null(control$Q.phi.start))
      {
        Q.phi.start<-q_start 
      }else{
        Q.phi.start<-control$Q.phi.start
      }
       active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,n%*%s))
      Z_aktuell<-Z_alles[,active]
      lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
      
      if(rnd.len==1)
      {
        if(s==1)
        {
          Q.phi.start<-diag(Q.phi.start,s)
          P_akt<-c(rep(0,lin_akt),rep((Q.phi.start^(-1)),n*s))
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
        }else{
          Q_inv.start<-chol2inv(chol(Q.phi.start))
          P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
          for(jf in 1:n)
            P_akt[(lin_akt+(jf-1)*s+1):(lin_akt+jf*s),(lin_akt+(jf-1)*s+1):(lin_akt+jf*s)]<-Q_inv.start
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
        }
      }else{
        if(all(s==1))
        {
          Q.phi.start<-diag(diag(Q.phi.start),sum(s))
          P_akt<-c(rep(0,lin_akt),rep(diag(Q.phi.start)^(-1),n))
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
        }else{
          Q_inv.start<-list()
          Q_inv.start[[1]]<-chol2inv(chol(Q.phi.start[1:s[1],1:s[1]]))
          
          for (zu in 2:rnd.len)
          {
            Q_inv.start[[zu]]<-chol2inv(chol(Q.phi.start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
          }
          
          P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
          for(jf in 1:n[1])
            P_akt[(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1]),(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1])]<-Q_inv.start[[1]]
          
          for (zu in 2:rnd.len)
          {
            for(jf in 1:n[zu])
              P_akt[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                    (lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.start[[zu]]
          }
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
        }
      }
      InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
      if(class(InvFisher2)=="try-error")
        InvFisher2<-try(solve(F_gross),silent=T)
      if(class(InvFisher2)=="try-error")
        stop("Fisher matrix not invertible")  
      
      if(all(Mu==0) &family$family=="gaussian")
      {
        phi<-1  
      }else{
        Hat<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))#E-Uu
        phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(Hat)))
      }  
      Sigma<-Sigma*phi
    }  
  }else{
    phi<-1
  }
  
  
  Q<-list()
  Q[[1]]<-Q_start
  
  l=1
  
  
  score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)
  
  if (BLOCK)
  {
    grad.1<-gradient.lasso.block(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
  }else{
    grad.1<-gradient.lasso(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
  }
  
  score_vec<-c(score_vec[1:q],grad.1,score_vec[(lin+1):(lin+n%*%s)])
  
  if(rnd.len==1)
  {
    if(s==1)
    {
      P1<-c(rep(0,lin),rep((Q_start^(-1)),n*s))
    }else{
      Q_inv.start<-chol2inv(chol(Q_start))
      P1<-matrix(0,lin+n%*%s,lin+n%*%s)
      for(j in 1:n)
        P1[(lin+(j-1)*s+1):(lin+j*s),(lin+(j-1)*s+1):(lin+j*s)]<-Q_inv.start
    }
  }else{
    if(all(s==1))
    {
      P1<-c(rep(0,lin),rep(diag(Q_start)^(-1),n))
    }else{
      Q_inv.start<-list()
      Q_inv.start[[1]]<-chol2inv(chol(Q_start[1:s[1],1:s[1]]))
      P1<-matrix(0,lin+n%*%s,lin+n%*%s)
      for(jf in 1:n[1])
        P1[(lin+(jf-1)*s[1]+1):(lin+jf*s[1]),(lin+(jf-1)*s[1]+1):(lin+jf*s[1])]<-Q_inv.start[[1]]
      
      for (zu in 2:rnd.len)
      {
        Q_inv.start[[zu]]<-chol2inv(chol(Q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
        for(jf in 1:n[zu])
          P1[(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
             (lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.start[[zu]]
      }
    }
  }
  
  if(all(s==1))
  {
    F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+diag(P1)
  }else{
    F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
  }
  
  
  crit.obj<-t.change(grad=score_vec[(q+1):lin],b=Delta[1,(q+1):lin])
  t_edge<-crit.obj$min.rate
  
  grad.2<-t(score_vec)%*%F_gross%*%score_vec
  
  t_opt<-l2norm(score_vec)$length/grad.2
  nue<-control$nue

  
  half.index<-0
  solve.test<-FALSE
  ######### big while loop for testing if the update leads to Fisher matrix which can be inverted
  while(!solve.test)
  {  
    
    solve.test2<-FALSE  
    while(!solve.test2)
    {  
      
      if(half.index>50)
      {
        stop("Fisher matrix not invertible")
      }
      
  Delta[1,]<-Delta_start+min(t_opt,t_edge)*nue*(0.5^half.index)*score_vec
  if(t_opt>t_edge & half.index==0)
  Delta[1,crit.obj$whichmin+q]<-0  
  Eta<-Z_alles%*%Delta[1,]
  Mu<-as.vector(family$linkinv(Eta))
  Sigma<-as.vector(family$variance(Mu))
  D<-as.vector(family$mu.eta(Eta))
  
  active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,n%*%s))
  Z_aktuell<-Z_alles[,active]
  lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
  
  if (control$method=="EM" || control$overdispersion)
  {  
  if(rnd.len==1)
  {
    if(s==1)
    {
      P_akt<-c(rep(0,lin_akt),rep((Q_start^(-1)),n*s))
      F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
    }else{
      P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
      for(jf in 1:n)
        P_akt[(lin_akt+(jf-1)*s+1):(lin_akt+jf*s),(lin_akt+(jf-1)*s+1):(lin_akt+jf*s)]<-Q_inv.start
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
        P_akt[(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1]),(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1])]<-Q_inv.start[[1]]
      
      for (zu in 2:rnd.len)
      {
        for(jf in 1:n[zu])
          P_akt[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                (lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.start[[zu]]
      }
      F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
    }
  }
  InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
  if(class(InvFisher2)=="try-error")
    InvFisher2<-try(solve(F_gross),silent=T)
  if(class(InvFisher2)=="try-error")
  {
    #stop("Fisher matrix not invertible")  
    half.index<-half.index+1  
  }else{
    solve.test2<-TRUE 
  }}else{
    solve.test2<-TRUE
  }
 }
  
  betaact<-Delta[1,active]
  
  if (control$method=="EM")
  {   
    ############################# Q update ################
    if(rnd.len==1)
    {
      Q1<-InvFisher2[(lin_akt+1):(lin_akt+s),(lin_akt+1):(lin_akt+s)]+Delta[1,(lin+1):(lin+s)]%*%t(Delta[1,(lin+1):(lin+s)])
      for (i in 2:n)
        Q1<-Q1+InvFisher2[(lin_akt+(i-1)*s+1):(lin_akt+i*s),(lin_akt+(i-1)*s+1):(lin_akt+i*s)]+Delta[1,(lin+(i-1)*s+1):(lin+i*s)]%*%t(Delta[1,(lin+(i-1)*s+1):(lin+i*s)])
      Q1<-1/n*Q1
    }else{
      Q1<-matrix(0,sum(s),sum(s))
      Q1[1:s[1],1:s[1]]<-InvFisher2[(lin_akt+1):(lin_akt+s[1]),(lin_akt+1):(lin_akt+s[1])]+Delta[1,(lin+1):(lin+s[1])]%*%t(Delta[1,(lin+1):(lin+s[1])])
      for (i in 2:n[1])
        Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher2[(lin_akt+(i-1)*s[1]+1):(lin_akt+i*s[1]),(lin_akt+(i-1)*s[1]+1):(lin_akt+i*s[1])]+Delta[1,(lin+(i-1)*s[1]+1):(lin+i*s[1])]%*%t(Delta[1,(lin+(i-1)*s[1]+1):(lin+i*s[1])])
      Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]
      
      for (zu in 2:rnd.len)
      {
        Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher2[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
        for (i in 2:n[zu])
          Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher2[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
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
        optim.obj<-nlminb(sqrt(Q_start),likelihood_nlminb,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp)
        Q1<-as.matrix(optim.obj$par)^2
      }else{
        q_start_vec<-c(diag(q_start),q_start[lower.tri(q_start)])
        up1<-min(20,50*max(q_start_vec))
        upp<-rep(up1,length(q_start_vec))
        low<-c(rep(0,s),rep(-up1,0.5*(s^2-s)))
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
         
        for (zu in 2:rnd.len)
        {
          q_start_vec<-c(q_start_vec,c(diag(q_start)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])],q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]))
          up1<-min(20,50*max(q_start_vec))
          low<-c(low,c(rep(0,s[zu]),rep(-up1,0.5*(s[zu]^2-s[zu]))))
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
    FinalHat<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))#E-Uu
    phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(FinalHat)))
    Sigma<-Sigma*phi
  }
  

 Eta.old<-Eta
 
 
  vorz<-F

  score_vec2<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)
 
  score_old2<-score_vec2
 
 if (BLOCK)
 {
   grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
 }else{
   grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
 }
 score_vec2<-c(score_vec2[1:q],grad.1,score_vec2[(lin+1):(lin+n%*%s)])
 

 if(rnd.len==1)
 {
   if(s==1)
   {
     P1<-c(rep(0,lin),rep((Q1^(-1)),n*s))
   }else{
     Q_inv<-solve(Q1)
     P1<-matrix(0,lin+n%*%s,lin+n%*%s)
     for(j in 1:n)
       P1[(lin+(j-1)*s+1):(lin+j*s),(lin+(j-1)*s+1):(lin+j*s)]<-Q_inv
   }
 }else{
   if(all(s==1))
   {
     P1<-c(rep(0,lin),rep(diag(Q1)^(-1),n))
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
   }
 }
 
 
 if(all(s==1))
 {
   F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+diag(P1)
 }else{
   F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
 }
 
 crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin])
 t_edge<-crit.obj$min.rate
 
 grad.2<-t(score_vec2)%*%F_gross%*%score_vec2
 
 t_opt<-l2norm(score_vec2)$length/grad.2
 tryNR<- (t_opt<t_edge) #&& !(all(active_old==active)  && !NRstep)  
 
 if(tryNR) 
 {
   lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
   if(rnd.len==1)
   {
     if(s==1)
     {
       P_akt<-c(rep(0,lin_akt),rep((Q1^(-1)),n*s))
       F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
     }else{
       P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
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
   

   InvFisher<-try(chol2inv(chol(F_gross)),silent=T)
   if(class(InvFisher)=="try-error")
     InvFisher<-try(solve(F_gross),silent=T)
   if(class(InvFisher)=="try-error")
   {
     half.index<-half.index+1  
   }else{
     solve.test<-TRUE 
     Delta.test<-Delta[1,active]+nue*InvFisher%*%score_old2[active]
     if(lin_akt>q)
     {
       vorz<-all(sign(Delta.test[(q+1):lin_akt])==sign(betaact[(q+1):lin_akt]))  
     }else{
       vorz<-T  
     }
   }
 }else{
   solve.test<-TRUE  
 }
}   
 
Eta.ma[2,]<-Eta

score_old<-score_old2
 score_vec<-score_vec2
 Q1.old<-Q1
 Q_inv.old<-Q_inv    

 Q1.very.old<-Q_start
 Q_inv.very.old<-Q_inv.start

###############################################################################################################################################
################################################################### Main Iteration ###################################################################
if(control$steps!=1)
  {
    for (l in 2:control$steps)
    {
      if(control$print.iter)
        print(paste("Iteration ", l,sep=""))
      
      if(!vorz)
        tryNR<-F
      
       
      half.index<-0
      solve.test<-FALSE
      ######### big while loop for testing if the update leads to Fisher matrix which can be inverted
      while(!solve.test)
      {  
        
        solve.test2<-FALSE  
        while(!solve.test2)
        {  
          if(half.index>50)
          {
            half.index<-Inf;Q1.old<-Q1.very.old;Q_inv.old<-Q_inv.very.old
          }
                    
       if(tryNR)
      {
        Delta[l,active]<-Delta[l-1,active]+nue*(0.5^half.index)*InvFisher%*%score_old[active]
        NRstep<-T
 #       print("NR-step!!!")
       }else{
        Delta[l,]<-Delta[l-1,]+min(t_opt,t_edge)*nue*(0.5^half.index)*nue*score_vec
        if(t_opt>t_edge & half.index==0)
          Delta[l,crit.obj$whichmin+q]<-0  
        NRstep<-F
      }
      
      Eta<-Z_alles%*%Delta[l,]
      Mu<-as.vector(family$linkinv(Eta))
      Sigma<-as.vector(family$variance(Mu))
      D<-as.vector(family$mu.eta(Eta))
      
      active_old<-active
      active<-c(rep(T,q),!is.element(Delta[l,(q+1):lin],0),rep(T,n%*%s))
      Z_aktuell<-Z_alles[,active]
      lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))
      
      if (control$method=="EM" || control$overdispersion)
      {  
      if(rnd.len==1)
      {
        if(s==1)
        {
          P_akt<-c(rep(0,lin_akt),rep((Q1.old^(-1)),n*s))
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
        }else{
          P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
          for(jf in 1:n)
            P_akt[(lin_akt+(jf-1)*s+1):(lin_akt+jf*s),(lin_akt+(jf-1)*s+1):(lin_akt+jf*s)]<-Q_inv.old
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
        }
      }else{
        if(all(s==1))
        {
          P_akt<-c(rep(0,lin_akt),rep(diag(Q1.old)^(-1),n))
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
        }else{
          P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
          for(jf in 1:n[1])
            P_akt[(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1]),(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1])]<-Q_inv.old[[1]]
          
          for (zu in 2:rnd.len)
          {
            for(jf in 1:n[zu])
              P_akt[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                    (lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.old[[zu]]
          }
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
        }
      }
      InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
      if(class(InvFisher2)=="try-error")
        InvFisher2<-try(solve(F_gross),silent=T)
      if(class(InvFisher2)=="try-error")
      {
         half.index<-half.index+1  
      }else{
        solve.test2<-TRUE 
      }}else{
        solve.test2<-TRUE 
      }
    }

      betaact<-Delta[l,active]
      
      if (control$method=="EM")
      {        
        ############################# Q update ################
        if(rnd.len==1)
        {
          Q1<-InvFisher2[(lin_akt+1):(lin_akt+s),(lin_akt+1):(lin_akt+s)]+Delta[l,(lin+1):(lin+s)]%*%t(Delta[l,(lin+1):(lin+s)])
          for (i in 2:n)
            Q1<-Q1+InvFisher2[(lin_akt+(i-1)*s+1):(lin_akt+i*s),(lin_akt+(i-1)*s+1):(lin_akt+i*s)]+Delta[l,(lin+(i-1)*s+1):(lin+i*s)]%*%t(Delta[l,(lin+(i-1)*s+1):(lin+i*s)])
          Q1<-1/n*Q1
        }else{
          Q1<-matrix(0,sum(s),sum(s))
          Q1[1:s[1],1:s[1]]<-InvFisher2[(lin_akt+1):(lin_akt+s[1]),(lin_akt+1):(lin_akt+s[1])]+Delta[l,(lin+1):(lin+s[1])]%*%t(Delta[l,(lin+1):(lin+s[1])])
          for (i in 2:n[1])
            Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher2[(lin_akt+(i-1)*s[1]+1):(lin_akt+i*s[1]),(lin_akt+(i-1)*s[1]+1):(lin_akt+i*s[1])]+Delta[l,(lin+(i-1)*s[1]+1):(lin+i*s[1])]%*%t(Delta[l,(lin+(i-1)*s[1]+1):(lin+i*s[1])])
          Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]
          
          for (zu in 2:rnd.len)
          {
            Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher2[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
            for (i in 2:n[zu])
              Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher2[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
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
            
            optim.obj<-nlminb(sqrt(Q1),likelihood_nlminb,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp)
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
        FinalHat<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))#E-Uu
        phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(FinalHat)))
        Sigma<-Sigma*phi
      }
      
      Q[[l+1]]<-Q1
            
    score_vec2<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)
    
   score_old2<-score_vec2

 if (BLOCK)
 {
   grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda,block=block)
 }else{
   grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda)
 }
 score_vec2<-c(score_vec2[1:q],grad.1,score_vec2[(lin+1):(lin+n%*%s)])
 
 if(rnd.len==1)
 {
   if(s==1)
   {
     P1<-c(rep(0,lin),rep((Q1^(-1)),n*s))
   }else{
     Q_inv<-solve(Q1)
     P1<-matrix(0,lin+n%*%s,lin+n%*%s)
     for(j in 1:n)
       P1[(lin+(j-1)*s+1):(lin+j*s),(lin+(j-1)*s+1):(lin+j*s)]<-Q_inv
   }
 }else{
   if(all(s==1))
   {
     P1<-c(rep(0,lin),rep(diag(Q1)^(-1),n))
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
   }
 }
 
 
 if(all(s==1))
 {
   F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+diag(P1)
 }else{
   F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
 }
 
 crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin])
 t_edge<-crit.obj$min.rate

 grad.2<-t(score_vec2)%*%F_gross%*%score_vec2
 
 t_opt<-l2norm(score_vec2)$length/grad.2
 tryNR<- (t_opt<t_edge) && !(all(active_old==active)  && !NRstep)  
 
 if(tryNR) 
 {
   lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))
   if(rnd.len==1)
   {
     if(s==1)
     {
       P_akt<-c(rep(0,lin_akt),rep((Q1^(-1)),n*s))
       F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
     }else{
       P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
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
   InvFisher3<-try(chol2inv(chol(F_gross)),silent=T)
   if(class(InvFisher3)=="try-error")
     InvFisher3<-try(solve(F_gross),silent=T)
   if(class(InvFisher3)=="try-error")
   {
     half.index<-half.index+1  
   }else{
     
     solve.test<-TRUE 
     Delta.test<-Delta[l,active]+nue*InvFisher3%*%score_old2[active]
     if(lin_akt>q)
     {
       vorz<-all(sign(Delta.test[(q+1):lin_akt])==sign(betaact[(q+1):lin_akt]))  
     }else{
       vorz<-T  
     }
   }
 }else{
   solve.test<-TRUE  
 }
}  
 
  if(tryNR) 
  InvFisher<-InvFisher3

   score_old<-score_old2
   score_vec<-score_vec2
   Q1.very.old<-Q1.old
   Q_inv.very.old<-Q_inv.old
   Q1.old<-Q1
   Q_inv.old<-Q_inv    

   Eta.ma[l+1,]<-Eta


   finish<-(sqrt(sum((Eta.ma[l,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l,])^2))<control$epsilon)
   finish2<-(sqrt(sum((Eta.ma[l-1,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l-1,])^2))<control$epsilon)
   if(finish ||  finish2) #|| (all(grad.1 == 0) ))
      break
    Eta.old<-Eta
}}

FinalHat.df<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))
df<-sum(diag(FinalHat.df))

  
conv.step<-l
phi.med<-phi

if(conv.step==control$steps)
{
  cat("Warning:\n")
  cat("Algorithm did not converge!\n")
}


  Delta_neu<-Delta[l,]
  Eta_opt<-Z_alles%*%Delta_neu
  Mu_opt<-as.vector(family$linkinv(Eta_opt))
  Sigma_opt<-as.vector(family$variance(Mu_opt))    
  D_opt<-as.vector(family$mu.eta(Eta_opt))
  Qfinal<-Q[[l+1]]
  
if(rnd.len==1)
{
  if(s==1)
  {
    P1.ran<-rep((Qfinal^(-1)),n*s)
    P1.ran<-diag(P1.ran)
  }else{
    P1.ran<-matrix(0,n*s,n*s)
    for(jf in 1:n)
      P1.ran[((jf-1)*s+1):(jf*s),((jf-1)*s+1):(jf*s)]<-chol2inv(chol(Qfinal))
  }
}else{
  if(all(s==1))
  {
    P1.ran<-rep(diag(Qfinal)^(-1),n)
    P1.ran<-diag(P1.ran)
  }else{
    P1.ran<-matrix(0,n%*%s,n%*%s)
    inv.act<-chol2inv(chol(Qfinal[1:s[1],1:s[1]]))
    for(jf in 1:n[1])
      P1.ran[( (jf-1)*s[1]+1):( jf*s[1]),( (jf-1)*s[1]+1):( jf*s[1])]<-inv.act
  
    for (zu in 2:rnd.len)
    {
      inv.act<-chol2inv(chol(Qfinal[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
      for(jf in 1:n[zu])
        P1.ran[( n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):( n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
               ( n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):( n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-inv.act
    }
  }
}
ranef.logLik<--0.5*t(Delta_neu[(lin+1):(lin+n%*%s)])%*%P1.ran%*%Delta_neu[(lin+1):(lin+n%*%s)]


if(final.re)
{    
############ final re-estimation
  aaa<-!is.element(Delta_neu[1:(lin)],0)
  
  if(rnd.len==1 && s==1)
  {  
  Q.max<-max(sqrt(unlist(Q)))
  Q.min<-min(sqrt(unlist(Q)))
  }else{
  Q.max<-max(Qfinal)+1
  Q.min<-min(Qfinal)-1e-10
  }

  if(rnd.len==1)
  {
    glmm_fin<-try(glmm_final(y,Z_fastalles[,aaa],W,k,n,q_start=Qfinal,
                             Delta_start=Delta_neu[c(aaa,rep(T,n%*%s))],s,steps=control$maxIter,
                             family=family,method=control$method.final,overdispersion=control$overdispersion,
                             phi=phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final,
                             Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac))#,silent = TRUE)
    if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
    {  
    glmm_fin2<-try(glmm_final(y,Z_fastalles[,aaa],W,k,n,q_start=q_start,
                             Delta_start=Delta_start[c(aaa,rep(T,n%*%s))],s,steps=control$maxIter,
                             family=family,method=control$method.final,overdispersion=control$overdispersion,
                             phi=control$phi,print.iter.final=control$print.iter.final,
                             eps.final=control$eps.final,Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac))
        if(class(glmm_fin2)!="try-error")
        {    
          if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
            glmm_fin<-glmm_fin2 
        }
    }
  }else{
    glmm_fin<-try(glmm_final_multi_random(y,Z_fastalles[,aaa],W,k,q_start=Qfinal,
                              Delta_start=Delta_neu[c(aaa,rep(T,n%*%s))],s,n,steps=control$maxIter,
                              family=family,method=control$method.final,overdispersion=control$overdispersion,
                              phi=phi,rnd.len=rnd.len,print.iter.final=control$print.iter.final,
                              eps.final=control$eps.final,Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac))#,silent = TRUE)
    if(class(glmm_fin)=="try-error"|| glmm_fin$opt>control$maxIter-10)
    {  
    glmm_fin2<-try(glmm_final_multi_random(y,Z_fastalles[,aaa],W,k,q_start=q_start,
                              Delta_start=Delta_start[c(aaa,rep(T,n%*%s))],s,n,steps=control$maxIter,
                              family=family,method=control$method.final,overdispersion=control$overdispersion,
                              phi=control$phi,rnd.len=rnd.len,print.iter.final=control$print.iter.final,
                              eps.final=control$eps.final,Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac))    
        if(class(glmm_fin2)!="try-error")
        {    
          if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
            glmm_fin<-glmm_fin2 
        }
    }
  }
  
  if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
  {
    cat("Warning:\n")
    cat("Final Fisher scoring reestimation did not converge!\n")
  }
#######
}else{
glmm_fin<-NA  
class(glmm_fin)<-"try-error"
}  
  Delta_neu2<-Delta_neu
  Standard_errors<-rep(NA,length(Delta_neu))
  
  if(class(glmm_fin)!="try-error")
  {
    Delta_neu2[c(aaa,rep(T,n%*%s))]<-glmm_fin$Delta
    Standard_errors[c(aaa,rep(T,n%*%s))]<-glmm_fin$Standard_errors
    Qfinal<-glmm_fin$Q
    phi<-glmm_fin$phi
    complexity<-glmm_fin$complexity
  }else{
    glmm_fin<-list()
    glmm_fin$ranef.logLik<-ranef.logLik
    complexity<-df
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
    
    if(s[1]==1)
      Qfinal[[1]]<-sqrt(Qfinal[[1]])
    
    for (zu in 2:rnd.len)
    {
      Qfinal[[zu]]<-as.matrix(Qfinal_old[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])
      colnames(Qfinal[[zu]])<-random.labels[[zu]]
      rownames(Qfinal[[zu]])<-random.labels[[zu]]

      if(s[zu]==1)
        Qfinal[[zu]]<-sqrt(Qfinal[[zu]])
      
    }
  }
  
  names(Delta_neu)[1:dim(X)[2]]<-colnames(X)
  names(Standard_errors)[1:dim(X)[2]]<-colnames(X)
  
  if(lin>1)
  {
    names(Delta_neu)[(dim(X)[2]+1):lin]<-colnames(U)
    names(Standard_errors)[(dim(X)[2]+1):lin]<-colnames(U)
  }
  
  names(Delta_neu)[(lin+1):(lin+n%*%s)]<-colnames(W)
  names(Standard_errors)[(lin+1):(lin+n%*%s)]<-colnames(W)
  colnames(Delta)<-c(final.names,colnames(W))
  
  aic<-NaN
  bic<-NaN
  
if (is.element(family$family,c("gaussian", "binomial", "poisson"))) 
{

  loglik<-logLik.glmmLasso(y=y,mu=Mu_opt,beta=Delta_neu[(q+1):(lin)],ranef.logLik=glmm_fin$ranef.logLik,
                           lambda=lambda,family=family,penal=FALSE)
  

if(control$complexity!="hat.matrix")  
{  
  if(rnd.len==1)
    {
      complexity<-0.5*(s*(s+1))
    }else{
      complexity<-0.5*(s[1]*(s[1]+1))
      for(zu in 2:rnd.len)
        complexity<-complexity+0.5*(s[zu]*(s[zu]+1))
    }
    complexity<-complexity+sum(Delta_neu[1:(lin)]!=0)
}    
      aic<--2*loglik+2*complexity
      bic<--2*loglik+log(N)*complexity
}else{
  warning("For the specified family (so far) no AIC and BIC are available!")  
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
  ret.obj$fix<-fix.old
  ret.obj$conv.step<-conv.step
  ret.obj$newrndfrml<-newrndfrml
  ret.obj$subject<-subject.names
  ret.obj$data<-data
  ret.obj$rnd.len<-rnd.len
  ret.obj$phi.med<-phi.med
  ret.obj$y <- y
  ret.obj$df<-df
  ret.obj$loglik<-loglik
  return(ret.obj)
##############################################################  
######################## 2. Smooth ###########################  
##############################################################  
}else{
  
  smooth<-control$smooth
  
  if(attr(terms(smooth$formula), "intercept")==0)
  {
    variables<-attr(terms(smooth$formula),"term.labels")
    smooth$formula<- "~ +1"
    for (ir in 1:length(variables))
      smooth$formula <- paste(smooth$formula, variables[ir], sep="+")
    smooth$formula <-formula(smooth$formula)
  }
  
  B <- model.matrix(smooth$formula, data)
  B.names<-attr(B,"dimnames")[[2]]
  
  B<-as.matrix(B[,-1])
  attr(B,"dimnames")[[2]]<-B.names[-1]
  
  nbasis<-smooth$nbasis
  diff.ord<-smooth$diff.ord
  spline.degree<-smooth$spline.degree
  knots.no<-nbasis-1
  if(spline.degree<3 && (spline.degree-diff.ord)<2)
  knots.no<-knots.no+1  
  penal<-smooth$penal
    
  if(!(diff.ord<spline.degree))
  stop("Order of differences must be lower than degree of B-spline polynomials!")
  
  m<-dim(B)[2]
  
  Phi<-numeric()
  
  for (r in 1:m)
  {
    Basis<-bs.design(B[,r],diff.ord=diff.ord,spline.degree=spline.degree,knots.no=knots.no)
    Phi_temp<-cbind(Basis$X[,-1],Basis$Z)
    colnames(Phi_temp)<-paste(colnames(B)[r],rep(1:dim(Phi_temp)[2],each=1), sep=".")
    Phi<-cbind(Phi,Phi_temp)
  }

  dim.smooth<-dim(Phi)[2]
  
  if(is.null(smooth$start))
  smooth$start<-rep(0,dim.smooth)  
  
  smooth_null<-smooth$start
  
  if(lin>1)
  {
    Eta_start<-X%*%beta_null+Phi%*%smooth_null+W%*%ranef_null
  }else{
    Eta_start<-rep(beta_null,N)+Phi%*%smooth_null+W%*%ranef_null
  }
  
  D<-as.vector(family$mu.eta(Eta_start))
  Mu<-as.vector(family$linkinv(Eta_start))
  Sigma<-as.vector(family$variance(Mu))
  
  if(rnd.len==1)
  {
    if(s==1)
    {
      Q_start<-diag(q_start,s)
    }else{
      Q_start<-q_start
    }
  }else{
    if(all(s==1))
    {
      Q_start<-diag(diag(q_start),sum(s))
    }else{
      Q_start<-matrix(0,sum(s),sum(s))
      Q_start[1:s[1],1:s[1]]<-q_start[1:s[1],1:s[1]]
      
      for (zu in 2:rnd.len)
      {
        Q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
      }
    }
  }
  
  U<-X[,!is.na(index.new),drop=FALSE]
  X<-X[,is.na(index.new),drop=FALSE]
  
  final.names<-c(colnames(X),colnames(U))
  
  q<-dim(X)[2]
  
  Z_alles<-cbind(X,U,Phi,W)
  ########################################################## some definitions ################################################
  Delta<-matrix(0,control$steps,(lin+dim.smooth+n%*%s))
  Delta[1,1:lin]<-beta_null[final.names]
  Delta[1,(lin+1):(lin+dim.smooth)]<-smooth_null
  Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]<-t(ranef_null)
  
  control$epsilon<-control$epsilon*sqrt(dim(Delta)[2])
  
  Delta_start<-Delta[1,]
  
  active_old<-!is.element(Delta[1,],0)

  Eta.ma<-matrix(0,control$steps+1,N)
  Eta.ma[1,]<-Eta_start
  
  
  Q_inv<-NULL
  Q_inv.old.temp<-NULL
  Q_inv.start<-NULL
    
  ## start value for overdispersion
  if(control$overdispersion)
  {
    if(!is.null(control$phi_start))
    {  
      phi<-control$phi_start
    }else{
      if(is.null(control$Q.phi.start))
      {
        Q.phi.start<-q_start 
      }else{
        Q.phi.start<-control$Q.phi.start
      }
      active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,dim.smooth+n%*%s))
      Z_aktuell<-Z_alles[,active]
      lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
      
      if(rnd.len==1)
      {
        if(s==1)
        {
          Q.phi.start<-diag(Q.phi.start,s)
          P_akt<-c(rep(0,lin_akt+dim.smooth),rep((Q.phi.start^(-1)),n*s))
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
        }else{
          Q_inv.start<-chol2inv(chol(Q.phi.start))
          P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
          for(jf in 1:n)
            P_akt[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),
                  (lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv.start
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
        }
      }else{
        if(all(s==1))
        {
          Q.phi.start<-diag(diag(Q.phi.start),sum(s))
          P_akt<-c(rep(0,lin_akt+dim.smooth),rep(diag(Q.phi.start)^(-1),dim.smooth+n))
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
        }else{
          Q_inv.start<-list()
          Q_inv.start[[1]]<-chol2inv(chol(Q.phi.start[1:s[1],1:s[1]]))
          
          for (zu in 2:rnd.len)
          {
            Q_inv.start[[zu]]<-chol2inv(chol(Q.phi.start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
          }
          
          P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
          for(jf in 1:n[1])
            P_akt[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),
                  (lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv.start[[1]]
          
          for (zu in 2:rnd.len)
          {
            for(jf in 1:n[zu])
              P_akt[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                    (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.start[[zu]]
          }
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
        }
      }
      InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
      if(class(InvFisher2)=="try-error")
        InvFisher2<-try(solve(F_gross),silent=T)
      if(class(InvFisher2)=="try-error")
        stop("Fisher matrix not invertible")  
      
      if(all(Mu==0) &family$family=="gaussian")
      {
        phi<-1  
      }else{
        Hat<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))#E-Uu
      phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(Hat)))
      }
      Sigma<-Sigma*phi
    }  
  }else{
    phi<-1
  }
  
  
  
  
  
  Q<-list()
  Q[[1]]<-Q_start
  
  l=1
  
  if(diff.ord>1)
  {
    k2<-c(rep(0,diff.ord-1),rep(1,nbasis-diff.ord+1))
  }else{
    k2<-rep(1,nbasis)
  }
  k22<-rep(k2,m)
  penal.vec<-penal*k22
    
  Q_inv<-NULL
  Q_inv.old.temp<-NULL
  Q_inv.start<-NULL
  
  #if(rnd.len==1)
  #{
  #  if(s==1)
  #  {
  #    P.smooth<-c(rep(0,lin),penal.vec,rep(0,n*s))
  #  }else{
  #    P.smooth<-c(rep(0,lin),penal.vec,rep(0,n%*%s))
  #  }
  #}else{
  #  if(all(s==1))
  #  {
  #    P.smooth<-c(rep(0,lin),penal.vec,rep(0,rnd.len*n))
   # }else{
  #    P.smooth<-c(rep(0,lin),penal.vec,rep(0,n%*%s))
  #  }
  #}

  score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)#-t(t(P.smooth)*Delta[1,])
  lambda.max<-max(abs(score_vec[(q+1):lin]))
  
  
  if (BLOCK)
  {
    grad.1<-gradient.lasso.block(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
  }else{
    grad.1<-gradient.lasso(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
  }
  
  score_vec<-c(score_vec[1:q],grad.1,score_vec[(lin+1):(lin+dim.smooth+n%*%s)])
  
  if(rnd.len==1)
  {
    if(s==1)
    {
      P1<-c(rep(0,lin),penal.vec,rep((Q_start^(-1)),n*s))
    }else{
      Q_inv.start<-chol2inv(chol(Q_start))
      P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
      diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
      for(j in 1:n)
        P1[(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s),(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s)]<-Q_inv
    }
  }else{
    if(all(s==1))
    {
      P1<-c(rep(0,lin),penal.vec,rep(diag(Q_start)^(-1),n))
    }else{
      Q_inv.start<-list()
      Q_inv.start[[1]]<-chol2inv(chol(Q_start[1:s[1],1:s[1]]))
      P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
      diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
      for(jf in 1:n[1])
        P1[(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1]),(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1])]<-Q_inv[[1]]
      
      for (zu in 2:rnd.len)
      {
        Q_inv.start[[zu]]<-chol2inv(chol(Q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
        for(jf in 1:n[zu])
          P1[(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
             (lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
      }
    }
  }

  if(all(s==1))
  {
    F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+diag(P1)
  }else{
    F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
  }
  
  crit.obj<-t.change(grad=score_vec[(q+1):lin],b=Delta[1,(q+1):lin])
  t_edge<-crit.obj$min.rate
  
  grad.2<-t(score_vec)%*%F_gross%*%score_vec
  
  t_opt<-l2norm(score_vec)$length/grad.2
  nue<-control$nue

  half.index<-0
  solve.test<-FALSE
  ######### big while loop for testing if the update leads to Fisher matrix which can be inverted
  while(!solve.test)
  {  
    
    solve.test2<-FALSE  
    while(!solve.test2)
    {  
      
  if(half.index>50)
  {
        stop("Fisher matrix not invertible")
  }
      
  Delta[1,]<-Delta_start+min(t_opt,t_edge)*nue*(0.5^half.index)*score_vec
  if(t_opt>t_edge & half.index==0)
    Delta[1,crit.obj$whichmin+q]<-0  
  Eta<-Z_alles%*%Delta[1,]
  Mu<-as.vector(family$linkinv(Eta))
  Sigma<-as.vector(family$variance(Mu))
  D<-as.vector(family$mu.eta(Eta))
  
  active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,dim.smooth+n%*%s))
  Z_aktuell<-Z_alles[,active]
  lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
  
  if (control$method=="EM" || control$overdispersion)
  {  
  if(rnd.len==1)
  {
    if(s==1)
    {
      P_akt<-c(rep(0,lin_akt),penal.vec,rep((Q_start^(-1)),n*s))
      F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
    }else{
      P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
      diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
      for(jf in 1:n)
        P_akt[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv.start
      F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
    }
  }else{
    if(all(s==1))
    {
      P_akt<-c(rep(0,lin_akt),penal.vec,rep(diag(Q_start)^(-1),n))
      F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
    }else{
      P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
      diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
      for(jf in 1:n[1])
        P_akt[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv.start[[1]]
      
      for (zu in 2:rnd.len)
      {
        for(jf in 1:n[zu])
          P_akt[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.start[[zu]]
      }
      F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
    }
  }
  InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
  if(class(InvFisher2)=="try-error")
    InvFisher2<-try(solve(F_gross),silent=T)
  if(class(InvFisher2)=="try-error")
  {
    #stop("Fisher matrix not invertible")  
    half.index<-half.index+1  
  }else{
    solve.test2<-TRUE 
  }}else{
    solve.test2<-TRUE
  }
  }
  
  betaact<-Delta[1,active]
  
  if (control$method=="EM")
  {   
    ############################# Q update ################
    if(rnd.len==1)
    {
      Q1<-InvFisher2[(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s),(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s)]+Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+s)]%*%t(Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+s)])
      for (i in 2:n)
        Q1<-Q1+InvFisher2[(lin_akt+dim.smooth+(i-1)*s+1):(lin_akt+dim.smooth+i*s),(lin_akt+dim.smooth+(i-1)*s+1):(lin_akt+dim.smooth+i*s)]+Delta[1,(lin+dim.smooth+(i-1)*s+1):(lin+dim.smooth+i*s)]%*%t(Delta[1,(lin+dim.smooth+(i-1)*s+1):(lin+dim.smooth+i*s)])
      Q1<-1/n*Q1
    }else{
      Q1<-matrix(0,sum(s),sum(s))
      Q1[1:s[1],1:s[1]]<-InvFisher2[(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s[1]),(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s[1])]+Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+s[1])]%*%t(Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+s[1])])
      for (i in 2:n[1])
        Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher2[(lin_akt+dim.smooth+(i-1)*s[1]+1):(lin_akt+dim.smooth+i*s[1]),(lin_akt+dim.smooth+(i-1)*s[1]+1):(lin_akt+dim.smooth+i*s[1])]+Delta[1,(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])]%*%t(Delta[1,(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])])
      Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]
      
      for (zu in 2:rnd.len)
      {
        Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher2[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
        for (i in 2:n[zu])
          Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher2[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
        Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-1/n[zu]*Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
      }
    }
  }else{
    Eta_tilde<-Eta+(y-Mu)*1/D
    Betadach<-Delta[1,1:(lin+dim.smooth)]     
    aktuell_vec<-!is.element(Delta[1,1:(lin)],0)
    X_aktuell<-Z_fastalles[,aktuell_vec]
    
    if(rnd.len==1)
    {
      
      if(s==1)
      {
        upp<-min(20,50*Q_start)
        low<-1e-14
        optim.obj<-try(nlminb(sqrt(Q_start),likelihood_nlminb,D=D,Sigma=Sigma,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),
                          Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp))
        if(class(optim.obj)=="try-error")
        optim.obj<-try(bobyqa(sqrt(Q_start),likelihood_nlminb,D=D,Sigma=Sigma,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),
                          Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp))
        Q1<-as.matrix(optim.obj$par)^2
      }else{
        q_start_vec<-c(diag(q_start),q_start[lower.tri(q_start)])
        up1<-min(20,50*max(q_start_vec))#[(s+1):(s*(s+1)*0.5)]))
        upp<-rep(up1,length(q_start_vec))
        low<-c(rep(0,s),rep(-up1,0.5*(s^2-s)))
        #   kkk_vec<-c(rep(-1,s),rep(0.5,0.5*(s^2-s)))
        optim.obj<-try(bobyqa(q_start_vec,likelihood,D=D,Sigma=Sigma,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),
                              Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp))
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
        optim.obj<-try(bobyqa(sqrt(q_start_vec),likelihood_diag,D=D,Sigma=Sigma,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
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
        optim.obj<-try(bobyqa(q_start_vec,likelihood_block,D=D,Sigma=Sigma,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
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
    FinalHat<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))#E-Uu
    phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(FinalHat)))
    Sigma<-Sigma*phi
  }
      
  Eta.old<-Eta
    
  vorz<-F
    
  score_vec2<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)#-t(t(P.smooth)*Delta[1,])
  lambda.max<-max(abs(score_vec2[(q+1):lin]))
  
  score_old2<-score_vec2
  
  if (BLOCK)
  {
    grad.1<-gradient.lasso.block(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
  }else{
    grad.1<-gradient.lasso(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
  }
  score_vec2<-c(score_vec2[1:q],grad.1,score_vec2[(lin+1):(lin+dim.smooth+n%*%s)])
  
  if(rnd.len==1)
  {
    if(s==1)
    {
      P1<-c(rep(0,lin),penal.vec,rep((Q1^(-1)),n*s))
    }else{
      Q_inv<-chol2inv(chol(Q1))
      Q_inv.old.temp<-Q_inv
      P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
      diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
      for(j in 1:n)
        P1[(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s),(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s)]<-Q_inv
    }
  }else{
    if(all(s==1))
    {
      P1<-c(rep(0,lin),penal.vec,rep(diag(Q1)^(-1),n))
    }else{
      Q_inv<-list()
      Q_inv[[1]]<-chol2inv(chol(Q1[1:s[1],1:s[1]]))
      P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
      diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
      for(jf in 1:n[1])
        P1[(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1]),(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1])]<-Q_inv[[1]]
      
      for (zu in 2:rnd.len)
      {
        Q_inv[[zu]]<-chol2inv(chol(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
        for(jf in 1:n[zu])
          P1[(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
             (lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
      }
    }
  }
  
  if(all(s==1))
  {
    F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+diag(P1)
  }else{
    F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
  }
  
  crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin])
  t_edge<-crit.obj$min.rate
  
  grad.2<-t(score_vec2)%*%F_gross%*%score_vec2
  
  t_opt<-l2norm(score_vec2)$length/grad.2
  tryNR<- (t_opt<t_edge) #&& !(all(active_old==active)  && !NRstep)
  
  if(tryNR) 
  {
    lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
    if(rnd.len==1)
    {
      if(s==1)
      {
        P_akt<-c(rep(0,lin_akt),penal.vec,rep((Q1^(-1)),n*s))
        F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
      }else{
        P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
        diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
        for(jf in 1:n)
          P_akt[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),
                (lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv
        F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
      }
    }else{
      if(all(s==1))
      {
        P_akt<-c(rep(0,lin_akt),penal.vec,rep(diag(Q1)^(-1),n))
        F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
      }else{
        P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
        diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
        for(jf in 1:n[1])
          P_akt[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),
                (lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv[[1]]
        
        for (zu in 2:rnd.len)
        {
          for(jf in 1:n[zu])
            P_akt[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                  (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
        }
        F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
      }
    }
    InvFisher<-try(chol2inv(chol(F_gross)),silent=T)
    if(class(InvFisher)=="try-error")
      InvFisher<-try(solve(F_gross),silent=T)
    if(class(InvFisher)=="try-error")
      {
        half.index<-half.index+1  
      }else{
        solve.test<-TRUE 
        Delta.test<-Delta[1,active]+nue*InvFisher%*%score_old2[active]
        if(lin_akt>q)
        {
          vorz<-all(sign(Delta.test[(q+1):lin_akt])==sign(betaact[(q+1):lin_akt]))  
        }else{
          vorz<-T  
        }
      }
}else{
  solve.test<-TRUE  
}
}   
  
Eta.ma[2,]<-Eta

score_old<-score_old2
score_vec<-score_vec2
Q1.old<-Q1
Q_inv.old<-Q_inv    

Q1.very.old<-Q_start
Q_inv.very.old<-Q_inv.start

###############################################################################################################################################
################################################################### Main Iteration ###################################################################
if(control$steps!=1)
{
  for (l in 2:control$steps)
  {
    if(control$print.iter)
      print(paste("Iteration ", l,sep=""))
      
    if(!vorz)
      tryNR<-F
    
half.index<-0
solve.test<-FALSE
######### big while loop for testing if the update leads to Fisher matrix which can be inverted
while(!solve.test)
{  
        
    solve.test2<-FALSE  
    while(!solve.test2)
    {  
     
    if(half.index>50)
    {
      half.index<-Inf;Q1.old<-Q1.very.old;Q_inv.old<-Q_inv.very.old;
    }
    
      if(tryNR)
      {
        Delta[l,active]<- Delta[l-1,active]+nue*(0.5^half.index)*InvFisher%*%score_old[active]
        NRstep<-T
      }else{
        Delta[l,]<-Delta[l-1,]+min(t_opt,t_edge)*nue*(0.5^half.index)*score_vec
        if(t_opt>t_edge & half.index==0)
          Delta[l,crit.obj$whichmin+q]<-0  
        NRstep<-F
      }
      
      Eta<-Z_alles%*%Delta[l,]
      Mu<-as.vector(family$linkinv(Eta))
      Sigma<-as.vector(family$variance(Mu))
      D<-as.vector(family$mu.eta(Eta))
      
      active_old<-active
      active<-c(rep(T,q),!is.element(Delta[l,(q+1):lin],0),rep(T,dim.smooth+n%*%s))
      Z_aktuell<-Z_alles[,active]
      lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))
      
      if (control$method=="EM" || control$overdispersion)
      {  
      if(rnd.len==1)
      {
        if(s==1)
        {
          P_akt<-c(rep(0,lin_akt),penal.vec,rep((Q1.old^(-1)),n*s))
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
        }else{
          P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
          diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
          for(jf in 1:n)
            P_akt[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),
                  (lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv.old
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
        }
      }else{
        if(all(s==1))
        {
          P_akt<-c(rep(0,lin_akt),penal.vec,rep(diag(Q1.old)^(-1),n))
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
        }else{
          P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
          diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
          for(jf in 1:n[1])
            P_akt[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),
                  (lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv.old[[1]]
          
          for (zu in 2:rnd.len)
          {
            for(jf in 1:n[zu])
              P_akt[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                    (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.old[[zu]]
          }
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
        }
      }
      InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
      if(class(InvFisher2)=="try-error")
        InvFisher2<-try(solve(F_gross),silent=T)
      if(class(InvFisher2)=="try-error")
      {
        half.index<-half.index+1  
      }else{
        solve.test2<-TRUE 
      }}else{
        solve.test2<-TRUE 
      }
    }
    
      betaact<-Delta[l,active]
      
      if (control$method=="EM")
      {        
        ############################# Q update ################
        if(rnd.len==1)
        {
          Q1<-InvFisher2[(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s),(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s)]+Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+s)]%*%t(Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+s)])
          for (i in 2:n)
          Q1<-Q1+InvFisher2[(lin_akt+dim.smooth+(i-1)*s+1):(lin_akt+dim.smooth+i*s),(lin_akt+dim.smooth+(i-1)*s+1):(lin_akt+dim.smooth+i*s)]+Delta[l,(lin+dim.smooth+(i-1)*s+1):(lin+dim.smooth+i*s)]%*%t(Delta[l,(lin+dim.smooth+(i-1)*s+1):(lin+dim.smooth+i*s)])
          Q1<-1/n*Q1
        }else{
          Q1<-matrix(0,sum(s),sum(s))
          Q1[1:s[1],1:s[1]]<-InvFisher2[(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s[1]),(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s[1])]+Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+s[1])]%*%t(Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+s[1])])
          for (i in 2:n[1])
            Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher2[(lin_akt+dim.smooth+(i-1)*s[1]+1):(lin_akt+dim.smooth+i*s[1]),(lin_akt+dim.smooth+(i-1)*s[1]+1):(lin_akt+dim.smooth+i*s[1])]+Delta[l,(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])]%*%t(Delta[l,(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])])
          Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]
          
          for (zu in 2:rnd.len)
          {
            Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher2[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[l,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[l,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
            for (i in 2:n[zu])
              Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher2[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[l,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[l,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
            Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-1/n[zu]*Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
          }
        }  
      }else{
        Eta_tilde<-Eta+(y-Mu)*1/D
        
        Betadach<-Delta[l,1:(lin+dim.smooth)]
        
        aktuell_vec<-!is.element(Delta[l,1:(lin)],0)
        X_aktuell<-Z_fastalles[,aktuell_vec]
        
        if(rnd.len==1)
        {
          
          if(s==1)
          {
            if(Q1<1e-14)
              low<-0
            
            optim.obj<-try(nlminb(sqrt(Q1),likelihood_nlminb,D=D,Sigma=Sigma,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp))
            if(class(optim.obj)=="try-error")
            optim.obj<-bobyqa(sqrt(Q1),likelihood_nlminb,D=D,Sigma=Sigma,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp)

            Q1<-as.matrix(optim.obj$par)^2
          }else{
            Q1_vec<-c(diag(Q1),Q1[lower.tri(Q1)])
            optim.obj<-try(bobyqa(Q1_vec,likelihood,D=D,Sigma=Sigma,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp))
            
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
            optim.obj<-try(bobyqa(sqrt(Q1_vec),likelihood_diag,D=D,Sigma=Sigma,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
            Q1<-diag(optim.obj$par)^2
          }else{
            Q1_vec<-c(diag(Q1)[1:s[1]],Q1[1:s[1],1:s[1]][lower.tri(Q1[1:s[1],1:s[1]])])
            
            for (zu in 2:rnd.len)
              Q1_vec<-c(Q1_vec,c(diag(Q1)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])],Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]))
            
            optim.obj<-try(bobyqa(Q1_vec,likelihood_block,D=D,Sigma=Sigma,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
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
        FinalHat<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))#E-Uu
        phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(FinalHat)))
        Sigma<-Sigma*phi
      }

    Q[[l+1]]<-Q1
    
    score_vec2<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)#-t(t(P.smooth)*Delta[l,])
    lambda.max<-max(abs(score_vec2[(q+1):lin]))
    
    score_old2<-score_vec2
    
    if (BLOCK)
    {
      grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda,block=block)
    }else{
      grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda)
    }
    score_vec2<-c(score_vec2[1:q],grad.1,score_vec2[(lin+1):(lin+dim.smooth+n%*%s)])
    
    if(rnd.len==1)
    {
      if(s==1)
      {
        P1<-c(rep(0,lin),penal.vec,rep((Q1^(-1)),n*s))
      }else{
        Q_inv<-chol2inv(chol(Q1))
        P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
        diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
        for(j in 1:n)
          P1[(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s),(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s)]<-Q_inv
      }
    }else{
      if(all(s==1))
      {
        P1<-c(rep(0,lin),penal.vec,rep(diag(Q1)^(-1),n))
      }else{
        Q_inv<-list()
        Q_inv[[1]]<-chol2inv(chol(Q1[1:s[1],1:s[1]]))
        P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
        diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
        for(jf in 1:n[1])
          P1[(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1]),(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1])]<-Q_inv[[1]]
        
        for (zu in 2:rnd.len)
        {
          Q_inv[[zu]]<-chol2inv(chol(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
          for(jf in 1:n[zu])
            P1[(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
               (lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
        }
      }
    }
    
    if(all(s==1))
    {
      F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+diag(P1)
    }else{
      F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
    }
    
    crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin])
    t_edge<-crit.obj$min.rate
    
    grad.2<-t(score_vec2)%*%F_gross%*%score_vec2
    
    t_opt<-l2norm(score_vec2)$length/grad.2
    tryNR<- (t_opt<t_edge) && !(all(active_old==active)  && !NRstep)
    
    if(tryNR) 
    {
      lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))
      if(rnd.len==1)
      {
        if(s==1)
        {
          P_akt<-c(rep(0,lin_akt),penal.vec,rep((Q1^(-1)),n*s))
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
        }else{
          P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
          diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
          for(jf in 1:n)
            P_akt[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),
                  (lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
        }
      }else{
        if(all(s==1))
        {
          P_akt<-c(rep(0,lin_akt),penal.vec,rep(diag(Q1)^(-1),n))
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
        }else{
          P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
          diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
          for(jf in 1:n[1])
            P_akt[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),
                  (lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv[[1]]
          
          for (zu in 2:rnd.len)
          {
            for(jf in 1:n[zu])
              P_akt[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                    (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
          }
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
        }
      }
      InvFisher3<-try(chol2inv(chol(F_gross)),silent=T)
      if(class(InvFisher3)=="try-error")
        InvFisher3<-try(solve(F_gross),silent=T)
      if(class(InvFisher3)=="try-error")
      {
        half.index<-half.index+1  
      }else{
        
                
        solve.test<-TRUE 
        Delta.test<-Delta[l,active]+nue*InvFisher3%*%score_old2[active]
        if(lin_akt>q)
        {
          vorz<-all(sign(Delta.test[(q+1):lin_akt])==sign(betaact[(q+1):lin_akt]))  
        }else{
          vorz<-T  
        }
      }
    }else{
      solve.test<-TRUE  
    }
    
}  

   if(tryNR) 
     InvFisher<-InvFisher3

    score_old<-score_old2
    score_vec<-score_vec2
    Q1.very.old<-Q1.old
    Q_inv.very.old<-Q_inv.old
    Q1.old<-Q1
    Q_inv.old<-Q_inv    

    Eta.ma[l+1,]<-Eta


    finish<-(sqrt(sum((Eta.ma[l,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l,])^2))<control$epsilon)
    finish2<-(sqrt(sum((Eta.ma[l-1,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l-1,])^2))<control$epsilon)
      if(finish ||  finish2) #|| (all(grad.1 == 0) ))
       break
    Eta.old<-Eta
  }}
 
  FinalHat.df<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))
  df<-sum(diag(FinalHat.df))

  ######## Final calculation
  conv.step<-l
  phi.med<-phi
  
  if(conv.step==control$steps)
  {
    cat("Warning:\n")
    cat("Algorithm did not converge!\n")
  }
  
  Delta_neu<-Delta[l,]
  Eta_opt<-Z_alles%*%Delta_neu
  Mu_opt<-as.vector(family$linkinv(Eta_opt))
  Sigma_opt<-as.vector(family$variance(Mu_opt))    
  D_opt<-as.vector(family$mu.eta(Eta_opt))
  Qfinal<-Q[[l+1]]
  
  aaa<-!is.element(Delta_neu[1:(lin)],0)

if(rnd.len==1)
{
  if(s==1)
  {
    P1.ran<-rep((Qfinal^(-1)),n*s)
    P1.ran<-diag(P1.ran)
  }else{
    P1.ran<-matrix(0,n*s,n*s)
    for(jf in 1:n)
      P1.ran[((jf-1)*s+1):(jf*s),((jf-1)*s+1):(jf*s)]<-chol2inv(chol(Qfinal))
  }
}else{
  if(all(s==1))
  {
    P1.ran<-rep(diag(Qfinal)^(-1),n)
    P1.ran<-diag(P1.ran)
  }else{
    P1.ran<-matrix(0,n%*%s,n%*%s)
    inv.act<-chol2inv(chol(Qfinal[1:s[1],1:s[1]]))
    for(jf in 1:n[1])
      P1.ran[( (jf-1)*s[1]+1):( jf*s[1]),( (jf-1)*s[1]+1):( jf*s[1])]<-inv.act
    
    for (zu in 2:rnd.len)
    {
      inv.act<-chol2inv(chol(Qfinal[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
      for(jf in 1:n[zu])
        P1.ran[( n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):( n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
               ( n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):( n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-inv.act
    }
  }  
}  

ranef.logLik<--0.5*t(Delta_neu[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)])%*%P1.ran%*%Delta_neu[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]


if(final.re)
{    
############ final re-estimation
  
  if(s==1)
  {  
    Q.max<-max(sqrt(unlist(Q)))
    Q.min<-min(sqrt(unlist(Q)))
  }else{
    Q.max<-max(Qfinal)+1
    Q.min<-min(Qfinal)-1e-10
  }
  

  if(rnd.len==1)
  {
    glmm_fin<-try(glmm_final_smooth(y,Z_fastalles[,aaa],Phi,W,k,penal.vec,q_start=Qfinal,
                      Delta_start=Delta_neu[c(aaa,rep(T,dim.smooth+n%*%s))],
                      s,steps=control$maxIter,family=family,method=control$method.final,overdispersion=control$overdispersion,
                      phi=phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final,
                      Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac))#,silent = TRUE)
    if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
    {  
      glmm_fin2<-try(glmm_final_smooth(y,Z_fastalles[,aaa],Phi,W,k,penal.vec,q_start=q_start,
                      Delta_start=Delta_start[c(aaa,rep(T,dim.smooth+n%*%s))],
                      s,steps=control$maxIter,family=family,method=control$method.final,overdispersion=control$overdispersion,
                      phi=control$phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final,
                      Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac))
        if(class(glmm_fin2)!="try-error")
        {    
          if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
          glmm_fin<-glmm_fin2 
        }
    }
  }else{
    glmm_fin<-try(glmm_final_multi_random_smooth(y,Z_fastalles[,aaa],Phi,W,k,penal.vec,q_start=Qfinal,
                    Delta_start=Delta_neu[c(aaa,rep(T,dim.smooth+n%*%s))],
                    s,n,steps=control$maxIter,family=family,method=control$method.final,overdispersion=control$overdispersion,
                    phi=phi,rnd.len=rnd.len,print.iter.final=control$print.iter.final,eps.final=control$eps.final,
                    Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac))#,silent = TRUE)
    if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
    {  
      glmm_fin2<-try(glmm_final_multi_random_smooth(y,Z_fastalles[,aaa],Phi,W,k,penal.vec,q_start=q_start,
                    Delta_start=Delta_start[c(aaa,rep(T,dim.smooth+n%*%s))],
                    s,n,steps=control$maxIter,family=family,method=control$method.final,overdispersion=control$overdispersion,
                    phi=control$phi,rnd.len=rnd.len,print.iter.final=control$print.iter.final,
                    eps.final=control$eps.final,Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac))
        if(class(glmm_fin2)!="try-error")
        {    
          if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
          glmm_fin<-glmm_fin2 
        }
    }
  }

  if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
  {
    cat("Warning:\n")
    cat("Final Fisher scoring reestimation did not converge!\n")
  }
  
#######
}else{
  glmm_fin<-NA  
  class(glmm_fin)<-"try-error"
}  
  
  Delta_neu2<-Delta_neu
  Standard_errors<-rep(NA,length(Delta_neu))
  

  if(class(glmm_fin)!="try-error")
  {
    EDF.matrix<-glmm_fin$EDF.matrix
    complexity.smooth<-sum(diag(EDF.matrix)[c(T,rep(F,sum(aaa)-1),rep(T,dim.smooth),rep(F,n%*%s))])  
    Delta_neu2[c(aaa,rep(T,dim.smooth+n%*%s))]<-glmm_fin$Delta
    Standard_errors[c(aaa,rep(T,dim.smooth+n%*%s))]<-glmm_fin$Standard_errors
    Qfinal<-glmm_fin$Q
    phi<-glmm_fin$phi
    complexity<-glmm_fin$complexity
  }else{
  glmm_fin<-list()
  glmm_fin$ranef.logLik<-ranef.logLik
  complexity<-df
  
  lin_akt<-q+sum(!is.element(Delta_neu[(q+1):lin],0))
  if(rnd.len==1)
  {
    if(s==1)
    {
      P1<-c(rep(0,lin_akt),penal.vec,rep((Q1^(-1)),n*s))
      P1a<-c(rep(0,lin_akt+dim.smooth),rep((Q1^(-1)),n*s))
      F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P1)
    }else{
      P1<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
      P1a<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
      diag(P1)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
      for(jf in 1:n)
      {
        P1[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),
           (lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv
        P1a[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),
            (lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv
      }  
      F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P1
    }
  }else{
    if(all(s==1))
    {
      P1<-c(rep(0,lin_akt),penal.vec,rep(diag(Q1)^(-1),n))
      P1a<-c(rep(0,lin_akt+dim.smooth),rep(diag(Q1)^(-1),n))
      F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P1)
    }else{
      P1<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
      P1a<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
      diag(P1)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
      for(jf in 1:n[1])
      {
        P1[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),
           (lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv[[1]]
        P1a[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),
            (lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv[[1]]
      }
      for (zu in 2:rnd.len)
      {
        for(jf in 1:n[zu])
        {
          P1[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
             (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
          P1a[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
              (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
        }}
      F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P1
    }
  }
  InvFisher<-try(chol2inv(chol(F_gross)),silent=T)
  if(class(InvFisher)=="try-error")
    InvFisher<-try(solve(F_gross),silent=T)
  if(class(InvFisher)=="try-error")
  {
    warning("No EDF's for smooth functions available, as Fisher matrix not invertible!")
    complexity.smooth<-dim.smooth
    }else{  
    ###### EDF of spline; compare Wood's Book on page 167
    EDF.matrix<-InvFisher%*%(t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P1a)
    complexity.smooth<-sum(diag(EDF.matrix)[c(T,rep(F,sum(aaa)-1),rep(T,dim.smooth),rep(F,n%*%s))])
    }
  }  
  
  if(!(complexity.smooth>=1 && complexity.smooth<=dim.smooth))
  complexity.smooth<-dim.smooth
     
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
  
  names(Delta_neu)[1:dim(X)[2]]<-colnames(X)
  names(Standard_errors)[1:dim(X)[2]]<-colnames(X)

  if(lin>1)
  {
    names(Delta_neu)[(dim(X)[2]+1):lin]<-colnames(U)
    names(Standard_errors)[(dim(X)[2]+1):lin]<-colnames(U)
  }

  names(Delta_neu)[(lin+1):(lin+dim.smooth)]<-colnames(Phi)
  names(Standard_errors)[(lin+1):(lin+dim.smooth)]<-colnames(Phi)
  names(Delta_neu)[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]<-colnames(W)
  names(Standard_errors)[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]<-colnames(W)
  colnames(Delta)<-c(old.names,colnames(Phi),colnames(W))
  
  aic<-NaN
  bic<-NaN
  
if (is.element(family$family,c("gaussian", "binomial", "poisson"))) 
{
    loglik<-logLik.glmmLasso(y=y,mu=Mu_opt,beta=Delta_neu[(q+1):(lin)],ranef.logLik=glmm_fin$ranef.logLik,
                             lambda=lambda,family=family,penal=FALSE)

if(control$complexity!="hat.matrix")  
{  
    if(rnd.len==1)
    {
      complexity<-0.5*(s*(s+1))
    }else{
      complexity<-0.5*(s[1]*(s[1]+1))
      for(zu in 2:rnd.len)
        complexity<-complexity+0.5*(s[zu]*(s[zu]+1))
    }
    complexity<-complexity+sum(Delta_neu[1:(lin)]!=0)+complexity.smooth
}
    aic<--2*loglik+2*complexity
    bic<--2*loglik+log(N)*complexity
}else{
  warning("For the specified family (so far) no AIC and BIC are available!")  
}
  
  ret.obj=list()
  ret.obj$aic<-aic
  ret.obj$bic<-bic
  ret.obj$Deltamatrix<-Delta
  ret.obj$smooth<-Delta_neu[(lin+1):(lin+dim.smooth)]
  ret.obj$ranef<-Delta_neu[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]
  ret.obj$coefficients<-Delta_neu[1:(lin)]
  ret.obj$fixerror<-Standard_errors[1:(lin)]
  ret.obj$ranerror<-Standard_errors[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]
  ret.obj$Q_long<-Q
  ret.obj$Q<-Qfinal
  ret.obj$y_hat<-Mu_opt
  ret.obj$phi<-phi
  ret.obj$family<-family
  ret.obj$fix<-fix.old
  ret.obj$newrndfrml<-newrndfrml
  ret.obj$subject<-names(rnd)
  ret.obj$data<-data
  ret.obj$rnd.len<-rnd.len
  ret.obj$B<-B
  ret.obj$nbasis<-nbasis
  ret.obj$spline.degree<-spline.degree
  ret.obj$diff.ord<-diff.ord
  ret.obj$knots.no<-knots.no
  ret.obj$conv.step<-conv.step
  ret.obj$phi.med<-phi.med
  ret.obj$complexity.smooth<-complexity.smooth
  ret.obj$y <- y
  ret.obj$df<-df
  ret.obj$loglik<-loglik
  ret.obj$lambda.max<-lambda.max
  return(ret.obj)
}  

#######################################################################  
######################## no switch to Newton Raphson ###############
#######################################################################  
}else{
  #######################################################################  
  ###########################  1. No Smooth #############################  
  #######################################################################  
  if(is.null(control$smooth))
  {  
    # browser()
    
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
      }else{
        Q_start<-q_start
      }
    }else{
      lin0<-sum(beta_null!=0)
      if(all(s==1))
      {
        Q_start<-diag(diag(q_start),sum(s))
      }else{
        Q_start<-matrix(0,sum(s),sum(s))
        Q_start[1:s[1],1:s[1]]<-q_start[1:s[1],1:s[1]]
        
        for (zu in 2:rnd.len)
        {
          Q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
         }
      }
    }
    
    U<-X[,!is.na(index.new),drop=FALSE]
    X<-X[,is.na(index.new),drop=FALSE]
    
    final.names<-c(colnames(X),colnames(U))
    
    q<-dim(X)[2]
    
    Z_alles<-cbind(X,U,W)
    ########################################################## some definitions ################################################
    Delta<-matrix(0,control$steps,(lin+n%*%s))
    Delta[1,1:lin]<-beta_null[final.names]
    Delta[1,(lin+1):(lin+n%*%s)]<-t(ranef_null)

    Eta.ma<-matrix(0,control$steps+1,N)
    Eta.ma[1,]<-Eta_start
    
    control$epsilon<-control$epsilon*sqrt(dim(Delta)[2])
    
    Delta_start<-Delta[1,]
    
    active_old<-!is.element(Delta[1,],0)

    Q_inv<-NULL
    Q_inv.old.temp<-NULL
    Q_inv.start<-NULL
    
    ## start value for overdispersion
    if(control$overdispersion)
    {
      if(!is.null(control$phi_start))
      {  
        phi<-control$phi_start
      }else{
        if(is.null(control$Q.phi.start))
        {
          Q.phi.start<-q_start 
        }else{
          Q.phi.start<-control$Q.phi.start
        }
        active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,n%*%s))
        Z_aktuell<-Z_alles[,active]
        lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
        
        if(rnd.len==1)
        {
          if(s==1)
          {
            Q.phi.start<-diag(Q.phi.start,s)
            P_akt<-c(rep(0,lin_akt),rep((Q.phi.start^(-1)),n*s))
            F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
          }else{
            Q_inv.start<-chol2inv(chol(Q.phi.start))
            P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
            for(jf in 1:n)
              P_akt[(lin_akt+(jf-1)*s+1):(lin_akt+jf*s),(lin_akt+(jf-1)*s+1):(lin_akt+jf*s)]<-Q_inv.start
            F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
          }
        }else{
          if(all(s==1))
          {
            Q.phi.start<-diag(diag(Q.phi.start),sum(s))
            P_akt<-c(rep(0,lin_akt),rep(diag(Q.phi.start)^(-1),n))
            F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
          }else{
            Q_inv.start<-list()
            Q_inv.start[[1]]<-chol2inv(chol(Q.phi.start[1:s[1],1:s[1]]))
            
            for (zu in 2:rnd.len)
            {
              Q_inv.start[[zu]]<-chol2inv(chol(Q.phi.start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
            }
            
            P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
            for(jf in 1:n[1])
              P_akt[(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1]),(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1])]<-Q_inv.start[[1]]
            
            for (zu in 2:rnd.len)
            {
              for(jf in 1:n[zu])
                P_akt[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                      (lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.start[[zu]]
            }
            F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
          }
        }
        InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
        if(class(InvFisher2)=="try-error")
          InvFisher2<-try(solve(F_gross),silent=T)
        if(class(InvFisher2)=="try-error")
          stop("Fisher matrix not invertible")  
        
        if(all(Mu==0) &family$family=="gaussian")
        {
          phi<-1  
        }else{
          Hat<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))#E-Uu
        phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(Hat)))
        }
        Sigma<-Sigma*phi
      }  
    }else{
      phi<-1
    }
    
    Q<-list()
    Q[[1]]<-Q_start
    
    l=1
        
    score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)
    

    if (BLOCK)
    {
      grad.1<-gradient.lasso.block(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
    }else{
      grad.1<-gradient.lasso(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
    }
    
    score_vec<-c(score_vec[1:q],grad.1,score_vec[(lin+1):(lin+n%*%s)])
    
    if(rnd.len==1)
    {
      if(s==1)
      {
        P1<-c(rep(0,lin),rep((Q_start^(-1)),n*s))
      }else{
        Q_inv.start<-chol2inv(chol(Q_start))
        P1<-matrix(0,lin+n%*%s,lin+n%*%s)
        for(j in 1:n)
          P1[(lin+(j-1)*s+1):(lin+j*s),(lin+(j-1)*s+1):(lin+j*s)]<-Q_inv.start
      }
    }else{
      if(all(s==1))
      {
        P1<-c(rep(0,lin),rep(diag(Q_start)^(-1),n))
      }else{
        Q_inv.start<-list()
        Q_inv.start[[1]]<-chol2inv(chol(Q_start[1:s[1],1:s[1]]))
        P1<-matrix(0,lin+n%*%s,lin+n%*%s)
        for(jf in 1:n[1])
          P1[(lin+(jf-1)*s[1]+1):(lin+jf*s[1]),(lin+(jf-1)*s[1]+1):(lin+jf*s[1])]<-Q_inv.start[[1]]
        
        for (zu in 2:rnd.len)
        {
          Q_inv.start[[zu]]<-chol2inv(chol(Q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
          for(jf in 1:n[zu])
            P1[(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
               (lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.start[[zu]]
        }
      }
    }

    if(all(s==1))
    {
      F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+diag(P1)
    }else{
      F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
    }
        
    crit.obj<-t.change(grad=score_vec[(q+1):lin],b=Delta[1,(q+1):lin])
    t_edge<-crit.obj$min.rate
    
    grad.2<-t(score_vec)%*%F_gross%*%score_vec
    
    t_opt<-l2norm(score_vec)$length/grad.2
    nue<-control$nue
    
    
    half.index<-0

      solve.test2<-FALSE  
      while(!solve.test2)
      {  
        
        if(half.index>50)
        {
          stop("Fisher matrix not invertible")
        }
        
        Delta[1,]<-Delta_start+min(t_opt,t_edge)*nue*(0.5^half.index)*score_vec
        if(t_opt>t_edge & half.index==0)
          Delta[1,crit.obj$whichmin+q]<-0  
        Eta<-Z_alles%*%Delta[1,]
        Mu<-as.vector(family$linkinv(Eta))
        Sigma<-as.vector(family$variance(Mu))
        D<-as.vector(family$mu.eta(Eta))

        active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,n%*%s))
        Z_aktuell<-Z_alles[,active]
        lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
       
        
        
        if (control$method=="EM" || control$overdispersion)
        {  
          if(rnd.len==1)
          {
            if(s==1)
            {
              P_akt<-c(rep(0,lin_akt),rep((Q_start^(-1)),n*s))
              F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
            }else{
              P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
              for(jf in 1:n)
                P_akt[(lin_akt+(jf-1)*s+1):(lin_akt+jf*s),(lin_akt+(jf-1)*s+1):(lin_akt+jf*s)]<-Q_inv.start
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
                P_akt[(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1]),(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1])]<-Q_inv.start[[1]]
              
              for (zu in 2:rnd.len)
              {
                for(jf in 1:n[zu])
                  P_akt[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                        (lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.start[[zu]]
              }
              F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
            }
          }
          InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
          if(class(InvFisher2)=="try-error")
            InvFisher2<-try(solve(F_gross),silent=T)
          if(class(InvFisher2)=="try-error")
          {
            #stop("Fisher matrix not invertible")  
            half.index<-half.index+1  
          }else{
            solve.test2<-TRUE 
          }}else{
            solve.test2<-TRUE
          }
      }
      
      betaact<-Delta[1,active]
          
    if (control$method=="EM")
      {   
        ############################# Q update ################
        if(rnd.len==1)
        {
          Q1<-InvFisher2[(lin_akt+1):(lin_akt+s),(lin_akt+1):(lin_akt+s)]+Delta[1,(lin+1):(lin+s)]%*%t(Delta[1,(lin+1):(lin+s)])
          for (i in 2:n)
            Q1<-Q1+InvFisher2[(lin_akt+(i-1)*s+1):(lin_akt+i*s),(lin_akt+(i-1)*s+1):(lin_akt+i*s)]+Delta[1,(lin+(i-1)*s+1):(lin+i*s)]%*%t(Delta[1,(lin+(i-1)*s+1):(lin+i*s)])
          Q1<-1/n*Q1
        }else{
          Q1<-matrix(0,sum(s),sum(s))
          Q1[1:s[1],1:s[1]]<-InvFisher2[(lin_akt+1):(lin_akt+s[1]),(lin_akt+1):(lin_akt+s[1])]+Delta[1,(lin+1):(lin+s[1])]%*%t(Delta[1,(lin+1):(lin+s[1])])
          for (i in 2:n[1])
            Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher2[(lin_akt+(i-1)*s[1]+1):(lin_akt+i*s[1]),(lin_akt+(i-1)*s[1]+1):(lin_akt+i*s[1])]+Delta[1,(lin+(i-1)*s[1]+1):(lin+i*s[1])]%*%t(Delta[1,(lin+(i-1)*s[1]+1):(lin+i*s[1])])
          Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]
          
          for (zu in 2:rnd.len)
          {
            Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher2[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
            for (i in 2:n[zu])
              Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher2[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[1,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
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
            optim.obj<-nlminb(sqrt(Q_start),likelihood_nlminb,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp)
            Q1<-as.matrix(optim.obj$par)^2
          }else{
            q_start_vec<-c(diag(q_start),q_start[lower.tri(q_start)])
            up1<-min(20,50*max(q_start_vec))#[(s+1):(s*(s+1)*0.5)]))
            upp<-rep(up1,length(q_start_vec))
            low<-c(rep(0,s),rep(-up1,0.5*(s^2-s)))
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
            
            for (zu in 2:rnd.len)
            {
              q_start_vec<-c(q_start_vec,c(diag(q_start)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])],q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]))
              up1<-min(20,50*max(q_start_vec))
              low<-c(low,c(rep(0,s[zu]),rep(-up1,0.5*(s[zu]^2-s[zu]))))
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
        FinalHat<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))#E-Uu
        phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(FinalHat)))
        Sigma<-Sigma*phi
      }
      
     Eta.old<-Eta
            
    score_vec2<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)

      if (BLOCK)
      {
        grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
      }else{
        grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
      }
      score_vec2<-c(score_vec2[1:q],grad.1,score_vec2[(lin+1):(lin+n%*%s)])
      
    if(rnd.len==1)
    {
      if(s==1)
      {
        P1<-c(rep(0,lin),rep((Q1^(-1)),n*s))
      }else{
        Q_inv<-solve(Q1)
        P1<-matrix(0,lin+n%*%s,lin+n%*%s)
        for(j in 1:n)
          P1[(lin+(j-1)*s+1):(lin+j*s),(lin+(j-1)*s+1):(lin+j*s)]<-Q_inv
      }
    }else{
      if(all(s==1))
      {
        P1<-c(rep(0,lin),rep(diag(Q1)^(-1),n))
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
      }
    }
    
    if(all(s==1))
      {
        F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+diag(P1)
      }else{
        F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
      }
      
    crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin])
    t_edge<-crit.obj$min.rate
    
      grad.2<-t(score_vec2)%*%F_gross%*%score_vec2
      
      t_opt<-l2norm(score_vec2)$length/grad.2
    
    Eta.ma[2,]<-Eta

    score_vec<-score_vec2
    Q1.old<-Q1
    Q_inv.old<-Q_inv    
    
    Q1.very.old<-Q_start
    Q_inv.very.old<-Q_inv.start
    
    ###############################################################################################################################################
    ################################################################### Main Iteration ###################################################################
    if(control$steps!=1)
    {
      for (l in 2:control$steps)
      {
         if(control$print.iter)
          print(paste("Iteration ", l,sep=""))
         
        half.index<-0

           solve.test2<-FALSE  
          while(!solve.test2)
          {  
            
            if(half.index>50)
            {
              half.index<-Inf;Q1.old<-Q1.very.old;Q_inv.old<-Q_inv.very.old
            }
            
              Delta[l,]<-Delta[l-1,]+min(t_opt,t_edge)*nue*(0.5^half.index)*score_vec
            if(t_opt>t_edge & half.index==0)
              Delta[l,crit.obj$whichmin+q]<-0  
          
            Eta<-Z_alles%*%Delta[l,]
            Mu<-as.vector(family$linkinv(Eta))
            Sigma<-as.vector(family$variance(Mu))
            D<-as.vector(family$mu.eta(Eta))
            
            active_old<-active
            active<-c(rep(T,q),!is.element(Delta[l,(q+1):lin],0),rep(T,n%*%s))
            Z_aktuell<-Z_alles[,active]
            lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))
            
            if (control$method=="EM" || control$overdispersion)
            {  
              if(rnd.len==1)
              {
                if(s==1)
                {
                  P_akt<-c(rep(0,lin_akt),rep((Q1.old^(-1)),n*s))
                  F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
                }else{
                  P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
                  for(jf in 1:n)
                    P_akt[(lin_akt+(jf-1)*s+1):(lin_akt+jf*s),(lin_akt+(jf-1)*s+1):(lin_akt+jf*s)]<-Q_inv.old
                  F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
                }
              }else{
                if(all(s==1))
                {
                  P_akt<-c(rep(0,lin_akt),rep(diag(Q1.old)^(-1),n))
                  F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
                }else{
                  P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
                  for(jf in 1:n[1])
                    P_akt[(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1]),(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1])]<-Q_inv.old[[1]]
                  
                  for (zu in 2:rnd.len)
                  {
                    for(jf in 1:n[zu])
                      P_akt[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                            (lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.old[[zu]]
                  }
                  F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
                }
              }
              InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
              if(class(InvFisher2)=="try-error")
                InvFisher2<-try(solve(F_gross),silent=T)
              if(class(InvFisher2)=="try-error")
              {
                half.index<-half.index+1  
              }else{
                solve.test2<-TRUE 
              }}else{
                solve.test2<-TRUE 
              }
          }
          
          betaact<-Delta[l,active]
          
          if (control$method=="EM")
          {        
            ############################# Q update ################
            if(rnd.len==1)
            {
              Q1<-InvFisher2[(lin_akt+1):(lin_akt+s),(lin_akt+1):(lin_akt+s)]+Delta[l,(lin+1):(lin+s)]%*%t(Delta[l,(lin+1):(lin+s)])
              for (i in 2:n)
                Q1<-Q1+InvFisher2[(lin_akt+(i-1)*s+1):(lin_akt+i*s),(lin_akt+(i-1)*s+1):(lin_akt+i*s)]+Delta[l,(lin+(i-1)*s+1):(lin+i*s)]%*%t(Delta[l,(lin+(i-1)*s+1):(lin+i*s)])
              Q1<-1/n*Q1
            }else{
              Q1<-matrix(0,sum(s),sum(s))
              Q1[1:s[1],1:s[1]]<-InvFisher2[(lin_akt+1):(lin_akt+s[1]),(lin_akt+1):(lin_akt+s[1])]+Delta[l,(lin+1):(lin+s[1])]%*%t(Delta[l,(lin+1):(lin+s[1])])
              for (i in 2:n[1])
                Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher2[(lin_akt+(i-1)*s[1]+1):(lin_akt+i*s[1]),(lin_akt+(i-1)*s[1]+1):(lin_akt+i*s[1])]+Delta[l,(lin+(i-1)*s[1]+1):(lin+i*s[1])]%*%t(Delta[l,(lin+(i-1)*s[1]+1):(lin+i*s[1])])
              Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]
              
              for (zu in 2:rnd.len)
              {
                Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher2[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
                for (i in 2:n[zu])
                  Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher2[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[l,(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
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
                
                optim.obj<-nlminb(sqrt(Q1),likelihood_nlminb,D=D,Sigma=Sigma,X=Z_fastalles,X_aktuell=X_aktuell,Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp)
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
            FinalHat<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))#E-Uu
            phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(FinalHat)))
            Sigma<-Sigma*phi
          }
          
          Q[[l+1]]<-Q1
          
        score_vec2<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)
        
          if (BLOCK)
          {
            grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda,block=block)
          }else{
            grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda)
          }
        
          score_vec2<-c(score_vec2[1:q],grad.1,score_vec2[(lin+1):(lin+n%*%s)])
        
        if(rnd.len==1)
        {
          if(s==1)
          {
            P1<-c(rep(0,lin),rep((Q1^(-1)),n*s))
          }else{
            Q_inv<-solve(Q1)
            P1<-matrix(0,lin+n%*%s,lin+n%*%s)
            for(j in 1:n)
              P1[(lin+(j-1)*s+1):(lin+j*s),(lin+(j-1)*s+1):(lin+j*s)]<-Q_inv
          }
        }else{
          if(all(s==1))
          {
            P1<-c(rep(0,lin),rep(diag(Q1)^(-1),n))
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
          }
        }
        
          if(all(s==1))
          {
            F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+diag(P1)
          }else{
            F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
          }
          
        crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin])
        t_edge<-crit.obj$min.rate
        
          grad.2<-t(score_vec2)%*%F_gross%*%score_vec2
          
          t_opt<-l2norm(score_vec2)$length/grad.2
        
        score_vec<-score_vec2
        Q1.very.old<-Q1.old
        Q_inv.very.old<-Q_inv.old
        Q1.old<-Q1
        Q_inv.old<-Q_inv    
        
        Eta.ma[l+1,]<-Eta
        
        finish<-(sqrt(sum((Eta.ma[l,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l,])^2))<control$epsilon)
        finish2<-(sqrt(sum((Eta.ma[l-1,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l-1,])^2))<control$epsilon)
        if(finish ||  finish2) #|| (all(grad.1 == 0) ))
          break
        Eta.old<-Eta
      }}
    

    if(control$method=="REML")
    {
      if(rnd.len==1)
      {
        if(s==1)
        {
          P_akt<-c(rep(0,lin_akt),rep((Q1.old^(-1)),n*s))
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
        }else{
          P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
          for(jf in 1:n)
            P_akt[(lin_akt+(jf-1)*s+1):(lin_akt+jf*s),(lin_akt+(jf-1)*s+1):(lin_akt+jf*s)]<-Q_inv.old
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
        }
      }else{
        if(all(s==1))
        {
          P_akt<-c(rep(0,lin_akt),rep(diag(Q1.old)^(-1),n))
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
        }else{
          P_akt<-matrix(0,lin_akt+n%*%s,lin_akt+n%*%s)
          for(jf in 1:n[1])
            P_akt[(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1]),(lin_akt+(jf-1)*s[1]+1):(lin_akt+jf*s[1])]<-Q_inv.old[[1]]
          
          for (zu in 2:rnd.len)
          {
            for(jf in 1:n[zu])
              P_akt[(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                    (lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.old[[zu]]
          }
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
        }
      }
      InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
      if(class(InvFisher2)=="try-error")
        InvFisher2<-try(solve(F_gross),silent=T)
    }

    FinalHat.df<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))
    df<-sum(diag(FinalHat.df))

    conv.step<-l
    phi.med<-phi
    
    if(conv.step==control$steps)
    {
      cat("Warning:\n")
      cat("Algorithm did not converge!\n")
    }
    
    
    Delta_neu<-Delta[l,]
    Eta_opt<-Z_alles%*%Delta_neu
    Mu_opt<-as.vector(family$linkinv(Eta_opt))

    Sigma_opt<-as.vector(family$variance(Mu_opt))    
    D_opt<-as.vector(family$mu.eta(Eta_opt))
    Qfinal<-Q[[l+1]]
    
    
    
    if(rnd.len==1)
    {
      if(s==1)
      {
        P1.ran<-rep((Qfinal^(-1)),n*s)
        P1.ran<-diag(P1.ran)
      }else{
        P1.ran<-matrix(0,n*s,n*s)
        for(jf in 1:n)
          P1.ran[((jf-1)*s+1):(jf*s),((jf-1)*s+1):(jf*s)]<-chol2inv(chol(Qfinal))
      }
    }else{
      if(all(s==1))
      {
        P1.ran<-rep(diag(Qfinal)^(-1),n)
        P1.ran<-diag(P1.ran)
      }else{
        P1.ran<-matrix(0,n%*%s,n%*%s)
        inv.act<-chol2inv(chol(Qfinal[1:s[1],1:s[1]]))
        for(jf in 1:n[1])
          P1.ran[( (jf-1)*s[1]+1):( jf*s[1]),( (jf-1)*s[1]+1):( jf*s[1])]<-inv.act
        
        for (zu in 2:rnd.len)
        {
          inv.act<-chol2inv(chol(Qfinal[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
          for(jf in 1:n[zu])
            P1.ran[( n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):( n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                   ( n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):( n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-inv.act
        }
      }
    }
    ranef.logLik<--0.5*t(Delta_neu[(lin+1):(lin+n%*%s)])%*%P1.ran%*%Delta_neu[(lin+1):(lin+n%*%s)]
    
    
    if(final.re)
    {    
      ############ final re-estimation
      aaa<-!is.element(Delta_neu[1:(lin)],0)
      
      if(rnd.len==1 && s==1)
      {  
        Q.max<-max(sqrt(unlist(Q)))
        Q.min<-min(sqrt(unlist(Q)))
      }else{
        Q.max<-max(Qfinal)+1
        Q.min<-min(Qfinal)-1e-10
      }
      
      if(rnd.len==1)
      {
        glmm_fin<-try(glmm_final(y,Z_fastalles[,aaa],W,k,n,q_start=Qfinal,
                                 Delta_start=Delta_neu[c(aaa,rep(T,n%*%s))],s,steps=control$maxIter,
                                 family=family,method=control$method.final,overdispersion=control$overdispersion,
                                 phi=phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final,
                                 Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac))#,silent = TRUE)
        if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
        {  
          glmm_fin2<-try(glmm_final(y,Z_fastalles[,aaa],W,k,n,q_start=q_start,
                                    Delta_start=Delta_start[c(aaa,rep(T,n%*%s))],s,steps=control$maxIter,
                                    family=family,method=control$method.final,overdispersion=control$overdispersion,
                                    phi=control$phi,print.iter.final=control$print.iter.final,
                                    eps.final=control$eps.final,Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac))
          if(class(glmm_fin2)!="try-error")
          {    
            if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
              glmm_fin<-glmm_fin2 
          }
        }
      }else{
        glmm_fin<-try(glmm_final_multi_random(y,Z_fastalles[,aaa],W,k,q_start=Qfinal,
                                              Delta_start=Delta_neu[c(aaa,rep(T,n%*%s))],s,n,steps=control$maxIter,
                                              family=family,method=control$method.final,overdispersion=control$overdispersion,
                                              phi=phi,rnd.len=rnd.len,print.iter.final=control$print.iter.final,
                                              eps.final=control$eps.final,Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac))#,silent = TRUE)
        if(class(glmm_fin)=="try-error"|| glmm_fin$opt>control$maxIter-10)
        {  
          glmm_fin2<-try(glmm_final_multi_random(y,Z_fastalles[,aaa],W,k,q_start=q_start,
                                                 Delta_start=Delta_start[c(aaa,rep(T,n%*%s))],s,n,steps=control$maxIter,
                                                 family=family,method=control$method.final,overdispersion=control$overdispersion,
                                                 phi=control$phi,rnd.len=rnd.len,print.iter.final=control$print.iter.final,
                                                 eps.final=control$eps.final,Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac))    
          if(class(glmm_fin2)!="try-error")
          {    
            if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
              glmm_fin<-glmm_fin2 
          }
        }
      }
      
      if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
      {
        cat("Warning:\n")
        cat("Final Fisher scoring reestimation did not converge!\n")
      }
      #######
    }else{
      glmm_fin<-NA  
      class(glmm_fin)<-"try-error"
    }  
    Delta_neu2<-Delta_neu
    Standard_errors<-rep(NA,length(Delta_neu))
    
   if(class(glmm_fin)!="try-error")
    {
      Delta_neu2[c(aaa,rep(T,n%*%s))]<-glmm_fin$Delta
      Standard_errors[c(aaa,rep(T,n%*%s))]<-glmm_fin$Standard_errors
      Qfinal<-glmm_fin$Q
      phi<-glmm_fin$phi
      complexity<-glmm_fin$complexity
    }else{
      glmm_fin<-list()
      glmm_fin$ranef.logLik<-ranef.logLik
      complexity<-df
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
      
      if(s[1]==1)
        Qfinal[[1]]<-sqrt(Qfinal[[1]])
      
      for (zu in 2:rnd.len)
      {
        Qfinal[[zu]]<-as.matrix(Qfinal_old[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])
        colnames(Qfinal[[zu]])<-random.labels[[zu]]
        rownames(Qfinal[[zu]])<-random.labels[[zu]]
        
        if(s[zu]==1)
          Qfinal[[zu]]<-sqrt(Qfinal[[zu]])
        
      }
    }
    
    names(Delta_neu)[1:dim(X)[2]]<-colnames(X)
    names(Standard_errors)[1:dim(X)[2]]<-colnames(X)
    
    if(lin>1)
    {
      names(Delta_neu)[(dim(X)[2]+1):lin]<-colnames(U)
      names(Standard_errors)[(dim(X)[2]+1):lin]<-colnames(U)
    }
    
    names(Delta_neu)[(lin+1):(lin+n%*%s)]<-colnames(W)
    names(Standard_errors)[(lin+1):(lin+n%*%s)]<-colnames(W)
    colnames(Delta)<-c(final.names,colnames(W))
    
    aic<-NaN
    bic<-NaN
    
    if (is.element(family$family,c("gaussian", "binomial", "poisson"))) 
    {
      
      loglik<-logLik.glmmLasso(y=y,mu=Mu_opt,beta=Delta_neu[(q+1):(lin)],ranef.logLik=glmm_fin$ranef.logLik,
                               lambda=lambda,family=family,penal=FALSE)
      
if(control$complexity!="hat.matrix")  
{  
        if(rnd.len==1)
      {
        complexity<-0.5*(s*(s+1))
      }else{
        complexity<-0.5*(s[1]*(s[1]+1))
        for(zu in 2:rnd.len)
          complexity<-complexity+0.5*(s[zu]*(s[zu]+1))
      }
      complexity<-complexity+sum(Delta_neu[1:(lin)]!=0)
}      
      aic<--2*loglik+2*complexity
      bic<--2*loglik+log(N)*complexity
    }else{
      warning("For the specified family (so far) no AIC and BIC are available!")  
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
    ret.obj$fix<-fix.old
    ret.obj$conv.step<-conv.step
    ret.obj$newrndfrml<-newrndfrml
    ret.obj$subject<-subject.names
    ret.obj$data<-data
    ret.obj$rnd.len<-rnd.len
    ret.obj$phi.med<-phi.med
    ret.obj$y <- y
    ret.obj$df<-df
    ret.obj$loglik<-loglik
    return(ret.obj)
    ##############################################################  
    ######################## 2. Smooth ###########################  
    ##############################################################  
  }else{
    
   # browser()
    
    smooth<-control$smooth
    
    if(attr(terms(smooth$formula), "intercept")==0)
    {
      variables<-attr(terms(smooth$formula),"term.labels")
      smooth$formula<- "~ +1"
      for (ir in 1:length(variables))
        smooth$formula <- paste(smooth$formula, variables[ir], sep="+")
      smooth$formula <-formula(smooth$formula)
    }
    
    B <- model.matrix(smooth$formula, data)
    B.names<-attr(B,"dimnames")[[2]]
    
    B<-as.matrix(B[,-1])
    attr(B,"dimnames")[[2]]<-B.names[-1]
    
    nbasis<-smooth$nbasis
    diff.ord<-smooth$diff.ord
    spline.degree<-smooth$spline.degree
    knots.no<-nbasis-1
    if(spline.degree<3 && (spline.degree-diff.ord)<2)
      knots.no<-knots.no+1  
    penal<-smooth$penal
    
    if(!(diff.ord<spline.degree))
      stop("Order of differences must be lower than degree of B-spline polynomials!")
    
    m<-dim(B)[2]
    
    Phi<-numeric()
    
    for (r in 1:m)
    {
      Basis<-bs.design(B[,r],diff.ord=diff.ord,spline.degree=spline.degree,knots.no=knots.no)
      Phi_temp<-cbind(Basis$X[,-1],Basis$Z)
      colnames(Phi_temp)<-paste(colnames(B)[r],rep(1:dim(Phi_temp)[2],each=1), sep=".")
      Phi<-cbind(Phi,Phi_temp)
    }
    
    dim.smooth<-dim(Phi)[2]
    
    if(is.null(smooth$start))
      smooth$start<-rep(0,dim.smooth)  
    
    smooth_null<-smooth$start
    
    if(lin>1)
    {
      Eta_start<-X%*%beta_null+Phi%*%smooth_null+W%*%ranef_null
    }else{
      Eta_start<-rep(beta_null,N)+Phi%*%smooth_null+W%*%ranef_null
    }
    
    D<-as.vector(family$mu.eta(Eta_start))
    Mu<-as.vector(family$linkinv(Eta_start))
    Sigma<-as.vector(family$variance(Mu))
     
    if(rnd.len==1)
    {
      if(s==1)
      {
        Q_start<-diag(q_start,s)
      }else{
        Q_start<-q_start
      }
    }else{
      if(all(s==1))
      {
        Q_start<-diag(diag(q_start),sum(s))
      }else{
        Q_start<-matrix(0,sum(s),sum(s))
        Q_start[1:s[1],1:s[1]]<-q_start[1:s[1],1:s[1]]
        
        for (zu in 2:rnd.len)
        {
          Q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
        }
      }
    }
    
    U<-X[,!is.na(index.new),drop=FALSE]
    X<-X[,is.na(index.new),drop=FALSE]
    
    final.names<-c(colnames(X),colnames(U))
    
    q<-dim(X)[2]
    
    
    Z_alles<-cbind(X,U,Phi,W)
    ########################################################## some definitions ################################################
    Delta<-matrix(0,control$steps,(lin+dim.smooth+n%*%s))
    Delta[1,1:lin]<-beta_null[final.names]
    Delta[1,(lin+1):(lin+dim.smooth)]<-smooth_null
    Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]<-t(ranef_null)
    
    control$epsilon<-control$epsilon*sqrt(dim(Delta)[2])
    
    Delta_start<-Delta[1,]
    
    active_old<-!is.element(Delta[1,],0)
    
    
    Eta.ma<-matrix(0,control$steps+1,N)
    Eta.ma[1,]<-Eta_start
    
    
    Q_inv<-NULL
    Q_inv.old.temp<-NULL
    Q_inv.start<-NULL
    
    ## start value for overdispersion
    if(control$overdispersion)
    {
      if(!is.null(control$phi_start))
      {  
        phi<-control$phi_start
      }else{
        if(is.null(control$Q.phi.start))
        {
          Q.phi.start<-q_start 
        }else{
          Q.phi.start<-control$Q.phi.start
        }
        active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,dim.smooth+n%*%s))
        Z_aktuell<-Z_alles[,active]
        lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
        
        if(rnd.len==1)
        {
          if(s==1)
          {
            Q.phi.start<-diag(Q.phi.start,s)
            P_akt<-c(rep(0,lin_akt+dim.smooth),rep((Q.phi.start^(-1)),n*s))
            F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
          }else{
            Q_inv.start<-chol2inv(chol(Q.phi.start))
            P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
            for(jf in 1:n)
              P_akt[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),
                    (lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv.start
            F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
          }
        }else{
          if(all(s==1))
          {
            Q.phi.start<-diag(diag(Q.phi.start),sum(s))
            P_akt<-c(rep(0,lin_akt+dim.smooth),rep(diag(Q.phi.start)^(-1),n))
            F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
          }else{
            Q_inv.start<-list()
            Q_inv.start[[1]]<-chol2inv(chol(Q.phi.start[1:s[1],1:s[1]]))
            
            for (zu in 2:rnd.len)
            {
              Q_inv.start[[zu]]<-chol2inv(chol(Q.phi.start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
            }
            
            P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
            for(jf in 1:n[1])
              P_akt[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),
                    (lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv.start[[1]]
            
            for (zu in 2:rnd.len)
            {
              for(jf in 1:n[zu])
                P_akt[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                      (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.start[[zu]]
            }
            F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
          }
        }
        InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
        if(class(InvFisher2)=="try-error")
          InvFisher2<-try(solve(F_gross),silent=T)
        if(class(InvFisher2)=="try-error")
          stop("Fisher matrix not invertible")  
        
        if(all(Mu==0) &family$family=="gaussian")
        {
          phi<-1  
        }else{
          Hat<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))#E-Uu
        phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(Hat)))
        }
        Sigma<-Sigma*phi
      }  
    }else{
      phi<-1
    }
    
    
    Q<-list()
    Q[[1]]<-Q_start
    
    l=1
    
    if(diff.ord>1)
    {
      k2<-c(rep(0,diff.ord-1),rep(1,nbasis-diff.ord+1))
    }else{
      k2<-rep(1,nbasis)
    }
    k22<-rep(k2,m)
    penal.vec<-penal*k22
    
     #if(rnd.len==1)
    #{
    #  if(s==1)
    #  {
    #    P.smooth<-c(rep(0,lin),penal.vec,rep(0,n*s))
    #  }else{
    #    P.smooth<-c(rep(0,lin),penal.vec,rep(0,n%*%s))
    #  }
    #}else{
    #  if(all(s==1))
    #  {
    #    P.smooth<-c(rep(0,lin),penal.vec,rep(0,rnd.len*n))
    # }else{
    #    P.smooth<-c(rep(0,lin),penal.vec,rep(0,n%*%s))
    #  }
    #}
    
   #browser()
    
    score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)#-t(t(P.smooth)*Delta[1,])
   
    lambda.max<-max(abs(score_vec[q:lin]))
   
    if (BLOCK)
    {
      grad.1<-gradient.lasso.block(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
    }else{
      grad.1<-gradient.lasso(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
    }
    
    score_vec<-c(score_vec[1:q],grad.1,score_vec[(lin+1):(lin+dim.smooth+n%*%s)])
    
    if(rnd.len==1)
    {
      if(s==1)
      {
        P1<-c(rep(0,lin),penal.vec,rep((Q_start^(-1)),n*s))
      }else{
        Q_inv.start<-chol2inv(chol(Q_start))
        P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
        diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
        for(j in 1:n)
          P1[(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s),(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s)]<-Q_inv
      }
    }else{
      if(all(s==1))
      {
        P1<-c(rep(0,lin),penal.vec,rep(diag(Q_start)^(-1),n))
      }else{
        Q_inv.start<-list()
        Q_inv.start[[1]]<-chol2inv(chol(Q_start[1:s[1],1:s[1]]))
        P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
        diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
        for(jf in 1:n[1])
          P1[(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1]),(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1])]<-Q_inv[[1]]
        
        for (zu in 2:rnd.len)
        {
          Q_inv.start[[zu]]<-chol2inv(chol(Q_start[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
          for(jf in 1:n[zu])
            P1[(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
               (lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
        }
      }
    }

    if(all(s==1))
    {
      F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+diag(P1)
    }else{
      F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
    }
    
    crit.obj<-t.change(grad=score_vec[(q+1):lin],b=Delta[1,(q+1):lin])
    t_edge<-crit.obj$min.rate
    
    grad.2<-t(score_vec)%*%F_gross%*%score_vec
    
    t_opt<-l2norm(score_vec)$length/grad.2
    nue<-control$nue
    
    half.index<-0

      solve.test2<-FALSE  
      while(!solve.test2)
      {  
        
        if(half.index>50)
        {
          stop("Fisher matrix not invertible")
        }
        
        Delta[1,]<-Delta_start+min(t_opt,t_edge)*nue*(0.5^half.index)*score_vec
        if(t_opt>t_edge & half.index==0)
          Delta[1,crit.obj$whichmin+q]<-0  
        
        Eta<-Z_alles%*%Delta[1,]
        Mu<-as.vector(family$linkinv(Eta))
        Sigma<-as.vector(family$variance(Mu))
        D<-as.vector(family$mu.eta(Eta))
        
        active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,dim.smooth+n%*%s))
        Z_aktuell<-Z_alles[,active]
        lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
        
        if (control$method=="EM" || control$overdispersion)
        {  
          if(rnd.len==1)
          {
            if(s==1)
            {
              P_akt<-c(rep(0,lin_akt),penal.vec,rep((Q_start^(-1)),n*s))
              F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
            }else{
              P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
              diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
              for(jf in 1:n)
                P_akt[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv.start
              F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
            }
          }else{
            if(all(s==1))
            {
              P_akt<-c(rep(0,lin_akt),penal.vec,rep(diag(Q_start)^(-1),n))
              F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
            }else{
              P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
              diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
              for(jf in 1:n[1])
                P_akt[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv.start[[1]]
              
              for (zu in 2:rnd.len)
              {
                for(jf in 1:n[zu])
                  P_akt[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                        (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.start[[zu]]
              }
              F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
            }
          }
          InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
          if(class(InvFisher2)=="try-error")
            InvFisher2<-try(solve(F_gross),silent=T)
          if(class(InvFisher2)=="try-error")
          {
            #stop("Fisher matrix not invertible")  
            half.index<-half.index+1  
          }else{
            solve.test2<-TRUE 
          }}else{
            solve.test2<-TRUE
          }
      }
      
      betaact<-Delta[1,active]
      
      if (control$method=="EM")
      {   
        ############################# Q update ################
        if(rnd.len==1)
        {
          Q1<-InvFisher2[(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s),(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s)]+Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+s)]%*%t(Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+s)])
          for (i in 2:n)
            Q1<-Q1+InvFisher2[(lin_akt+dim.smooth+(i-1)*s+1):(lin_akt+dim.smooth+i*s),(lin_akt+dim.smooth+(i-1)*s+1):(lin_akt+dim.smooth+i*s)]+Delta[1,(lin+dim.smooth+(i-1)*s+1):(lin+dim.smooth+i*s)]%*%t(Delta[1,(lin+dim.smooth+(i-1)*s+1):(lin+dim.smooth+i*s)])
          Q1<-1/n*Q1
        }else{
          Q1<-matrix(0,sum(s),sum(s))
          Q1[1:s[1],1:s[1]]<-InvFisher2[(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s[1]),(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s[1])]+Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+s[1])]%*%t(Delta[1,(lin+dim.smooth+1):(lin+dim.smooth+s[1])])
          for (i in 2:n[1])
            Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher2[(lin_akt+dim.smooth+(i-1)*s[1]+1):(lin_akt+dim.smooth+i*s[1]),(lin_akt+dim.smooth+(i-1)*s[1]+1):(lin_akt+dim.smooth+i*s[1])]+Delta[1,(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])]%*%t(Delta[1,(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])])
          Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]
          
          for (zu in 2:rnd.len)
          {
            Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher2[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
            for (i in 2:n[zu])
              Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher2[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[1,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
            Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-1/n[zu]*Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
          }
        }
      }else{
        Eta_tilde<-Eta+(y-Mu)*1/D
        Betadach<-Delta[1,1:(lin+dim.smooth)]     
        aktuell_vec<-!is.element(Delta[1,1:(lin)],0)
        X_aktuell<-Z_fastalles[,aktuell_vec]
        
        if(rnd.len==1)
        {
          
          if(s==1)
          {
            upp<-min(20,50*Q_start)
            low<-1e-14
            optim.obj<-try(nlminb(sqrt(Q_start),likelihood_nlminb,D=D,Sigma=Sigma,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),
                                  Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp))
            if(class(optim.obj)=="try-error")
              optim.obj<-try(bobyqa(sqrt(Q_start),likelihood_nlminb,D=D,Sigma=Sigma,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),
                                    Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp))
            Q1<-as.matrix(optim.obj$par)^2
          }else{
            q_start_vec<-c(diag(q_start),q_start[lower.tri(q_start)])
            up1<-min(20,50*max(q_start_vec))#[(s+1):(s*(s+1)*0.5)]))
            upp<-rep(up1,length(q_start_vec))
            low<-c(rep(0,s),rep(-up1,0.5*(s^2-s)))
            #   kkk_vec<-c(rep(-1,s),rep(0.5,0.5*(s^2-s)))
            optim.obj<-try(bobyqa(q_start_vec,likelihood,D=D,Sigma=Sigma,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),
                                  Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp))
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
            optim.obj<-try(bobyqa(sqrt(q_start_vec),likelihood_diag,D=D,Sigma=Sigma,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
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
            optim.obj<-try(bobyqa(q_start_vec,likelihood_block,D=D,Sigma=Sigma,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
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
        FinalHat<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))#E-Uu
        phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(FinalHat)))
        Sigma<-Sigma*phi
      }
      
      Eta.old<-Eta
      
      vorz<-F
    
   # browser()
    score_vec2<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)#-t(t(P.smooth)*Delta[1,])
    lambda.max<-max(abs(score_vec2[q:lin]))
   
      if (BLOCK)
      {
        grad.1<-gradient.lasso.block(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
      }else{
        grad.1<-gradient.lasso(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
      }
      score_vec2<-c(score_vec2[1:q],grad.1,score_vec2[(lin+1):(lin+dim.smooth+n%*%s)])
      
    if(rnd.len==1)
    {
      if(s==1)
      {
        P1<-c(rep(0,lin),penal.vec,rep((Q1^(-1)),n*s))
      }else{
        Q_inv<-chol2inv(chol(Q1))
        Q_inv.old.temp<-Q_inv
        P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
        diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
        for(j in 1:n)
          P1[(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s),(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s)]<-Q_inv
      }
    }else{
      if(all(s==1))
      {
        P1<-c(rep(0,lin),penal.vec,rep(diag(Q1)^(-1),n))
      }else{
        Q_inv<-list()
        Q_inv[[1]]<-chol2inv(chol(Q1[1:s[1],1:s[1]]))
        P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
        diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
        for(jf in 1:n[1])
          P1[(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1]),(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1])]<-Q_inv[[1]]
        
        for (zu in 2:rnd.len)
        {
          Q_inv[[zu]]<-chol2inv(chol(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
          for(jf in 1:n[zu])
            P1[(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
               (lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
        }
      }
    }
    
    if(all(s==1))
      {
        F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+diag(P1)
      }else{
        F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
      }
      
    crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin])
    t_edge<-crit.obj$min.rate
    
      grad.2<-t(score_vec2)%*%F_gross%*%score_vec2
      
      t_opt<-l2norm(score_vec2)$length/grad.2
       
    Eta.ma[2,]<-Eta
    
    score_vec<-score_vec2
    Q1.old<-Q1
    Q_inv.old<-Q_inv    
    
    Q1.very.old<-Q_start
    Q_inv.very.old<-Q_inv.start
    
    ###############################################################################################################################################
    ################################################################### Main Iteration ###################################################################
    if(control$steps!=1)
    {
      for (l in 2:control$steps)
      {
        if(control$print.iter)
          print(paste("Iteration ", l,sep=""))
        
 
        
        
        half.index<-0

        solve.test2<-FALSE  
          while(!solve.test2)
          {  
            
            if(half.index>50)
            {
              half.index<-Inf;Q1.old<-Q1.very.old;Q_inv.old<-Q_inv.very.old;
            }
            
              Delta[l,]<-Delta[l-1,]+min(t_opt,t_edge)*nue*(0.5^half.index)*score_vec
              if(t_opt>t_edge & half.index==0)
                Delta[l,crit.obj$whichmin+q]<-0
            
            
            Eta<-Z_alles%*%Delta[l,]
            Mu<-as.vector(family$linkinv(Eta))
            Sigma<-as.vector(family$variance(Mu))
            D<-as.vector(family$mu.eta(Eta))
            
            active_old<-active
            active<-c(rep(T,q),!is.element(Delta[l,(q+1):lin],0),rep(T,dim.smooth+n%*%s))
            Z_aktuell<-Z_alles[,active]
            lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))
            
            if (control$method=="EM" || control$overdispersion)
            {  
              if(rnd.len==1)
              {
                if(s==1)
                {
                  P_akt<-c(rep(0,lin_akt),penal.vec,rep((Q1.old^(-1)),n*s))
                  F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
                }else{
                  P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
                  diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
                  for(jf in 1:n)
                    P_akt[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),
                          (lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv.old
                  F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
                }
              }else{
                if(all(s==1))
                {
                  P_akt<-c(rep(0,lin_akt),penal.vec,rep(diag(Q1.old)^(-1),n))
                  F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
                }else{
                  P_akt<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
                  diag(P_akt)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
                  for(jf in 1:n[1])
                    P_akt[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),
                          (lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv.old[[1]]
                  
                  for (zu in 2:rnd.len)
                  {
                    for(jf in 1:n[zu])
                      P_akt[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                            (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv.old[[zu]]
                  }
                  F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P_akt
                }
              }
              InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
              if(class(InvFisher2)=="try-error")
                InvFisher2<-try(solve(F_gross),silent=T)
              if(class(InvFisher2)=="try-error")
              {
                half.index<-half.index+1  
              }else{
                solve.test2<-TRUE 
              }}else{
                solve.test2<-TRUE 
              }
          }
          
          betaact<-Delta[l,active]
          
          if (control$method=="EM")
          {        
            ############################# Q update ################
            if(rnd.len==1)
            {
              Q1<-InvFisher2[(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s),(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s)]+Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+s)]%*%t(Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+s)])
              for (i in 2:n)
                Q1<-Q1+InvFisher2[(lin_akt+dim.smooth+(i-1)*s+1):(lin_akt+dim.smooth+i*s),(lin_akt+dim.smooth+(i-1)*s+1):(lin_akt+dim.smooth+i*s)]+Delta[l,(lin+dim.smooth+(i-1)*s+1):(lin+dim.smooth+i*s)]%*%t(Delta[l,(lin+dim.smooth+(i-1)*s+1):(lin+dim.smooth+i*s)])
              Q1<-1/n*Q1
            }else{
              Q1<-matrix(0,sum(s),sum(s))
              Q1[1:s[1],1:s[1]]<-InvFisher2[(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s[1]),(lin_akt+dim.smooth+1):(lin_akt+dim.smooth+s[1])]+Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+s[1])]%*%t(Delta[l,(lin+dim.smooth+1):(lin+dim.smooth+s[1])])
              for (i in 2:n[1])
                Q1[1:s[1],1:s[1]]<-Q1[1:s[1],1:s[1]]+InvFisher2[(lin_akt+dim.smooth+(i-1)*s[1]+1):(lin_akt+dim.smooth+i*s[1]),(lin_akt+dim.smooth+(i-1)*s[1]+1):(lin_akt+dim.smooth+i*s[1])]+Delta[l,(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])]%*%t(Delta[l,(lin+dim.smooth+(i-1)*s[1]+1):(lin+dim.smooth+i*s[1])])
              Q1[1:s[1],1:s[1]]<-1/n[1]*Q1[1:s[1],1:s[1]]
              
              for (zu in 2:rnd.len)
              {
                Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-InvFisher2[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu]),(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]+Delta[l,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])]%*%t(Delta[l,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+s[zu])])
                for (i in 2:n[zu])
                  Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]+InvFisher2[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu]),(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]+Delta[l,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])]%*%t(Delta[l,(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(i-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+i*s[zu])])
                Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]<-1/n[zu]*Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]
              }
            }  
          }else{
            Eta_tilde<-Eta+(y-Mu)*1/D
            
            Betadach<-Delta[l,1:(lin+dim.smooth)]
            
            aktuell_vec<-!is.element(Delta[l,1:(lin)],0)
            X_aktuell<-Z_fastalles[,aktuell_vec]
            
            if(rnd.len==1)
            {
              
              if(s==1)
              {
                if(Q1<1e-14)
                  low<-0
                
                optim.obj<-try(nlminb(sqrt(Q1),likelihood_nlminb,D=D,Sigma=Sigma,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp))
                if(class(optim.obj)=="try-error")
                  optim.obj<-bobyqa(sqrt(Q1),likelihood_nlminb,D=D,Sigma=Sigma,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,Betadach=Betadach,W=W, lower = low, upper = upp)
                
                Q1<-as.matrix(optim.obj$par)^2
              }else{
                Q1_vec<-c(diag(Q1),Q1[lower.tri(Q1)])
                optim.obj<-try(bobyqa(Q1_vec,likelihood,D=D,Sigma=Sigma,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp))
                
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
                optim.obj<-try(bobyqa(sqrt(Q1_vec),likelihood_diag,D=D,Sigma=Sigma,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
                Q1<-diag(optim.obj$par)^2
              }else{
                Q1_vec<-c(diag(Q1)[1:s[1]],Q1[1:s[1],1:s[1]][lower.tri(Q1[1:s[1],1:s[1]])])
                
                for (zu in 2:rnd.len)
                  Q1_vec<-c(Q1_vec,c(diag(Q1)[(sum(s[1:(zu-1)])+1):sum(s[1:zu])],Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])][lower.tri(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])])]))
                
                optim.obj<-try(bobyqa(Q1_vec,likelihood_block,D=D,Sigma=Sigma,X=cbind(Z_fastalles,Phi),X_aktuell=cbind(X_aktuell,Phi),Eta_tilde=Eta_tilde,n=n,s=s,k=k,Betadach=Betadach,W=W, lower=low,upper=upp,rnd.len=rnd.len))
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
            FinalHat<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))#E-Uu
            phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(FinalHat)))
            Sigma<-Sigma*phi
          }
          
          Q[[l+1]]<-Q1
          
        score_vec2<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)#-t(t(P.smooth)*Delta[l,])
        score.pure<-score_vec2
        lambda.max<-max(abs(score_vec2[(q+1):lin]))
        
        
          if (BLOCK)
          {
            grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda,block=block)
          }else{
            grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda)
          }
          score_vec2<-c(score_vec2[1:q],grad.1,score_vec2[(lin+1):(lin+dim.smooth+n%*%s)])
          
        if(rnd.len==1)
        {
          if(s==1)
          {
            P1<-c(rep(0,lin),penal.vec,rep((Q1^(-1)),n*s))
          }else{
            Q_inv<-chol2inv(chol(Q1))
            P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
            diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
            for(j in 1:n)
              P1[(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s),(lin+dim.smooth+(j-1)*s+1):(lin+dim.smooth+j*s)]<-Q_inv
          }
        }else{
          if(all(s==1))
          {
            P1<-c(rep(0,lin),penal.vec,rep(diag(Q1)^(-1),n))
          }else{
            Q_inv<-list()
            Q_inv[[1]]<-chol2inv(chol(Q1[1:s[1],1:s[1]]))
            P1<-matrix(0,lin+dim.smooth+n%*%s,lin+dim.smooth+n%*%s)
            diag(P1)[(lin+1):(lin+dim.smooth)]<-penal.vec
            for(jf in 1:n[1])
              P1[(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1]),(lin+dim.smooth+(jf-1)*s[1]+1):(lin+dim.smooth+jf*s[1])]<-Q_inv[[1]]
            
            for (zu in 2:rnd.len)
            {
              Q_inv[[zu]]<-chol2inv(chol(Q1[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
              for(jf in 1:n[zu])
                P1[(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                   (lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
            }
          }
        }
        
        if(all(s==1))
          {
            F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+diag(P1)
          }else{
            F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
          }
          
        crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin])
        t_edge<-crit.obj$min.rate
        
          grad.2<-t(score_vec2)%*%F_gross%*%score_vec2
          
          t_opt<-l2norm(score_vec2)$length/grad.2
          
        score_vec<-score_vec2
        Q1.very.old<-Q1.old
        Q_inv.very.old<-Q_inv.old
        Q1.old<-Q1
        Q_inv.old<-Q_inv    
                
        Eta.ma[l+1,]<-Eta
        
        finish<-(sqrt(sum((Eta.ma[l,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l,])^2))<control$epsilon)
        finish2<-(sqrt(sum((Eta.ma[l-1,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l-1,])^2))<control$epsilon)
        if(finish ||  finish2) #|| (all(grad.1 == 0) ))
          break
        Eta.old<-Eta
      }}
    
    
    ######## Final calculation
    FinalHat.df<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))
    df<-sum(diag(FinalHat.df))
    
    conv.step<-l
    phi.med<-phi

   #browser()
   
   
    if(conv.step==control$steps)
    {
      cat("Warning:\n")
      cat("Algorithm did not converge!\n")
    }
    
    Delta_neu<-Delta[l,]
    Eta_opt<-Z_alles%*%Delta_neu
    Mu_opt<-as.vector(family$linkinv(Eta_opt))
    Sigma_opt<-as.vector(family$variance(Mu_opt))    
    D_opt<-as.vector(family$mu.eta(Eta_opt))
    Qfinal<-Q[[l+1]]
    
    aaa<-!is.element(Delta_neu[1:(lin)],0)
    
    if(rnd.len==1)
    {
      if(s==1)
      {
        P1.ran<-rep((Qfinal^(-1)),n*s)
        P1.ran<-diag(P1.ran)
      }else{
        P1.ran<-matrix(0,n*s,n*s)
        for(jf in 1:n)
          P1.ran[((jf-1)*s+1):(jf*s),((jf-1)*s+1):(jf*s)]<-chol2inv(chol(Qfinal))
      }
    }else{
      if(all(s==1))
      {
        P1.ran<-rep(diag(Qfinal)^(-1),n)
        P1.ran<-diag(P1.ran)
      }else{
        P1.ran<-matrix(0,n%*%s,n%*%s)
        inv.act<-chol2inv(chol(Qfinal[1:s[1],1:s[1]]))
        for(jf in 1:n[1])
          P1.ran[( (jf-1)*s[1]+1):( jf*s[1]),( (jf-1)*s[1]+1):( jf*s[1])]<-inv.act
        
        for (zu in 2:rnd.len)
        {
          inv.act<-chol2inv(chol(Qfinal[(sum(s[1:(zu-1)])+1):sum(s[1:zu]),(sum(s[1:(zu-1)])+1):sum(s[1:zu])]))
          for(jf in 1:n[zu])
            P1.ran[( n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):( n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                   ( n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):( n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-inv.act
        }
      }  
    }  
    
    ranef.logLik<--0.5*t(Delta_neu[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)])%*%P1.ran%*%Delta_neu[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]
    
    
    if(final.re)
    {    
      ############ final re-estimation
      
      if(s==1)
      {  
        Q.max<-max(sqrt(unlist(Q)))
        Q.min<-min(sqrt(unlist(Q)))
      }else{
        Q.max<-max(Qfinal)+1
        Q.min<-min(Qfinal)-1e-10
      }
      
      
      if(rnd.len==1)
      {
        glmm_fin<-try(glmm_final_smooth(y,Z_fastalles[,aaa],Phi,W,k,penal.vec,q_start=Qfinal,
                                        Delta_start=Delta_neu[c(aaa,rep(T,dim.smooth+n%*%s))],
                                        s,steps=control$maxIter,family=family,method=control$method.final,overdispersion=control$overdispersion,
                                        phi=phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final,
                                        Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac))#,silent = TRUE)
        if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
        {  
          glmm_fin2<-try(glmm_final_smooth(y,Z_fastalles[,aaa],Phi,W,k,penal.vec,q_start=q_start,
                                           Delta_start=Delta_start[c(aaa,rep(T,dim.smooth+n%*%s))],
                                           s,steps=control$maxIter,family=family,method=control$method.final,overdispersion=control$overdispersion,
                                           phi=control$phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final,
                                           Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac))
          if(class(glmm_fin2)!="try-error")
          {    
            if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
              glmm_fin<-glmm_fin2 
          }
        }
      }else{
        glmm_fin<-try(glmm_final_multi_random_smooth(y,Z_fastalles[,aaa],Phi,W,k,penal.vec,q_start=Qfinal,
                                                     Delta_start=Delta_neu[c(aaa,rep(T,dim.smooth+n%*%s))],
                                                     s,n,steps=control$maxIter,family=family,method=control$method.final,overdispersion=control$overdispersion,
                                                     phi=phi,rnd.len=rnd.len,print.iter.final=control$print.iter.final,eps.final=control$eps.final,
                                                     Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac))#,silent = TRUE)
        if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
        {  
          glmm_fin2<-try(glmm_final_multi_random_smooth(y,Z_fastalles[,aaa],Phi,W,k,penal.vec,q_start=q_start,
                                                        Delta_start=Delta_start[c(aaa,rep(T,dim.smooth+n%*%s))],
                                                        s,n,steps=control$maxIter,family=family,method=control$method.final,overdispersion=control$overdispersion,
                                                        phi=control$phi,rnd.len=rnd.len,print.iter.final=control$print.iter.final,
                                                        eps.final=control$eps.final,Q.max=Q.max,Q.min=Q.min,Q.fac=control$Q.fac))
          if(class(glmm_fin2)!="try-error")
          {    
            if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
              glmm_fin<-glmm_fin2 
          }
        }
      }
      
      if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
      {
        cat("Warning:\n")
        cat("Final Fisher scoring reestimation did not converge!\n")
      }
      
      #######
    }else{
      glmm_fin<-NA  
      class(glmm_fin)<-"try-error"
    }  
    
    Delta_neu2<-Delta_neu
    Standard_errors<-rep(NA,length(Delta_neu))
    
    
    if(class(glmm_fin)!="try-error")
    {
      EDF.matrix<-glmm_fin$EDF.matrix
      complexity.smooth<-sum(diag(EDF.matrix)[c(T,rep(F,sum(aaa)-1),rep(T,dim.smooth),rep(F,n%*%s))])  
      Delta_neu2[c(aaa,rep(T,dim.smooth+n%*%s))]<-glmm_fin$Delta
      Standard_errors[c(aaa,rep(T,dim.smooth+n%*%s))]<-glmm_fin$Standard_errors
      Qfinal<-glmm_fin$Q
      phi<-glmm_fin$phi
      complexity<-glmm_fin$complexity
    }else{
      glmm_fin<-list()
      glmm_fin$ranef.logLik<-ranef.logLik
      complexity<-df
      
      lin_akt<-q+sum(!is.element(Delta_neu[(q+1):lin],0))
      if(rnd.len==1)
      {
        if(s==1)
        {
          P1<-c(rep(0,lin_akt),penal.vec,rep((Q1^(-1)),n*s))
          P1a<-c(rep(0,lin_akt+dim.smooth),rep((Q1^(-1)),n*s))
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P1)
        }else{
          P1<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
          P1a<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
          diag(P1)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
          for(jf in 1:n)
          {
            P1[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),
               (lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv
            P1a[(lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s),
                (lin_akt+dim.smooth+(jf-1)*s+1):(lin_akt+dim.smooth+jf*s)]<-Q_inv
          }  
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P1
        }
      }else{
        if(all(s==1))
        {
          P1<-c(rep(0,lin_akt),penal.vec,rep(diag(Q1)^(-1),n))
          P1a<-c(rep(0,lin_akt+dim.smooth),rep(diag(Q1)^(-1),n))
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P1)
        }else{
          P1<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
          P1a<-matrix(0,lin_akt+dim.smooth+n%*%s,lin_akt+dim.smooth+n%*%s)
          diag(P1)[(lin_akt+1):(lin_akt+dim.smooth)]<-penal.vec
          for(jf in 1:n[1])
          {
            P1[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),
               (lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv[[1]]
            P1a[(lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1]),
                (lin_akt+dim.smooth+(jf-1)*s[1]+1):(lin_akt+dim.smooth+jf*s[1])]<-Q_inv[[1]]
          }
          for (zu in 2:rnd.len)
          {
            for(jf in 1:n[zu])
            {
              P1[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                 (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
              P1a[(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu]),
                  (lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+(jf-1)*s[zu]+1):(lin_akt+dim.smooth+n[1:(zu-1)]%*%s[1:(zu-1)]+jf*s[zu])]<-Q_inv[[zu]]
            }}
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P1
        }
      }
      InvFisher<-try(chol2inv(chol(F_gross)),silent=T)
      if(class(InvFisher)=="try-error")
        InvFisher<-try(solve(F_gross),silent=T)
      if(class(InvFisher)=="try-error")
      {
        warning("No EDF's for smooth functions available, as Fisher matrix not invertible!")
        complexity.smooth<-dim.smooth
      }else{  
        ###### EDF of spline; compare Wood's Book on page 167
        EDF.matrix<-InvFisher%*%(t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+P1a)
        complexity.smooth<-sum(diag(EDF.matrix)[c(T,rep(F,sum(aaa)-1),rep(T,dim.smooth),rep(F,n%*%s))])
      }
    }  
    
    if(!(complexity.smooth>=1 && complexity.smooth<=dim.smooth))
      complexity.smooth<-dim.smooth
    
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
    
    names(Delta_neu)[1:dim(X)[2]]<-colnames(X)
    names(Standard_errors)[1:dim(X)[2]]<-colnames(X)
    
    if(lin>1)
    {
      names(Delta_neu)[(dim(X)[2]+1):lin]<-colnames(U)
      names(Standard_errors)[(dim(X)[2]+1):lin]<-colnames(U)
    }
    
    names(Delta_neu)[(lin+1):(lin+dim.smooth)]<-colnames(Phi)
    names(Standard_errors)[(lin+1):(lin+dim.smooth)]<-colnames(Phi)
    names(Delta_neu)[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]<-colnames(W)
    names(Standard_errors)[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]<-colnames(W)
    colnames(Delta)<-c(old.names,colnames(Phi),colnames(W))
    
    aic<-NaN
    bic<-NaN
    
    if (is.element(family$family,c("gaussian", "binomial", "poisson"))) 
    {
      loglik<-logLik.glmmLasso(y=y,mu=Mu_opt,beta=Delta_neu[(q+1):(lin)],ranef.logLik=glmm_fin$ranef.logLik,
                               lambda=lambda,family=family,penal=FALSE)
if(control$complexity!="hat.matrix")  
{  
        if(rnd.len==1)
      {
        complexity<-0.5*(s*(s+1))
      }else{
        complexity<-0.5*(s[1]*(s[1]+1))
        for(zu in 2:rnd.len)
          complexity<-complexity+0.5*(s[zu]*(s[zu]+1))
      }
      complexity<-complexity+sum(Delta_neu[1:(lin)]!=0)+complexity.smooth
}      
      aic<--2*loglik+2*complexity
      bic<--2*loglik+log(N)*complexity
    }else{
      warning("For the specified family (so far) no AIC and BIC are available!")  
    }
    
    ret.obj=list()
    ret.obj$aic<-aic
    ret.obj$bic<-bic
    ret.obj$Deltamatrix<-Delta
    ret.obj$smooth<-Delta_neu[(lin+1):(lin+dim.smooth)]
    ret.obj$ranef<-Delta_neu[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]
    ret.obj$coefficients<-Delta_neu[1:(lin)]
    ret.obj$fixerror<-Standard_errors[1:(lin)]
    ret.obj$ranerror<-Standard_errors[(lin+dim.smooth+1):(lin+dim.smooth+n%*%s)]
    ret.obj$Q_long<-Q
    ret.obj$Q<-Qfinal
    ret.obj$y_hat<-Mu_opt
    ret.obj$phi<-phi
    ret.obj$family<-family
    ret.obj$fix<-fix.old
    ret.obj$newrndfrml<-newrndfrml
    ret.obj$subject<-names(rnd)
    ret.obj$data<-data
    ret.obj$rnd.len<-rnd.len
    ret.obj$B<-B
    ret.obj$nbasis<-nbasis
    ret.obj$spline.degree<-spline.degree
    ret.obj$diff.ord<-diff.ord
    ret.obj$knots.no<-knots.no
    ret.obj$conv.step<-conv.step
    ret.obj$phi.med<-phi.med
    ret.obj$complexity.smooth<-complexity.smooth
    ret.obj$y <- y
    ret.obj$df<-df
    ret.obj$loglik<-loglik
    ret.obj$lambda.max<-lambda.max
ret.obj$score.pure<-score.pure
    return(ret.obj)
  }  
}
##################################################################
##################################################################
}
##################################################################
##################################################################
##################################################################
##################################################################



glmmLasso <- function(fix=formula, rnd=formula, data, lambda, family=NULL, switch.NR = TRUE, final.re=FALSE, 
                      control=list()) UseMethod("glmmLasso")

glmmLasso.formula <- function(fix, rnd, ...)
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

  if(!is.null(x$smooth))
  {
    cat("\nSmooth Effects:\n")
    print(colnames(x$B))
  }  
  
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
              coefficients=TAB,smooth.eff=colnames(object$B),StdDev=object$StdDev)
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
  
  if(!is.null(x$smooth))
  {
    cat("\nSmooth Effects:\n")
    print(x$smooth.eff)
  }  
  
  cat("\nRandom Effects:\n")
  cat("\nStdDev:\n")
  print(x$StdDev)
}


predict.glmmLasso <- function(object,newdata=NULL,new.random.design=NULL,...)
{
  if(is.null(newdata))
  {
    y<-fitted(object)
  }else{  
    rnd.len<-object$rnd.len                 
    family<-object$family
    
    X <- model.matrix(formula(object$fix), newdata)
    
   if(is.null(new.random.design))
   {   
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
        if(nrow(X)!=1)
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
      {
        if(krit.random[zu])
          dim.W.single[zu+1]<-dim(W.single[[zu]])[2]
      }
      
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
  }else{
    Design<-cbind(X,new.random.design)
    Delta<-c(object$coef,object$ranef)
    if(ncol(Design)!=length(Delta))
       stop("Wrong dimension of random effects design matrix!")
    y<- as.vector(family$linkinv(Design%*%Delta))   
  }}
  y
}


plot.glmmLasso <- function(x,which=NULL,plot.data=TRUE,include.icept=FALSE,ylab=NULL,main=NULL,...)
{
  if(is.null(ylab))
    ylab<-""
  
  if(is.null(main))
    main<-""
  
    
  if(is.null(x$B))
  stop("No smooth terms to plot!")
  
  Phi<-x$B
  m<-dim(Phi)[2]
  
  if(is.null(which))
    which<-1:m
  
  p<-length(which)
  if(p>9)
    stop("Too many smooth functions! Please specify at maximum nine.")
  
  a<-ceiling(sqrt(p))
  b<-floor(sqrt(p))
  if(b==0)
    b<-1
  
  nbasis<-x$nbasis
  diff.ord<-x$diff.ord
  spline.degree<-x$spline.degree
  
  knots.no<-nbasis-1
  if(spline.degree<3 && (spline.degree-diff.ord)<2)
    knots.no<-knots.no+1  
  
  spline.ma<-list()
  Design<-list()
  smooth.ma<-matrix(0,m,dim(Phi)[1])
  
  for(i in which)
  {
    if(plot.data)
    {      
    spline.ma[[i]]<-bs.design(sort(Phi[,i]), diff.ord=diff.ord, spline.degree=spline.degree, knots.no=knots.no)
    }else{
    smooth.ma<-matrix(0,m,1000)
    data.seq<-seq(min(Phi[,i]),max(Phi[,i]),length.out=1000)  
    spline.ma[[i]]<-bs.design(data.seq, diff.ord=diff.ord, spline.degree=spline.degree, knots.no=knots.no)
    }
    Design[[i]]<-cbind(spline.ma[[i]]$X[,-1],spline.ma[[i]]$Z)
    smooth.ma[i,]<-Design[[i]]%*%x$smooth[((i-1)*nbasis+1):(i*nbasis)]
  }
  
  
  par(mfrow=c(a,b))
  
  for(i in which)
  {
    if(include.icept  && is.element("(Intercept)",names(x$coef)))
    {
      
      if(plot.data)
      {      
      plot(sort(Phi[,i]), x$coef[match("(Intercept)",names(x$coef))]+smooth.ma[i,], 
           type = "l", lwd=2, xlab=paste(colnames(Phi)[i]), 
           ylab=ylab, main=main,cex.lab=2,cex.axis=2,...)
      }else{
      data.seq<-seq(min(Phi[,i]),max(Phi[,i]),length.out=1000)    
      plot(data.seq, x$coef[match("(Intercept)",names(x$coef))]+smooth.ma[i,], type = "l", lwd=2, 
           xlab=paste(colnames(Phi)[i]), ylab=ylab,main=main,cex.lab=2,cex.axis=2,...)
      }
       rug(jitter(Phi[,i]))
    }else{
      if(plot.data)
      {      
      plot(sort(Phi[,i]), smooth.ma[i,], type = "l", lwd=2, xlab=paste(colnames(Phi)[i]),
           ylab=ylab,main=main,cex.lab=2,cex.axis=2,...)
      }else{
      data.seq<-seq(min(Phi[,i]),max(Phi[,i]),length.out=1000)  
      plot(data.seq, smooth.ma[i,], type = "l", lwd=2, xlab=paste(colnames(Phi)[i]),
           ylab=ylab,main=main,cex.lab=2,cex.axis=2,...)
      }
       rug(jitter(Phi[,i]))
    }   
  }
}
