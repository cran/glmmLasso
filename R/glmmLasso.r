taylor.opt<-function(t_opt,y,X,fixef,ranef,Grad,family,P)
{
  delta<-c(fixef,ranef)+t_opt*Grad
  Eta<- X%*%delta
  mu<-as.vector(family$linkinv(Eta))
  
  loglik <- logLik.glmmLasso(y=y,mu=mu,family=family,ranef.logLik=NULL,penal=FALSE) 
    
  loglik<- -loglik + 0.5*t(delta[(length(fixef)+1):length(delta)]) %*% P %*% delta[(length(fixef)+1):length(delta)]
  return(loglik)
}

taylor.opt.noRE<-function(t_opt,y,X,fixef,Grad,family)
{
  delta<-fixef+t_opt*Grad
  Eta<- X%*%delta
  mu<-as.vector(family$linkinv(Eta))
  
  loglik <- -logLik.glmmLasso(y=y,mu=mu,family=family,ranef.logLik=NULL,penal=FALSE) 
  
  return(loglik)
}


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

#############################################################################################################################################               
###################################################### main Lasso function #################################################################

est.glmmLasso<-function(fix,rnd,data,lambda,family=gaussian(link = "identity"),
                        final.re=FALSE,switch.NR=T,control=list())
{  
  if(!is.null(rnd))
  {
    return.obj <- est.glmmLasso.RE(fix=fix,rnd=rnd,data=data,lambda=lambda,family=family,
                                   final.re=final.re,switch.NR=switch.NR,control=control)
  }else{
    return.obj <- est.glmmLasso.noRE(fix=fix,data=data,lambda=lambda,family=family,
                                     final.re=final.re,switch.NR=switch.NR,control=control)
  }
return(return.obj)
}
    
#############################################################################################################################################               
#############################################################################################################################################               
  
  
  
  
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
  
  if(!is.null(x$rnd))
  {  
  cat("\nRandom Effects:\n")
  cat("\nStdDev:\n")
  print(x$StdDev)
  }else{
    cat("\nNo random effects included!\n")
  }
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
  
  if(!is.null(x$rnd))
  {  
    cat("\nRandom Effects:\n")
    cat("\nStdDev:\n")
    print(x$StdDev)
  }else{
    cat("\nNo random effects included!\n")
  }
}


predict.glmmLasso <- function(object,newdata=NULL,new.random.design=NULL,...)
{
  if(is.null(newdata))
  {
    y<-fitted(object)
  }else{  
    X <- model.matrix(formula(object$fix), newdata)
    family<-object$family
    
    if(!is.null(object$rnd))
    {  
    rnd.len<-object$rnd.len                 
    
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
  }}else{
    y <- X%*%object$coef
  }
  }
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
