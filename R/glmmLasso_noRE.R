est.glmmLasso.noRE<-function(fix,data,lambda,family=gaussian(link = "identity"),
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
  
  if(ncol(X)==1)
  {
    if(colnames(X)=="(Intercept)")
      stop("No terms to select! Use glmer, glmmPQL or glmmML!")
  }

  old.names<-attr(X,"dimnames")[[2]]

  if(control$print.iter)
    message("Iteration 1")

  very.old.names<-very.old.names[!is.na(control$index)]
  
  block<-as.numeric(table(index.new[!is.na(index.new)]))
  
  BLOCK<-FALSE
  if(!all(block==1))
    BLOCK<-TRUE
  
  lin<-ncol(X)
  
  if(is.null(control$start))
    control$start<-c(rep(0,lin))
  
  N<-length(y)
  
  beta_null<-control$start[1:lin]
  if(is.null(attr(beta_null,"names")))
    attr(beta_null,"names")<-orig.names
  beta_null<-beta_null[colnames(X)]
  
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
        Eta_start<-X%*%beta_null
      }else{
        Eta_start<-rep(beta_null,N)
      }
      
      D<-as.vector(family$mu.eta(Eta_start))
      Mu<-as.vector(family$linkinv(Eta_start))
      Sigma<-as.vector(family$variance(Mu))
      
      U<-X[,!is.na(index.new),drop=FALSE]
      X<-X[,is.na(index.new),drop=FALSE]
      
      final.names<-c(colnames(X),colnames(U))
      
      q<-ncol(X)
      
      Z_alles<-cbind(X,U)
      ########################################################## some definitions ################################################
      Delta<-matrix(0,control$steps,lin)
      Delta[1,1:lin]<-beta_null[final.names]

      Eta.ma<-matrix(0,control$steps+1,N)
      Eta.ma[1,]<-Eta_start
      
      logLik.vec<-c()
      
      control$epsilon<-control$epsilon*sqrt(dim(Delta)[2])
      
      Delta_start<-Delta[1,]
      
      active_old<-!is.element(Delta[1,],0)
      
     ## start value for overdispersion
      if(control$overdispersion)
      {
        if(!is.null(control$phi_start))
        {  
          phi<-control$phi_start
        }else{
          active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0))
          Z_aktuell<-Z_alles[,active]
          lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
          
          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)

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
      
      l=1

      score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)
      lambda.max<-max(abs(score_vec[(q+1):lin]))
      
      if (BLOCK)
      {
        grad.1<-gradient.lasso.block(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
      }else{
        grad.1<-gradient.lasso(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
      }
      
      score_vec<-c(score_vec[1:q],grad.1)
      
      F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)
      
      crit.obj<-t.change(grad=score_vec[(q+1):lin],b=Delta[1,(q+1):lin])
      t_edge<-crit.obj$min.rate
      
      grad.2<-t(score_vec)%*%F_gross%*%score_vec
      
      optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt.noRE,y=y,X=Z_alles,fixef=Delta_start[1:lin],
                                         Grad=score_vec,family=family,lower = 0, upper = Inf))
      
      t_opt<-optim.obj$par
      
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
          
          logLik.vec[1]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=NULL,family=family,penal=F)
          
          active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0))
          Z_aktuell<-Z_alles[,active]
          lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
          
          if (control$overdispersion)
          {  
                F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)

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
        
        if(control$overdispersion)
        {    
          FinalHat<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))#E-Uu
          phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(FinalHat)))
          Sigma<-Sigma*phi
        }
        
         Eta.old<-Eta
        
        vorz<-F
        
        score_vec2<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)
        lambda.max<-max(abs(score_vec2[(q+1):lin]))
        
        score_old2<-score_vec2
        
        if (BLOCK)
        {
          grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
        }else{
          grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
        }
        score_vec2<-c(score_vec2[1:q],grad.1)
        
        F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)
        
        crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin])
        t_edge<-crit.obj$min.rate
        
        grad.2<-t(score_vec2)%*%F_gross%*%score_vec2
        
        optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt.noRE,y=y,X=Z_alles,fixef=Delta[1,1:lin],
                                           Grad=score_vec2,family=family,lower = 0, upper = Inf))
        
        t_opt<-optim.obj$par
        
        tryNR<- (t_opt<t_edge) #&& !(all(active_old==active)  && !NRstep)  
        
        if(tryNR) 
        {
          lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))

          F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)

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
      ###############################################################################################################################################
      ################################################################### Main Iteration ###################################################################
      
      if(control$steps!=1)
      {
        for (l in 2:control$steps)
        {
          if(control$print.iter)
            message("Iteration ",l)
          #print(paste("Iteration ", l,sep=""))
          
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
                  half.index<-Inf

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
              
              logLik.vec[l]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=NULL,family=family,penal=F)
              
              active_old<-active
              active<-c(rep(T,q),!is.element(Delta[l,(q+1):lin],0))
              Z_aktuell<-Z_alles[,active]
              lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))
              
              if (control$overdispersion)
              {  
               F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)
               
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
            
            if(control$overdispersion)
            {
              FinalHat<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))#E-Uu
              phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(FinalHat)))
              Sigma<-Sigma*phi
            }
            
            score_vec2<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)
            lambda.max<-max(abs(score_vec2[(q+1):lin]))
            
            score_old2<-score_vec2
            
            if (BLOCK)
            {
              grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda,block=block)
            }else{
              grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda)
            }
            score_vec2<-c(score_vec2[1:q],grad.1)
            
            F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)
            
            crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin])
            t_edge<-crit.obj$min.rate
            
            grad.2<-t(score_vec2)%*%F_gross%*%score_vec2
            
            optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt.noRE,y=y,X=Z_alles,fixef=Delta[l,1:lin],
                                               Grad=score_vec2,family=family,lower = 0, upper = Inf))
            
            t_opt<-optim.obj$par
            
            tryNR<- (t_opt<t_edge) #&& !(all(active_old==active)  && !NRstep)  
            
            if(tryNR) 
            {
              lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))
              F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)

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

          Eta.ma[l+1,]<-Eta
          
          finish<-(sqrt(sum((Eta.ma[l,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l,])^2))<control$epsilon)
          finish2<-(sqrt(sum((Eta.ma[l-1,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l-1,])^2))<control$epsilon)
          if(finish ||  finish2) #|| (all(grad.1 == 0) ))
            break
          Eta.old<-Eta
        }}
      
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

      if(final.re)
      {    
        ############ final re-estimation
        aaa<-!is.element(Delta_neu[1:lin],0)
        
          glmm_fin<-try(glmm_final_noRE(y,Z_fastalles[,aaa],
                                   Delta_start=Delta_neu[aaa],steps=control$maxIter,
                                   family=family,overdispersion=control$overdispersion,
                                   phi=phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final),silent = TRUE)
          
          if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
          {  
            glmm_fin2<-try(glmm_final_noRE(y,Z_fastalles[,aaa],
                                      Delta_start=Delta_start[aaa],steps=control$maxIter,
                                      family=family,overdispersion=control$overdispersion,
                                      phi=control$phi,print.iter.final=control$print.iter.final,
                                      eps.final=control$eps.final),silent = TRUE)
            if(class(glmm_fin2)!="try-error")
            {    
              if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
                glmm_fin<-glmm_fin2 
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
        Delta_neu2[aaa]<-glmm_fin$Delta
        Standard_errors[aaa]<-glmm_fin$Standard_errors
        phi<-glmm_fin$phi
        complexity<-glmm_fin$complexity
      }else{
        glmm_fin<-list()
        F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)
        InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
        if(class(InvFisher2)=="try-error")
          InvFisher2<-try(solve(F_gross),silent=T)
        if(class(InvFisher2)!="try-error")
        {  
        FinalHat.df<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))
        df<-sum(diag(FinalHat.df))
        }else{
        df <- NA
        }
        complexity<-df
      }
      
      Delta_neu<-Delta_neu2
      
      Eta_opt<-Z_alles%*%Delta_neu
      Mu_opt<-as.vector(family$linkinv(Eta_opt))
      
      names(Delta_neu)[1:dim(X)[2]]<-colnames(X)
      names(Standard_errors)[1:dim(X)[2]]<-colnames(X)
      
      if(lin>1)
      {
        names(Delta_neu)[(dim(X)[2]+1):lin]<-colnames(U)
        names(Standard_errors)[(dim(X)[2]+1):lin]<-colnames(U)
      }
      
      colnames(Delta)<-final.names
      
      aic<-NaN
      bic<-NaN
      
      if(is.element(family$family,c("gaussian", "binomial", "poisson"))) 
      {
        
        loglik<-logLik.glmmLasso(y=y,mu=Mu_opt,ranef.logLik=NULL,family=family,penal=F)
        
        
        if(control$complexity!="hat.matrix")  
           complexity<-sum(Delta_neu[1:(lin)]!=0)
 
        aic<--2*loglik+2*complexity
        bic<--2*loglik+log(N)*complexity
      }else{
        warning("For the specified family (so far) no AIC and BIC are available!")  
      }

      ret.obj=list()
      ret.obj$aic<-aic
      ret.obj$bic<-bic
      ret.obj$Deltamatrix<-Delta
      ret.obj$coefficients<-Delta_neu[1:(lin)]
      ret.obj$fixerror<-Standard_errors[1:(lin)]
      ret.obj$y_hat<-Mu_opt
      ret.obj$phi<-phi
      ret.obj$family<-family
      ret.obj$fix<-fix.old
      ret.obj$conv.step<-conv.step
      ret.obj$data<-data
      ret.obj$phi.med<-phi.med
      ret.obj$y <- y
      ret.obj$df<-df
      ret.obj$loglik<-loglik
      ret.obj$lambda.max<-lambda.max
      ret.obj$logLik.vec<-logLik.vec
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
        Eta_start<-X%*%beta_null+Phi%*%smooth_null
      }else{
        Eta_start<-rep(beta_null,N)+Phi%*%smooth_null
      }
      
      D<-as.vector(family$mu.eta(Eta_start))
      Mu<-as.vector(family$linkinv(Eta_start))
      Sigma<-as.vector(family$variance(Mu))
      
      U<-X[,!is.na(index.new),drop=FALSE]
      X<-X[,is.na(index.new),drop=FALSE]
      
      final.names<-c(colnames(X),colnames(U))
      
      q<-dim(X)[2]
      
      Z_alles<-cbind(X,U,Phi)
      ########################################################## some definitions ################################################
      Delta<-matrix(0,control$steps,lin+dim.smooth)
      Delta[1,1:lin]<-beta_null[final.names]
      Delta[1,(lin+1):(lin+dim.smooth)]<-smooth_null

      control$epsilon<-control$epsilon*sqrt(dim(Delta)[2])
      
      Delta_start<-Delta[1,]
      
      logLik.vec<-c()
      
      active_old<-!is.element(Delta[1,],0)
      
      Eta.ma<-matrix(0,control$steps+1,N)
      Eta.ma[1,]<-Eta_start
      
      ## start value for overdispersion
      if(control$overdispersion)
      {
        if(!is.null(control$phi_start))
        {  
          phi<-control$phi_start
        }else{
          active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,dim.smooth))
          Z_aktuell<-Z_alles[,active]
          lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
          
              P_akt<-rep(0,lin_akt+dim.smooth)
              F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)
          
          InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
          if(class(InvFisher2)=="try-error")
            InvFisher2<-try(solve(F_gross),silent=T)
          if(class(InvFisher2)=="try-error")
            stop("Fisher matrix not invertible")  
          
          if(all(Mu==0) & family$family=="gaussian")
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

      l=1
      
      if(diff.ord>1)
      {
        k2<-c(rep(0,diff.ord-1),rep(1,nbasis-diff.ord+1))
      }else{
        k2<-rep(1,nbasis)
      }
      k22<-rep(k2,m)
      penal.vec<-penal*k22
      
          P1<-c(rep(0,lin),penal.vec)
          P1<-diag(P1)

      score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-P1%*%Delta[1,]
      lambda.max<-max(abs(score_vec[(q+1):lin]))
      
      
      if (BLOCK)
      {
        grad.1<-gradient.lasso.block(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
      }else{
        grad.1<-gradient.lasso(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
      }
      
      score_vec<-c(score_vec[1:q],grad.1,score_vec[(lin+1):(lin+dim.smooth)])
      
      F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
      
      crit.obj<-t.change(grad=score_vec[(q+1):lin],b=Delta[1,(q+1):lin])
      t_edge<-crit.obj$min.rate
      
      grad.2<-t(score_vec)%*%F_gross%*%score_vec
      
      optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt.noRE,y=y,X=Z_alles,fixef=Delta_start[1:(lin+dim.smooth)],
                                         Grad=score_vec,family=family,lower = 0, upper = Inf))
      
      t_opt<-optim.obj$par
      
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
          
          logLik.vec[1]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=NULL,family=family,penal=F)
          
          active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,dim.smooth))
          Z_aktuell<-Z_alles[,active]
          lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
          
          if (control$overdispersion)
          {  

                P_akt<-c(rep(0,lin_akt),penal.vec)
                F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)

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
        
        if(control$overdispersion)
        {
          FinalHat<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))#E-Uu
          phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(FinalHat)))
          Sigma<-Sigma*phi
        }
        
        Eta.old<-Eta
        
        vorz<-F
        
            P1<-c(rep(0,lin),penal.vec)
            P1<-diag(P1)

        score_vec2<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-P1%*%Delta[1,]
        lambda.max<-max(abs(score_vec2[(q+1):lin]))
        
        score_old2<-score_vec2
        
        if (BLOCK)
        {
          grad.1<-gradient.lasso.block(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
        }else{
          grad.1<-gradient.lasso(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
        }
        score_vec2<-c(score_vec2[1:q],grad.1,score_vec2[(lin+1):(lin+dim.smooth)])
        
        F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
        
        crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin])
        t_edge<-crit.obj$min.rate
        
        grad.2<-t(score_vec2)%*%F_gross%*%score_vec2
        
        optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt.noRE,y=y,X=Z_alles,fixef=Delta[1,1:(lin+dim.smooth)],
                                           Grad=score_vec2,family=family,lower = 0, upper = Inf))
        
        t_opt<-optim.obj$par
        
        tryNR<- (t_opt<t_edge) #&& !(all(active_old==active)  && !NRstep)
        
        if(tryNR) 
        {
          lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))

              P_akt<-c(rep(0,lin_akt),penal.vec)
              F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)

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

      ###############################################################################################################################################
      ################################################################### Main Iteration ###################################################################
      if(control$steps!=1)
      {
        for (l in 2:control$steps)
        {
          
          
          if(control$print.iter)
            message("Iteration ",l)
          #print(paste("Iteration ", l,sep=""))
          
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
                half.index<-Inf;
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
              
              logLik.vec[l]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=NULL,family=family,penal=F)
              
              active_old<-active
              active<-c(rep(T,q),!is.element(Delta[l,(q+1):lin],0),rep(T,dim.smooth))
              Z_aktuell<-Z_alles[,active]
              lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))
              
              if (control$overdispersion)
              {  

                    P_akt<-c(rep(0,lin_akt),penal.vec)
                    F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)

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
            
                P1<-c(rep(0,lin),penal.vec)
                P1<-diag(P1)
                           
            score_vec2<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-P1%*%Delta[l,]
            lambda.max<-max(abs(score_vec2[(q+1):lin]))
            
            score_old2<-score_vec2
            
            if (BLOCK)
            {
              grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda,block=block)
            }else{
              grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda)
            }
            score_vec2<-c(score_vec2[1:q],grad.1,score_vec2[(lin+1):(lin+dim.smooth)])
            
            F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
            
            crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin])
            t_edge<-crit.obj$min.rate
            
            grad.2<-t(score_vec2)%*%F_gross%*%score_vec2
            
            optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt.noRE,y=y,X=Z_alles,fixef=Delta[l,1:(lin+dim.smooth)],
                                               Grad=score_vec2,family=family,lower = 0, upper = Inf))
            
            t_opt<-optim.obj$par
            
            tryNR<- (t_opt<t_edge) && !(all(active_old==active)  && !NRstep)
            
            if(tryNR) 
            {
              lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))

                  P_akt<-c(rep(0,lin_akt),penal.vec)
                  F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)

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

          Eta.ma[l+1,]<-Eta
          
          finish<-(sqrt(sum((Eta.ma[l,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l,])^2))<control$epsilon)
          finish2<-(sqrt(sum((Eta.ma[l-1,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l-1,])^2))<control$epsilon)
          if(finish ||  finish2) #|| (all(grad.1 == 0) ))
            break
          Eta.old<-Eta
        }}
      
      
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

      aaa<-!is.element(Delta_neu[1:(lin)],0)
      
      
      if(final.re)
      {    
        ############ final re-estimation
          glmm_fin<-try(glmm_final_smooth_noRE(y,Z_fastalles[,aaa],Phi,penal.vec,
                                          Delta_start=Delta_neu[aaa],steps=control$maxIter,
                                          family=family,overdispersion=control$overdispersion,
                                          phi=phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final),silent = TRUE)
          if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
          {  
            glmm_fin2<-try(glmm_final_smooth_noRE(y,Z_fastalles[,aaa],Phi,penal.vec,
                                             Delta_start=Delta_start[aaa],steps=control$maxIter,
                                             family=family,overdispersion=control$overdispersion,
                                             phi=control$phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final),silent = TRUE)
            if(class(glmm_fin2)!="try-error")
            {    
              if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
                glmm_fin<-glmm_fin2 
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
        complexity.smooth<-sum(diag(EDF.matrix)[c(T,rep(F,sum(aaa)-1),rep(T,dim.smooth))])  
        Delta_neu2[c(aaa,rep(T,dim.smooth))]<-glmm_fin$Delta
        Standard_errors[c(aaa,rep(T,dim.smooth))]<-glmm_fin$Standard_errors
        phi<-glmm_fin$phi
        complexity<-glmm_fin$complexity
      }else{
        glmm_fin<-list()
        P_akt<-c(rep(0,lin_akt),penal.vec)
        F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
        InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
        if(class(InvFisher2)=="try-error")
          InvFisher2<-try(solve(F_gross),silent=T)
        if(class(InvFisher2)!="try-error")
        {  
          FinalHat.df<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))
          df<-sum(diag(FinalHat.df))
        }else{
          df <- NA
        }
        complexity<-df
        
        lin_akt<-q+sum(!is.element(Delta_neu[(q+1):lin],0))

            P1<-c(rep(0,lin_akt),penal.vec)
            P1a<-c(rep(0,lin_akt+dim.smooth))
            F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P1)

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
          complexity.smooth<-sum(diag(EDF.matrix)[c(T,rep(F,sum(aaa)-1),rep(T,dim.smooth))])
        }
      }  
      
      if(!(complexity.smooth>=1 && complexity.smooth<=dim.smooth))
        complexity.smooth<-dim.smooth
      
      Delta_neu<-Delta_neu2
      
      Eta_opt<-Z_alles%*%Delta_neu
      Mu_opt<-as.vector(family$linkinv(Eta_opt))

      names(Delta_neu)[1:dim(X)[2]]<-colnames(X)
      names(Standard_errors)[1:dim(X)[2]]<-colnames(X)
      
      if(lin>1)
      {
        names(Delta_neu)[(dim(X)[2]+1):lin]<-colnames(U)
        names(Standard_errors)[(dim(X)[2]+1):lin]<-colnames(U)
      }
      
      names(Delta_neu)[(lin+1):(lin+dim.smooth)]<-colnames(Phi)
      names(Standard_errors)[(lin+1):(lin+dim.smooth)]<-colnames(Phi)
      colnames(Delta)<-c(old.names,colnames(Phi))
      
      aic<-NaN
      bic<-NaN
      
      if (is.element(family$family,c("gaussian", "binomial", "poisson"))) 
      {
        loglik<-logLik.glmmLasso(y=y,mu=Mu_opt,ranef.logLik=NULL,family=family,penal=F)
        
        if(control$complexity!="hat.matrix")  
        {  
          complexity<-sum(Delta_neu[1:(lin)]!=0)+complexity.smooth
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
      ret.obj$coefficients<-Delta_neu[1:(lin)]
      ret.obj$fixerror<-Standard_errors[1:(lin)]
      ret.obj$y_hat<-Mu_opt
      ret.obj$phi<-phi
      ret.obj$family<-family
      ret.obj$fix<-fix.old
      ret.obj$data<-data
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
      ret.obj$logLik.vec<-logLik.vec
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
      
       if(lin>1)
      {
        Eta_start<-X%*%beta_null
      }else{
        Eta_start<-rep(beta_null,N)
      }
      
      D<-as.vector(family$mu.eta(Eta_start))
      Mu<-as.vector(family$linkinv(Eta_start))
      Sigma<-as.vector(family$variance(Mu))
      
        lin0<-sum(beta_null!=0)

      U<-X[,!is.na(index.new),drop=FALSE]
      X<-X[,is.na(index.new),drop=FALSE]
      
      final.names<-c(colnames(X),colnames(U))
      
      q<-dim(X)[2]
      
      Z_alles<-cbind(X,U)
      
      ########################################################## some definitions ################################################
      Delta<-matrix(0,control$steps,lin)
      Delta[1,1:lin]<-beta_null[final.names]

      Eta.ma<-matrix(0,control$steps+1,N)
      Eta.ma[1,]<-Eta_start
      
      logLik.vec<-c()

      control$epsilon<-control$epsilon*sqrt(dim(Delta)[2])
      
      Delta_start<-Delta[1,]
      
      active_old<-!is.element(Delta[1,],0)
      
      ## start value for overdispersion
      if(control$overdispersion)
      {
        if(!is.null(control$phi_start))
        {  
          phi<-control$phi_start
        }else{
          active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0))
          Z_aktuell<-Z_alles[,active]
          lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
          
              F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)

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
      
      l=1
      
      score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)
      lambda.max<-max(abs(score_vec[(q+1):lin]))
      
      #browser()
      if (BLOCK)
      {
        grad.1<-gradient.lasso.block(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
      }else{
        grad.1<-gradient.lasso(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
      }
      
      score_vec<-c(score_vec[1:q],grad.1)
      
      F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)
      
      crit.obj<-t.change(grad=score_vec[(q+1):lin],b=Delta[1,(q+1):lin])
      t_edge<-crit.obj$min.rate
      
      grad.2<-t(score_vec)%*%F_gross%*%score_vec
      
      optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt.noRE,y=y,X=Z_alles,fixef=Delta_start[1:lin],
                                         Grad=score_vec,family=family,lower = 0, upper = Inf))
      
      t_opt<-optim.obj$par
      
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
        
        logLik.vec[1]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=NULL,family=family,penal=F)
        
        active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0))
        Z_aktuell<-Z_alles[,active]
        lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
        
        if (control$overdispersion)
        {  
        F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)

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

      if(control$overdispersion)
      {    
        FinalHat<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))#E-Uu
        phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(FinalHat)))
        Sigma<-Sigma*phi
      }
      
      Eta.old<-Eta

      score_vec2<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)
      lambda.max<-max(abs(score_vec2[(q+1):lin]))
      
      if (BLOCK)
      {
        grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
      }else{
        grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
      }
      score_vec2<-c(score_vec2[1:q],grad.1)
      
      F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)
      
      crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin])
      t_edge<-crit.obj$min.rate
      
      grad.2<-t(score_vec2)%*%F_gross%*%score_vec2
      
      optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt.noRE,y=y,X=Z_alles,fixef=Delta[1,1:lin],
                                         Grad=score_vec2,family=family,lower = 0, upper = Inf))
      
      t_opt<-optim.obj$par
      
      Eta.ma[2,]<-Eta
      
      score_vec<-score_vec2

      ###############################################################################################################################################
      ################################################################### Main Iteration ###################################################################
      if(control$steps!=1)
      {
        for (l in 2:control$steps)
        {
          #browser() 
          if(control$print.iter)
            message("Iteration ",l)
          #print(paste("Iteration ", l,sep=""))
          
          half.index<-0
          
          solve.test2<-FALSE  
          while(!solve.test2)
          {  
            
            if(half.index>50)
            {
              half.index<-Inf;
            }
            
            Delta[l,]<-Delta[l-1,]+min(t_opt,t_edge)*nue*(0.5^half.index)*score_vec
            if(t_opt>t_edge & half.index==0)
              Delta[l,crit.obj$whichmin+q]<-0  
            
            Eta<-Z_alles%*%Delta[l,]
            Mu<-as.vector(family$linkinv(Eta))
            Sigma<-as.vector(family$variance(Mu))
            D<-as.vector(family$mu.eta(Eta))
            
            logLik.vec[l]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=NULL,family=family,penal=F)
            
            active_old<-active
            active<-c(rep(T,q),!is.element(Delta[l,(q+1):lin],0))
            Z_aktuell<-Z_alles[,active]
            lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))
            
            if (control$overdispersion)
            {  
                  F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)

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
          
          score_vec.unpen<-score_vec2<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)
          lambda.max<-max(abs(score_vec2[(q+1):lin]))
          
          if (BLOCK)
          {
            grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda,block=block)
          }else{
            grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda)
          }
          
          score_vec2<-c(score_vec2[1:q],grad.1)
          
          F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)
          
          crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin])
          t_edge<-crit.obj$min.rate
          
          grad.2<-t(score_vec2)%*%F_gross%*%score_vec2
          

          optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt.noRE,y=y,X=Z_alles,fixef=Delta[l,1:lin],
                                             Grad=score_vec2,family=family,lower = 0, upper = Inf))
          
          t_opt<-optim.obj$par
          
          score_vec<-score_vec2

          Eta.ma[l+1,]<-Eta
          finish<-(sqrt(sum((Eta.ma[l,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l,])^2))<control$epsilon)
          finish2<-(sqrt(sum((Eta.ma[l-1,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l-1,])^2))<control$epsilon)
          if(finish ||  finish2) #|| (all(grad.1 == 0) ))
            break
          Eta.old<-Eta
        }}
      
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

      if(final.re)
      {    
        ############ final re-estimation
        aaa<-!is.element(Delta_neu[1:(lin)],0)
        

        glmm_fin<-try(glmm_final_noRE(y,Z_fastalles[,aaa],
                                      Delta_start=Delta_neu[aaa],steps=control$maxIter,
                                      family=family,overdispersion=control$overdispersion,
                                      phi=phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final),silent = TRUE)
 

        if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
        {  
          glmm_fin2<-try(glmm_final_noRE(y,Z_fastalles[,aaa],
                                         Delta_start=Delta_start[aaa],steps=control$maxIter,
                                         family=family,overdispersion=control$overdispersion,
                                         phi=control$phi,print.iter.final=control$print.iter.final,
                                         eps.final=control$eps.final),silent = TRUE)
            if(class(glmm_fin2)!="try-error")
            {    
              if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
                glmm_fin<-glmm_fin2 
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
        Delta_neu2[aaa]<-glmm_fin$Delta
        Standard_errors[aaa]<-glmm_fin$Standard_errors
        phi<-glmm_fin$phi
        complexity<-glmm_fin$complexity
      }else{
        glmm_fin<-list()
        F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)
        InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
        if(class(InvFisher2)=="try-error")
          InvFisher2<-try(solve(F_gross),silent=T)
        if(class(InvFisher2)!="try-error")
        {  
          FinalHat.df<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))
          df<-sum(diag(FinalHat.df))
        }else{
          df <- NA
        }
        complexity<-df
      }
      
      Delta_neu<-Delta_neu2
      
      Eta_opt<-Z_alles%*%Delta_neu
      Mu_opt<-as.vector(family$linkinv(Eta_opt))
      
      names(Delta_neu)[1:dim(X)[2]]<-colnames(X)
      names(Standard_errors)[1:dim(X)[2]]<-colnames(X)
      
      if(lin>1)
      {
        names(Delta_neu)[(dim(X)[2]+1):lin]<-colnames(U)
        names(Standard_errors)[(dim(X)[2]+1):lin]<-colnames(U)
      }
      
      colnames(Delta)<-final.names
      
      aic<-NaN
      bic<-NaN
      
      if (is.element(family$family,c("gaussian", "binomial", "poisson"))) 
      {
        
        loglik<-logLik.glmmLasso(y=y,mu=Mu_opt,ranef.logLik=NULL,family=family,penal=F)
        

        if(control$complexity!="hat.matrix")  
          complexity<-sum(Delta_neu[1:(lin)]!=0)

        aic<--2*loglik+2*complexity
        bic<--2*loglik+log(N)*complexity
      }else{
        warning("For the specified family (so far) no AIC and BIC are available!")  
      }
      
      ret.obj=list()
      ret.obj$aic<-aic
      ret.obj$bic<-bic
      ret.obj$Deltamatrix<-Delta
      ret.obj$coefficients<-Delta_neu[1:(lin)]
      ret.obj$fixerror<-Standard_errors[1:(lin)]
      ret.obj$y_hat<-Mu_opt
      ret.obj$phi<-phi
      ret.obj$family<-family
      ret.obj$fix<-fix.old
      ret.obj$conv.step<-conv.step
      ret.obj$data<-data
      ret.obj$phi.med<-phi.med
      ret.obj$y <- y
      ret.obj$df<-df
      ret.obj$loglik<-loglik
      ret.obj$lambda.max<-lambda.max
      ret.obj$logLik.vec<-logLik.vec
      ret.obj$score_vec.unpen<- score_vec.unpen
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
        Eta_start<-X%*%beta_null+Phi%*%smooth_null
      }else{
        Eta_start<-rep(beta_null,N)+Phi%*%smooth_null
      }
      
      D<-as.vector(family$mu.eta(Eta_start))
      Mu<-as.vector(family$linkinv(Eta_start))
      Sigma<-as.vector(family$variance(Mu))
      
      U<-X[,!is.na(index.new),drop=FALSE]
      X<-X[,is.na(index.new),drop=FALSE]
      
      final.names<-c(colnames(X),colnames(U))
      
      q<-dim(X)[2]
      
      Z_alles<-cbind(X,U,Phi)
      ########################################################## some definitions ################################################
      Delta<-matrix(0,control$steps,(lin+dim.smooth))
      Delta[1,1:lin]<-beta_null[final.names]
      Delta[1,(lin+1):(lin+dim.smooth)]<-smooth_null

      control$epsilon<-control$epsilon*sqrt(dim(Delta)[2])
      
      Delta_start<-Delta[1,]
      
      active_old<-!is.element(Delta[1,],0)
      logLik.vec<-c()
      
      Eta.ma<-matrix(0,control$steps+1,N)
      Eta.ma[1,]<-Eta_start
      
      ## start value for overdispersion
      if(control$overdispersion)
      {
        if(!is.null(control$phi_start))
        {  
          phi<-control$phi_start
        }else{
          active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0))
          Z_aktuell<-Z_alles[,active]
          lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
          
             P_akt<-c(rep(0,lin_akt+dim.smooth))
              F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)

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

      l=1
      
      if(diff.ord>1)
      {
        k2<-c(rep(0,diff.ord-1),rep(1,nbasis-diff.ord+1))
      }else{
        k2<-rep(1,nbasis)
      }
      k22<-rep(k2,m)
      penal.vec<-penal*k22
      
          P1<-c(rep(0,lin),penal.vec)
          P1<-diag(P1)

      score_vec<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-P1%*%Delta[1,]
      
      lambda.max<-max(abs(score_vec[(q+1):lin]))
      
      if (BLOCK)
      {
        grad.1<-gradient.lasso.block(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
      }else{
        grad.1<-gradient.lasso(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
      }
      
      score_vec<-c(score_vec[1:q],grad.1,score_vec[(lin+1):(lin+dim.smooth)])
      
      F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
      
      crit.obj<-t.change(grad=score_vec[(q+1):lin],b=Delta[1,(q+1):lin])
      t_edge<-crit.obj$min.rate
      
      grad.2<-t(score_vec)%*%F_gross%*%score_vec
      
      optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt.noRE,y=y,X=Z_alles,fixef=Delta_start[1:(lin+dim.smooth)],
                                         Grad=score_vec,family=family,lower = 0, upper = Inf))
      
      t_opt<-optim.obj$par
      
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
        
        logLik.vec[1]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=NULL,family=family,penal=F)
        
        active<-c(rep(T,q),!is.element(Delta[1,(q+1):lin],0),rep(T,dim.smooth))
        Z_aktuell<-Z_alles[,active]
        lin_akt<-q+sum(!is.element(Delta[1,(q+1):lin],0))
        
        if (control$overdispersion)
        {  
              P_akt<-c(rep(0,lin_akt),penal.vec)
              F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)

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
      

            if(control$overdispersion)
      {
        FinalHat<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))#E-Uu
        phi<-(sum((y-Mu)^2/Mu))/(N-sum(diag(FinalHat)))
        Sigma<-Sigma*phi
      }
      
      Eta.old<-Eta
      
      vorz<-F
      
          P1<-c(rep(0,lin),penal.vec)
          P1<-diag(P1)

      score_vec2<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-P1%*%Delta[1,]
      
      lambda.max<-max(abs(score_vec2[(q+1):lin]))
      
      if (BLOCK)
      {
        grad.1<-gradient.lasso.block(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda,block=block)
      }else{
        grad.1<-gradient.lasso(score.beta=score_vec[(q+1):lin],b=Delta[1,(q+1):lin],lambda.b=lambda)
      }
      score_vec2<-c(score_vec2[1:q],grad.1,score_vec2[(lin+1):(lin+dim.smooth)])
      
      F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
      
      crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[1,(q+1):lin])
      t_edge<-crit.obj$min.rate
      
      grad.2<-t(score_vec2)%*%F_gross%*%score_vec2
      
      optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt.noRE,y=y,X=Z_alles,fixef=Delta[1,1:(lin+dim.smooth)],
                                         Grad=score_vec2,family=family,lower = 0, upper = Inf))
      
      t_opt<-optim.obj$par
      
      Eta.ma[2,]<-Eta
      
      score_vec<-score_vec2
      ###############################################################################################################################################
      ################################################################### Main Iteration ###################################################################
      if(control$steps!=1)
      {
        for (l in 2:control$steps)
        {
          if(control$print.iter)
            message("Iteration ",l)
          #print(paste("Iteration ", l,sep=""))
          
          half.index<-0
          
          solve.test2<-FALSE  
          while(!solve.test2)
          {  
            
            if(half.index>50)
            {
              half.index<-Inf;
            }
            
            Delta[l,]<-Delta[l-1,]+min(t_opt,t_edge)*nue*(0.5^half.index)*score_vec
            if(t_opt>t_edge & half.index==0)
              Delta[l,crit.obj$whichmin+q]<-0
            
            
            Eta<-Z_alles%*%Delta[l,]
            Mu<-as.vector(family$linkinv(Eta))
            Sigma<-as.vector(family$variance(Mu))
            D<-as.vector(family$mu.eta(Eta))
            
            logLik.vec[l]<-logLik.glmmLasso(y=y,mu=Mu,ranef.logLik=NULL,family=family,penal=F)
            
            active_old<-active
            active<-c(rep(T,q),!is.element(Delta[l,(q+1):lin],0),rep(T,dim.smooth))
            Z_aktuell<-Z_alles[,active]
            lin_akt<-q+sum(!is.element(Delta[l,(q+1):lin],0))
            
            if (control$overdispersion)
            {  
                  P_akt<-c(rep(0,lin_akt),penal.vec)
                  F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)

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
                  
              P1<-c(rep(0,lin),penal.vec)
              P1<-diag(P1)

          score_vec2<-t(Z_alles)%*%((y-Mu)*D*1/Sigma)-P1%*%Delta[l,]
          score.pure<-score_vec2
          lambda.max<-max(abs(score_vec2[(q+1):lin]))
          
          
          if (BLOCK)
          {
            grad.1<-gradient.lasso.block(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda,block=block)
          }else{
            grad.1<-gradient.lasso(score.beta=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin],lambda.b=lambda)
          }
          score_vec2<-c(score_vec2[1:q],grad.1,score_vec2[(lin+1):(lin+dim.smooth)])
          
          F_gross<-t(Z_alles)%*%(Z_alles*D*1/Sigma*D)+P1
          
          crit.obj<-t.change(grad=score_vec2[(q+1):lin],b=Delta[l,(q+1):lin])
          t_edge<-crit.obj$min.rate
          
          grad.2<-t(score_vec2)%*%F_gross%*%score_vec2
          
          optim.obj<-suppressWarnings(nlminb(1e-16,taylor.opt.noRE,y=y,X=Z_alles,fixef=Delta[l,1:(lin+dim.smooth)],
                                             Grad=score_vec2,family=family,lower = 0, upper = Inf))
          
          t_opt<-optim.obj$par
          
          score_vec<-score_vec2
         
          Eta.ma[l+1,]<-Eta
          
          finish<-(sqrt(sum((Eta.ma[l,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l,])^2))<control$epsilon)
          finish2<-(sqrt(sum((Eta.ma[l-1,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l-1,])^2))<control$epsilon)
          if(finish ||  finish2) #|| (all(grad.1 == 0) ))
            break
          Eta.old<-Eta
        }}
      
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

      aaa<-!is.element(Delta_neu[1:(lin)],0)

      if(final.re)
      {    
        ############ final re-estimation
        
        glmm_fin<-try(glmm_final_smooth_noRE(y,Z_fastalles[,aaa],Phi,penal.vec,
                                             Delta_start=Delta_neu[c(aaa,rep(T,dim.smooth))],steps=control$maxIter,
                                             family=family,overdispersion=control$overdispersion,
                                             phi=phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final),silent = TRUE)
        if(class(glmm_fin)=="try-error" || glmm_fin$opt>control$maxIter-10)
        {  
          glmm_fin2<-try(glmm_final_smooth_noRE(y,Z_fastalles[,aaa],Phi,penal.vec,
                                                Delta_start=Delta_start[c(aaa,rep(T,dim.smooth))],steps=control$maxIter,
                                                family=family,overdispersion=control$overdispersion,
                                                phi=control$phi,print.iter.final=control$print.iter.final,eps.final=control$eps.final),silent = TRUE)
          if(class(glmm_fin2)!="try-error")
            {    
              if(class(glmm_fin)=="try-error" || (glmm_fin$opt>control$maxIter-10 && glmm_fin2$opt<control$maxIter-10))
                glmm_fin<-glmm_fin2 
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
        complexity.smooth<-sum(diag(EDF.matrix)[c(T,rep(F,sum(aaa)-1),rep(T,dim.smooth))])  
        Delta_neu2[c(aaa,rep(T,dim.smooth))]<-glmm_fin$Delta
        Standard_errors[c(aaa,rep(T,dim.smooth))]<-glmm_fin$Standard_errors
        phi<-glmm_fin$phi
        complexity<-glmm_fin$complexity
      }else{
        glmm_fin<-list()
        P_akt<-c(rep(0,lin_akt),penal.vec)
        F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P_akt)
        InvFisher2<-try(chol2inv(chol(F_gross)),silent=T)
        if(class(InvFisher2)=="try-error")
          InvFisher2<-try(solve(F_gross),silent=T)
        if(class(InvFisher2)!="try-error")
        {  
          FinalHat.df<-(Z_aktuell*sqrt(Sigma*D*1/Sigma*D))%*%InvFisher2%*%t(Z_aktuell*sqrt(D*1/Sigma*D*1/Sigma))
          df<-sum(diag(FinalHat.df))
        }else{
          df <- NA
        }
        complexity<-df
        
        lin_akt<-q+sum(!is.element(Delta_neu[(q+1):lin],0))

            P1<-c(rep(0,lin_akt),penal.vec)
            P1a<-c(rep(0,lin_akt+dim.smooth))
            F_gross<-t(Z_aktuell)%*%(Z_aktuell*D*1/Sigma*D)+diag(P1)

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
          complexity.smooth<-sum(diag(EDF.matrix)[c(T,rep(F,sum(aaa)-1),rep(T,dim.smooth))])
        }
      }  
      
      if(!(complexity.smooth>=1 && complexity.smooth<=dim.smooth))
        complexity.smooth<-dim.smooth
      
      Delta_neu<-Delta_neu2
      
      Eta_opt<-Z_alles%*%Delta_neu
      Mu_opt<-as.vector(family$linkinv(Eta_opt))
      
      names(Delta_neu)[1:dim(X)[2]]<-colnames(X)
      names(Standard_errors)[1:dim(X)[2]]<-colnames(X)
      
      if(lin>1)
      {
        names(Delta_neu)[(dim(X)[2]+1):lin]<-colnames(U)
        names(Standard_errors)[(dim(X)[2]+1):lin]<-colnames(U)
      }
      
      names(Delta_neu)[(lin+1):(lin+dim.smooth)]<-colnames(Phi)
      names(Standard_errors)[(lin+1):(lin+dim.smooth)]<-colnames(Phi)
      colnames(Delta)<-c(old.names,colnames(Phi))
      
      aic<-NaN
      bic<-NaN
      
      
      if (is.element(family$family,c("gaussian", "binomial", "poisson"))) 
      {
        loglik<-logLik.glmmLasso(y=y,mu=Mu_opt,ranef.logLik=NULL,family=family,penal=F)
        
        if(control$complexity!="hat.matrix")  
        {  
           complexity<-sum(Delta_neu[1:(lin)]!=0)+complexity.smooth
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
      ret.obj$coefficients<-Delta_neu[1:(lin)]
      ret.obj$fixerror<-Standard_errors[1:(lin)]
      ret.obj$y_hat<-Mu_opt
      ret.obj$phi<-phi
      ret.obj$family<-family
      ret.obj$fix<-fix.old
      ret.obj$data<-data
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
      ret.obj$logLik.vec<-logLik.vec
      return(ret.obj)
    }  
  }
  ##################################################################
  ##################################################################
}