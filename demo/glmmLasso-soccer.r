library(glmmLasso)
data("soccer")
## generalized additive mixed model
## grid for the smoothing parameter

## center all metric variables so that also the 
## starting values with glmmPQL are in the correct scaling

soccer[,c(4,5,9:16)]<-scale(soccer[,c(4,5,9:16)],center=T,scale=T)
soccer<-data.frame(soccer)


lambda <- seq(500,0,by=-5)

family <- poisson(link = log)


################## First Simple Method ############################################
## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda

BIC_vec<-rep(Inf,length(lambda))

## first fit good starting model
library(MASS);library(nlme)
PQL<-glmmPQL(points~1,random = ~1|team,family=family,data=soccer)
Delta.start<-c(as.numeric(PQL$coef$fixed),rep(0,6),as.numeric(t(PQL$coef$random$team)))
Q.start<-as.numeric(VarCorr(PQL)[1,1])

for(j in 1:length(lambda))
{
print(paste("Iteration ", j,sep=""))
  
glm1 <- try(glmmLasso(points~transfer.spendings  
        + ave.unfair.score + ball.possession
        + tackles + ave.attend + sold.out, rnd = list(team=~1),  
        family = family, data = soccer, lambda=lambda[j],switch.NR=FALSE,final.re=FALSE,
        control=list(start=Delta.start,q_start=Q.start)), silent=TRUE)  


if(!inherits(glm1, "try-error"))
{  
BIC_vec[j]<-glm1$bic
}
        
}
    
opt<-which.min(BIC_vec)
        
glm1_final <- glmmLasso(points~transfer.spendings  
        + ave.unfair.score + ball.possession
        + tackles + ave.attend + sold.out, rnd = list(team=~1),  
        family = family, data = soccer, lambda=lambda[opt],switch.NR=FALSE,final.re=FALSE,
        control=list(start=Delta.start,q_start=Q.start))
         
  
        
summary(glm1_final)





################## Second Simple Method ###########################
## Using 5-fold CV to determine the optimal tuning parameter lambda

### set seed
set.seed(123)
N<-dim(soccer)[1]
ind<-sample(N,N)
lambda <- seq(500,0,by=-5)

## set number of folds
kk<-5
nk <- floor(N/kk)

Devianz_ma<-matrix(Inf,ncol=kk,nrow=length(lambda))

## first fit good starting model
library(MASS);library(nlme)
PQL<-glmmPQL(points~1,random = ~1|team,family=family,data=soccer)
Delta.start<-c(as.numeric(PQL$coef$fixed),rep(0,6),as.numeric(t(PQL$coef$random$team)))
Q.start<-as.numeric(VarCorr(PQL)[1,1])


for(j in 1:length(lambda))
{
print(paste("Iteration ", j,sep=""))
  
  for (i in 1:kk)
  {
    if (i < kk)
    {
    indi <- ind[(i-1)*nk+(1:nk)]
    }else{
    indi <- ind[((i-1)*nk+1):N]
    }
  
soccer.train<-soccer[-indi,]
soccer.test<-soccer[indi,]
  
glm2 <- try(glmmLasso(points~transfer.spendings  
        + ave.unfair.score  + ball.possession
        + tackles + ave.attend + sold.out, rnd = list(team=~1),  
        family = family, data = soccer.train, lambda=lambda[j],switch.NR=FALSE,final.re=FALSE,
        control=list(start=Delta.start,q_start=Q.start))
        ,silent=TRUE) 
        
    if(!inherits(glm2, "try-error"))
    {  
    y.hat<-predict(glm2,soccer.test)    

    Devianz_ma[j,i]<-sum(family$dev.resids(soccer.test$points,y.hat,wt=rep(1,length(y.hat))))
    }
}
print(sum(Devianz_ma[j,]))
}
    
Devianz_vec<-apply(Devianz_ma,1,sum)
opt2<-which.min(Devianz_vec)
       
       
glm2_final <- glmmLasso(points~transfer.spendings  
                        + ave.unfair.score + ball.possession
                        + tackles + ave.attend + sold.out, rnd = list(team=~1),  
                        family = family, data = soccer, lambda=lambda[opt2],switch.NR=FALSE,final.re=FALSE,
                        control=list(start=Delta.start,q_start=Q.start))



summary(glm2_final)


################## More Elegant Method ############################################
## Idea: start with big lambda and use the estimates of the previous fit (BUT: before
## the final re-estimation Fisher scoring is performed!) as starting values for the next fit;
## make sure, that your lambda sequence starts at a value big enough such that all covariates are
## shrunk to zero;

## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda
## (here, a log-grid is used)
lambda <- exp(seq(log(460),log(1e-4),length.out = 50))
#lambda <- seq(460,0,length.out = 51)


BIC_vec<-rep(Inf,length(lambda))
family <- poisson(link = log)

# specify starting values for the very first fit; pay attention that Delta.start has suitable length! 
Delta.start <- as.matrix(t(rep(0,7+23)))
# set reasonable starting value for intercept (depending on GLM family)
Delta.start[1,1] <- family$linkfun(mean(soccer$points))
Q.start<-0.1  

for(j in 1:length(lambda))
{
  print(paste("Iteration ", j,sep=""))
  
  glm3 <- glmmLasso(points~1 +transfer.spendings
                    + ave.unfair.score 
                    + tackles  
                    + sold.out 
                    + ball.possession
                    + ave.attend
                    ,rnd = list(team=~1),  
                    family = family, data = soccer, 
                    lambda=lambda[j], switch.NR=FALSE,final.re=FALSE,
                    control = list(start=Delta.start[j,],q_start=Q.start[j]))  
  
  print(colnames(glm3$Deltamatrix)[2:7][glm3$Deltamatrix[glm3$conv.step,2:7]!=0])
  BIC_vec[j]<-glm3$bic
  Delta.start<-rbind(Delta.start,glm3$Deltamatrix[glm3$conv.step,])
  Q.start<-c(Q.start,glm3$Q_long[[glm3$conv.step+1]])
}

opt3<-which.min(BIC_vec)

glm3_final <- glmmLasso(points~transfer.spendings + ave.unfair.score 
                        + ball.possession
                        + tackles + ave.attend + sold.out, rnd = list(team=~1),  
                        family = family, data = soccer, lambda=lambda[opt3],
                        switch.NR=FALSE,final.re=FALSE,
                        control = list(start=Delta.start[opt3,],q_start=Q.start[opt3]))  


summary(glm3_final)

## plot coefficient paths
# quartz()
par(mar=c(6,6,4,4))
plot(lambda,Delta.start[2:(length(lambda)+1),2],type="l",lwd=2,ylim=c(-0.05,0.05),ylab=expression(hat(beta[j])))
abline(h=0,lty=3,lwd=2,col="grey")
for(i in 3:7){
  lines(lambda[1:length(lambda)],Delta.start[2:(length(lambda)+1),i], lwd=2)
}
axis(3,at=lambda[opt3],label=substitute(paste(lambda["opt"]," = ",nn), list(nn=lambda[opt3])))
abline(v=lambda[opt3],lty=2,lwd=2,col=2)



################## Elegant Cross-Validation ###########################
## Using 5-fold CV to determine the optimal tuning parameter lambda
## Idea: on each training data, similar to the previous method, start 
## with big lambda and use the estimates of the previous fit (BUT: before
## the final re-estimation Fisher scoring is performed!) as starting values for the next fit;
## make sure, that your lambda sequence starts at a value big enough such that all
## covariates are shrunk to zero;


### set seed
set.seed(1909)
N<-dim(soccer)[1]
ind<-sample(N,N)
lambda <- seq(500,0,by=-5)

kk<-5
nk <- floor(N/kk)

Devianz_ma<-matrix(Inf,ncol=kk,nrow=length(lambda))

## first fit good starting model
library(MASS);library(nlme)
PQL <- glmmPQL(points~1,random = ~1|team,family=family,data=soccer)

Delta.start <- as.matrix(t(c(as.numeric(PQL$coef$fixed),rep(0,6),as.numeric(t(PQL$coef$random$team)))))
Q.start <- as.numeric(VarCorr(PQL)[1,1])


## loop over the folds  
for (i in 1:kk)
{
  print(paste("CV Loop ", i,sep=""))
  
    if (i < kk)
    {
      indi <- ind[(i-1)*nk+(1:nk)]
    }else{
      indi <- ind[((i-1)*nk+1):N]
    }
    
    soccer.train<-soccer[-indi,]
    soccer.test<-soccer[indi,]
    
    Delta.temp <- Delta.start
    Q.temp <- Q.start
    ## loop over lambda grid
    for(j in 1:length(lambda))
    {
    #print(paste("Lambda Iteration ", j,sep=""))
      
    glm4 <- try(glmmLasso(points~transfer.spendings  
                          + ave.unfair.score  + ball.possession
                          + tackles + ave.attend + sold.out, rnd = list(team=~1),  
                          family = family, data = soccer.train, lambda=lambda[j],switch.NR=FALSE,final.re=FALSE,
                          control=list(start=Delta.temp[j,],q_start=Q.temp[j]))
                ,silent=TRUE) 
    
    if(!inherits(glm4, "try-error"))
    {  
      y.hat<-predict(glm4,soccer.test)    
      Delta.temp<-rbind(Delta.temp,glm4$Deltamatrix[glm4$conv.step,])
      Q.temp<-c(Q.temp,glm4$Q_long[[glm4$conv.step+1]])
      
      Devianz_ma[j,i]<-sum(family$dev.resids(soccer.test$points,y.hat,wt=rep(1,length(y.hat))))
    }
  }
}

Devianz_vec<-apply(Devianz_ma,1,sum)
opt4<-which.min(Devianz_vec)

## now fit full model until optimnal lambda (which is at opt4)
for(j in 1:opt4)
{  
glm4.big <- glmmLasso(points~transfer.spendings  
                        + ave.unfair.score + ball.possession
                        + tackles + ave.attend + sold.out, rnd = list(team=~1),  
                        family = family, data = soccer, lambda=lambda[j],switch.NR=FALSE,final.re=FALSE,
                        control=list(start=Delta.start[j,],q_start=Q.start[j]))
Delta.start<-rbind(Delta.start,glm4.big$Deltamatrix[glm4.big$conv.step,])
Q.start<-c(Q.start,glm4.big$Q_long[[glm4.big$conv.step+1]])
}

glm4_final <- glm4.big

summary(glm4_final)
