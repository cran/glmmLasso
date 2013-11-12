library(glmmLasso)
data("soccer")
## generalized additive mixed model
## grid for the smoothing parameter

## center all metric variables
soccer[,c(4,5,9:16)]<-scale(soccer[,c(4,5,9:16)],center=T,scale=T)
soccer<-data.frame(soccer)


lambda <- seq(0,200,by=1)

family = poisson(link = log)


################## First Simple Method ############################################
## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda

BIC_vec<-rep(Inf,length(lambda))


for(j in 1:length(lambda))
{
print(paste("Iteration ", j,sep=""))
  
glm1 <- try(glmmLasso(points~transfer.spendings  
        + ave.unfair.score + transfer.receits + ball.possession
        + tackles + ave.attend + sold.out, rnd = list(team=~1),  
        family = family, data = soccer, lambda=lambda[j],
        control = list(overdispersion=TRUE)),silent=TRUE)  

if(class(glm1)!="try-error")
{  
BIC_vec[j]<-glm1$bic
}
        
}
    
opt<-match(min(BIC_vec),BIC_vec)
        
glm1_final <- glmmLasso(points~transfer.spendings  
        + ave.unfair.score + transfer.receits + ball.possession
        + tackles + ave.attend + sold.out, rnd = list(team=~1),  
        family = family, data = soccer, lambda=lambda[opt],
        control = list(overdispersion=TRUE)) 
  
        
summary(glm1_final)





################## Second Simple Method ###########################
## Using 5-fold CV to determine the optimal tuning parameter lambda

### set seed
set.seed(0815)
N<-dim(soccer)[1]
ind<-sample(N,N)

kk<-5
nk <- floor(N/kk)

Devianz_ma<-matrix(Inf,ncol=kk,nrow=length(lambda))


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
        + ave.unfair.score + transfer.receits + ball.possession
        + tackles + ave.attend + sold.out, rnd = list(team=~1),  
        family = family, data = soccer.train, lambda=lambda[j],
        control = list(overdispersion=TRUE)),silent=TRUE) 
        
    if(class(glm2)!="try-error")
    {  
    y.hat<-predict(glm2,soccer.test)    

    Devianz_ma[j,i]<-sum(family$dev.resids(soccer.test$points,y.hat,wt=rep(1,length(y.hat))))
    }
}}
    
Devianz_vec<-apply(Devianz_ma,1,sum)
opt2<-match(min(Devianz_vec),Devianz_vec)
       
       
glm2_final <- glmmLasso(points~transfer.spendings  
        + ave.unfair.score + transfer.receits + ball.possession
        + tackles + ave.attend + sold.out, rnd = list(team=~1),  
        family = family, data = soccer, lambda=lambda[opt2],
        control = list(overdispersion=TRUE))  
    
        
summary(glm2_final)


################## More Elegant Method ############################################
## Idea: start with big lambda and use the estimates of the previous fit (BUT: before
## the final re-estimation Fisher scoring is performed!) as starting values for the next fit;
## make sure, that your lambda sequence starts at a value big enough such that all covariates are
## shrinked to zero;

## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda
lambda <- seq(200,0,by=-1)

BIC_vec<-rep(Inf,length(lambda))
family = poisson(link = log)

# specify starting values for the very first fit; pay attention that Delta.start has suitable length! 
Delta.start<-as.matrix(t(rep(0,8+23)))   
Q.start<-0.1  
phi.start<-1

for(j in 1:length(lambda))
{
  print(paste("Iteration ", j,sep=""))
  
  glm3 <- glmmLasso(points~transfer.spendings  
                    + ave.unfair.score + transfer.receits + tackles
                    + sold.out 
                    + ball.possession+ ave.attend,
                    rnd = list(team=~1),  
                    family = family, data = soccer, lambda=lambda[j],
                    control = list(overdispersion=TRUE,start=Delta.start[j,],
                    q_start=Q.start[j],phi_start=phi.start[j]))  
  
  BIC_vec[j]<-glm3$bic
  Delta.start<-rbind(Delta.start,glm3$Deltamatrix[glm3$conv.step,])
  Q.start<-c(Q.start,glm3$Q_long[[glm3$conv.step+1]])
  phi.start<-c(phi.start,glm3$phi.med)
}

opt3<-match(min(BIC_vec),BIC_vec)

glm3_final <- glmmLasso(points~transfer.spendings  
                        + ave.unfair.score + transfer.receits + ball.possession
                        + tackles + ave.attend + sold.out, rnd = list(team=~1),  
                        family = family, data = soccer, lambda=lambda[opt3],
                        control = list(overdispersion=TRUE,start=Delta.start[opt3,],
                        q_start=Q.start[opt3],phi_start=phi.start[opt3]))  


summary(glm3_final)


