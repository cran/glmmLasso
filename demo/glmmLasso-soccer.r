library(glmmLasso)
data("soccer")
## generalized additive mixed model
## grid for the smoothing parameter

## center all metric variables
soccer[,c(4,5,9:16)]<-scale(soccer[,c(4,5,9:16)],center=T,scale=T)
soccer<-data.frame(soccer)


lambda <- seq(0,1,by=0.1)

### set seed
set.seed(0815)

N<-dim(soccer)[1]
ind<-sample(N,N)


## 10-fold CV

kk<-5
nk <- floor(N/kk)

Devianz_ma<-matrix(Inf,ncol=kk,nrow=length(lambda))

family = poisson(link = log)

for(j in 1:length(lambda))
{
print("j")
print(j)

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
  
glm1 <- glmmLasso(points~transfer.spendings  + I(transfer.spendings^2)
        + ave.unfair.score + transfer.receits + ball.possession
        + tackles + ave.attend + sold.out, rnd = list(team=~1),  
        family = family, data = soccer.train, lambda=lambda[j],
        control = list(start=c(1,rep(0,30)),overdispersion=TRUE))
        
y.hat<-predict(glm1,soccer.test)    

Devianz_ma[j,i]<-sum(family$dev.resids(soccer.test$points,y.hat,wt=rep(1,length(y.hat))))
}}
    
Devianz_vec<-apply(Devianz_ma,1,sum)
opt<-match(min(Devianz_vec),Devianz_vec)
        
glm1_final <- glmmLasso(points~transfer.spendings  + 
        + as.factor(red.card) + transfer.receits + ball.possession
        + tackles + ave.attend + sold.out, rnd = list(team=~1),  
        family = family, data = soccer.train, lambda=lambda[opt],
        control = list(start=c(1,rep(0,30)),overdispersion=TRUE))
    
        
summary(glm1_final)
