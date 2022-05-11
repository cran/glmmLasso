###############
# preliminaries
library(glmmLasso)
data("soccer")

soccer[,c(4,5,9:16)]<-scale(soccer[,c(4,5,9:16)],center=TRUE,scale=TRUE)
soccer<-data.frame(soccer)

###############
## linear mixed model
lm1 <- glmmLasso(points ~ transfer.spendings + ave.unfair.score 
                 + ball.possession + tackles 
                 + ave.attend + sold.out, rnd = list(team=~1), 
                 lambda=50, data = soccer)

summary(lm1)

## similar linear model without random effects
lm1b <- glmmLasso(points ~ transfer.spendings + ave.unfair.score 
                  + ball.possession + tackles 
                  + ave.attend + sold.out, rnd = NULL, 
                  lambda=50, data = soccer)

summary(lm1b)


## linear mixed model with slope on ave.attend;  
## the coefficient of ave.attend is not penalized;
lm2 <- glmmLasso(points~transfer.spendings + ave.unfair.score 
                 + ball.possession + tackles + ave.attend 
                 + sold.out, rnd = list(team=~1 + ave.attend), lambda=10, 
                 data = soccer, control = list(index=c(1,2,3,4,NA,5), 
                                               method="REML",print.iter=TRUE))

summary(lm2)

## linear mixed model with categorical covariates
## and final Fisher scoring re-estimation step
lm3 <- glmmLasso(points ~ transfer.spendings + as.factor(red.card)  
                 + as.factor(yellow.red.card) + ball.possession 
                 + tackles + ave.attend + sold.out, rnd = list(team=~1), 
                 data = soccer, lambda=100, final.re=TRUE,
                 control = list(print.iter=TRUE,print.iter.final=TRUE))

summary(lm3)

## generalized linear mixed model
## with starting values
glm1 <- glmmLasso(points~transfer.spendings  
                  + ave.unfair.score + sold.out 
                  + tackles + ave.attend + ball.possession, rnd = list(team=~1),  
                  family = poisson(link = log), data = soccer, lambda=100, 
                  control = list(print.iter=TRUE,start=c(1,rep(0,29)),q_start=0.7)) 

