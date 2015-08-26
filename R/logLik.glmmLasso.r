logLik.glmmLasso<-function(y,mu,ranef.logLik=NULL,family,penal=FALSE) 
{
  fam <- family$family

  if(fam=="poisson")
  loglik <- sum(y*log(mu)-mu)
  
  if(fam=="binomial")
  loglik <- sum(log(mu[y==1]))+sum(log((1-mu)[y==0]))
  
  if(fam=="gaussian")
  loglik <- sum(y*mu-0.5*(mu^2))
  
  if(penal)
  loglik <- loglik + ranef.logLik 

return(loglik)
}
