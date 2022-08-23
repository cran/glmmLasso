cumulative <- function()
{
family <- "cumulative"


linkinv <- function(eta_cat,K=NULL){
  eta <- c(eta_cat)
  linkin <- binomial()$linkinv(eta)
}
  
 
  #createSigmaInv <- function(mu){
    #Sigma <- mu%*%t(1-mu)
    #Sigma[lower.tri(Sigma)]<-t(Sigma)[lower.tri(Sigma)]
    #SigmaInv <- try(solve(Sigma),silent=T); k=1
    #while(class(SigmaInv)=="try-error")
    #{  
    #Sigma.new <- Sigma + diag(10^(-7+k),length(mu))
    #SigmaInv <- try(solve(Sigma.new),silent=T);k<-k+1
    #}
    #return(SigmaInv)
  #}
  
mulist <- function(mu, K){
  split(mu, as.integer((seq_along(mu) -1) / K))
}
  
  SigmaInv <- function(mu, K){
    #sig.list <- lapply(cumulative()$mulist(mu),cumulative()$createSigmaInv)
    #SigmaInv <- bdiag_m(sig.list)   
    sig.list <- lapply(cumulative()$mulist(mu, K),RcppEigenSigmaInv)
    SigmaInv <- bdiag_m(sig.list)
    }
  
  deriv.mat <- function(Eta,K=NULL){
    eta <- c(Eta)
    RcppEigenDiagSp(binomial()$mu.eta(eta)) 
    }
  
  multivariate <- TRUE
  
  ret.list <- list(linkinv = linkinv, SigmaInv = SigmaInv,
              ##createSigmaInv=#createSigmaInv,     
              deriv.mat = deriv.mat,
                   mulist = mulist, multivariate = multivariate, family = family)
  return(ret.list)

}
  

