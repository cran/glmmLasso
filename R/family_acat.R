acat <- function()
{
family <- "acat"

  responseHelp <- function(eta, K){
    eta.help <- matrix(rep(c(0,eta),each=K+1),ncol=K+1)
    eta.help[upper.tri(eta.help)] <- 0
    pi.temp <- cumprod(c(1,exp(eta[-K])))/sum(apply(exp(eta.help),1,prod))
    pi.temp
  }
  
  responseFun <- function(eta, K){
    eta.temp <- matrix(eta, byrow = TRUE, ncol = K)
    pi <- c(apply(eta.temp,1,responseHelp, K=K))
    pi
  } 
    
  linkinv <- function(eta,K){
    linkin <- acat()$responseFun(eta,K)
  }
  
  createSigmaInv <- function(mu){
    Sigma <- diag(mu) - mu %*% t(mu)
    RcppEigenInvMa(Sigma)
  }
  
  mulist <- function(mu,K){
    mu.temp <- matrix(mu,ncol=K)
    mu.list <- split(mu.temp, rep(1:nrow(mu.temp), ncol(mu.temp)))
    mu.list
  }
  
  SigmaInv <- function(mu,K){
    SigmaInv <- as(as(as(as.matrix(bdiag(lapply(acat()$mulist(mu,K),acat()$createSigmaInv))), "dMatrix"), "generalMatrix"), "CsparseMatrix")
    SigmaInv
  }
  
  createD <- function(mu, K){
#    browser()
    D2 <- matrix(0,K,K)
    diag(D2) <- -(1/mu)
    
    if(K==2){
      D2[2,1] <- 1/mu[-1]
    }else{
      diag(D2[2:K,1:(K-1)]) <- 1/mu[-1]
    }
    
    D2[,K] <- -1/(1-sum(mu))
    D2[K,K] <- -(1-sum(mu[-K]))/((1-sum(mu))*mu[K])
    
    D <- solve(D2)
    D
  }
  
  deriv.mat <- function(eta,K){
    mu <- linkinv(eta, K = K)
    d.temp <- as(as(as(as.matrix(bdiag(lapply(acat()$mulist(mu,K),acat()$createD, K = K))), "dMatrix"), "generalMatrix"), "CsparseMatrix")
    d.temp
  }
  
  multivariate <- TRUE
  
  ret.list <- list(responseFun = responseFun, linkinv = linkinv, SigmaInv = SigmaInv,
                   createSigmaInv = createSigmaInv, createD = createD, deriv.mat = deriv.mat,
                   mulist = mulist,
                   multivariate = multivariate, family = family)
  return(ret.list)

}
  
