glmmLassoControl<-function(nue=1,lin="(Intercept)",start=NULL,q_start=NULL, steps=2000,method="EM",overdispersion=FALSE,epsilon=1e-7,maxIter=1000,print.iter=FALSE)
{                       
if (lin[1]!="(Intercept)")
lin<-c("(Intercept)",lin)

list(nue = nue, lin = lin, start = start, q_start = q_start, steps = steps, method = method,         
        overdispersion = overdispersion, epsilon = epsilon, maxIter = maxIter,print.iter = print.iter)
}
