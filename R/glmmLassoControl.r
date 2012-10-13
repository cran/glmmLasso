glmmLassoControl<-function(nue=1,lin="(Intercept)",start=NULL,q_start=NULL, phi_start = 1, steps=2000,method="EM", overdispersion=FALSE,epsilon=1e-7,maxIter=1000,print.iter=FALSE,print.iter.final=FALSE,method.final="REML")
{                       
if (lin[1]!="(Intercept)")
lin<-c("(Intercept)",lin)

list(nue = nue, lin = lin, start = start, q_start = q_start, phi_start = phi_start, steps = steps, method = method,         
        overdispersion = overdispersion, epsilon = epsilon, maxIter = maxIter,print.iter = print.iter, print.iter.final=print.iter.final, method.final = method.final)
}
