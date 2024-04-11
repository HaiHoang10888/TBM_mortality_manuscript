simp.reg.boot<-function(x,y,b,k,alpha)
{
n<-length(x)
a1<-data.frame(x,y)
a<-data.matrix(a1)

inta<-n*sum(x*y)
intb<-sum(x)*sum(y)
intc<-n*sum(x^2)
intd<-(sum(x))^2
inte<-inta-intb
intf<-intc-intd
beta1<-inte/intf
ybar<-sum(y)/n
xbar<-sum(x)/n
beta0<-ybar-(beta1*xbar)

cat("Least Square Method","\n")
cat("beta zero =",beta0,"\n")
cat("beta one =",beta1,"\n")

yhat<-beta0+(beta1*x)
mse<-sum((y-yhat)^2)/(n-2)
varbeta1<-mse/sum((x-xbar)^2)
sebeta1<-sqrt(varbeta1)
sxx<-sum((x-xbar)^2)
varbeta0<-mse*(sum(x^2)/(n*sxx))
sebeta0<-sqrt(varbeta0)

cat("Classical Method","\n")
cat("The accuracy of beta zero","\n")
cat("bias of beta zero =",0,"\n")
cat("standard error of beta zero =",sebeta0,"\n")
cat(((1-alpha)*100),"% confidence interval for beta zero","\n")
lslobound0=beta0-(qt(1-(alpha/2),(n-2))*sebeta0)
lsupbound0=beta0+(qt(1-(alpha/2),(n-2))*sebeta0)

cat("lower bound =",lslobound0,"\n")
cat("upper bound =",lsupbound0,"\n")

cat("\n")

cat("The accuracy of beta one","\n")
cat("bias of beta one =",0,"\n")
cat("standard error of beta one =",sebeta1,"\n")
cat(((1-alpha)*100),"% confidence interval for beta one","\n")
lslobound1=beta1-(qt(1-(alpha/2),(n-2))*sebeta1)
lsupbound1=beta1+(qt(1-(alpha/2),(n-2))*sebeta1)

cat("lower bound =",lslobound1,"\n")
cat("upper bound =",lsupbound1,"\n")
cat("\n")


bootbeta0=NULL
bootbeta1=NULL
for(i in 1:b)
{
v<-sample(1:n,k,replace=TRUE)
indep<-a[v]
dep<-a[v,2]
int1<-k*sum(indep*dep)
int2<-sum(indep)*sum(dep)
int3<-k*sum(indep^2)
int4<-(sum(indep))^2
int5<-int1-int2
int6<-int3-int4
betaboot1<-int5/int6
depbar<-sum(dep)/k
indepbar<-sum(indep)/k
betaboot0<-depbar-(betaboot1*indepbar)
{
bootbeta0[i]=betaboot0
bootbeta1[i]=betaboot1
}
}

cat("\n")
cat("Bootstrap for Regression, Correlation Model","\n")

meanbetaboot0=mean(bootbeta0)
cat("Bootstrap's beta zero :",meanbetaboot0,"\n")
cat("The Bootstrap's accuracy measures of beta zero","\n")
biasbeta0=meanbetaboot0-beta0
cat("bootstrap's bias for beta zero=",biasbeta0,"\n")

varbetaboot0=var(bootbeta0)
sebetaboot0=sqrt(varbetaboot0)
cat("bootstrap's standard error for beta zero=",sebetaboot0,"\n")

cat(((1-alpha)*100),"% confidence interval for beta zero","\n")
lobound0=quantile(bootbeta0,(alpha/2))
upbound0=quantile(bootbeta0,(1-(alpha/2)))

cat("lower bound =",lobound0,"\n")
cat("upper bound =",upbound0,"\n")

cat("\n")

meanbetaboot1=mean(bootbeta1)
cat("Bootstrap's beta one :",meanbetaboot1,"\n")
cat("The Bootstrap's accuracy measures of beta one","\n")
biasbeta1=meanbetaboot1-beta1
cat("bootstrap's bias for beta one=",biasbeta1,"\n")

varbetaboot1=var(bootbeta1)
sebetaboot1=sqrt(varbetaboot1)
cat("bootstrap's standard error for beta one=",sebetaboot1,"\n")

cat(((1-alpha)*100),"% confidence interval for beta one","\n")
lobound1=quantile(bootbeta1,(alpha/2))
upbound1=quantile(bootbeta1,(1-(alpha/2)))

cat("lower bound =",lobound1,"\n")
cat("upper bound =",upbound1,"\n")

par(mfrow=c(1,2))
hist(bootbeta0)
hist(bootbeta1)
}
