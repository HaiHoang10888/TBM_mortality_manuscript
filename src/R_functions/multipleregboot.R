multiple.reg.boot<-function(x,y,p,b,k,alpha)
{
n<-length(y)
int1<-t(x)
int2<-int1%*%x
int3<-solve(int2)
int4<-int1%*%y
koef.beta<-int3%*%int4
index.beta<-matrix(0:(p-1))
out<-cbind(index.beta,koef.beta)
cat("Least Square Estimator","\n")
print(out)

cat("\n")
cat("Standard Error for Beta","\n")
int5<-t(y)
int6<-t(koef.beta)
int7<-int5%*%y
int8<-int6%*%int1
int9<-int8%*%y
sse<-int7-int9
mse<-sse/(n-p)
msedata<-mse[1]
var.all<-msedata*int3
var.koef.beta<-diag(var.all)
se.koef.beta<-sqrt(var.koef.beta)
print(se.koef.beta)
cat("\n")

cat(((1-(alpha))*100),"% Confidence Interval for Beta","\n")
low<-koef.beta-qt((1-(alpha/2)),(n-p))*se.koef.beta
up<-koef.beta+qt((1-(alpha/2)),(n-p))*se.koef.beta
cat("Lower Bound","\n")
print(low)
cat("\n")
cat("Upper Bound","\n")
print(up)
cat("\n")


matcoef<-matrix(0:(p-1))
beta.boot=NULL
for(i in 1:b)
{
v<-sample(1:n,k,replace=TRUE)
ystarpre<-y[v]
xstarpre<-x[v,1:p]
ystar<-data.matrix(ystarpre)
xstar<-data.matrix(xstarpre)
inta<-t(xstar)
intb<-inta%*%xstar
detb<-det(intb)

if(detb==0)
{
repeat
{
v<-sample(1:n,k,replace=TRUE)
ystarpre<-y[v]
xstarpre<-x[v,1:p]
ystar<-data.matrix(ystarpre)
xstar<-data.matrix(xstarpre)
inta<-t(xstar)
intb<-inta%*%xstar
detb<-det(intb)
if(detb!=0) next
}
}

intc<-solve(intb)
intd<-inta%*%ystar
koef.beta.boot<-intc%*%intd
{
beta.boot=koef.beta.boot
}
matcoef<-cbind(matcoef,beta.boot)
}

matcoef1<-matcoef[,2:(b+1)]
mean.beta.boot<-apply(matcoef1,1,mean)
bias.beta.boot<-mean.beta.boot-koef.beta
var.beta.boot<-apply(matcoef1,1,var)
se.beta.boot<-sqrt(var.beta.boot)

loboundpre<-matrix(0:(p-1))
upboundpre<-matrix(0:(p-1))
for (t in 1:p)
{
matquant<-matcoef1[t,]
loboundquant<-quantile(matquant,(alpha/2))
upboundquant<-quantile(matquant,(1-(alpha/2)))

loboundpre<-rbind(loboundpre,loboundquant)
upboundpre<-rbind(upboundpre,upboundquant)
}

lobound0<-loboundpre[(p+1):(2*p),]
upbound0<-upboundpre[(p+1):(2*p),]

lobound<-matrix(0:(p-1))
lobound<-cbind(lobound,lobound0)

upbound<-matrix(0:(p-1))
upbound<-cbind(upbound,upbound0)


cat("Bootstrap Correlation Model","\n")
cat("\n")
cat("Bootstrap estimator for beta","\n")
print(mean.beta.boot)
cat("\n")
cat("bias beta","\n")
print(bias.beta.boot)
cat("\n")
cat("standard error beta","\n")
print(se.beta.boot)
cat("\n")
cat((1-alpha)*100,"% Confidence Interval for beta","\n")
cat("lower bound =","\n")
print(lobound)
cat("upper bound =","\n")
print(upbound)
cat("\n")

}


