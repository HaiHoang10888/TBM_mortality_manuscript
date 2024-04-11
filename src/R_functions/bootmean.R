boot.mean<-function(x,n,k,alpha)
{
meanx=mean(x)
cat("Mean Original Sample :",meanx,"\n")

meanbootsamp=NULL
for(i in 1:k)
{
data=sample(x,n,replace=TRUE,prob=NULL)
meansamp=mean(data)
{
meanbootsamp[i]=meansamp
}
}

meanboot=mean(meanbootsamp)
cat("Mean Bootstrap Sample :",meanboot,"\n")

bootbias=meanboot-mean(x)
cat("bootstrap bias =",bootbias,"\n")
cat("\n")

varboot=var(meanbootsamp)
seboot=sqrt(varboot)
cat("bootstrap standard error =",seboot,"\n")
cat("\n")

lobound=quantile(meanbootsamp,(alpha/2))
upbound=quantile(meanbootsamp,(1-(alpha/2)))

cat(((1-alpha)*100),"% Confidence Interval for Mean Statistic:","\n")
cat("lower bound =",lobound,"\n")
cat("upper bound =",upbound,"\n")

hist(meanbootsamp)
}
