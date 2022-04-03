library(multtest)
data(golub)
gol.fac=factor(golub.cl, levels=0:1, labels=c("ALL","AML"))
var.test(golub[1, ] ~ gol.fac)
n1=sum(gol.fac == "ALL")
n2=sum(gol.fac == "AML")
f=function(x) {  # F density function
    df(x, df1 = n1-1, df2 = n2-1) 
}
curve(f, 0, 3)
abline(v=var.test(golub[1, ]~gol.fac)$stat,col=2)  #statistics is F value here
f.test.fun=function(x) {
  var.test(x ~ gol.fac)$p.value 
}
f.test.fun(golub[1, ])
pvals=apply(golub, 1, f.test.fun)
pvals

#???????????????????????????:Bonferroni?????????:????????????10000???,??????????????????,??????N=5%/ 10000=0.000005
m=nrow(golub)
sum(pvals < 0.05/m)
sum(pvals < 1 - 0.95^(1/m))
max(which(sort(pvals) <= 0.05*seq(1, m)/m))
plot(sort(pvals)[1:360], xlab="Genes sorted by p-value", ylab="p-values")
abline(0, 0.05 / m, lwd=2, lty=2, col=2) # BH
abline(h=0.05 / m, lwd=2, lty=3, col=3) # Bonferroni
abline(h=1-(1-0.05)^(1/m), lwd=2, lty=4, col=4) # Sidek
