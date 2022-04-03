library(multtest)
data(golub)



y = factor(golub.cl, levels = c(0, 1), labels = c("ALL", "AML"))
x = golub[1, ] # pick the first gene
fit1 = glm(y ~ x, family = "binomial")
summary(fit1)

f = function(x) {
   exp(0.7894 + 1.5203 * x) / (1 + exp(0.7894 + 1.5203 * x))
   }
curve(f, -2, 2, ylab="P(Y=1|X=x)", xlab="x", ylim=c(0,1))
points(x,as.numeric(y)-1,col=2) # add data to plot

est = coef(summary(fit1))
est
est[2, 4] # p-value for the coefficient of X

f(0.7)
predict(fit1, newdata = data.frame(x = 0.7), type = "response")


newdata1 = data.frame(x = c(-0.3, 0.7, 1.2))
predict(fit1, newdata = newdata1, type = "response")
f(newdata1$x)

pvals = 1:3051
for (j in 1:3051) {
   fit.j = glm(y ~ golub[j, ], family = "binomial")
   pvals[j] = coef(summary(fit.j))[2, 4]
   }

pos = sort(pvals, decreasing = FALSE, index.return = TRUE)$ix
pos = pos[1:200]
pos

install.packages('glmnet')
library(glmnet)
X = t(golub[pos, ]) # t() transposes the input matrix
lasso.fit = glmnet(X, y, family = "binomial", alpha = 1)
plot(lasso.fit, xvar = "lambda", lwd=2)
