path <- "http://stat-gen.org/book.e1/data/FMS_data.txt"
fms = read.delim(path, header = TRUE, sep = "\t")
attach(fms)
gtype = as.factor(actn3_r577x)
lm.fit = lm(NDRM.CH ~ gtype + Gender)
summary(lm.fit)
head(fms)
View(fms)
View(NDRM.CH)
-23.472 + c(-1, 1) * qt(0.975, df = 599) * 2.565
confint(lm.fit)

lm.fit2 = lm(NDRM.CH ~ gtype * Gender)
anova(lm.fit, lm.fit2)
plot(fitted(lm.fit2), residuals(lm.fit2))
abline(h = 0, lty = 2)
qqnorm(residuals(lm.fit2))
qqline(residuals(lm.fit2))
plot(fitted(lm.fit2))


library(multtest)
data(golub)
gol.fac = factor(golub.cl, levels = 0:1, labels = c("ALL", "AML"))
n1 = summary(gol.fac)[1] # number of ALL
n2 = summary(gol.fac)[2] # number of AML
n = n1 + n2
m = 10
X = golub[1:m, ]
muALL = rowMeans(X[, gol.fac == "ALL"])
muAML = rowMeans(X[, gol.fac == "AML"])
Y <- X
Y[, gol.fac == "ALL"] <- Y[, gol.fac == "ALL"] -
   matrix(muALL, m, n1, byrow = FALSE)
Y[, gol.fac == "AML"] <- Y[, gol.fac == "AML"] -
   matrix(muAML, m, n2, byrow = FALSE)
B = 250 # number of bootstrap samples
F.mat = matrix(nrow = B, ncol = m)
for (b in 1:B) {
   # Choose the n individuals assigned to group 1 and to group 2 by
     # drawing with replacement from the total sample. Since the newly
     # sampled individuals all come from the same distribution, we have
     # enforced the null hypothesis to be true.
     indices = sample(1:n, n, replace = TRUE)
    
       # Compute the F statistics. Note that, due to how we have created
       # the bootstrap samples, these are simulations from the null
       # distribution.
       f.test.fun = function(x) {
         var.test(x ~ gol.fac, conf.level = 0.95)$statistic
         }
       F.mat[b, ] = apply(Y[, indices], 1, f.test.fun)
      }