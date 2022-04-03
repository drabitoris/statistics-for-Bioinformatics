#1a
setwd("C:/Users/iefad/Desktop/course2021/bioinfosta/assignment/data-assignment1")
gt1 <- read.csv('genotypes1.csv')
gt2 <- read.csv('genotypes2.csv')
gt3 <- read.csv('genotypes3.csv')
pt1 <- read.csv('phenotype1.csv')
pt2 <- read.csv('phenotype2.csv')

tinytex::install_tinytex()

View(gt1)
dim(gt1)
dim(pt2)

as.numeric(pt1[1:5,2])
pt2[1:5,2]


#1b anova(mao) or just summary???
#test
model.1 <- lm(pt1$BMI ~ as.numeric(gt1$snp001))
summary(model.1)$coefficients[2,4]
summary(model.1)$coefficients
anova(model.1)[1,5]


model.2 <- lm(pt1$BMI ~ gt1$snp002)
summary(model.2)$coefficients[2,4]
summary(model.2)$coefficients
summary(model.2)
#answer
#genotype       pvalue
#1     snp001 2.788045e-01
#2     snp002 4.213277e-02
#3     snp003 1.106544e-02
#4     snp004 4.598057e-01
#5     snp005 4.773127e-01
recordpvalue <- data.frame(genotype = colnames(gt1)[2:201], pvalue = rep(0, 200))
for (i in 1:200) {
  model <-  lm(pt1$BMI ~ gt1[,i+1])
  recordpvalue[i,2] = anova(model)[1,5]
}
recordpvalue
#or anova????

#1c??????ggplot2
ggplot(data = recordpvalue) + geom_point(aes(x = 1:200, y = -log10(recordpvalue$pvalue)))
#1d
recordpvalue[order(recordpvalue$pvalue), ]
#answer: snp052 2.280280e-07

#2a
colnames(pt1)[2] = 'overweight'
pt1[,2] = as.logical(pt1[,2] > 25)
pt1[1,2]

pt = rbind(pt1,pt2)
pt

gt = rbind(gt1,gt2)
gt

#2b summary or anova?
combinepvalue <- data.frame(genotype = colnames(gt)[2:201], pvalue = rep(0, 200))
for (i in 1:200) {
  model <-  lm(as.numeric(pt$overweight) ~ as.numeric(gt[,i+1]))
  combinepvalue[i,2] = summary(model)$coefficients[2,4]
}

combinepvalue
# modellol <-  lm(as.numeric(pt$overweight) ~ as.numeric(gt[,2]))
# modellol$conf.int
#  genotype       pvalue
# 1     snp001 7.347895e-01
# 2     snp002 2.754467e-02
# 3     snp003 7.606861e-01
# 4     snp004 6.202933e-01
# 5     snp005 6.447473e-02


#2c manhatton plot? scatter plot?
plot(x=1:200,y=-log10(recordpvalue$pvalue),xlab = 'index', ylab = 'pvalue',xlim = c(0,14))
plot(x=1:200,y=-log10(combinepvalue$pvalue),xlab = 'index', ylab = 'pvalue' )

help(plot)
#2d
combinepvalue[order(combinepvalue$pvalue), ]
sum(combinepvalue$pvalue > 0.01)
sum(recordpvalue$pvalue > 0.01)
#answer: snp052 2.853124e-14

#2e(i)
m = length(combinepvalue$pvalue)
sum(combinepvalue$pvalue < 0.05/m)
#answer:16

#2e(ii)
max(which(sort(combinepvalue$pvalue) <= 0.05*seq(1, m)/m))
#answer: 29

#2e(iii)
hist(combinepvalue$pvalue,n=20)
nsig = sum(combinepvalue$pvalue < 0.001)
nullfrac = mean(combinepvalue$pvalue>0.1)/0.9
falsep = nullfrac*m*0.001
FDR = falsep/nsig
print(list(nsig,FDR))
# [[1]]
# [1] 21
# 
# [[2]]
# [1] 0.007407407

#2f???????


#3a
combinepvalue[order(combinepvalue$pvalue),][1:8, ]
# genotype       pvalue
# 52    snp052 2.853124e-14
# 51    snp051 1.090388e-13
# 53    snp053 5.688590e-11
# 107   snp107 7.394690e-11
# 56    snp056 1.073751e-10
# 105   snp105 1.986769e-10
# 106   snp106 5.010754e-10
# 54    snp054 7.754202e-10

#3b??????too long
small8 <- combinepvalue[order(combinepvalue$pvalue),][1:8, ]
cormatrix <- matrix(ncol = 8, nrow = 8)

colnames(cormatrix) <- small8$genotype
rownames(cormatrix) <- small8$genotype

index <- as.numeric(row.names(small8))

for (i in 1:8) {
  for (j in 1:8) {
    cormatrix[i,j] <- cor(as.numeric(gt[,index[i]+1]),as.numeric(gt[,index[j]+1]))
  }
}
cormatrix
nrow(cormatrix)

#3c
i <- 1
j <- 1
while (i <= 2){
  while (j <= nrow(cormatrix)){
    if (cormatrix[i,j] > 0.5 & (i != j)){
      cormatrix <- cormatrix[-j,-j]
    }
    else {
      j <- j+1
    }
  }
  j <- 1
  i <- i+1
}
cormatrix


left_index
index[!index %in% (index[3])]

#4a
View(gt)
overweightstate <- pt$overweight
has_052 <- as.numeric(gt[,53])
has_107 <- as.numeric(gt[,108])

fit1 = glm(overweightstate ~ has_052 + has_107,family="binomial")
summary(fit1)
coef(summary(fit1))[2:3,1:2]
#Estimate Std. Error
# has_052 0.5029238 0.07573746
# has_107 0.4025903 0.07375765

#4b
summary(fit1)
coef(fit1)
beta052 <- as.numeric(coef(fit1)[2])
beta107 <- as.numeric(coef(fit1)[3])

beta052
beta107

od_052 <- exp(beta052)
od_107 <- exp(beta107)

od_052
od_107

#4c
se1 <- coef(summary(fit1))[2,2]
se2 <- coef(summary(fit1))[3,2]
se1
se2
beta052
beta107
#?????
C.I_052 <- exp(beta052 + c(-1,1) * qnorm(0.975) * se1)
C.I_107 <- exp(beta107 + c(-1,1) * qnorm(0.975) * se2)

C.I_052 <- exp(beta052 + c(-1,1) * 1.96 * se1)
C.I_107 <- exp(beta107 + c(-1,1) * 1.96 * se2)

C.I_052
C.I_107
#5a
View(gt3)
new <- data.frame(has_052 = as.numeric(gt3[,53]), has_107 = as.numeric(gt3[,108]))
prerisk <- predict.glm(fit1,new,type = "response")
prerisk

help("predict.glm")
#5b
premodel <- data.frame(individuals = gt3$X, risk = prerisk)
premodel[order(premodel$risk), ]
plot(premodel[order(premodel$risk), ]$risk, xlab = "sorted individuals", ylab = "risk", main = "risks with the individuals sorted in order of increasing risk")

#5c
premodel[premodel$individuals=='indiv2001',]

#5d
maxrisk <- max(prerisk)
minrisk <- min(prerisk)
maxrisk
minrisk
ODmax <- maxrisk/(1 - maxrisk)
ODmin <- minrisk/(1 - minrisk)
OR <- ODmax/ODmin

## Q1c
```{r pressure, echo=FALSE}
#From the plot we can conclude that most p-values are larger than 0.01, which indicates a less significant association between the corresponding SNPs and BMI. The peak y value is higher than 6, which stands for less than e-6 p-value and the strongest association among these SNPs. 
install.packages("ggplot2")
library("ggplot2")
snps <- 1:200
pval1 <- recordpvalue$pvalue
ggplot(recordpvalue, aes(snps, -log10(pval1))) + geom_point()

```

tinytex::install_tinytex()
Sys.setenv(R_GSCMD="D:/gs9.54.0/bin/gswin64.exe")
tinytex::tlmgr_install("pdfcrop")

## Q2c
```{}
#From the plot we can conclude that most p-values are still larger than 0.01, which means weaker association with the overweight state. However, compared with the plot before combining data from study1 and study2, peak y values are higher and less y values are below 2, which indicates an overall stronger association with the overweight state(trait) here in this plot.
pval2 <- combinepvalue$pvalue
ggplot(combinepvalue, aes(x=snps, y=-log10(pval2))) + geom_point()
```
