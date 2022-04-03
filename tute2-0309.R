path <- "http://stat-gen.org/book.e1/data/FMS_data.txt"
fms <- read.delim(path, header = TRUE, sep = "\t")
attach(fms)
?attach # without assigning a new name and u can access it

install.packages=('genetics')
library(genetics)
mySNP <- genotype(akt1_t10726c_t12868c, sep = "")
head(fms)
View(mySNP)
HWE.exact(mySNP)
HWE.chisq(mySNP) # hypergeometric 
table(as.factor(fms$Race))

out <- HWE.exact(mySNP)
names(out)
out$p.value

mySNP <- genotype(akt1_t10726c_t12868c[Race=="Caucasian"], sep = "")
HWE.exact(mySNP)

m <- 25 # number of SNPs
pval.seq1 <- 1:m
pval.seq2 <- 1:m
for (i in 1:m) {
     mySNP1 <- genotype(fms[, i + 3], sep = "")
     mySNP2 <- genotype(fms[Race=="Caucasian", i + 3], sep = "")
     pval.seq1[i] <- HWE.exact(mySNP1)$p.value
     pval.seq2[i] <- HWE.exact(mySNP2)$p.value
   }
hist(pval.seq1,n=20)
hist(pval.seq2,n=20)
