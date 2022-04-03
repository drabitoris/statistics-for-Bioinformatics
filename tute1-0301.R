#tute1
install.packages('SNPassoc')
library(SNPassoc)
data(SNPs)
dim(SNPs)   #dimension columns and then rows
colnames(SNPs)
snp10001 <- SNPs$snp10001 # extracts a SNP
is.numeric(snp10001) # is it numeric? No.
class(snp10001) # so what is then? It's a factor.
levels(snp10001) # what are the levels of the factor?
length(snp10001) # sample size
freq <- table(snp10001) # frequency table
plot(freq)
View(snp10001)
x <-table(SNPs$casco,SNPs$snp10001)
View(SNPs$casco)
x
myChisq <- function(Z) {
    n <- sum(Z) # table total = sample size
    out <- 0 # initialize the output
    for (r in 1:2) { # loop across rows and cols
        for (c in 1:3) {
            o <- Z[r, c] # observed
            e <- sum(Z[r, ]) * sum(Z[, c]) / n # expected, "blanc means all"
            out <- out + (o - e)^2 / e # increment the output
        }
    }
    return(out)
}
myChisq(x)
qchisq(0.95, df = 2)
chisq.test(x)
BiocManager::install("multtest")
library(multtest)
data(golub)
ncol(golub) # sample size
nrow(golub) # number of genes
View(golub)
x <- golub[1042, ] # extracts expression values
golub.gnames[1042, 3] # determines label of gene
boxplot(x[golub.cl == 0], x[golub.cl == 1])
pos = 1080:1085
M = golub[pos,]
View(M)
rownames(M) = golub.gnames[pos, 3]
View(rownames(M))
heatmap(M)
levelplot(M)
head(golub.cl) # to see the label 'classification'?
?golub # to dig more if not familiar
# Cochran-armitage trend test
M <- data.frame(t(M))
plot(M)
library(lattice)
levelplot(cor(M))
?data.frame
?levelplot
Calculate
P100


# Excises
#2
install.packages("xfun")
library(SNPassoc)
data(SNPs)
a <- table(SNPs$snp10001)
a
(a[1] + a[2]/2)/(12+53+92)
# Here's a neater way:
alleleFreq.C <- sum(a * c(2, 1, 0)) / sum(a * 2)
alleleFreq.T <- sum(a * c(0, 1, 2)) / sum(a * 2)
alleleFreq <- c(C = alleleFreq.C, T = alleleFreq.T)
alleleFreq

#3. For gene 1042 in the golub dataset:
#(a) Explore the distribution of expression values using basic summary statistics, a his-
#    togram and a QQ-plot.
library(multtest)
data(golub)
x <- golub[1080, ]

View(x)
y <- golub[1085, ]
boxplot(x ~ golub.cl, col = 'blue', xlab = "x", ylab = "golub.cl")
boxplot(y ~ golub.cl, col = 'pink', xlab = "y", ylab = "golub.cl")

??multtest
View(golub)
x <- golub[1042, ]
length(x)
summary(x)
sd(x)
hist(x,col = "pink")
qqnorm(x)
qqline(x)
boxplot(x, horizontal = TRUE, col = "lightblue")

#(b) Describe the shape of the distribution, and produce summary statistics to quantify
#its centre and spread. Does it appear reasonable to assume that the data are from
#a normal distribution? Are there any outliers? Justify your answers.

#(c) When comparing the gene expression distributions between the two leukemia types,
#comment on the relative location, spread and shapes of the distributions. Are there
#any outliers? Name a measure of location and a measure of spread that you would
#recommend for these data and motivate your choices.
x0 <- x[golub.cl == 0]
x1 <- x[golub.cl == 1]
summary(x0)
summary(x1)
#interquantile IQR = Q3-Q1
boxplot(x ~ golub.cl, horizontal = TRUE, col = "lightblue",xlab = "x", ylab = "golub.cl")  



#4
dim(SNPs)
colnames(SNPs)

# golub: expression level; golub.cl: 0 or 1 (two genotypes?????????????????????) 
                                                         