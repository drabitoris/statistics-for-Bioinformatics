BiocManager::install("ALL")
library(ALL)
data(ALL)
View(ALL)
probe.names <- c("1389_at", "35991_at", "40440_at")
ALLB123 <- ALL[ , ALL$BT %in% c("B1", "B2", "B3")]
X <- t(exprs(ALLB123)[probe.names, ])
G <- factor(ALLB123$BT, levels = c("B1", "B2", "B3"))
D <- data.frame(G, X)

library(nnet)
multinomial.fit <- multinom(G ~ ., data = D)
summary(multinomial.fit)
