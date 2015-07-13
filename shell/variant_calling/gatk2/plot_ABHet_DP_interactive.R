
#############
### Task 1 - plot histograms and boxplots of ABHet and DP for variants in dbSNP129 (known) and variants not in dbSN19 (novel).
#############

# data file is in tab delimited format with a header line
# It needs to contain columns ABHet, DP and dbSNP129.RS
# missing values are NA

#define variables
data.file <- "aitman_taadx.ABHet.DP.tsv"
project <- "aitman_taadx"
dbSNP <- "dbSNP129"

#script

data <- read.table(file = data.file, header = TRUE, stringsAsFactors = FALSE)

known <- subset( x = data, subset = !is.na(dbSNP129.RS)) 
novel <- subset( x = data, subset = is.na(dbSNP129.RS))

pdf( paste(project, "_ABHet_DP_plots_", dbSNP, ".pdf", sep = "") )
boxplot(known$ABHet, novel$ABHet, names = c(paste("known (", dbSNP, ")", sep=""), "novel"), main = paste(project, "ABHet values", sep = " "), ylab = "ABHet")

hist(known$ABHet, breaks = 40, col=rgb(1,0,0,0.5), xlim=c(0,1), freq = FALSE, main = paste(project, "ABHet density histogram", sep = " "), xlab = "ABHet" )
hist(novel$ABHet, breaks = 40, col=rgb(0,0,1,0.5), freq = FALSE, add = TRUE)
legend( x = "topright", legend = c(paste("known (", dbSNP, ")", sep = ""), "novel"), fill = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)) )

boxplot(known$DP, novel$DP, names = c(paste("known (", dbSNP, ")", sep = ""), "novel"), main = paste(project, "DP values", sep = " "), ylab = "DP")

hist(known$DP, breaks = 40, col=rgb(1,0,0,0.5), freq = FALSE, main = paste(project, "DP density histogram", sep = " "), xlab = "DP" )
hist(novel$DP, breaks = 40, col=rgb(0,0,1,0.5), freq = FALSE, add = TRUE)
legend( x = "topright", legend = c(paste("known (", dbSNP, ")", sep = ""), "novel"), fill = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)) )

plot(known$ABHet ~ known$DP, main = paste(project, ", known (", dbSNP, ")", sep = ""), xlab = "DP", ylab = "ABHet")
plot(novel$ABHet ~ novel$DP, main = paste(project, ", novel (", dbSNP, ")", sep = ""), xlab = "DP", ylab = "ABHet")

dev.off()

summary(known)

rm(list=ls(all=TRUE))

#############
### Task 2 - fit mixture distribution and select filtering threshold
#############

# data file is in tab delimited format with a header line
# It needs to contain columns ABHet, DP and dbSNP129.RS
# missing values are NA

#define variables
data.file <- "aitman_taadx.ABHet.DP.tsv"
project <- "aitman_taadx"
dbSNP <- "dbSNP129"

#script

data <- read.table(file = data.file, header = TRUE, stringsAsFactors = FALSE)

#remove missing ABHet values and transform data to get right skewed distribution

ABHet <- data$ABHet[!is.na(data$ABHet)]
x <- 1-ABHet

#inspect the data (if running interactively)
hist(x, freq=FALSE, breaks = 40)

#define mixture distribution functions
dmixture.p <- function(x,p)
    p[1]*dnorm(x,p[2],p[3]) + (1-p[1])*dgamma(x,p[4],p[5])
pmixture.p <- function(x,p)
    p[1]*pnorm(x,p[2],p[3]) + (1-p[1])*pgamma(x,p[4],p[5])

f.neglog <- function(p) -sum(log(dmixture.p(x,p)))

#define starting parameters

start.params <- c(0.5, 0.5, 0.08, 11, 59) # best for mixture with gamma (aitman_taadx project)

#run MLE
mle.fit <- optim(start.params, f.neglog)
mle.fit$par

#Kolmogorov-Smirnoff test
ks.test(x, "pmixture.p", mle.fit$par)

#inspect interactively
hist(x, freq=FALSE, breaks = 40, main = paste(project, " ABHet mixture distributions", sep = ""), xlab = "1 - ABHet")
t <- seq(0,1,len=100)
lines(t,dmixture.p(t,mle.fit$par), col = "red")
lines(t,mle.fit$par[1]*dnorm(t,mle.fit$par[2],mle.fit$par[3]), col = "blue")
lines(t,(1-mle.fit$par[1])*dgamma(t,mle.fit$par[4],mle.fit$par[5]), col = "green")
legend( x = "topright", legend = c("mixture", "negatives (gamma distribution)", "positives (Gaussian distribution)"), lty = c(1,1,1), col = c("red", "green", "blue"))


#find threshold z for chosen FDR
z <- 0.2805
fdr <- (1 - mle.fit$par[1])*(1-pgamma(z,mle.fit$par[4], mle.fit$par[5])) / ( mle.fit$par[1]*(1-pnorm(z,mle.fit$par[2],mle.fit$par[3])) +  (1 - mle.fit$par[1])*(1-pgamma(z,mle.fit$par[4], mle.fit$par[5])))
fdr

#sensitivity at chosen z
sens <- mle.fit$par[1]*(1-pnorm(z,mle.fit$par[2],mle.fit$par[3]))/mle.fit$par[1]
sens

# sensitivity = TP/P = TP/(TP+FN)
# specificity = TN/N = TN/(FP+TN)
# roc = sensitivity ~ (1-specificity)

## plot to a file
pdf(paste(project, "_mixture_plots.pdf", sep = ""))
#plot histogram and fitted distributions
hist(x, freq=FALSE, breaks = 40, main = paste(project, " ABHet mixture distributions", sep = ""), xlab = "1 - ABHet")
t <- seq(0,1,len=100)
lines(t,dmixture.p(t,mle.fit$par), col = "red")
lines(t,mle.fit$par[1]*dnorm(t,mle.fit$par[2],mle.fit$par[3]), col = "blue")
lines(t,(1-mle.fit$par[1])*dgamma(t,mle.fit$par[4],mle.fit$par[5]), col = "green")
legend( x = "topright", legend = c("mixture", "negatives (gamma distribution)", "positives (Gaussian distribution)"), lty = c(1,1,1), col = c("red", "green", "blue"))

#plot FDR against ABHet threshold
z <- seq(0,1,len=100)
fdr <- (1 - mle.fit$par[1])*(1-pgamma(z,mle.fit$par[4], mle.fit$par[5])) / ( mle.fit$par[1]*(1-pnorm(z,mle.fit$par[2],mle.fit$par[3])) +  (1 - mle.fit$par[1])*(1-pgamma(z,mle.fit$par[4], mle.fit$par[5])))

plot(fdr~z, type="l", main = paste(project, " false discovery rate", sep = ""), xlab = "1 - ABHet", ylab = "FDR")

#plot ROC curves
# ROC based on mixture model
z <- seq(0,1,len=100)
sens <- mle.fit$par[1]*(1-pnorm(z,mle.fit$par[2],mle.fit$par[3]))/mle.fit$par[1]
spec <- (1 - mle.fit$par[1])*pgamma(z,mle.fit$par[4], mle.fit$par[5])/(1-mle.fit$par[1])
spec.roc <- 1 - spec

plot(sens ~ spec.roc , type = "l", main = paste(project, " ROC based on mixture model", sep = ""), xlab = "1 - specificity", ylab = "sensitivity")

# ROC assuming all novel variants are false positives (negatives)
known <- subset( x = data, subset = !is.na(dbSNP129.RS)) 
novel <- subset( x = data, subset = is.na(dbSNP129.RS))

z <- seq(0,1,len=100)
pos <- known$ABHet[!is.na(known$ABHet)]
n.true <- length(pos)
n.tp <- numeric()
neg <- novel$ABHet[!is.na(known$ABHet)]
n.neg <- length(neg)
n.tn <- numeric()

for ( i in 1:100 ) {
    n.tp[i] <- length(pos[pos < z[i]])
    n.tn[i] <- length(neg[neg > z[i]])
}

sens.counts <- n.tp/n.true
spec.counts <- n.tn/n.neg
spec.counts.roc <- 1-spec.counts

plot(sens.counts ~ spec.counts.roc , type = "l", main = paste(project, " ROC \nassuming all variants not in dbSNP129 are false positives", sep = ""), xlab = "1 - specificity", ylab = "sensitivity")

dev.off()

rm(list=ls(all=TRUE))
