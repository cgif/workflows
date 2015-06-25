
#############
### Task 1 - plot histograms and boxplots of ABHet and DP for variants in dbSNP129 (known) and variants not in dbSN19 (novel).
#############

# data files are in format "variant\tinfo.ABHet\tinfo.DP", with a header line

#define variables
known.file <- "aitman_taadz_dbSNP129_ABHet_DP.txt"
novel.file <- "aitman_taadz_not_dbSNP129_ABHet_DP.txt"
project <- "aitman_taadz"
dbSNP <- "dbSNP129"

#script

known <- read.table(file = known.file, header = TRUE, stringsAsFactors = FALSE)
novel <- read.table(file = novel.file, header = TRUE, stringsAsFactors = FALSE)

pdf( paste(project, "_ABHet_DP_plots_", dbSNP, ".pdf", sep = "") )
boxplot(known$info.ABHet, novel$info.ABHet, names = c(paste("known", dbSNP, sep=" "), "novel"), main = paste(project, dbSNP, "ABHet values", sep = " "), ylab = "ABHet")

hist(known$info.ABHet, breaks = 40, col=rgb(1,0,0,0.5), xlim=c(0,1), freq = FALSE, main = paste(project, dbSNP, "ABHet density histogram", sep = " "), xlab = "ABHet" )
hist(novel$info.ABHet, breaks = 40, col=rgb(0,0,1,0.5), freq = FALSE, add = TRUE)
legend( x = "topright", legend = c("known", "novel"), fill = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)) )

boxplot(known$info.DP, novel$info.DP, names = c(paste("known", dbSNP, sep = " "), "novel"), main = paste(project, dbSNP, "DP values", sep = " "), ylab = "DP")

hist(known$info.DP, breaks = 40, col=rgb(1,0,0,0.5), freq = FALSE, main = paste(project, dbSNP, "DP density histogram", sep = " "), xlab = "DP" )
hist(novel$info.DP, breaks = 40, col=rgb(0,0,1,0.5), freq = FALSE, add = TRUE)
legend( x = "topright", legend = c("known", "novel"), fill = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)) )

plot(known$info.ABHet ~ known$info.DP, main = paste(project, ", known ", dbSNP, sep = ""))
plot(novel$info.ABHet ~ novel$info.DP, main = paste(project, ", novel ", dbSNP, sep = ""))

dev.off()

summary(known)

rm(list=ls(all=TRUE))

#############
### Task 2 - fit mixture distribution and select filtering threshold
#############

# data files are in format "variant\tinfo.ABHet\tinfo.DP", with header

#define variables
known.file <- "aitman_taadz_dbSNP129_ABHet_DP.txt"
novel.file <- "aitman_taadz_not_dbSNP129_ABHet_DP.txt"
project <- "aitman_taadz"
dbSNP <- "dbSNP129"

#script

known <- read.table(file = known.file, header = TRUE, stringsAsFactors = FALSE)
novel <- read.table(file = novel.file, header = TRUE, stringsAsFactors = FALSE)

#combine ABHet values from both sets and remove missing values
x <- c(novel$info.ABHet, known$info.ABHet)
x <- x[!is.na(x)]
x <- 1-x

#inspect the data (if running interactively)
hist(x, freq=FALSE, breaks = 40)

#define mixture distribution functions
dmixture.p <- function(x,p)
    p[1]*dnorm(x,p[2],p[3]) + (1-p[1])*dgamma(x,p[4],p[5])
pmixture.p <- function(x,p)
    p[1]*pnorm(x,p[2],p[3]) + (1-p[1])*pgamma(x,p[4],p[5])

f.neglog <- function(p) -sum(log(dmixture.p(x,p)))

#define starting parameters
start.params <- c(0.5, 0.5, 0.03, 12, 57) # best for mixture with gamma

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
z <- 0.3026
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
z <- seq(0,1,len=100)
pos <- known$info.ABHet[!is.na(known$info.ABHet)]
n.true <- length(pos)
n.tp <- numeric()
neg <- novel$info.ABHet[!is.na(known$info.ABHet)]
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
