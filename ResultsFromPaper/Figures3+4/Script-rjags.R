##
## Load the library and the problematic data
##
library(rjags)
load("ProblemData.rda")

##
## We are fitting a quartic polynomial (K=2, q=4) to the data that is
##  monotone increasing (alpha=1) on the interval [0,1] (a=0, b=1)
##
N <- length(x)
K <- 2
q <- 2*K
alpha <- 1
a <- 0
b <- 1

##
## Create Figure 3 of the paper
##
pdf("problem.pdf", height=6, width=9)
plot(x, y, ylab="", main=expression("p(x) = "*2 *x^4 + x^2), cex=1.5)
dev.off()

##
## Create initial values for running two chains
##
inits <- list()


##
## Initial values based on fitting a quadratic polynomial to the data
##
fm <- lm(y~x+I(x^2))
b0 <- as.numeric(coef(fm)[1])
bn <- matrix(0, ncol=2, nrow=K+1)
bn[1,1] <- sqrt(abs(coef(fm)[2]+2*a*coef(fm)[3]))
bn[1,2] <- sqrt(2*abs(coef(fm)[3]))
tauy <- 1/mean(resid(fm)^2)
inits[[1]] <- list(b0=b0, tauy=tauy, bn=bn,
                   .RNG.name="base::Mersenne-Twister",
                   .RNG.seed=170117)

##
## Initial values based on fitting a constant polynomial to the data
##
b0 <- mean(y)
bn <- matrix(0, ncol=2, nrow=K+1)
tauy <- 1/var(y)
inits[[2]] <- list(b0=b0, tauy=tauy, bn=bn, 
                   .RNG.name="base::Mersenne-Twister",
                   .RNG.seed=170117)

##
## Compile JAGS model
##
m <- jags.model("PolMonCompactEven.bug",
                data=list(N=N, x=x, y=y, d=q+1, dm1=q,
                          K=K, Kp1=K+1, alpha=alpha, a=a, b=b),
                inits=inits,
                n.chains=length(inits))
##
## burn-in
##
update(m, n.iter=10000)

##
## set monitors on some parameters for the next iterations
##
out <-  coda.samples(m, c("beta", "sigy"),
                     n.iter=20000)

##
## create Figure 4 of the paper
##
## Old Version before final revision
##
## pdf("ProblemTrace.pdf", width=36, height=18)
## par(mfcol=c(3,4))
## par(cex=2*par()$cex)
## plot(out[,1:3,drop=FALSE], density=FALSE, auto.layout=FALSE, 
##      col=c("darkgrey", "black"), lty=c("solid", "dashed"))
## plot(out[,1:3,drop=FALSE], trace=FALSE, auto.layout=FALSE)
## plot(out[,4:6,drop=FALSE], density=FALSE, auto.layout=FALSE, 
##      col=c("darkgrey", "black"), lty=c("solid", "dashed"))
## plot(out[,4:6,drop=FALSE], trace=FALSE, auto.layout=FALSE)
## dev.off()

pdf("ProblemTrace.pdf", width=36, height=72)

gcol <- gray(c(0.5, 0))
par(mfcol=c(6,2), mar=c(4.1, 2.1, 2.1, 0.1), cex=3)
plot(out[,1,drop=FALSE], density=FALSE, auto.layout=FALSE, 
     main = expression("Trace of " * beta[0]),
     col=gcol, lty=c("solid", "dashed"))
plot(out[,2,drop=FALSE], density=FALSE, auto.layout=FALSE, 
     main = expression("Trace of " * beta[1]),
     col=gcol, lty=c("solid", "dashed"))
plot(out[,3,drop=FALSE], density=FALSE, auto.layout=FALSE, 
     main = expression("Trace of " * beta[2]),
     col=gcol, lty=c("solid", "dashed"))
plot(out[,4,drop=FALSE], density=FALSE, auto.layout=FALSE, 
     main = expression("Trace of " * beta[3]),
     col=gcol, lty=c("solid", "dashed"))
plot(out[,5,drop=FALSE], density=FALSE, auto.layout=FALSE, 
     main = expression("Trace of " * beta[4]),
     col=gcol, lty=c("solid", "dashed"))
plot(out[,6,drop=FALSE], density=FALSE, auto.layout=FALSE, 
     main = expression("Trace of " * sigma[epsilon]),
     col=gcol, lty=c("solid", "dashed"))
plot(out[,1,drop=FALSE], trace=FALSE, auto.layout=FALSE, 
     main = expression("Density of " * beta[0]))
plot(out[,2,drop=FALSE], trace=FALSE, auto.layout=FALSE, 
     main = expression("Density of " * beta[1]))
plot(out[,3,drop=FALSE], trace=FALSE, auto.layout=FALSE, 
     main = expression("Density of " * beta[2]))
plot(out[,4,drop=FALSE], trace=FALSE, auto.layout=FALSE, 
     main = expression("Density of " * beta[3]))
plot(out[,5,drop=FALSE], trace=FALSE, auto.layout=FALSE, 
     main = expression("Density of " * beta[4]))
plot(out[,6,drop=FALSE], trace=FALSE, auto.layout=FALSE, 
     main = expression("Density of " * sigma[epsilon]))

dev.off()

sessionInfo()
