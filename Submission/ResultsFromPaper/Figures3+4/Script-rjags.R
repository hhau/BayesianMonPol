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
## Create Figure 2 of the paper
##
pdf("problem.pdf", height=8, width=12)
plot(x, y, ylab="", main=expression("p(x) = "*2 *x^4 + x^2))
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
## create Figure 3 of the paper
##
pdf("ProblemTrace.pdf", width=36, height=18)
par(mfrow=c(3,4))
plot(out, auto.layout=FALSE)
dev.off()

sessionInfo()
