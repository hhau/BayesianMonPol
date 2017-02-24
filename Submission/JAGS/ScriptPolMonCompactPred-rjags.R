##
## R script to drive the files PolMonCompactEvenPred.bug and
## PolMonCompactOddPred.bug, i.e. to fit polynomials that are
## constrained to be montone over a compact interval
##

## 
## Load some packages that we need and set the seed for reproducibility
##
library(MonoPoly)
library(rjags)
set.seed(13116)

##
## Generate some data
##
N <- 41
x <- seq(from=-1, to=1, length=N)
y <- x^2 + 2*x^4 + rnorm(x, sd=0.1)

bt <- c(0,0,1,0,2) ## The vector of coefficients of the underlying polynomial

##
## Say we want to fit a quartic or quintic (K=2) polynomial that is
## monotone over the interval [0,1] to these data
##
K <- 2
q <- 2*K         ## choose whether the polynomial is supposed to be even
## q <- 2*K + 1  ## or odd by commenting out one line

alpha <- 1       ## we want the polynomial to be monotone increasing on [0,1]
a <- 0
b <- 1

##
## Grid of values on which we want to have predictions
##
Nnew <- 101
xnew <- seq(from=-1,to=1,length=Nnew)

##
## Initial plot of the data
##
plot(x,y)

##
## Create initial values for running several chains
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
inits[[1]] <- list(beta0=b0, tauy=tauy, bn=bn)

##
## Initial values based on fitting a constant polynomial to the data
##
b0 <- mean(y)
bn <- matrix(0, ncol=2, nrow=K+1)
tauy <- 1/var(y)
inits[[2]] <- list(beta0=b0, tauy=tauy, bn=bn)

##
## Compile the appropriate model
## This also loads the necessary data and the initial values
##
## Note that the number of chains is automatically decided based on
## the length of the list with initial values.  Thus, adding to that
## list will increase the numbers of chains that are run
##
if( q%%2 == 0 ){
  m <- jags.model("PolMonCompactEvenPred.bug",
                  data=list(N=N, x=x, y=y, q=q, K=K, 
                      alpha=alpha, a=a, b=b,
                      xnew=xnew, Nnew=Nnew),
                  inits=inits,
                  n.chains=length(inits))
}else{
  m <- jags.model("PolMonCompactOddPred.bug",
                  data=list(N=N, x=x, y=y, q=q, K=K, 
                      alpha=alpha, a=a, b=b,
                      xnew=xnew, Nnew=Nnew),
                  inits=inits,
                  n.chains=length(inits))
}

##
## Update the model for burn-in
##
update(m, n.iter=10000)

##
## Further updates (MCMC interations) while monitoring some parameters.
## 
out <-  coda.samples(m, c("beta", "sigy"),
                     n.iter=20000)

##
## The summary and trace plots offered by the rjags package.
##
summary(out)
par(ask=TRUE)
plot(out)
par(ask=FALSE)

##
## Further updates (MCMC interations) while monitoring some parameters.
## 
out <-  coda.samples(m, c("beta", "mupred", "ypred"),
                     n.iter=10000)

## The posterior mean of the beta parameters
beta <- colMeans(do.call("rbind",lapply(out, function(x) x[,1:(q+1)])))

##
## Create a plot that shows the data, the true polynomial and the
## fitted line.
##
plot(x,y)
lines(xnew, evalPol(xnew, bt), lty="dashed")
lines(xnew, evalPol(xnew, beta))

##
## Add credible intervals and prediction intervals using the generated
## values for mupred and ypred
##
mupred <- do.call("rbind", lapply(out, function(x) x[,(q+1)+1:Nnew]))
lines(xnew, apply(mupred, 2, function(x) quantile(x, 0.025)), col="red")
lines(xnew, apply(mupred, 2, function(x) quantile(x, 0.975)), col="red")

ypred <- do.call("rbind", lapply(out, function(x) x[,(q+1)+Nnew+1:Nnew]))
lines(xnew, apply(ypred, 2, function(x) quantile(x, 0.025)), col="green")
lines(xnew, apply(ypred, 2, function(x) quantile(x, 0.975)), col="green")
