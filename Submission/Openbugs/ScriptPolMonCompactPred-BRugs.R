##
## R script to drive the files PolMonCompactEvenPred.bug and
## PolMonCompactOddPred.bug, i.e. to fit polynomials that are
## constrained to be montone over a compact interval
##

## 
## Load some packages that we need and set the seed for reproducibility
##
library(BRugs)
library(MonoPoly)
set.seed(13116)

##
## Generate some data
##
N <- 41
x <- seq(from=-1, to=1, length=N)
y <- x^2 + 2*x^4 + rnorm(x, sd=0.1)

bt <- c(0,0,1,0,2) ## The vector of coefficients of the underlying polynomial

##
## The *.bug files are currently set up for K=2, if you change K here,
## you must also modify the *.bug files using the R snipped in the
## article!
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
##
if( q%%2 == 0 ){
  modelCheck("PolMonCompactEvenPred.bug")
}else{
  modelCheck("PolMonCompactOddPred.bug")
}

##
## Load the necessary data
##
bugsData(list("N", "x", "y", "q", "K", "alpha", "a", "b",
              "xnew", "Nnew"), file="data.txt")
modelData("data.txt")

##
## Compile the chains, load the inital values and generate default
## initals for all nodes that need to be initialised but were not.
##
## Note that the number of chains is automatically decided based on
## the length of the list with initial values.  Thus, adding to that
## list will increase the numbers of chains that are run
##
nc <- length(inits)
modelCompile(numChains=nc)
bugsInits(inits, numChain=nc, fileName=paste("init", 1:nc, ".txt", sep=""))
modelInits(paste("init", 1:nc, ".txt", sep=""))
modelGenInits()

##
## Update the model, set monitors on some parameters and update the
## model again to get MCMC values for those parameters that are now
## monitored.
##
modelUpdate(10000)
samplesSet(c("beta", "sigy"))
modelUpdate(50000)

##
## The various diagnostics offered by the BRugs package.
##
samplesStats("*")
samplesHistory("*") 
samplesDensity("*")
samplesBgr("*")             
samplesAutoC("*", 1)

##
## Set monitors to track mupred and ypred too.
##
samplesSet(c("mupred", "ypred"))
modelUpdate(10000)

##
## Create a plot that shows the data, the true polynomial and the
## fitted line.
##
par(ask=FALSE, mfrow=c(1,1))
plot(x,y)
beta <- samplesStats("beta")[,1]  ## The posterior mean of the parameters
lines(xnew, evalPol(xnew, bt), lty="dashed")
lines(xnew, evalPol(xnew, beta))

##
## Add credible intervals and prediction intervals using the generated
## values for mupred and ypred
##
lines(xnew, samplesStats("mupred")[,4], col="red")
lines(xnew, samplesStats("mupred")[,6], col="red")
lines(xnew, samplesStats("ypred")[,4], col="green")
lines(xnew, samplesStats("ypred")[,6], col="green")
