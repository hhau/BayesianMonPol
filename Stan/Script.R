##
## Load required libraries and set options to run stan in parallel
##
library(MonoPoly)
library(fda)
library(rstan)
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())

##
## Data is the onechild data from package fda
## Rescale data so that x and y are both in [0, 1]
## Remember how y was scaled to undo the scaling in later plots
##

x <- with(onechild, (day-min(day))/diff(range(day)))
y <- onechild$height
miny <- min(y)
ry <- diff(range(y))
y <- (y-miny)/ry
N <- length(x)   

##
## Initial plot of the data
##

plot(x,y)


## 
## Choose degree of polynomial to fit,
## as well as the region of monotonicity
## Here we want to have a polynomial monotone on [0, infinity)
## we also specific the polynomial as monotone increasing,
## by choosing alpha to be one
##

q <- 9

lower.bound <- 0
upper.bound <- Inf

alpha <- 1

##
## Given order of polynomial and bounds, figure out K and op.mode
## This is just following Proposition 2.1
##

if (lower.bound == -Inf) {
  if (upper.bound == Inf) {
    op.mode <- 1 # whole real line
    k <- (q - 1)/2
    
  } else {
    op.mode <- 2 # (-inf, b]
    if(q %% 2 == 0) {
      k <- q / 2
      
    }  else {
      k <- (q - 1) / 2
      
    }
  }
} else {
  if (upper.bound == Inf) {
    op.mode <- 3 # [a, Inf)
    if(q %% 2 == 0) {
      k <- q / 2
      
    }  else {
      k <- (q - 1) / 2
      
    }
  } else {
    op.mode <- 4 # [a,b]
    if(q %% 2 == 0) {
      k <- q / 2
      
    }  else {
      k <- (q - 1) / 2
      
    }
  }
}

if ( (poly.degree %% 2 == 0) & (op.mode == 1) ) {
  stop("Fitting a monotonic polynomial to the whole real line requires q to be odd")
}

##
## Grid of values on which we want to have predictions
##

Nnew <- 101
xnew <- seq(from=min(x),to=max(x),length=Nnew)

##
## Put all the data to be passed to stan in a named list
##
data.in <- list(N = N, q = q, K = k,
                operation_mode = op.mode,
                y = y, x = x,
                a = lower.bound, b = upper.bound,
                alpha = alpha, Nnew = Nnew,
                xnew = xnew)

##
## Decide on number of iterations, this depends partly on the degree
## of the polynomial.
## See the 00README file
##

iter <- 5000

##
## Fit the model with stan
##

model.fit <- stan("MonPolyV0_0_4.stan", data = data.in, 
                  iter = iter, seed = 170117, chains = 2, cores = 2)


##
## beta final stuff, used for plotting

beta.final.samples <- extract(model.fit, "beta_final")
beta.final <- apply(beta.final.samples[[1]], 2, mean)

##
## plot the fitted polynomial
## 

plot(x,y)
lines(xnew, evalPol(xnew, beta.final))

##
## Add credible intervals and prediction intervals using the generated
## values for mupred and ypred
##
mupred_smp <- extract(model.fit, "mupred")[[1]]
cred.intervals <- apply(mupred_smp, 2, function(x) quantile(x=x, c(0.025, 0.975)))
# cred.intervals <- ry * cred.intervals + miny  ## undo the scaling of the y axis
lines(x = xnew, cred.intervals[1,], col = "black", lty = 3)
lines(x = xnew, cred.intervals[2,], col = "black", lty = 3)

ypred_smp <- extract(model.fit, "ypred")[[1]]
pred.intervals <- apply(ypred_smp, 2,  function(x) quantile(x=x, c(0.025, 0.975)))
# pred.intervals <- ry * pred.intervals + miny  ## undo the scaling of the y axis
lines(x = xnew, y = pred.intervals[1,], col = "cadetblue", lty = 1)
lines(x = xnew, y = pred.intervals[2,], col = "cadetblue", lty = 1)

