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
## Choose degree of polynomial to fit.
## Here we want to have a polynomial monotone on [0, infinity)
##
K <- 4
## q <- 2*K      ## choose whether the polynomial is supposed to be even
q <- 2*K + 1     ## or odd by commenting out one line

alpha <- 1       ## polynomial should be increasing on [0, infinity)
a <- 0

##
## Grid of values on which we want to have predictions
##
Nnew <- 101
xnew <- seq(from=min(x),to=max(x),length=Nnew)

##
## Put all the data to be passed to stan into a list
## The polynomial should be monotone on [0, infinity), so operation_mode = 3
##
data.in <- list(N = N, q = q, K = K,
                operation_mode = 3,
                y = y, x = x,
                a = 0, b = Inf,
                alpha = 1, Nnew = Nnew,
                xnew = xnew)
##
## Decide on number of iterations, this depends partly on the degree
## of the polynomial.
## See the 00README file
##
iter <- 50000

##
## Fit the model with stan
##
model.fit <- stan("MonPolyV0_0_4.stan", data = data.in, 
                  iter = iter, seed = 170117)


##
## Incidental plots of this script go into another pdf file
##
pdf(file=paste0("Scriptq=", q, ".pdf"), height=12, width=16)

##
## Plot the data, with y on the original scale
##
with(onechild, plot(x, height))

##
## Extract the beta parameters from the fitted model and calculate the
## posterior mean of the parameters, also adjust them to take the
## scaling of y into account.
##
beta_smp <- extract(model.fit, pars="beta_final")[[1]]
bhat <- colMeans(beta_smp)
bhat <- ry * bhat
bhat[1] <- bhat[1] + miny

##
## Add a line with the fitted regression line
##
poly.est <- evalPol(x, bhat)
lines(x = x, y = poly.est, col = "black", lty = 2)

##
## Add credible intervals and prediction intervals using the generated
## values for mupred and ypred
##
mupred_smp <- extract(model.fit, "mupred")[[1]]
cred.intervals <- apply(mupred_smp, 2, function(x) quantile(x=x, c(0.025, 0.975)))
cred.intervals <- ry * cred.intervals + miny  ## undo the scaling of the y axis
lines(x = xnew, cred.intervals[1,], col = "black", lty = 3)
lines(x = xnew, cred.intervals[2,], col = "black", lty = 3)

ypred_smp <- extract(model.fit, "ypred")[[1]]
pred.intervals <- apply(ypred_smp, 2,  function(x) quantile(x=x, c(0.025, 0.975)))
pred.intervals <- ry * pred.intervals + miny  ## undo the scaling of the y axis
lines(x = xnew, y = pred.intervals[1,], col = "cadetblue", lty = 1)
lines(x = xnew, y = pred.intervals[2,], col = "cadetblue", lty = 1)

##
## The following commands are producing the figures used in the paper
##
pdf(file=paste0("OneChildDeg", q, ".pdf"), width=6, height=6)
with(onechild, plot(x, height, type="n", ylim=c(123,131)))
polygon( c(xnew,rev(xnew)), c(pred.intervals[1,], rev(pred.intervals[2,])),
        col="lightgrey", border=NA)
polygon( c(xnew,rev(xnew)), c(cred.intervals[1,], rev(cred.intervals[2,])),
        col="darkgrey", border=NA)
with(onechild, points(x, height, ylim=c(120,133)))
lines(x = x, y = poly.est, col = "black")
dev.off()


##
## Diagnostic plots using the bayesplot package 
## Inspired by the vignettes of that package and refer to these
## vignettes for details
##
library(bayesplot)
posterior <- as.array(model.fit)

color_scheme_set("red")
mcmc_intervals(posterior, pars = "sd_y", regex_pars="beta_final*")

mcmc_areas(
  posterior, 
  pars = "sd_y",
  regex_pars = "beta_final*",
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)

##
## Trace plots
##
color_scheme_set("mix-blue-red")
mcmc_trace(posterior, regex_pars = "beta_final*")
mcmc_trace(posterior, pars=c("sd_y", "lp__"))

##
## ACF plots of chains
##
mcmc_acf(posterior, regex_pars="beta_final*", pars="sd_y")

##
## R hat statistics to assess convergence of chain
##
rhats <- rhat(model.fit)
ii <- grep("beta_final*", names(rhats))
ii <- c(ii, grep("sd_y", names(rhats)))
print(rhats[ii])
mcmc_rhat(as.vector(rhats[ii]))

##
##  Effective sample sizes
##
ratios <- neff_ratio(model.fit)
print(ratios[ii])
mcmc_neff(ratios[ii]) + yaxis_text()

dev.off()

##
## If not run as script but interactively, then one might want to
## explore the model fit using shinystan
## To do so, uncomment the following two lines
##
## library(shinystan)
## launch_shinystan(model.fit)
