library(MonoPoly)
library(fda)
library(rstan)
library(coda)
library(shinystan)
library(parallel)
library(RColorBrewer)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

K <- 4
q <- 2*K
##q <- 2*K + 1
alpha <- 1

x <- with(onechild, (day-min(day))/diff(range(day)))
y <- onechild$height
miny <- min(y)
ry <- diff(range(y))
y <- (y-miny)/ry

N <- length(x)
a <- 0

Nnew <- 101
xnew <- seq(from=min(x),to=max(x),length=Nnew)

plot(x,y)

data.in.temp <- list(N = N, q = q, K = K,
                  operation_mode = 3,
                  y = y, x = x,
                  a = 0, b = Inf,
                  alpha = 1, Nnew = Nnew,
                  xnew = xnew)

  # model.prefit <- stan("MonPolyV0_0_2.stan", chains = 0, data = data.in.temp)
  model.prefit <- stan("MonPolyV0_0_4.stan", chains = 0, data = data.in.temp)

  fit.mono.poly <- function(x.values, y.values, lower.bound, upper.bound, poly.degree,
                            alpha, n.new, multi, mc.num, chains, iter, prefit = model.prefit) {

    # figure out which of the cases it's going to be
    # 1 == whole real line
    # 2 == (-inf, b]
    # 3 == [a, Inf)
    # 4 == [a, b]

    op.mode <- 0
#  I really only need this section to figure out q and k
    if (lower.bound == -Inf) {
      if (upper.bound == Inf) {
        op.mode <- 1 # whole real line
        k <- (poly.degree - 1)/2

        } else {
          op.mode <- 2 # (-inf, b]
          if(poly.degree %% 2 == 0) {
            k <- poly.degree / 2

          }  else {
            k <- (poly.degree - 1) / 2

          }
      }
    } else {
        if (upper.bound == Inf) {
          op.mode <- 3 # [a, Inf)
          if(poly.degree %% 2 == 0) {
            k <- poly.degree / 2

          }  else {
            k <- (poly.degree - 1) / 2

          }
        } else {
            op.mode <- 4 # [a,b]
            if(poly.degree %% 2 == 0) {
              k <- poly.degree / 2

            }  else {
              k <- (poly.degree - 1) / 2

            }
        }
    }

    if ( (poly.degree %% 2 == 0) & (op.mode == 1) ) {
      stop("Fitting a monotonic polynomial to the whole real line requires q to be odd")
    }


    x.new <- seq(from = min(x.values), to = max(x.values),
                 length = n.new)
    data.in <- list(N = length(y.values), q = poly.degree, K = k,
                    operation_mode = op.mode,
                    y = y.values, x = x.values,
                    a = lower.bound, b = upper.bound,
                    alpha = alpha, Nnew = n.new,
                    xnew = x.new)

    # model.prefit <- stan(file = "MonPolyV0_0_2.stan", chains = 0, data = data.in)

    if(multi) {
      sflist <- mclapply(1:mc.num, mc.cores = mc.num, function(i) {stan(fit = prefit, data = data.in, chains = 1,
                                                              chain_id = i, iter = iter,
                                                              #control = list(max_treedepth = 30),
                                                              refresh = 100)})
      model.fit <- sflist2stanfit(sflist)

    }  else {

      model.fit <- stan(fit = prefit, data = data.in, chains = chains, iter = iter,
                        refresh = 100)
    }

    return(model.fit)

    # call stan with the prefit, drop the warmup please.
    # recall that extract() will grab all the samples, and i'll have to drop the first 1/2 or so for burn in
    # ^- this can be absolved with the save_warmup flag
    # computation time dependent.


  }

  model.fit <- fit.mono.poly(x.values = x, y.values = y,
                             lower.bound = 0, upper.bound =  Inf, 
                             poly.degree = q, alpha = 1, n.new = 101, 
                             multi = TRUE, mc.num = 2, chains = 2, 
                             iter = 500000) 


##  launch_shinystan(model.fit)
  samples <- extract(model.fit)

  #data plot
  with(onechild, plot(x, height))
  x.new <- xnew
  # make xmat
  x.mat <- cbind(1, x)
  data.order <- q
  for(i in 2:q){
    x.mat <- cbind(x.mat, x^i)
  }
  # true mean plot

  ####
  # this inherently assumes q > 1 btw ^
  ####



  # get beta est
  beta.estmates <- apply(samples$beta_final, 2, mean)
  beta.estmates <- ry*beta.estmates
  beta.estmates[1] <- beta.estmates[1]+miny

  # add estimated mean to plot
  poly.est <- x.mat %*% beta.estmates
  lines(x = x, y = poly.est, col = "black", lty = 2)

  # credible intervals
  # don't need to sort these because x is already in ascending order weewww
  fitted.vals <- x.mat %*% t(samples$beta_final)
  fitted.quants <- apply(fitted.vals, 1, function(x){quantile(x=x, c(0.025, 0.975))})
  fitted.quants <- ry*fitted.quants + miny
  lines(x = x, fitted.quants[1,], col = "black", lty = 3)
  lines(x = x, fitted.quants[2,], col = "black", lty = 3)

  # prediction intervals
  #x.new <- seq(from = min(simulated.data$x.values), to = max(simulated.data$x.values), length = 1000)
  pred.quants <- apply(samples$y_new, 2,  function(x){quantile(x=x, c(0.025, 0.975))})
  pred.quants <- pred.quants *ry + miny
  lines(x = x.new, y = pred.quants[1,], col = "cadetblue", lty = 1)
  lines(x = x.new, y = pred.quants[2,], col = "cadetblue", lty = 1)


pdf(file=paste0("OneChildDeg", q, ".pdf"), width=6, height=6)
with(onechild, plot(x, height, type="n", ylim=c(123,131)))
polygon( c(x.new,rev(x.new)), c(pred.quants[1,], rev(pred.quants[2,])),
        col="lightgrey", border=NA)
polygon( c(x,rev(x)), c(fitted.quants[1,], rev(fitted.quants[2,])),
        col="darkgrey", border=NA)
with(onechild, points(x, height, ylim=c(120,133)))
lines(x = x, y = poly.est, col = "black")
dev.off()




