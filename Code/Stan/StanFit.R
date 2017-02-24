# Fitting Monotone polynomials with STAN
# 16 / 2 / 2016

# 1st, get one sub set (of polynomials) working
  # this includes:
    # writing the stan file and convolution function
    # ensuring the associated R output matches that of the jags implementation
# 2nd hopefully by this point gettin the other sub set of polynomaills should be easier
# 3rd, automate process of interval selection (If a and/or b is specified to be Infty in R,
#  translate this to positive_infinity() or negative_inifity() in stan and change code accordingly)

#### WD SET
setwd("/Users/hilary/Dropbox/BayesianMonPol/Code/Stan/")
setwd("/home/andrew-local/Dropbox/BayesianMonPol/Code/Stan/")
setwd("C:/Users/Andrew/Dropbox/BayesianMonPol/Code/Stan/")
####

# Librarys

library(MonoPoly)
library(rstan)
library(coda)
library(shinystan)
library(parallel)
library(RColorBrewer)

# data generation
  # this should be a function that takes the following
  # upper and lower intervals (Even if we are going to fit a polynomial over the whole line)
  # order of polynomial from which to generate
  # type of polynomial (multiple different quadratics for example)
  # number of data points
  # noise value?

  GeneratePolyData <- function(data.lower, data.upper, poly.order, poly.type, n.data.points, sigma.noise,
                               noise.type, noise.df) {
    # poly.order corresponds to 0 = intercept only, 1  =  flat line, 2 = quadratic, etc etc
    # poly type should take values from 1 - ?
    x.values <- seq(from = data.lower, to = data.upper, length = n.data.points)
    if (noise.type == 1) {
      noise <- rnorm(n.data.points, sd = sigma.noise)
    } else {
      noise <- rt(n = n.data.points, df = noise.df)
    }

    if (poly.order < 1) {
      # this is a silly case
      print("Don't be silly")
      stop("incorrect poly order")

    } else if (poly.order == 1) {
      # this is a straight line with off set (for some reason?)
      beta.true <- c(1/4, 1.5)
      y.values <- 1/4 + 1.5 * x.values + noise

    } else if (poly.order == 2 & poly.type == 1) {
      # just x^2
      beta.true <- c(0, 0, 1)
      y.values <-  x.values^2 + noise

    } else if (poly.order == 2 & poly.type == 2) {
      # 1/4 + x + x^2
      beta.true <- c(1/4, 1, 1)
      y.values <- 1/4 + x.values + x.values^2 + noise

    } else if (poly.order == 3) {
      # the only cubic is just x^3
      beta.true <- c(0,0,0,1)
      y.values <-  x.values^3 + noise

    } else if (poly.order == 4) {
      # again only one type of quartic
      # x^2 + 2x^4
      beta.true <- c(0,0,1,0,2)
      y.values <- x.values^2 + 2 * x.values^4 + noise

    } else if (poly.order == 5 & poly.type == 1) {
      # only consider x^5
      beta.true <- c(0,0,0,0,0,1)
      y.values <- x.values^5 + noise

    } else if (poly.order == 5 & poly.type == 2) {
      beta.true <- c(0,3,2,0,0,1)
      y.values <- 3*x.values + 2*x.values^2 + x.values ^5 + noise

    } else if (poly.order == 7) {
      # seems like you need more noise for this one, leave up to adjusting sigma
      beta.true <- c(0,3,0,2,0,1,0,1)
      y.values <- 3 * x.values + 2 * x.values^3 + x.values^5 + x.values^7 + noise

    } else {
      print("Please choose appropriate polynomial order and type (or something else has goofed)")
    }

    output.list <- list(x.values = x.values, y.values = y.values, beta.true = beta.true)
    return(output.list)

  }

  simulated.data <- GeneratePolyData(data.lower = -1, data.upper =  1,
                                     poly.order = 5, poly.type = 2,
                                     n.data.points = 50, sigma.noise = 0.2,
                                     noise.type = 1, noise.df = 1)
  # this function should return a data frame of the following
  # X and Y values
  # coefficient values of the data generate polynomial

  # get data into right shape and form and things
  # probably needs truncation ? or can just do with upper and lower limits

# data prep

  # get into list format


  # get inital values ? can try with / without
  # either way I have to write the stan file first

# model fitting

  # dummy fit to skip compile time
  # some of this is going to inherentily be specific to this task

  x.new.temp <- seq(from = -1, to = 1,
               length = 2*length(simulated.data$x.values)+1)
  data.in.temp <- list(N = length(simulated.data$y.values), q = 1, K = 1000,
                  operation_mode = 1,
                  y = simulated.data$y.values, x = simulated.data$x.values,
                  a = -Inf, b = Inf,
                  alpha = 1, Nnew = 2*length(simulated.data$x.values)+1,
                  xnew = x.new.temp)

  # model.prefit <- stan("MonPolyV0_0_2.stan", chains = 0, data = data.in.temp)
    model.prefit <- stan("MonPolyV0_0_3.stan", chains = 0, data = data.in.temp)

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


  # converage plots
  # do at the 80% level.

  # this is a fairly specific function to generate the coverage plots, write the n_eff and r_hat values
  # and other quantities for the 6 specified cases of iterest.
  generate.coverage.plots <- function(plot.iter = 10, chain.iter = 100, data.lower = -1, data.upper = 1,  n.x = 100,
                                      poly.lower = -Inf, poly.upper = Inf, alpha = 1, prefit = model.prefit) {

    poly.cases <- list(c(poly.order = 5, poly.type = 1), c(poly.order = 5, poly.type = 2))
    x.values <- seq(from = data.lower, to = data.upper, length = n.x)
    for(case in poly.cases) {
      #print(case)


      order.to.fit <- c(3,5,7)

      beta.true <- (GeneratePolyData(data.lower = data.lower, data.upper = data.upper,
                                     poly.order = as.numeric(case[1]), poly.type = as.numeric(case[2]),
                                     n.data.points = n.x, sigma.noise = 0.1, noise.type = 1, noise.df = 1)) [3]

      x.mat <- cbind(1, x.values)
      data.order <- length(beta.true$beta.true) - 1


      for(i in 2:data.order){
        x.mat <- cbind(x.mat, x.values^i)
      }

      y.true <- x.mat %*% as.matrix(beta.true$beta.true)

      for(q in order.to.fit) {

        file.name <- paste("case", as.numeric(case[1]), "type", as.numeric(case[2]), "order", q, sep = "")

        n.eff.container <- c()
        r.hat.container <- c()

        coverage.mat.80 <- matrix(0, nrow = plot.iter, ncol = n.x)
        coverage.mat.95 <- matrix(0, nrow = plot.iter, ncol = n.x)

        x.mat.apparent <- cbind(1,x.values)
        for(i in 2:q) {
          x.mat.apparent <- cbind(x.mat.apparent, x.values^i)
        }

        for(i in i:plot.iter) {
          sim.data <- GeneratePolyData(data.lower = data.lower, data.upper = data.upper,
                                       poly.order = as.numeric(case[1]), poly.type = as.numeric(case[2]),
                                       n.data.points = n.x, sigma.noise = 0.1, noise.type = 1, noise.df = 1)

          model.fit <- fit.mono.poly(x.values = sim.data$x.values, y.values = sim.data$y.values,
                                     lower.bound = poly.lower, upper.bound =  poly.upper, poly.degree = q, alpha = alpha, n.new = 2*n.x + 1,
                                     multi = TRUE, mc.num = 3, chains = 3, iter = chain.iter, prefit = prefit)
          # now need to extract n_eff + r_rhat values
          samples <- extract(model.fit)
          temp.df <- as.data.frame(summary(model.fit))

          n.eff.container <- rbind(n.eff.container, temp.df$summary.n_eff[grep("beta_final", row.names(temp.df))])
          r.hat.container <- rbind(r.hat.container, temp.df$summary.Rhat[grep("beta_final", row.names(temp.df))])

          colnames(r.hat.container) <- row.names(temp.df)[grep("beta_final",row.names(temp.df))]
          colnames(n.eff.container) <- row.names(temp.df)[grep("beta_final",row.names(temp.df))]

          # now to coverage plots
          beta.estmates <- apply(samples$beta_final, 2, mean)
          poly.est <- x.mat.apparent %*% beta.estmates
          beta.quantiles <- apply(samples$beta_final, 2, function(x){quantile(x=x, c(0.025, 0.1, 0.9, 0.975))})
          # don't need to sort these because x is already in ascending order weewww
          fitted.vals <- x.mat.apparent %*% t(samples$beta_final)
          fitted.quants <- apply(fitted.vals, 1, function(x){quantile(x=x, c(0.025, 0.1, 0.9, 0.975))})

          # 80% coverage
          logical.coverage.80 <- (y.true > fitted.quants[2,]) & (y.true < fitted.quants[3,])
          coverage.mat.80[i,] <- t(logical.coverage.80)
          # 95% coverage
          logical.coverage.95 <- (y.true > fitted.quants[1,]) & (y.true < fitted.quants[4,])
          coverage.mat.95[i,] <- t(logical.coverage.95)
        }
        ####
        ## should probably take the mean of these
        ####
        write.csv(n.eff.container, paste(file.name, "_eff.csv", sep = ""))
        write.csv(r.hat.container, paste(file.name, "_rhat.csv", sep = ""))

        # produce coverage plots with dev and dev.off
        means.80 <- apply(coverage.mat.80, 2, mean)
        means.95 <- apply(coverage.mat.95, 2, mean)

        pdf(paste(file.name, "coverage80.pdf", sep = ""), width = 10, height = 5)
          plot(x = x.values, y = means.80, type = "l", ylim = c(0,1))
          abline(a = 0.8, b= 0, lty = 2)
        dev.off()

        pdf(paste(file.name, "coverage95.pdf", sep = ""), width = 10, height = 5)
          plot(x = x.values, y = means.95, type = "l", ylim = c(0,1))
          abline(a = 0.95, b= 0, lty = 2)
        dev.off()
      }
    }

  }

generate.coverage.plots()


  simulated.data <- GeneratePolyData(data.lower = -1, data.upper =  1,
                                     poly.order = 1, poly.type = 2,
                                     n.data.points = 50, sigma.noise = 0.2,
                                     noise.type = 1, noise.df = 1)

  model.fit <- fit.mono.poly(x.values = simulated.data$x.values, y.values = simulated.data$y.values,
                             lower.bound = -Inf, upper.bound =  Inf, poly.degree = 1, alpha = 1, n.new = 101,
                             multi = TRUE, mc.num = 2, chains = 2, iter = 1000)



  launch_shinystan(model.fit)
  samples <- extract(model.fit)

  #data plot
  plot(x = simulated.data$x.values, y = simulated.data$y.values)
  x.new <- seq(from = min(x.values), to = max(x.values),
               length = n.new)
  # make xmat
  x.mat <- cbind(1, simulated.data$x.values)
  data.order <- length(simulated.data$beta.true) - 1
  for(i in 2:data.order){
    x.mat <- cbind(x.mat, simulated.data$x.values^i)
  }
  # true mean plot

  ####
  # this inherently assumes q > 1 btw ^
  ####



  y.true <- apply(t(t(x.mat) * (simulated.data$beta.true)), 1, sum)

  lines(x = simulated.data$x.values, y = y.true,col = "blue" ,lty = 7, cex = 3.3)

  # get beta est
  beta.estmates <- apply(samples$beta_final, 2, mean)

  # add estimated mean to plot
  poly.est <- x.mat %*% beta.estmates
  lines(x = simulated.data$x.values, y = poly.est, col = "black", lty = 2)

  # credible intervals
  beta.quantiles <- apply(samples$beta_final, 2, function(x){quantile(x=x, c(0.025, 0.1, 0.9, 0.975))})
  # don't need to sort these because x is already in ascending order weewww
  fitted.vals <- x.mat %*% t(samples$beta_final)
  fitted.quants <- apply(fitted.vals, 1, function(x){quantile(x=x, c(0.025, 0.1, 0.9, 0.975))})
  lines(x = simulated.data$x.values, fitted.quants[1,], col = "black", lty = 3)
  lines(x = simulated.data$x.values, fitted.quants[2,], col = "black", lty = 5)
  lines(x = simulated.data$x.values, fitted.quants[3,], col = "black", lty = 5)
  lines(x = simulated.data$x.values, fitted.quants[4,], col = "black", lty = 3)

  # prediction intervals
  #x.new <- seq(from = min(simulated.data$x.values), to = max(simulated.data$x.values), length = 1000)
  pred.quants <- apply(samples$ypred, 2,  function(x){quantile(x=x, c(0.025, 0.1, 0.9, 0.975))})
  lines(x = x.new, y = pred.quants[1,], col = "cadetblue", lty = 1)
  lines(x = x.new, y = pred.quants[2,], col = "coral", lty = 1)
  lines(x = x.new, y = pred.quants[3,], col = "coral", lty = 1)
  lines(x = x.new, y = pred.quants[4,], col = "cadetblue", lty = 1)





