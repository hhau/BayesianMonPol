# Librarys

library(rstan)
library(parallel)

# data generation
# This function returns noisy data from a polynomial of requested order 


GeneratePolyData <- function(data.lower, data.upper, poly.order, poly.type, n.data.points, sigma.noise,
                               noise.type, noise.df) {
    # poly.order corresponds to 0 = intercept only, 1  =  flat line, 2 = quadratic, etc etc
    x.values <- seq(from = data.lower, to = data.upper, length = n.data.points)
    if (noise.type == 1) {
      noise <- rnorm(n.data.points, sd = sigma.noise)
    } else {
      noise <- rt(n = n.data.points, df = noise.df)
    }

    if (poly.order < 1) {
      stop("incorrect poly order")

    } else if (poly.order == 1) {
      beta.true <- c(1/4, 1.5)
      y.values <- 1/4 + 1.5 * x.values + noise

    } else if (poly.order == 2 & poly.type == 1) {
      beta.true <- c(0, 0, 1)
      y.values <-  x.values^2 + noise

    } else if (poly.order == 2 & poly.type == 2) {
      beta.true <- c(1/4, 1, 1)
      y.values <- 1/4 + x.values + x.values^2 + noise

    } else if (poly.order == 3) {
      beta.true <- c(0,0,0,1)
      y.values <-  x.values^3 + noise

    } else if (poly.order == 4) {
      beta.true <- c(0,0,1,0,2)
      y.values <- x.values^2 + 2 * x.values^4 + noise

    } else if (poly.order == 5 & poly.type == 1) {
      beta.true <- c(0,0,0,0,0,1)
      y.values <- x.values^5 + noise

    } else if (poly.order == 5 & poly.type == 2) {
      beta.true <- c(0,3,2,0,0,1)
      y.values <- 3*x.values + 2*x.values^2 + x.values ^5 + noise

    } else if (poly.order == 7) {
      beta.true <- c(0,3,0,2,0,1,0,1)
      y.values <- 3 * x.values + 2 * x.values^3 + x.values^5 + x.values^7 + noise

    } else {
      print("Please choose appropriate polynomial order and type")
    }

    output.list <- list(x.values = x.values, y.values = y.values, beta.true = beta.true)
    return(output.list)

  }



  fitMonoPoly <- function(x.values, y.values, lower.bound, upper.bound, poly.degree,
                            alpha, n.new, mc.num, chains, iter, prefit = model.prefit) {

    # figure out which of the cases it's going to be
    # 1 == whole real line
    # 2 == (-inf, b]
    # 3 == [a, Inf)
    # 4 == [a, b]
    
    # then assign the appropriate value to k

    op.mode <- 0
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

    # new X values for sampling the posterior predictive distribution
    x.new <- seq(from = min(x.values), to = max(x.values),
                 length = n.new)
    # put all data in named list to pass to Stan
    data.in <- list(N = length(y.values), q = poly.degree, K = k,
                    operation_mode = op.mode,
                    y = y.values, x = x.values,
                    a = lower.bound, b = upper.bound,
                    alpha = alpha, Nnew = n.new,
                    xnew = x.new)

   
    # fit the model in stan and return the stanfit object
    model.fit <- stan(fit = prefit, data = data.in, chains = chains, cores = mc.num,
                      iter = iter, refresh = 100)

    return(model.fit)

  }



  
  
  generateCoveragePlots <- function(plot.iter = 2000, chain.iter = 5000, data.lower = -1, data.upper = 1,  n.x = 100,
                                      poly.lower = -Inf, poly.upper = Inf, alpha = 1, prefit = model.prefit) {

    poly.cases <- list(c(poly.order = 5, poly.type = 1), c(poly.order = 5, poly.type = 2))
    x.values <- seq(from = data.lower, to = data.upper, length = n.x)
    file.name <- "twoby"
    pdf(paste(file.name, "coverage80.pdf", sep = ""))
    par(mfrow = c(2,2), mar = c(5.1, 2.1, 1.1, 1.1), oma = c(0,0,3,0))
  
      for(case in poly.cases) {
      
        order.to.fit <- c(5,7)

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

        
          coverage.mat.80 <- matrix(0, nrow = plot.iter, ncol = n.x)
        
          x.mat.apparent <- cbind(1,x.values)
        
          for(i in 2:q) {
            x.mat.apparent <- cbind(x.mat.apparent, x.values^i)
          }
          
        
          for(i in i:plot.iter) {
            sim.data <- GeneratePolyData(data.lower = data.lower, data.upper = data.upper,
                                         poly.order = as.numeric(case[1]), poly.type = as.numeric(case[2]),
                                         n.data.points = n.x, sigma.noise = 0.1, noise.type = 1, noise.df = 1)
  
            model.fit <- fitMonoPoly(x.values = sim.data$x.values, y.values = sim.data$y.values,
                                     lower.bound = poly.lower, upper.bound =  poly.upper, poly.degree = q, 
                                     alpha = alpha, n.new = 2*n.x + 1,
                                     mc.num = 3, chains = 3, iter = chain.iter, prefit = prefit)
            
            samples <- extract(model.fit)
  
            # now to coverage plots
            fitted.vals <- x.mat.apparent %*% t(samples$beta_final)
            fitted.quants <- apply(fitted.vals, 1, function(x){quantile(x=x, c(0.025, 0.1, 0.9, 0.975))})
  
            # 80% coverage
            logical.coverage.80 <- (y.true > fitted.quants[2,]) & (y.true < fitted.quants[3,])
            coverage.mat.80[i,] <- t(logical.coverage.80)
          }
  
          # produce coverage plots with dev and dev.off
          means.80 <- apply(coverage.mat.80, 2, mean)
          # means.95 <- apply(coverage.mat.95, 2, mean)
          
          # save the data for plot recreation at a later date
          write.csv(means.80, paste(file.name,"_80coverageprobs.csv",sep = ""))
          
          plot(x = x.values, y = means.80, type = "l", ylim = c(0.55,1), xlab = "x", ylab = "")
          abline(a = 0.8, b= 0, lty = 2, lwd = 2, col = "grey")
        }
      }
    mtext("Coverage Probabilites", outer = T, cex = 1.5)
    dev.off()
  }
  
  # Simulate some data so that we can compile a prefit object for Stan
  simulated.data <- GeneratePolyData(data.lower = -1, data.upper =  1,
                                     poly.order = 5, poly.type = 2,
                                     n.data.points = 50, sigma.noise = 0.2,
                                     noise.type = 1, noise.df = 1)
  
  x.new.temp <- seq(from = -1, to = 1,
                    length = 2*length(simulated.data$x.values)+1)
  data.in.temp <- list(N = length(simulated.data$y.values), q = 1, K = 1000,
                       operation_mode = 1,
                       y = simulated.data$y.values, x = simulated.data$x.values,
                       a = -Inf, b = Inf,
                       alpha = 1, Nnew = 2*length(simulated.data$x.values)+1,
                       xnew = x.new.temp)
  
  model.prefit <- stan("MonPolyV0_0_4.stan", chains = 0, data = data.in.temp)
  

  # don't necessarily need to pass any arguments if defaults are acceptible
  generateCoveragePlots(plot.iter = 800, chain.iter = 2000)
  
  ##########
  
  # there's probably some kind of bug up in there, here's some backup code to produce the plot just incase
  file.names <- Sys.glob("./*.csv")
  x.vals <- seq(from = -1, to = 1, length = 100)
  cov.probs <- list()
  
  for(ii in 1:4) {  
    cov.probs[[ii]] <- read.csv(file = file.names[ii], header = T)[,2]
  }

  
  cov.probs <-as.data.frame(cov.probs, col.names = c("1", "2", "3", "4"))
  # matplot(x = x.vals, y =  cov.probs, type = "l", ylim = c(0.55, 1))
  # abline (0.8, 0, lwd = 2, lty = 3, col = "grey")
  pdf(file = "backup.pdf")
  par(mfrow = c(2,2), mar = c(5.1, 2.1, 1.1,  1.1), oma = c(0, 0, 3, 0))
  
  for(ii in 1:4) {
    plot(x = x.vals, y = cov.probs[,ii], type = "l", ylim = c(0.55, 1), xlab = "x", ylab = "")
    abline(0.8, 0, lty = 2, lwd = 2, col = "grey")
  }
  mtext("Coverage Probabilites", outer = T, cex = 1.3)
  dev.off()
