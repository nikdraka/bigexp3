## Running all experiments for L1 regularisation, both unweighted and weighted

## Packages to run
library(tsutils)
library(forecast)
library(smooth) # the latest smooth 3.1.3.41006
library(MASS)
library(nloptr)
library(xtable)
library(greybox)
library(clusterGeneration)
library(parallel)
library(RColorBrewer)
library(devtools)

## Function to run
# generate time series
genseries <- function(dgp = c("ANN", "AAN", "AAdN","ANA", "AAA"), n = 100, freq = 7) {
  
  dgp <- dgp[1]
  n <- n[1]
  freq <- freq
  
  burnin <- 200
  
  # 1. Generating the time series, season initial values are generated randomly
  if (dgp == "ANN") {
    input.persist <- 0.4
    input.initial <- c(100)
    y <- sim.es(model = dgp, obs = burnin+n, nsim = 1, frequency = freq,
                bounds = "usual", persistence = input.persist,
                initial = input.initial)$data[(burnin+1):(burnin+n)]
  } else if (dgp == "AAN") {
    input.persist <- c(0.4, 0.3)
    input.initial <- c(100, 0.5)
    y <- sim.es(model = dgp, obs = burnin+n, nsim = 1, frequency = freq,
                bounds = "usual", persistence = input.persist,
                initial = input.initial)$data[(burnin+1):(burnin+n)]
  } else if (dgp == "ANA") {
    input.persist <- c(0.4, 0.2)
    input.initial <- c(100)
    y <- sim.es(model = dgp, obs = burnin+n, nsim = 1, frequency = freq,
                bounds = "usual", persistence = input.persist,
                initial = input.initial)$data[(burnin+1):(burnin+n)]
  } else if (dgp == "AAA") {
    input.persist <- c(0.4, 0.3, 0.2)
    input.initial <- c(100, 0.5)
    y <- sim.es(model = dgp, obs = burnin+n, nsim = 1, frequency = freq,
                bounds = "usual", persistence = input.persist,
                initial = input.initial)$data[(burnin+1):(burnin+n)]
  } else if (is.null(dgp)) {
    input.persist <- 0.4
    input.initial <- c(100)
    y <- sim.es(model = dgp, obs = burnin+n, nsim = 1, frequency = freq,
                bounds = "usual", persistence = input.persist,
                initial = input.initial)$data[(burnin+1):(burnin+n)]
  } else if (dgp == "AAdN") {
    input.persist <- c(0.4, 0.3)
    input.initial <- c(100, 0.5)
    y <- sim.es(model = dgp, obs = burnin+n, nsim = 1, frequency = freq,
                bounds = "usual", persistence = input.persist, phi = 0.9,
                initial = input.initial)$data[(burnin+1):(burnin+n)]
  } else if (dgp == "AAdA") {
    input.persist <- c(0.4, 0.3, 0.2)
    input.initial <- c(100, 0.5)
    y <- sim.es(model = dgp, obs = burnin+n, nsim = 1, frequency = freq,
                bounds = "usual", persistence = input.persist, phi = 0.9,
                initial = input.initial)$data[(burnin+1):(burnin+n)]
  } else if (dgp == "MAM") {
    input.persist <- c(0.3, 0.2, 0.1)
    input.initial <- c(200, 5)
    y <- sim.es(model = dgp, obs = n, nsim = 1, frequency = freq,
                bounds = "usual", persistence = input.persist,
                initial = input.initial, 
                randomizer = "rlnorm", mean = 0.0001, sd = 0.005)$data
  }
  
  return(ts(y, frequency = freq))
  
}
# loss rmse
loss.rmse <- function(actual, fitted, B) {
  error <- actual - fitted
  obsInSample <- length(actual)
  CFValue <- sqrt(sum(error^2)/obsInSample)
  
  return(CFValue)
}
# loss with weighted regularisation
loss.reg <- function(actual, fitted, B) {
  
  nB <- length(B)
  obsActual <- length(actual)
  yDenominator <- max(sd(diff(actual)),1)
  
  if (mdl == "ANN") {
    B[2] <- (B[2] - mean(actual[1:frequency(actual)]))/yDenominator # level
    
    error <- actual - fitted
    loss <- sqrt(sum((error/yDenominator)^2)/obsActual)
    
    if (type.loss == "l1") {
      CFValue <- (1-lambda)*loss + lambda*w1*abs(B[1])
    } else if (type.loss == "l2") {
      CFValue <- (1-lambda)*loss + lambda*w1*(B[1]^2)
    }
    
  } else if (mdl == "AAN") {
    B[3] <- (B[3] - mean(actual[1:frequency(actual)]))/yDenominator # level
    B[4] <- B[4]/ yDenominator # trend  
    
    error <- actual - fitted
    loss <- sqrt(sum((error/yDenominator)^2)/obsActual)
    
    if (type.loss == "l1") {
      CFValue <- (1-lambda)*loss + lambda*w1*abs(B[1]) + lambda*w2*abs(B[2])
    } else if (type.loss == "l2") {
      CFValue <- (1-lambda)*loss + lambda*w1*(B[1]^2) + lambda*w2*(B[2]^2)
    }
    
  } else if (mdl == "AAdN") {
    B[4] <- (B[4] - mean(actual[1:frequency(actual)]))/yDenominator # level
    B[5] <- B[5]/ yDenominator # trend  
    
    error <- actual - fitted
    loss <- sqrt(sum((error/yDenominator)^2)/obsActual)
    
    if (type.loss == "l1") {
      CFValue <- (1-lambda)*loss + lambda*w1*abs(B[1]) + lambda*w2*abs(B[2]) + lambda*w3*abs(1-B[3])
    } else if (type.loss == "l2") {
      CFValue <- (1-lambda)*loss + lambda*w1*(B[1]^2) + lambda*w2*(B[2]^2) + lambda*w3*(1-B[3])^2
    }
    
  } else if (mdl == "AAA") {
    B[4] <- (B[4] - mean(actual[1:frequency(actual)]))/yDenominator # level
    B[5] <- B[5]/ yDenominator # trend  
    B[6:nB] <- B[6:nB]/ yDenominator # season
    
    error <- actual - fitted
    loss <- sqrt(sum((error/yDenominator)^2)/obsActual)
    
    if (type.loss == "l1") {
      CFValue <- (1-lambda)*loss + lambda*w1*abs(B[1]) + lambda*w2*abs(B[2]) + lambda*w3*abs(B[3])
    } else if (type.loss == "l2") {
      CFValue <- (1-lambda)*loss + lambda*w1*(B[1]^2) + lambda*w2*(B[2]^2) + lambda*w3*(B[3]^2)
    }
    
  } else if (mdl == "AAdA") {
    B[5] <- (B[5] - mean(actual[1:frequency(actual)]))/yDenominator # level
    B[6] <- B[6]/ yDenominator # trend  
    B[7:nB] <- B[7:nB]/ yDenominator # season
    
    error <- actual - fitted
    loss <- sqrt(sum((error/yDenominator)^2)/obsActual)
    
    if (type.loss == "l1") {
      CFValue <- (1-lambda)*loss + lambda*w1*abs(B[1]) + lambda*w2*abs(B[2]) + lambda*w3*abs(B[3]) + lambda*w4*abs(1-B[4])
    } else if (type.loss == "l2") {
      CFValue <- (1-lambda)*loss + lambda*w1*(B[1]^2) + lambda*w2*(B[2]^2) + lambda*w3*(B[3]^2) + lambda*w4*(1-B[4])^2
    }
    
  } else if (mdl == "MAM") {
    B[4] <- (B[4] - mean(actual[1:frequency(actual)]))/yDenominator # level
    B[5] <- B[5]/ yDenominator # trend  
    B[6:nB] <- B[6:nB]/ yDenominator # season
    
    error <- actual - fitted
    loss <- sqrt(sum((error/yDenominator)^2)/obsActual)
    
    if (type.loss == "l1") {
      CFValue <- (1-lambda)*loss + lambda*w1*abs(B[1]) + lambda*w2*abs(B[2]) + lambda*w3*abs(B[3])
    } else if (type.loss == "l2") {
      CFValue <- (1-lambda)*loss + lambda*w1*(B[1]^2) + lambda*w2*(B[2]^2) + lambda*w3*(B[3]^2)
    }
    
  } else if (mdl == "ANA") {
    
    B[3] <- (B[3] - mean(actual[1:frequency(actual)]))/yDenominator # level
    B[4:nB] <- B[4:nB]/ yDenominator # season
    
    error <- actual - fitted
    loss <- sqrt(sum((error/yDenominator)^2)/obsActual)
    
    if (type.loss == "l1") {
      CFValue <- (1-lambda)*loss + lambda*w1*abs(B[1]) + lambda*w2*abs(B[2])
    } else if (type.loss == "l2") {
      CFValue <- (1-lambda)*loss + lambda*w1*(B[1]^2) + lambda*w2*(B[2]^2)
    }
    
  }
  
  return(CFValue)
}
# loss with unweighted regularisation
loss.reg.ur <- function(actual, fitted, B) {
  
  nB <- length(B)
  obsActual <- length(actual)
  yDenominator <- max(sd(diff(actual)),1)
  
  if (mdl == "ANN") {
    B[2] <- (B[2] - mean(actual[1:frequency(actual)]))/yDenominator # level
    
    error <- actual - fitted
    loss <- sqrt(sum((error/yDenominator)^2)/obsActual)
    
    if (type.loss == "l1") {
      CFValue <- (1-lambda)*loss + lambda*abs(B[1])
    } else if (type.loss == "l2") {
      CFValue <- (1-lambda)*loss + lambda*(B[1]^2)
    }
    
  } else if (mdl == "AAN") {
    B[3] <- (B[3] - mean(actual[1:frequency(actual)]))/yDenominator # level
    B[4] <- B[4]/ yDenominator # trend  
    
    error <- actual - fitted
    loss <- sqrt(sum((error/yDenominator)^2)/obsActual)
    
    if (type.loss == "l1") {
      CFValue <- (1-lambda)*loss + lambda*abs(B[1]) + lambda*abs(B[2])
    } else if (type.loss == "l2") {
      CFValue <- (1-lambda)*loss + lambda*(B[1]^2) + lambda*(B[2]^2)
    }
    
  } else if (mdl == "AAdN") {
    B[4] <- (B[4] - mean(actual[1:frequency(actual)]))/yDenominator # level
    B[5] <- B[5]/ yDenominator # trend  
    
    error <- actual - fitted
    loss <- sqrt(sum((error/yDenominator)^2)/obsActual)
    
    if (type.loss == "l1") {
      CFValue <- (1-lambda)*loss + lambda*abs(B[1]) + lambda*abs(B[2]) + lambda*abs(1-B[3])
    } else if (type.loss == "l2") {
      CFValue <- (1-lambda)*loss + lambda*(B[1]^2) + lambda*(B[2]^2) + lambda*(1-B[3])^2
    }
    
  } else if (mdl == "AAA") {
    B[4] <- (B[4] - mean(actual[1:frequency(actual)]))/yDenominator # level
    B[5] <- B[5]/ yDenominator # trend  
    B[6:nB] <- B[6:nB]/ yDenominator # season
    
    error <- actual - fitted
    loss <- sqrt(sum((error/yDenominator)^2)/obsActual)
    
    if (type.loss == "l1") {
      CFValue <- (1-lambda)*loss + lambda*abs(B[1]) + lambda*abs(B[2]) + lambda*abs(B[3])
    } else if (type.loss == "l2") {
      CFValue <- (1-lambda)*loss + lambda*(B[1]^2) + lambda*(B[2]^2) + lambda*(B[3]^2)
    }
    
  } else if (mdl == "AAdA") {
    B[5] <- (B[5] - mean(actual[1:frequency(actual)]))/yDenominator # level
    B[6] <- B[6]/ yDenominator # trend  
    B[7:nB] <- B[7:nB]/ yDenominator # season
    
    error <- actual - fitted
    loss <- sqrt(sum((error/yDenominator)^2)/obsActual)
    
    if (type.loss == "l1") {
      CFValue <- (1-lambda)*loss + lambda*abs(B[1]) + lambda*abs(B[2]) + lambda*abs(B[3]) + lambda*abs(1-B[4])
    } else if (type.loss == "l2") {
      CFValue <- (1-lambda)*loss + lambda*(B[1]^2) + lambda*(B[2]^2) + lambda*(B[3]^2) + lambda*(1-B[4])^2
    }
    
  } else if (mdl == "MAM") {
    B[4] <- (B[4] - mean(actual[1:frequency(actual)]))/yDenominator # level
    B[5] <- B[5]/ yDenominator # trend  
    B[6:nB] <- B[6:nB]/ yDenominator # season
    
    error <- actual - fitted
    loss <- sqrt(sum((error/yDenominator)^2)/obsActual)
    
    if (type.loss == "l1") {
      CFValue <- (1-lambda)*loss + lambda*abs(B[1]) + lambda*abs(B[2]) + lambda*abs(B[3])
    } else if (type.loss == "l2") {
      CFValue <- (1-lambda)*loss + lambda*(B[1]^2) + lambda*(B[2]^2) + lambda*(B[3]^2)
    }
    
  } else if (mdl == "ANA") {
    
    B[3] <- (B[3] - mean(actual[1:frequency(actual)]))/yDenominator # level
    B[4:nB] <- B[4:nB]/ yDenominator # season
    
    error <- actual - fitted
    loss <- sqrt(sum((error/yDenominator)^2)/obsActual)
    
    if (type.loss == "l1") {
      CFValue <- (1-lambda)*loss + lambda*abs(B[1]) + lambda*abs(B[2])
    } else if (type.loss == "l2") {
      CFValue <- (1-lambda)*loss + lambda*(B[1]^2) + lambda*(B[2]^2)
    }
    
  }
  
  return(CFValue)
}

# running model without optimised tuning parameter
simulStudy <- function(data = y, model = NULL, 
                       origin = 5, h = h, year.insample = 2,
                       type.loss = c("l1", "l2"), tuning = NULL, 
                       optimise = FALSE, weighted = TRUE) {
  y <<- data
  mdl <- model
  freq <- frequency(y)
  origin <- origin
  h <- h
  y.i <- year.insample
  type.loss <<- type.loss[1]
  
  if (weighted) {
    tuning <<- tuning
    lambda <<- tuning[1]
    w <<- tuning[2:length(tuning)]
    
    if (length(w) == 1) {
      w1 <<- 1
      
      w <<- w1
      
    } else if (length(w) == 2) {
      w1 <<- w[1]
      w2 <<- 1-w1
      
      if (w1 + w2 != 1 || w2 < 0) {
        return(c(1e300))
      }
      
      w <- c(w1, w2)
      
    } else if (length(w) == 3) {
      w1 <<- w[1]
      w2 <<- w[2]
      w3 <<- 1-(w1+w2)
      
      if (w1 + w2 + w3 != 1 || w3 < 0) {
        return(c(1e300))
      }
      
      w <<- c(w1, w2, w3)
      
    } else if (length(w) == 4) {
      w1 <<- w[1]
      w2 <<- w[2]
      w3 <<- w[3]
      w4 <<- 1-(w1+w2+w3)
      
      if (w1 + w2 + w3 + w4 != 1 || w4 < 0) {
        return(c(1e300))
      }
      
      w <- c(w1, w2, w3, w4)
      
    }
  } else if (!weighted) {
    
    lambda <<- tuning[1]
    w <<- NULL
    
  }
  
  
  list.fcst.reg <- vector("list", origin)
  list.error.reg <- array(NA, c(origin,h))
  list.scaled.error.reg <- array(NA, c(origin,h))
  
  list.fcst.bm <- vector("list", origin)
  list.error.bm <- array(NA, c(origin,h))
  list.scaled.error.bm <- array(NA, c(origin,h))
  
  data.train <- vector("list", origin)
  data.test <- array(NA, c(origin, h))
  
  list.param.reg <- vector("list", origin)
  list.param.bm <- vector("list", origin)
  
  # list.accpi.reg <- vector("list", origin)
  # list.accpi.bm <- vector("list", origin)
  # mat.param.reg <- array(NA, c(origin, h), dimnames = list(c(paste0("ori",1:origin)), c(paste0("t+",1:h))))
  # mat.param.bm <- array(NA, c(origin, h), dimnames = list(c(paste0("ori",1:origin)), c(paste0("t+",1:h))))
  
  if (optimise) {
    
    for (i in 1:origin) {
      
      yTrain <- window(y, start = c(1,1), end = c(y.i,freq+i))
      yTest <- window(y, start = c(y.i+1,1+i), end = c(y.i+1,freq+i))
      
      # Model with regularisation, column: DGP, row: horizon
      if (weighted) {
        mdl <<- model
        fit.reg <- adam(yTrain, model = mdl, lags = freq, loss = loss.reg)
      } else if (!weighted) {
        mdl <<- model
        fit.reg <- adam(yTrain, model = mdl, lags = freq, loss = loss.reg.ur)
      }
      
      fcst.reg <- forecast(fit.reg, interval = "empirical", level = c(0.90, 0.95, 0.99), h = h)
      
      error.reg <- yTest - fcst.reg$mean
      
      list.error.reg[i,] <- error.reg
      
    }
    
    return(mean(sqrt(list.error.reg[,1]^2)))
    
  } else {
    
    for (i in 1:origin) {
      
      yTrain <- window(y, start = c(1,1), end = c(y.i,freq+i))
      yTest <- window(y, start = c(y.i+1,1+i), end = c(y.i+1,freq+i))
      
      # Model with regularisation, column: DGP, row: horizon
      if (weighted) {
        mdl <<- model
        fit.reg <- adam(yTrain, model = mdl, lags = freq, loss = loss.reg)
      } else if (!weighted) {
        mdl <<- model
        fit.reg <- adam(yTrain, model = mdl, lags = freq, loss = loss.reg.ur)
      }
      
      fcst.reg <- forecast(fit.reg, interval = "empirical", level = c(0.90, 0.95, 0.99), h = h)
      
      error.reg <- yTest - fcst.reg$mean
      
      # Benchmark, column: DGP, row: horizon
      fit.bm <- adam(yTrain, model = mdl, lags = freq, loss = loss.rmse)
      
      fcst.bm <- forecast(fit.bm, interval = "empirical", level = c(0.90, 0.95, 0.99), h = h)
      
      error.bm <- yTest - fcst.bm$mean
      
      dy <- sum(abs(diff(yTrain)))
      
      list.fcst.reg[[i]] <- fcst.reg
      list.error.reg[i,] <- error.reg
      list.scaled.error.reg[i,] <- error.reg/dy
      list.fcst.bm[[i]] <- fcst.bm
      list.scaled.error.bm[i,] <- error.bm/dy
      list.error.bm[i,] <- error.bm
      data.train[[i]] <- yTrain
      data.test[i,] <- yTest
      list.param.reg[[i]] <- fit.reg$B
      list.param.bm[[i]] <- fit.bm$B
      
    }
    
    result.reg <- rbind(accuracy.etsreg(error = list.error.reg, horizon = h, scaled = FALSE),
                        accuracy.etsreg(error = list.scaled.error.reg, horizon = h, scaled = TRUE),
                        predinterval.etsreg(test = data.test, forecast = list.fcst.reg, origin = origin, h = h, level = 0.90),
                        predinterval.etsreg(test = data.test, forecast = list.fcst.reg, origin = origin, h = h, level = 0.95),
                        predinterval.etsreg(test = data.test, forecast = list.fcst.reg, origin = origin, h = h, level = 0.99))
    param.reg <- colMeans(t(sapply(list.param.reg, function(x) x)))
    
    result.bm <- rbind(accuracy.etsreg(error = list.error.bm, horizon = h, scaled = FALSE),
                       accuracy.etsreg(error = list.scaled.error.bm, horizon = h, scaled = TRUE),
                       predinterval.etsreg(test = data.test, forecast = list.fcst.bm, origin = origin, h = h, level = 0.90),
                       predinterval.etsreg(test = data.test, forecast = list.fcst.bm, origin = origin, h = h, level = 0.95),
                       predinterval.etsreg(test = data.test, forecast = list.fcst.bm, origin = origin, h = h, level = 0.99))
    param.bm <- colMeans(t(sapply(list.param.bm, function(x) x)))
    
    # return(list(fcst.bm = list.fcst.bm,
    #             fcst.reg = list.fcst.reg,
    #             error.bm = list.error.bm,
    #             error.reg = list.error.reg,
    #             data.train = t(sapply(data.train, function(x) x)),
    #             data.test = t(sapply(data.test, function(x) x)),
    #             param.bm = t(sapply(list.param.bm, function(x) x)),
    #             param.reg = t(sapply(list.param.reg, function(x) x)),
    #             tunpar = c(lambda, w)))
    
    return(list(result.reg = result.reg,
                result.bm = result.bm,
                param.reg = param.reg,
                param.bm = param.bm,
                tunpar = c(lambda, w),
                data = y))
    
  }
  
}

# optimise tuning parameter
optimal.tuning <- function(x0 = x0, data = y, model = NULL, 
                           origin = 5, h = h, year.insample = 2,
                           type.loss = c("l1", "l2"), weighted = TRUE) {
  
  opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD", "xtol_rel" = 1e-08, "maxeval" = 1000)
  
  y <- data
  mdl <- model
  origin <- origin
  h <- h
  y.i <- year.insample
  type.loss <- type.loss
  
  lb <- rep(0, length(x0))
  ub <- rep(1, length(x0))
  
  if (weighted) {
    x1 <- nloptr(x0, function(x) simulStudy(tuning = x, data = y, model = mdl,
                                            origin = origin, h = h, year.insample = y.i,
                                            type.loss = type.loss, optimise = TRUE, weighted = TRUE), 
                 lb = lb, ub = ub, opts = opts)
  } else if (!weighted) {
    x1 <- nloptr(x0, function(x) simulStudy(tuning = x, data = y, model = mdl,
                                            origin = origin, h = h, year.insample = y.i,
                                            type.loss = type.loss, optimise = TRUE, weighted = FALSE), 
                 lb = lb, ub = ub, opts = opts)
  }
  
  
  return(x1$solution)
}

# combine optimising the tuning parameter and estimate the model
compile.ets <- function(data = y, model = model, origin = origin, year.insample = 2,
                        h = h, type.loss = type.loss, 
                        tuning = tuning, weighted = FALSE) {
  
  y <- data
  mdl <- model
  origin <- origin
  freq <- frequency(y)
  h <- h
  type.loss <<- type.loss
  tuning <<- tuning
  weighted <- weighted
  y <- y
  y.ins <- year.insample
  # generating time series
  
  
  # small sample size
  tunpar1 <- optimal.tuning(x0 = tuning, data = y, model = mdl,
                            origin = origin, h = h, year.insample = y.ins, 
                            type.loss = type.loss, weighted = weighted)
  model.small <- simulStudy(data = y, model = mdl,
                            origin = origin, h = h, year.insample = y.ins,
                            type.loss = type.loss, tuning = tunpar1, 
                            optimise = FALSE, weighted = weighted)
  
  # # large sample size
  # tunpar2 <- optimal.tuning(x0 = tuning, y, model = mdl,
  #                           origin = origin, h = h, year.insample = 60, 
  #                           type.loss = "l1", weighted = weighted)
  # model.large <- simulStudy(y, model = mdl,
  #                           origin = origin, h = h, year.insample = 60,
  #                           type.loss = type.loss, tuning = tunpar2, 
  #                           optimise = FALSE, weighted = weighted)
  
  return(sample002 = model.small)
  
}

# generate time series, optimise the tuning parameter, and estimate the model
final.compile <- function(dgp = dgp, nobs = c(28, 420), freq = freq, model = model, h = h, 
                          origin = origin, type.loss = type.loss, 
                          tuning = NULL, weighted = FALSE) {
  
  dgp <- dgp
  nobs1 <- nobs[1]
  nobs2 <- nobs[2]
  mdl <<- model
  freq <- freq
  h <- h
  type.loss <- type.loss
  tuning <- tuning
  weighted <- weighted
  
  ## Small sample size
  y1 <-  ts(genseries(dgp = dgp, n = nobs1, freq = freq), frequency = freq)
  while (min(y1) < 0) {
    y1 <- ts(genseries(dgp = dgp, n = nobs1, freq = freq), frequency = freq)
  }
  
  model1 <- compile.ets(data = y1, model = mdl, origin = origin, h = h,
                        type.loss = type.loss, tuning = tuning, 
                        weighted = weighted, year.insample = 2)
  
  ## Large sample size
  y2 <-  ts(genseries(dgp = dgp, n = nobs2, freq = freq), frequency = freq)
  while (min(y2) < 0) {
    y2 <- ts(genseries(dgp = dgp, n = nobs2, freq = freq), frequency = freq)
  }
  
  model2 <- compile.ets(data = y2, model = mdl, origin = origin, h = h,
                        type.loss = type.loss, tuning = tuning, 
                        weighted = weighted, year.insample = 55)
  
  return(list(sample1 = model1,
              sample2 = model2))
}

# calculating accuracy
accuracy.etsreg <- function(error = error, horizon = horizon, scaled = FALSE) {
  
  error <- error
  h <- horizon
  
  if (!scaled) {
    acc <- array(NA, c(4, h), dimnames = list(c("RMSE", "MAE", "MdAE", "ME"), c(paste0("t+", 1:h))))
  } else {
    acc <- array(NA, c(4, h), dimnames = list(c("RMSSE", "MAsE", "MdAsE", "MsE"), c(paste0("t+", 1:h))))
  }
  
  for (i in 2:h) {
    
    acc[1,1] <- sqrt(mean((error[,1])^2))
    acc[2,1] <- mean(abs(error[,1]))
    acc[3,1] <- median(abs(error[,1]))
    acc[4,1] <- mean(error[,1])
    
    acc[1, i] <- mean(apply(error[,1:i], 1, function(x) sqrt(mean(x^2))))
    acc[2, i] <- mean(apply(error[,1:i], 1, function(x) mean(abs(x))))
    acc[3, i] <- mean(apply(error[,1:i], 1, function(x) median(abs(x))))
    acc[4, i] <- mean(apply(error[,1:i], 1, function(x) mean(x)))
    
  }
  return(acc)
}

# calculating the prediction intervals
predinterval.etsreg <- function(test = NULL, forecast = NULL, origin = origin, 
                                h = h, level = c(0.90, 0.95, 0.99)) {
  
  test <- test
  fcst <- forecast
  origin <- origin
  h <- h
  
  if (level == 0.90) {
    col.level <- 1
  } else if (level == 0.95) {
    col.level <- 2
  } else if (level == 0.99) {
    col.level <- 3
  }
  
  mat.mis <- array(NA, c(origin, h), dimnames = list(c(paste0("ori",1:origin)), c(paste0("t+1-",1:h))))
  mat.pinball <- array(NA, c(origin, h), dimnames = list(c(paste0("ori",1:origin)), c(paste0("t+1-",1:h))))
  mat.coverage <- array(NA, c(origin, h), dimnames = list(c(paste0("ori",1:origin)), c(paste0("t+1-",1:h))))
  for (i in 1:origin) {
    for (j in 1:h) {
      
      mat.mis[i,j] <- MIS(test[i,1:j], 
                          fcst[[i]]$lower[1:j,col.level], 
                          fcst[[i]]$upper[1:j,col.level], level)
      
      mat.pinball[i,j] <- pinball(test[i,1:j], fcst[[i]]$lower[1:j,col.level], (1-level)/2, loss = 1) +
        pinball(test[i,1:j], fcst[[i]]$upper[1:j,col.level], level+(1-level)/2, loss = 1)
      
      mat.coverage[i,j] <- mean(test[i,j] > fcst[[i]]$lower[1:j,col.level] & test[i,j] < fcst[[i]]$upper[1:j,col.level])
      
    }
  }
  
  mat.mis.mean <- apply(mat.mis, 2, mean)
  mat.pinball.mean <- apply(mat.pinball, 2, mean)
  mat.coverage.mean <- apply(mat.coverage, 2, mean)
  
  mat.interval <- rbind(mat.mis.mean, mat.pinball.mean, mat.coverage.mean)
  rownames(mat.interval) <- c(paste("MIS", level*100, sep = "_"), paste("Pinball", level*100, sep = "_"), 
                              paste("Coverage", level*100, sep = "_"))
  
  return(mat.interval)
}


# Running the parallel ----------------------------------------------------

set.seed(020193)

crs <- detectCores()
cl <- makeForkCluster(getOption("cl.cores", crs-1))
writeLines(paste("Running with", crs-1, 'cores'))
# Load packages to cluster
invisible(clusterCall(cl, function(pkgs) {
  library(clusterGeneration)
  library(smooth)
  library(nloptr)
}))

invisible(clusterExport(cl, "genseries"))
invisible(clusterExport(cl, "loss.rmse"))
invisible(clusterExport(cl, "loss.reg"))
invisible(clusterExport(cl, "loss.reg.ur"))
invisible(clusterExport(cl, "simulStudy"))
invisible(clusterExport(cl, "optimal.tuning"))
invisible(clusterExport(cl, "compile.ets"))
invisible(clusterExport(cl, "final.compile"))
invisible(clusterExport(cl, "accuracy.etsreg"))
invisible(clusterExport(cl, "predinterval.etsreg"))

runs <- 5

## DGP: ANN
# Without weights
start.timeUR <- Sys.time()
system.time({dgpANN_modelANN_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANN", nobs = c(28, 420), freq = 7, model = "ANN", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})
system.time({dgpANN_modelAAN_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANN", nobs = c(28, 420), freq = 7, model = "AAN", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})
system.time({dgpANN_modelANA_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANN", nobs = c(28, 420), freq = 7, model = "ANA", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})
system.time({dgpANN_modelAAA_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANN", nobs = c(28, 420), freq = 7, model = "AAA", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})
system.time({dgpANN_modelAAdA_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANN", nobs = c(28, 420), freq = 7, model = "AAdA", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1), weighted = FALSE))})
system.time({dgpANN_modelAAdN_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANN", nobs = c(28, 420), freq = 7, model = "AAdN", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1), weighted = FALSE))})
system.time({dgpANN_modelMAM_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANN", nobs = c(28, 420), freq = 7, model = "MAM", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})

# With weights
system.time({dgpANN_modelANN_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANN", nobs = c(28, 420), freq = 7, model = "ANN", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.1), weighted = TRUE))})
system.time({dgpANN_modelAAN_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANN", nobs = c(28, 420), freq = 7, model = "AAN", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.5, 0.5), weighted = TRUE))})
system.time({dgpANN_modelANA_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANN", nobs = c(28, 420), freq = 7, model = "ANA", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.5, 0.5), weighted = TRUE))})
system.time({dgpANN_modelAAA_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANN", nobs = c(28, 420), freq = 7, model = "AAA", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.33, 0.33, 0.33), weighted = TRUE))})
system.time({dgpANN_modelAAdA_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANN", nobs = c(28, 420), freq = 7, model = "AAdA", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1, 0.25, 0.25, 0.25, 0.25), weighted = TRUE))})
system.time({dgpANN_modelAAdN_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANN", nobs = c(28, 420), freq = 7, model = "AAdN", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1, 0.33, 0.33, 0.33), weighted = TRUE))})
system.time({dgpANN_modelMAM_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANN", nobs = c(28, 420), freq = 7, model = "MAM", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.33, 0.33, 0.33), weighted = TRUE))})
end.timeWR <- Sys.time()
computetime1 <- end.timeWR - start.timeUR

save(dgpANN_modelANN_l1_UR, file = "dgpANN_modelANN_l1_UR.Rdata")
save(dgpANN_modelAAN_l1_UR, file = "dgpANN_modelAAN_l1_UR.Rdata")
save(dgpANN_modelANA_l1_UR, file = "dgpANN_modelANA_l1_UR.Rdata")
save(dgpANN_modelAAA_l1_UR, file = "dgpANN_modelAAA_l1_UR.Rdata")
save(dgpANN_modelAAdA_l1_UR, file = "dgpANN_modelAAdA_l1_UR.Rdata")
save(dgpANN_modelAAdN_l1_UR, file = "dgpANN_modelAAdN_l1_UR.Rdata")
save(dgpANN_modelMAM_l1_UR, file = "dgpANN_modelMAM_l1_UR.Rdata")

save(dgpANN_modelANN_l1_WR, file = "dgpANN_modelANN_l1_WR.Rdata")
save(dgpANN_modelAAN_l1_WR, file = "dgpANN_modelAAN_l1_WR.Rdata")
save(dgpANN_modelANA_l1_WR, file = "dgpANN_modelANA_l1_WR.Rdata")
save(dgpANN_modelAAA_l1_WR, file = "dgpANN_modelAAA_l1_WR.Rdata")
save(dgpANN_modelAAdA_l1_WR, file = "dgpANN_modelAAdA_l1_WR.Rdata")
save(dgpANN_modelAAdN_l1_WR, file = "dgpANN_modelAAdN_l1_WR.Rdata")
save(dgpANN_modelMAM_l1_WR, file = "dgpANN_modelMAM_l1_WR.Rdata")

## DGP: AAN
# Without weights
start.timeUR <- Sys.time()
system.time({dgpAAN_modelANN_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAN", nobs = c(28, 420), freq = 7, model = "ANN", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAN_modelAAN_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAN", nobs = c(28, 420), freq = 7, model = "AAN", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAN_modelANA_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAN", nobs = c(28, 420), freq = 7, model = "ANA", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAN_modelAAA_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAN", nobs = c(28, 420), freq = 7, model = "AAA", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAN_modelAAdA_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAN", nobs = c(28, 420), freq = 7, model = "AAdA", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAN_modelAAdN_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAN", nobs = c(28, 420), freq = 7, model = "AAdN", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAN_modelMAM_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAN", nobs = c(28, 420), freq = 7, model = "MAM", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})

# With weights
system.time({dgpAAN_modelANN_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAN", nobs = c(28, 420), freq = 7, model = "ANN", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.1), weighted = TRUE))})
system.time({dgpAAN_modelAAN_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAN", nobs = c(28, 420), freq = 7, model = "AAN", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.5, 0.5), weighted = TRUE))})
system.time({dgpAAN_modelANA_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAN", nobs = c(28, 420), freq = 7, model = "ANA", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.5, 0.5), weighted = TRUE))})
system.time({dgpAAN_modelAAA_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAN", nobs = c(28, 420), freq = 7, model = "AAA", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.33, 0.33, 0.33), weighted = TRUE))})
system.time({dgpAAN_modelAAdA_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAN", nobs = c(28, 420), freq = 7, model = "AAdA", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1, 0.25, 0.25, 0.25, 0.25), weighted = TRUE))})
system.time({dgpAAN_modelAAdN_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAN", nobs = c(28, 420), freq = 7, model = "AAdN", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1, 0.33, 0.33, 0.33), weighted = TRUE))})
system.time({dgpAAN_modelMAM_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAN", nobs = c(28, 420), freq = 7, model = "MAM", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.33, 0.33, 0.33), weighted = TRUE))})
end.timeWR <- Sys.time()
computetime2 <- end.timeWR - start.timeUR

save(dgpAAN_modelANN_l1_UR, file = "dgpAAN_modelANN_l1_UR.Rdata")
save(dgpAAN_modelAAN_l1_UR, file = "dgpAAN_modelAAN_l1_UR.Rdata")
save(dgpAAN_modelANA_l1_UR, file = "dgpAAN_modelANA_l1_UR.Rdata")
save(dgpAAN_modelAAA_l1_UR, file = "dgpAAN_modelAAA_l1_UR.Rdata")
save(dgpAAN_modelAAdA_l1_UR, file = "dgpAAN_modelAAdA_l1_UR.Rdata")
save(dgpAAN_modelAAdN_l1_UR, file = "dgpAAN_modelAAdN_l1_UR.Rdata")
save(dgpAAN_modelMAM_l1_UR, file = "dgpAAN_modelMAM_l1_UR.Rdata")

save(dgpAAN_modelANN_l1_WR, file = "dgpAAN_modelANN_l1_WR.Rdata")
save(dgpAAN_modelAAN_l1_WR, file = "dgpAAN_modelAAN_l1_WR.Rdata")
save(dgpAAN_modelANA_l1_WR, file = "dgpAAN_modelANA_l1_WR.Rdata")
save(dgpAAN_modelAAA_l1_WR, file = "dgpAAN_modelAAA_l1_WR.Rdata")
save(dgpAAN_modelAAdA_l1_WR, file = "dgpAAN_modelAAdA_l1_WR.Rdata")
save(dgpAAN_modelAAdN_l1_WR, file = "dgpAAN_modelAAdN_l1_WR.Rdata")
save(dgpAAN_modelMAM_l1_WR, file = "dgpAAN_modelMAM_l1_WR.Rdata")

## DGP: ANA
# Without weights
start.timeUR <- Sys.time()
system.time({dgpANA_modelANN_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANA", nobs = c(28, 420), freq = 7, model = "ANN", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})
system.time({dgpANA_modelAAN_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANA", nobs = c(28, 420), freq = 7, model = "AAN", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})
system.time({dgpANA_modelANA_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANA", nobs = c(28, 420), freq = 7, model = "ANA", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})
system.time({dgpANA_modelAAA_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANA", nobs = c(28, 420), freq = 7, model = "AAA", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})
system.time({dgpANA_modelAAdA_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANA", nobs = c(28, 420), freq = 7, model = "AAdA", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1), weighted = FALSE))})
system.time({dgpANA_modelAAdN_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANA", nobs = c(28, 420), freq = 7, model = "AAdN", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1), weighted = FALSE))})
system.time({dgpANA_modelMAM_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANA", nobs = c(28, 420), freq = 7, model = "MAM", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})

# With weights
system.time({dgpANA_modelANN_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANA", nobs = c(28, 420), freq = 7, model = "ANN", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.1), weighted = TRUE))})
system.time({dgpANA_modelAAN_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANA", nobs = c(28, 420), freq = 7, model = "AAN", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.5, 0.5), weighted = TRUE))})
system.time({dgpANA_modelANA_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANA", nobs = c(28, 420), freq = 7, model = "ANA", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.5, 0.5), weighted = TRUE))})
system.time({dgpANA_modelAAA_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANA", nobs = c(28, 420), freq = 7, model = "AAA", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.33, 0.33, 0.33), weighted = TRUE))})
system.time({dgpANA_modelAAdA_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANA", nobs = c(28, 420), freq = 7, model = "AAdA", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1, 0.25, 0.25, 0.25, 0.25), weighted = TRUE))})
system.time({dgpANA_modelAAdN_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANA", nobs = c(28, 420), freq = 7, model = "AAdN", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1, 0.33, 0.33, 0.33), weighted = TRUE))})
system.time({dgpANA_modelMAM_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "ANA", nobs = c(28, 420), freq = 7, model = "MAM", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.33, 0.33, 0.33), weighted = TRUE))})
end.timeWR <- Sys.time()
computetime3 <- end.timeWR - start.timeUR

save(dgpANA_modelANN_l1_UR, file = "dgpANA_modelANN_l1_UR.Rdata")
save(dgpANA_modelAAN_l1_UR, file = "dgpANA_modelAAN_l1_UR.Rdata")
save(dgpANA_modelANA_l1_UR, file = "dgpANA_modelANA_l1_UR.Rdata")
save(dgpANA_modelAAA_l1_UR, file = "dgpANA_modelAAA_l1_UR.Rdata")
save(dgpANA_modelAAdA_l1_UR, file = "dgpANA_modelAAdA_l1_UR.Rdata")
save(dgpANA_modelAAdN_l1_UR, file = "dgpANA_modelAAdN_l1_UR.Rdata")
save(dgpANA_modelMAM_l1_UR, file = "dgpANA_modelMAM_l1_UR.Rdata")

save(dgpANA_modelANN_l1_WR, file = "dgpANA_modelANN_l1_WR.Rdata")
save(dgpANA_modelAAN_l1_WR, file = "dgpANA_modelAAN_l1_WR.Rdata")
save(dgpANA_modelANA_l1_WR, file = "dgpANA_modelANA_l1_WR.Rdata")
save(dgpANA_modelAAA_l1_WR, file = "dgpANA_modelAAA_l1_WR.Rdata")
save(dgpANA_modelAAdA_l1_WR, file = "dgpANA_modelAAdA_l1_WR.Rdata")
save(dgpANA_modelAAdN_l1_WR, file = "dgpANA_modelAAdN_l1_WR.Rdata")
save(dgpANA_modelMAM_l1_WR, file = "dgpANA_modelMAM_l1_WR.Rdata")


## DGP: AAA
# Without weights
start.timeUR <- Sys.time()
system.time({dgpAAA_modelANN_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAA", nobs = c(28, 420), freq = 7, model = "ANN", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAA_modelAAN_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAA", nobs = c(28, 420), freq = 7, model = "AAN", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAA_modelANA_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAA", nobs = c(28, 420), freq = 7, model = "ANA", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAA_modelAAA_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAA", nobs = c(28, 420), freq = 7, model = "AAA", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAA_modelAAdA_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAA", nobs = c(28, 420), freq = 7, model = "AAdA", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAA_modelAAdN_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAA", nobs = c(28, 420), freq = 7, model = "AAdN", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAA_modelMAM_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAA", nobs = c(28, 420), freq = 7, model = "MAM", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})

# With weights
system.time({dgpAAA_modelANN_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAA", nobs = c(28, 420), freq = 7, model = "ANN", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.1), weighted = TRUE))})
system.time({dgpAAA_modelAAN_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAA", nobs = c(28, 420), freq = 7, model = "AAN", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.5, 0.5), weighted = TRUE))})
system.time({dgpAAA_modelANA_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAA", nobs = c(28, 420), freq = 7, model = "ANA", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.5, 0.5), weighted = TRUE))})
system.time({dgpAAA_modelAAA_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAA", nobs = c(28, 420), freq = 7, model = "AAA", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.33, 0.33, 0.33), weighted = TRUE))})
system.time({dgpAAA_modelAAdA_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAA", nobs = c(28, 420), freq = 7, model = "AAdA", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1, 0.25, 0.25, 0.25, 0.25), weighted = TRUE))})
system.time({dgpAAA_modelAAdN_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAA", nobs = c(28, 420), freq = 7, model = "AAdN", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1, 0.33, 0.33, 0.33), weighted = TRUE))})
system.time({dgpAAA_modelMAM_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAA", nobs = c(28, 420), freq = 7, model = "MAM", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.33, 0.33, 0.33), weighted = TRUE))})
end.timeWR <- Sys.time()
computetime4 <- end.timeWR - start.timeUR

save(dgpAAA_modelANN_l1_UR, file = "dgpAAA_modelANN_l1_UR.Rdata")
save(dgpAAA_modelAAN_l1_UR, file = "dgpAAA_modelAAN_l1_UR.Rdata")
save(dgpAAA_modelANA_l1_UR, file = "dgpAAA_modelANA_l1_UR.Rdata")
save(dgpAAA_modelAAA_l1_UR, file = "dgpAAA_modelAAA_l1_UR.Rdata")
save(dgpAAA_modelAAdA_l1_UR, file = "dgpAAA_modelAAdA_l1_UR.Rdata")
save(dgpAAA_modelAAdN_l1_UR, file = "dgpAAA_modelAAdN_l1_UR.Rdata")
save(dgpAAA_modelMAM_l1_UR, file = "dgpAAA_modelMAM_l1_UR.Rdata")

save(dgpAAA_modelANN_l1_WR, file = "dgpAAA_modelANN_l1_WR.Rdata")
save(dgpAAA_modelAAN_l1_WR, file = "dgpAAA_modelAAN_l1_WR.Rdata")
save(dgpAAA_modelANA_l1_WR, file = "dgpAAA_modelANA_l1_WR.Rdata")
save(dgpAAA_modelAAA_l1_WR, file = "dgpAAA_modelAAA_l1_WR.Rdata")
save(dgpAAA_modelAAdA_l1_WR, file = "dgpAAA_modelAAdA_l1_WR.Rdata")
save(dgpAAA_modelAAdN_l1_WR, file = "dgpAAA_modelAAdN_l1_WR.Rdata")
save(dgpAAA_modelMAM_l1_WR, file = "dgpAAA_modelMAM_l1_WR.Rdata")

## DGP: AAdA
# Without weights
start.timeUR <- Sys.time()
system.time({dgpAAdA_modelANN_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdA", nobs = c(28, 420), freq = 7, model = "ANN", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAdA_modelAAN_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdA", nobs = c(28, 420), freq = 7, model = "AAN", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAdA_modelANA_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdA", nobs = c(28, 420), freq = 7, model = "ANA", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAdA_modelAAA_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdA", nobs = c(28, 420), freq = 7, model = "AAA", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAdA_modelAAdA_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdA", nobs = c(28, 420), freq = 7, model = "AAdA", 
                                                                                             h = 7, origin = 5, type.loss = "l1", 
                                                                                             tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAdA_modelAAdN_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdA", nobs = c(28, 420), freq = 7, model = "AAdN", 
                                                                                             h = 7, origin = 5, type.loss = "l1", 
                                                                                             tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAdA_modelMAM_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdA", nobs = c(28, 420), freq = 7, model = "MAM", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1), weighted = FALSE))})

# With weights
system.time({dgpAAdA_modelANN_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdA", nobs = c(28, 420), freq = 7, model = "ANN", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1, 0.1), weighted = TRUE))})
system.time({dgpAAdA_modelAAN_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdA", nobs = c(28, 420), freq = 7, model = "AAN", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1, 0.5, 0.5), weighted = TRUE))})
system.time({dgpAAdA_modelANA_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdA", nobs = c(28, 420), freq = 7, model = "ANA", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1, 0.5, 0.5), weighted = TRUE))})
system.time({dgpAAdA_modelAAA_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdA", nobs = c(28, 420), freq = 7, model = "AAA", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1, 0.33, 0.33, 0.33), weighted = TRUE))})
system.time({dgpAAdA_modelAAdA_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdA", nobs = c(28, 420), freq = 7, model = "AAdA", 
                                                                                             h = 7, origin = 5, type.loss = "l1", 
                                                                                             tuning = c(0.1, 0.25, 0.25, 0.25, 0.25), weighted = TRUE))})
system.time({dgpAAdA_modelAAdN_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdA", nobs = c(28, 420), freq = 7, model = "AAdN", 
                                                                                             h = 7, origin = 5, type.loss = "l1", 
                                                                                             tuning = c(0.1, 0.33, 0.33, 0.33), weighted = TRUE))})
system.time({dgpAAdA_modelMAM_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdA", nobs = c(28, 420), freq = 7, model = "MAM", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1, 0.33, 0.33, 0.33), weighted = TRUE))})
end.timeWR <- Sys.time()
computetime5 <- end.timeWR - start.timeUR

save(dgpAAdA_modelANN_l1_UR, file = "dgpAAdA_modelANN_l1_UR.Rdata")
save(dgpAAdA_modelAAN_l1_UR, file = "dgpAAdA_modelAAN_l1_UR.Rdata")
save(dgpAAdA_modelANA_l1_UR, file = "dgpAAdA_modelANA_l1_UR.Rdata")
save(dgpAAdA_modelAAA_l1_UR, file = "dgpAAdA_modelAAA_l1_UR.Rdata")
save(dgpAAdA_modelAAdA_l1_UR, file = "dgpAAdA_modelAAdA_l1_UR.Rdata")
save(dgpAAdA_modelAAdN_l1_UR, file = "dgpAAdA_modelAAdN_l1_UR.Rdata")
save(dgpAAdA_modelMAM_l1_UR, file = "dgpAAdA_modelMAM_l1_UR.Rdata")

save(dgpAAdA_modelANN_l1_WR, file = "dgpAAdA_modelANN_l1_WR.Rdata")
save(dgpAAdA_modelAAN_l1_WR, file = "dgpAAdA_modelAAN_l1_WR.Rdata")
save(dgpAAdA_modelANA_l1_WR, file = "dgpAAdA_modelANA_l1_WR.Rdata")
save(dgpAAdA_modelAAA_l1_WR, file = "dgpAAdA_modelAAA_l1_WR.Rdata")
save(dgpAAdA_modelAAdA_l1_WR, file = "dgpAAdA_modelAAdA_l1_WR.Rdata")
save(dgpAAdA_modelAAdN_l1_WR, file = "dgpAAdA_modelAAdN_l1_WR.Rdata")
save(dgpAAdA_modelMAM_l1_WR, file = "dgpAAdA_modelMAM_l1_WR.Rdata")

## DGP: AAdN
# Without weights
start.timeUR <- Sys.time()
system.time({dgpAAdN_modelANN_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdN", nobs = c(28, 420), freq = 7, model = "ANN", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAdN_modelAAN_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdN", nobs = c(28, 420), freq = 7, model = "AAN", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAdN_modelANA_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdN", nobs = c(28, 420), freq = 7, model = "ANA", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAdN_modelAAA_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdN", nobs = c(28, 420), freq = 7, model = "AAA", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAdN_modelAAdA_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdN", nobs = c(28, 420), freq = 7, model = "AAdA", 
                                                                                             h = 7, origin = 5, type.loss = "l1", 
                                                                                             tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAdN_modelAAdN_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdN", nobs = c(28, 420), freq = 7, model = "AAdN", 
                                                                                             h = 7, origin = 5, type.loss = "l1", 
                                                                                             tuning = c(0.1), weighted = FALSE))})
system.time({dgpAAdN_modelMAM_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdN", nobs = c(28, 420), freq = 7, model = "MAM", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1), weighted = FALSE))})

# With weights
system.time({dgpAAdN_modelANN_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdN", nobs = c(28, 420), freq = 7, model = "ANN", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1, 0.1), weighted = TRUE))})
system.time({dgpAAdN_modelAAN_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdN", nobs = c(28, 420), freq = 7, model = "AAN", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1, 0.5, 0.5), weighted = TRUE))})
system.time({dgpAAdN_modelANA_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdN", nobs = c(28, 420), freq = 7, model = "ANA", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1, 0.5, 0.5), weighted = TRUE))})
system.time({dgpAAdN_modelAAA_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdN", nobs = c(28, 420), freq = 7, model = "AAA", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1, 0.33, 0.33, 0.33), weighted = TRUE))})
system.time({dgpAAdN_modelAAdA_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdN", nobs = c(28, 420), freq = 7, model = "AAdA", 
                                                                                             h = 7, origin = 5, type.loss = "l1", 
                                                                                             tuning = c(0.1, 0.25, 0.25, 0.25, 0.25), weighted = TRUE))})
system.time({dgpAAdN_modelAAdN_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdN", nobs = c(28, 420), freq = 7, model = "AAdN", 
                                                                                             h = 7, origin = 5, type.loss = "l1", 
                                                                                             tuning = c(0.1, 0.33, 0.33, 0.33), weighted = TRUE))})
system.time({dgpAAdN_modelMAM_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "AAdN", nobs = c(28, 420), freq = 7, model = "MAM", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1, 0.33, 0.33, 0.33), weighted = TRUE))})
end.timeWR <- Sys.time()
computetime6 <- end.timeWR - start.timeUR

save(dgpAAdN_modelANN_l1_UR, file = "dgpAAdN_modelANN_l1_UR.Rdata")
save(dgpAAdN_modelAAN_l1_UR, file = "dgpAAdN_modelAAN_l1_UR.Rdata")
save(dgpAAdN_modelANA_l1_UR, file = "dgpAAdN_modelANA_l1_UR.Rdata")
save(dgpAAdN_modelAAA_l1_UR, file = "dgpAAdN_modelAAA_l1_UR.Rdata")
save(dgpAAdN_modelAAdA_l1_UR, file = "dgpAAdN_modelAAdA_l1_UR.Rdata")
save(dgpAAdN_modelAAdN_l1_UR, file = "dgpAAdN_modelAAdN_l1_UR.Rdata")
save(dgpAAdN_modelMAM_l1_UR, file = "dgpAAdN_modelMAM_l1_UR.Rdata")

save(dgpAAdN_modelANN_l1_WR, file = "dgpAAdN_modelANN_l1_WR.Rdata")
save(dgpAAdN_modelAAN_l1_WR, file = "dgpAAdN_modelAAN_l1_WR.Rdata")
save(dgpAAdN_modelANA_l1_WR, file = "dgpAAdN_modelANA_l1_WR.Rdata")
save(dgpAAdN_modelAAA_l1_WR, file = "dgpAAdN_modelAAA_l1_WR.Rdata")
save(dgpAAdN_modelAAdA_l1_WR, file = "dgpAAdN_modelAAdA_l1_WR.Rdata")
save(dgpAAdN_modelAAdN_l1_WR, file = "dgpAAdN_modelAAdN_l1_WR.Rdata")
save(dgpAAdN_modelMAM_l1_WR, file = "dgpAAdN_modelMAM_l1_WR.Rdata")

## DGP: MAM
# Without weights
start.timeUR <- Sys.time()
system.time({dgpMAM_modelANN_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "ANN", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})
system.time({dgpMAM_modelAAN_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "AAN", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})
system.time({dgpMAM_modelANA_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "ANA", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})
system.time({dgpMAM_modelAAA_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "AAA", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})
system.time({dgpMAM_modelAAdA_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "AAdA", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1), weighted = FALSE))})
system.time({dgpMAM_modelAAdN_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "AAdN", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1), weighted = FALSE))})
system.time({dgpMAM_modelMAM_l1_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "MAM", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1), weighted = FALSE))})

# With weights
system.time({dgpMAM_modelANN_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "ANN", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.1), weighted = TRUE))})
system.time({dgpMAM_modelAAN_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "AAN", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.5, 0.5), weighted = TRUE))})
system.time({dgpMAM_modelANA_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "ANA", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.5, 0.5), weighted = TRUE))})
system.time({dgpMAM_modelAAA_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "AAA", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.33, 0.33, 0.33), weighted = TRUE))})
system.time({dgpMAM_modelAAdA_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "AAdA", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1, 0.25, 0.25, 0.25, 0.25), weighted = TRUE))})
system.time({dgpMAM_modelAAdN_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "AAdN", 
                                                                                            h = 7, origin = 5, type.loss = "l1", 
                                                                                            tuning = c(0.1, 0.33, 0.33, 0.33), weighted = TRUE))})
system.time({dgpMAM_modelMAM_l1_WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "MAM", 
                                                                                           h = 7, origin = 5, type.loss = "l1", 
                                                                                           tuning = c(0.1, 0.33, 0.33, 0.33), weighted = TRUE))})
end.timeWR <- Sys.time()
computetime7 <- end.timeWR - start.timeUR

save(dgpMAM_modelANN_l1_UR, file = "dgpMAM_modelANN_l1_UR.Rdata")
save(dgpMAM_modelAAN_l1_UR, file = "dgpMAM_modelAAN_l1_UR.Rdata")
save(dgpMAM_modelANA_l1_UR, file = "dgpMAM_modelANA_l1_UR.Rdata")
save(dgpMAM_modelAAA_l1_UR, file = "dgpMAM_modelAAA_l1_UR.Rdata")
save(dgpMAM_modelAAdA_l1_UR, file = "dgpMAM_modelAAdA_l1_UR.Rdata")
save(dgpMAM_modelAAdN_l1_UR, file = "dgpMAM_modelAAdN_l1_UR.Rdata")
save(dgpMAM_modelMAM_l1_UR, file = "dgpMAM_modelMAM_l1_UR.Rdata")

save(dgpMAM_modelANN_l1_WR, file = "dgpMAM_modelANN_l1_WR.Rdata")
save(dgpMAM_modelAAN_l1_WR, file = "dgpMAM_modelAAN_l1_WR.Rdata")
save(dgpMAM_modelANA_l1_WR, file = "dgpMAM_modelANA_l1_WR.Rdata")
save(dgpMAM_modelAAA_l1_WR, file = "dgpMAM_modelAAA_l1_WR.Rdata")
save(dgpMAM_modelAAdA_l1_WR, file = "dgpMAM_modelAAdA_l1_WR.Rdata")
save(dgpMAM_modelAAdN_l1_WR, file = "dgpMAM_modelAAdN_l1_WR.Rdata")
save(dgpMAM_modelMAM_l1_WR, file = "dgpMAM_modelMAM_l1_WR.Rdata")

save(list(computetime1, computetime2, computetime3, computetime4,
          computetime5, computetime6, computetime7), 
     file = "computetime_l1.RData")

stopCluster(cl)

