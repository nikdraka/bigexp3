## Running all experiments for L2 regularisation, both unweighted and weighted

## Packages to run
library(nloptr)
library(parallel)
library(snow)
library(clusterGeneration)
library(rlecuyer)
library(smooth)

## Function to run
## Function to run
# generate time series
genseries <- function(dgp = c("ANN", "AAN", "AAdN","ANA", "AAA"), n = 100, freq = 7) {
  
  dgp <- dgp[1]
  n <- n[1]
  freq <- freq
  
  burnin <- 0
  
  # 1. Generating the time series, season initial values are generated randomly
  if (dgp == "ANN") {
    input.persist <- 0.4
    input.initial <- c(200)
    y <- sim.es(model = dgp, obs = burnin+n, nsim = 1, frequency = freq,
                bounds = "usual", persistence = input.persist,
                initial = input.initial)$data[(burnin+1):(burnin+n)]
  } else if (dgp == "AAN") {
    input.persist <- c(0.4, 0.3)
    input.initial <- c(200, 0.5)
    y <- sim.es(model = dgp, obs = burnin+n, nsim = 1, frequency = freq,
                bounds = "usual", persistence = input.persist,
                initial = input.initial)$data[(burnin+1):(burnin+n)]
  } else if (dgp == "ANA") {
    input.persist <- c(0.4, 0.1)
    input.initial <- c(200)
    y <- sim.es(model = dgp, obs = burnin+n, nsim = 1, frequency = freq,
                bounds = "usual", persistence = input.persist,
                initial = input.initial)$data[(burnin+1):(burnin+n)]
  } else if (dgp == "AAA") {
    input.persist <- c(0.4, 0.3, 0.1)
    input.initial <- c(200, 0.5)
    y <- sim.es(model = dgp, obs = burnin+n, nsim = 1, frequency = freq,
                bounds = "usual", persistence = input.persist,
                initial = input.initial)$data[(burnin+1):(burnin+n)]
  } else if (is.null(dgp)) {
    input.persist <- 0.4
    input.initial <- c(200)
    y <- sim.es(model = dgp, obs = burnin+n, nsim = 1, frequency = freq,
                bounds = "usual", persistence = input.persist,
                initial = input.initial)$data[(burnin+1):(burnin+n)]
  } else if (dgp == "AAdN") {
    input.persist <- c(0.4, 0.3)
    input.initial <- c(200, 0.5)
    y <- sim.es(model = dgp, obs = burnin+n, nsim = 1, frequency = freq,
                bounds = "usual", persistence = input.persist, phi = 0.94,
                initial = input.initial)$data[(burnin+1):(burnin+n)]
  } else if (dgp == "AAdA") {
    input.persist <- c(0.4, 0.3, 0.1)
    input.initial <- c(100, 0.5)
    y <- sim.es(model = dgp, obs = burnin+n, nsim = 1, frequency = freq,
                bounds = "usual", persistence = input.persist, phi = 0.94,
                initial = input.initial)$data[(burnin+1):(burnin+n)]
  } else if (dgp == "MAM") {
    input.persist <- c(0.1, 0.075, 0.05)
    input.initial <- c(200, 5)
    y <- sim.es(model = dgp, obs = n, nsim = 1, frequency = freq,
                bounds = "usual", persistence = input.persist,
                initial = input.initial, 
                randomizer = "rnorm", mean = 1, sd = 0.0075)$data
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
    B[6:nB] <- log(B[6:nB]) # season
    
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
    B[6:nB] <- log(B[6:nB]) # season
    
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
  } else {
    
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
      } else {
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
      } else {
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
  } else {
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
  if (dgp == "MAM" || mdl == "MAM") {
    y1 <-  ts(genseries(dgp = dgp, n = nobs1, freq = freq), frequency = freq)
    while (sum((y1 > 0)) != length(y1)) {
      y1 <- ts(genseries(dgp = dgp, n = nobs1, freq = freq), frequency = freq)
    }
  } else {
    y1 <-  ts(genseries(dgp = dgp, n = nobs1, freq = freq), frequency = freq)
  }
  
  model1 <- compile.ets(data = y1, model = mdl, origin = origin, h = h,
                        type.loss = type.loss, tuning = tuning, 
                        weighted = weighted, year.insample = 2)
  
  ## Large sample size
  if (dgp == "MAM" || mdl == "MAM") {
    y2 <-  ts(genseries(dgp = dgp, n = nobs2, freq = freq), frequency = freq)
    while (sum((y2 > 0)) != length(y2)) {
      y2 <- ts(genseries(dgp = dgp, n = nobs2, freq = freq), frequency = freq)
    }
  } else {
    y2 <-  ts(genseries(dgp = dgp, n = nobs2, freq = freq), frequency = freq)
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
    
    acc[1,1] <- mean(sqrt((error[,1])^2))
    acc[2,1] <- mean(abs(error[,1]))
    acc[3,1] <- mean(abs(error[,1]))
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
      
      mat.pinball[i,j] <- pinball(test[i,1:j], fcst[[i]]$upper[1:j,col.level], level+(1-level)/2, loss = 1)
      
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

crs <- detectCores()
cl <- makeCluster(getOption("cl.cores", crs))
writeLines(paste("Running with", crs, 'cores'))
# Load packages to cluster
invisible(clusterCall(cl, function(pkgs) {
  library(clusterGeneration)
  library(smooth)
  library(nloptr)
  library(snow)
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

clusterSetupRNG(cl, seed = 020193)

runs <- 500

## DGP: MAM
# Without weights
dgpMAM_modelANN_l2_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "ANN", 
                                                                              h = 7, origin = 5, type.loss = "l2", 
                                                                              tuning = c(0.1), weighted = FALSE))
dgpMAM_modelAAN_l2_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "AAN", 
                                                                              h = 7, origin = 5, type.loss = "l2", 
                                                                              tuning = c(0.1), weighted = FALSE))
dgpMAM_modelANA_l2_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "ANA", 
                                                                              h = 7, origin = 5, type.loss = "l2", 
                                                                              tuning = c(0.1), weighted = FALSE))
dgpMAM_modelAAA_l2_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "AAA", 
                                                                              h = 7, origin = 5, type.loss = "l2", 
                                                                              tuning = c(0.1), weighted = FALSE))
dgpMAM_modelAAdA_l2_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "AAdA", 
                                                                              h = 7, origin = 5, type.loss = "l2", 
                                                                              tuning = c(0.1), weighted = FALSE))
dgpMAM_modelAAdN_l2_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "AAdN", 
                                                                              h = 7, origin = 5, type.loss = "l2", 
                                                                              tuning = c(0.1), weighted = FALSE))
dgpMAM_modelMAM_l2_UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "MAM", 
                                                                              h = 7, origin = 5, type.loss = "l2", 
                                                                              tuning = c(0.1), weighted = FALSE))

save(list = c("dgpMAM_modelANN_l2_UR", "dgpMAM_modelAAN_l2_UR", "dgpMAM_modelANA_l2_UR", "dgpMAM_modelAAA_l2_UR",
              "dgpMAM_modelAAdA_l2_UR", "dgpMAM_modelAAdN_l2_UR", "dgpMAM_modelMAM_l2_UR"),
     file = "L2_UR_dgpMAM.Rdata")

stopCluster(cl)

