set.seed(020193)
collect.ts <- array(NA, c(420, 500))
for (i in 1:500) {
  
  yts <- ts(genseries(dgp = "MAM", n = 420, freq = 7), frequency = 7)
  while (min(yts) < 0) {
    yts <- ts(genseries(dgp = "MAM", n = 420, freq = 7), frequency = 7)
  }
  
  collect.ts[,i] <- yts
  
}

colSums(apply(t(apply(collect.ts, 2, summary)), 2, function(x) x < 0))

y <- ts(collect.ts[,round(runif(1, min = 1, max = 500))], frequency = 7)
plot(y, type = "l", main = "60 weeks, daily data")

collect.fcst <- array(NA, c(7, 5))
collect.yTest <- array(NA, c(7, 5))
for (j in 1:5) {
  
  yTrain <- ts(y[1:(405+j)], frequency = 7)
  yTest <- ts(y[(406+j):(412+j)], frequency = 7)
  fit <<- adam(yTrain, model = "MAM", loss = loss.rmse, lags = 7)
  fcst <- forecast(fit, h = 7)$mean
  
  collect.fcst[,j] <- fcst
  collect.yTest[,j] <- yTest
  
  plot(y, type = "l", main = paste("Origin", j))
  lines(ts(fit$fitted, frequency = 7), col = "red", type = "l", lty = 2)
  lines(fcst, col = "blue", type = "l", lty = 2)
}

