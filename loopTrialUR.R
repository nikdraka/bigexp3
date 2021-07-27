set.seed(020193)

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
start.timeUR <- Sys.time()
system.time({trial1UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "ANN", 
                                                                            h = 7, origin = 5, type.loss = "l2", 
                                                                            tuning = c(0.33), weighted = FALSE))})
Sys.time()
Sys.time()
system.time({trial2UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "AAN", 
                                                                            h = 7, origin = 5, type.loss = "l2", 
                                                                            tuning = c(0.33), weighted = FALSE))})
Sys.time()
Sys.time()
system.time({trial3UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "ANA", 
                                                                            h = 7, origin = 5, type.loss = "l2", 
                                                                            tuning = c(0.33), weighted = FALSE))})
Sys.time()
Sys.time()
system.time({trial4UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "AAA", 
                                                                            h = 7, origin = 5, type.loss = "l2", 
                                                                            tuning = c(0.33), weighted = FALSE))})
Sys.time()
Sys.time()
system.time({trial5UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "AAdA", 
                                                                            h = 7, origin = 5, type.loss = "l2", 
                                                                            tuning = c(0.33), weighted = FALSE))})
Sys.time()
Sys.time()
system.time({trial6UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "AAdN", 
                                                                            h = 7, origin = 5, type.loss = "l2", 
                                                                            tuning = c(0.33), weighted = FALSE))})
Sys.time()
Sys.time()
system.time({trial7UR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "MAM", 
                                                                            h = 7, origin = 5, type.loss = "l2", 
                                                                            tuning = c(0.33), weighted = FALSE))})
end.timeUR <- Sys.time()
computetimeUR <- end.timeUR-start.timeUR
stopCluster(cl)
