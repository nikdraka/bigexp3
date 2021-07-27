set.seed(020193)

crs <- detectCores()
cl <- makeCluster(getOption("cl.cores", crs-1))
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
start.timeWR <- Sys.time()
system.time({trial1WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "ANN", 
                                                                            h = 7, origin = 5, type.loss = "l2", 
                                                                            tuning = c(0.1, 0.1), weighted = TRUE))})
Sys.time()
Sys.time()
system.time({trial2WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "AAN", 
                                                                            h = 7, origin = 5, type.loss = "l2", 
                                                                            tuning = c(0.1, 0.1, 0.1), weighted = TRUE))})
Sys.time()
Sys.time()
system.time({trial3WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "ANA", 
                                                                            h = 7, origin = 5, type.loss = "l2", 
                                                                            tuning = c(0.1, 0.1, 0.1), weighted = TRUE))})
Sys.time()
Sys.time()
system.time({trial4WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "AAA", 
                                                                            h = 7, origin = 5, type.loss = "l2", 
                                                                            tuning = c(0.1, 0.1, 0.1, 0.1), weighted = TRUE))})
Sys.time()
Sys.time()
system.time({trial5WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "AAdA", 
                                                                            h = 7, origin = 5, type.loss = "l2", 
                                                                            tuning = c(0.1, 0.1, 0.1, 0.1, 0.1), weighted = TRUE))})
Sys.time()
Sys.time()
system.time({trial6WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "AAdN", 
                                                                            h = 7, origin = 5, type.loss = "l2", 
                                                                            tuning = c(0.1, 0.1, 0.1, 0.1), weighted = TRUE))})
Sys.time()
Sys.time()
system.time({trial7WR <- clusterApplyLB(cl, 1:runs, function(x) final.compile(dgp = "MAM", nobs = c(28, 420), freq = 7, model = "MAM", 
                                                                            h = 7, origin = 5, type.loss = "l2", 
                                                                            tuning = c(0.1, 0.1, 0.1, 0.1), weighted = TRUE))})
end.timeWR <- Sys.time()
computetimeWR <- end.timeWR - start.timeWR
stopCluster(cl)
