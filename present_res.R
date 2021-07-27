dgpANN_modelANN_UR_l2[[1]]$sample002$result.reg


dgpANN_UR_l2 <- list(dgpANN_modelANN_UR_l2,
                     dgpANN_modelAAN_UR_l2,
                     dgpANN_modelANA_UR_l2,
                     dgpANN_modelAAA_UR_l2,
                     dgpANN_modelAAdN_UR_l2,
                     dgpANN_modelAAdA_UR_l2,
                     dgpANN_modelMAM_UR_l2)

dgpAAN_UR_l2 <- list(dgpAAN_modelANN_UR_l2,
                     dgpAAN_modelAAN_UR_l2,
                     dgpAAN_modelANA_UR_l2,
                     dgpAAN_modelAAA_UR_l2,
                     dgpAAN_modelAAdN_UR_l2,
                     dgpAAN_modelAAdA_UR_l2,
                     dgpAAN_modelMAM_UR_l2)

dgpANA_UR_l2 <- list(dgpANA_modelANN_UR_l2,
                     dgpANA_modelAAN_UR_l2,
                     dgpANA_modelANA_UR_l2,
                     dgpANA_modelAAA_UR_l2,
                     dgpANA_modelAAdN_UR_l2,
                     dgpANA_modelAAdA_UR_l2,
                     dgpANA_modelMAM_UR_l2)

dgpAAdA_UR_l2 <- list(dgpAAdA_modelANN_UR_l2,
                     dgpAAdA_modelAAN_UR_l2,
                     dgpAAdA_modelANA_UR_l2,
                     dgpAAdA_modelAAA_UR_l2,
                     dgpAAdA_modelAAdN_UR_l2,
                     dgpAAdA_modelAAdA_UR_l2,
                     dgpAAdA_modelMAM_UR_l2)

dgpMAM_UR_l2 <- list(dgpMAM_modelANN_UR_l2,
                     dgpMAM_modelAAN_UR_l2,
                     dgpMAM_modelANA_UR_l2,
                     dgpMAM_modelAAA_UR_l2,
                     dgpMAM_modelAAdN_UR_l2,
                     dgpMAM_modelAAdA_UR_l2,
                     dgpMAM_modelMAM_UR_l2)

## Small sample
# RMSE
rmse.dgpANN_UR_l2 <- sapply(dgpANN_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg[1,1]-x$sample002$result.bm[1,1])/x$sample002$result.bm[1,1]))
rmse.dgpAAN_UR_l2 <- sapply(dgpAAN_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg[1,1]-x$sample002$result.bm[1,1])/x$sample002$result.bm[1,1]))
rmse.dgpANA_UR_l2 <- sapply(dgpANA_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg[1,1]-x$sample002$result.bm[1,1])/x$sample002$result.bm[1,1]))
rmse.dgpAAdA_UR_l2 <- sapply(dgpAAdA_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg[1,1]-x$sample002$result.bm[1,1])/x$sample002$result.bm[1,1]))
rmse.dgpMAM_UR_l2 <- sapply(dgpMAM_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg[1,1]-x$sample002$result.bm[1,1])/x$sample002$result.bm[1,1]))

name.model <- c("ANN)", "AAN)", "ANA)", "AAA)", "AAdN)", "AAdA)", "MAM)")
matresult.rmse <- array(NA, c(7,7), dimnames = list(paste0("DGP: ETS(", name.model), paste0("MDL: ETS(", name.model)))
matresult.rmse["DGP: ETS(ANN)",] <- colMeans(rmse.dgpANN_UR_l2)*100
matresult.rmse["DGP: ETS(AAN)",] <- colMeans(rmse.dgpAAN_UR_l2)*100
matresult.rmse["DGP: ETS(ANA)",] <- colMeans(rmse.dgpANA_UR_l2)*100
matresult.rmse["DGP: ETS(AAdA)",] <- colMeans(rmse.dgpAAdA_UR_l2)*100
matresult.rmse["DGP: ETS(MAM)",] <- colMeans(rmse.dgpMAM_UR_l2, na.rm = TRUE)*100

# MdAE
mdae.dgpANN_UR_l2 <- sapply(dgpANN_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["MdAE",1]-x$sample002$result.bm["MdAE",1])/x$sample002$result.bm["MdAE",1]))
mdae.dgpAAN_UR_l2 <- sapply(dgpAAN_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["MdAE",1]-x$sample002$result.bm["MdAE",1])/x$sample002$result.bm["MdAE",1]))
mdae.dgpANA_UR_l2 <- sapply(dgpANA_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["MdAE",1]-x$sample002$result.bm["MdAE",1])/x$sample002$result.bm["MdAE",1]))
mdae.dgpAAdA_UR_l2 <- sapply(dgpAAdA_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["MdAE",1]-x$sample002$result.bm["MdAE",1])/x$sample002$result.bm["MdAE",1]))
mdae.dgpMAM_UR_l2 <- sapply(dgpMAM_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["MdAE",1]-x$sample002$result.bm["MdAE",1])/x$sample002$result.bm["MdAE",1]))

name.model <- c("ANN)", "AAN)", "ANA)", "AAA)", "AAdN)", "AAdA)", "MAM)")
matresult.mdae <- array(NA, c(7,7), dimnames = list(paste0("DGP: ETS(", name.model), paste0("MDL: ETS(", name.model)))
matresult.mdae["DGP: ETS(ANN)",] <- colMeans(mdae.dgpANN_UR_l2)*100
matresult.mdae["DGP: ETS(AAN)",] <- colMeans(mdae.dgpAAN_UR_l2)*100
matresult.mdae["DGP: ETS(ANA)",] <- colMeans(mdae.dgpANA_UR_l2)*100
matresult.mdae["DGP: ETS(AAdA)",] <- colMeans(mdae.dgpAAdA_UR_l2)*100
matresult.mdae["DGP: ETS(MAM)",] <- colMeans(mdae.dgpMAM_UR_l2)*100
matresult.mdae

# ME
me.dgpANN_UR_l2 <- sapply(dgpANN_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg[4,1]-x$sample002$result.bm[4,1])/x$sample002$result.bm[4,1]))
me.dgpAAN_UR_l2 <- sapply(dgpAAN_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg[4,1]-x$sample002$result.bm[4,1])/x$sample002$result.bm[4,1]))
me.dgpANA_UR_l2 <- sapply(dgpANA_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg[4,1]-x$sample002$result.bm[4,1])/x$sample002$result.bm[4,1]))
me.dgpAAdA_UR_l2 <- sapply(dgpAAdA_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg[4,1]-x$sample002$result.bm[4,1])/x$sample002$result.bm[4,1]))
me.dgpMAM_UR_l2 <- sapply(dgpMAM_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg[4,1]-x$sample002$result.bm[4,1])/x$sample002$result.bm[4,1]))

name.model <- c("ANN)", "AAN)", "ANA)", "AAA)", "AAdN)", "AAdA)", "MAM)")
matresult.me <- array(NA, c(7,7), dimnames = list(paste0("DGP: ETS(", name.model), paste0("MDL: ETS(", name.model)))
matresult.me["DGP: ETS(ANN)",] <- colMeans(me.dgpANN_UR_l2)*100
matresult.me["DGP: ETS(AAN)",] <- colMeans(me.dgpAAN_UR_l2)*100
matresult.me["DGP: ETS(ANA)",] <- colMeans(me.dgpANA_UR_l2)*100
matresult.me["DGP: ETS(AAdA)",] <- colMeans(me.dgpAAdA_UR_l2)*100
matresult.me["DGP: ETS(MAM)",] <- colMeans(me.dgpMAM_UR_l2)*100
matresult.me

# RMSSE
rmsse.dgpANN_UR_l2 <- sapply(dgpANN_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["RMSSE",1]-x$sample002$result.bm["RMSSE",1])/x$sample002$result.bm["RMSSE",1]))
rmsse.dgpAAN_UR_l2 <- sapply(dgpAAN_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["RMSSE",1]-x$sample002$result.bm["RMSSE",1])/x$sample002$result.bm["RMSSE",1]))
rmsse.dgpANA_UR_l2 <- sapply(dgpANA_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["RMSSE",1]-x$sample002$result.bm["RMSSE",1])/x$sample002$result.bm["RMSSE",1]))
rmsse.dgpAAdA_UR_l2 <- sapply(dgpAAdA_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["RMSSE",1]-x$sample002$result.bm["RMSSE",1])/x$sample002$result.bm["RMSSE",1]))
rmsse.dgpMAM_UR_l2 <- sapply(dgpMAM_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["RMSSE",1]-x$sample002$result.bm["RMSSE",1])/x$sample002$result.bm["RMSSE",1]))

name.model <- c("ANN)", "AAN)", "ANA)", "AAA)", "AAdN)", "AAdA)", "MAM)")
matresult.rmsse <- array(NA, c(7,7), dimnames = list(paste0("DGP: ETS(", name.model), paste0("MDL: ETS(", name.model)))
matresult.rmsse["DGP: ETS(ANN)",] <- colMeans(rmsse.dgpANN_UR_l2)*100
matresult.rmsse["DGP: ETS(AAN)",] <- colMeans(rmsse.dgpAAN_UR_l2)*100
matresult.rmsse["DGP: ETS(ANA)",] <- colMeans(rmsse.dgpANA_UR_l2)*100
matresult.rmsse["DGP: ETS(AAdA)",] <- colMeans(rmsse.dgpAAdA_UR_l2)*100
matresult.rmsse["DGP: ETS(MAM)",] <- colMeans(rmsse.dgpMAM_UR_l2, na.rm = TRUE)*100
matresult.rmsse

# MdAsE
mdase.dgpANN_UR_l2 <- sapply(dgpANN_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["MdAsE",1]-x$sample002$result.bm["MdAsE",1])/x$sample002$result.bm["MdAsE",1]))
mdase.dgpAAN_UR_l2 <- sapply(dgpAAN_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["MdAsE",1]-x$sample002$result.bm["MdAsE",1])/x$sample002$result.bm["MdAsE",1]))
mdase.dgpANA_UR_l2 <- sapply(dgpANA_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["MdAsE",1]-x$sample002$result.bm["MdAsE",1])/x$sample002$result.bm["MdAsE",1]))
mdase.dgpAAdA_UR_l2 <- sapply(dgpAAdA_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["MdAsE",1]-x$sample002$result.bm["MdAsE",1])/x$sample002$result.bm["MdAsE",1]))
mdase.dgpMAM_UR_l2 <- sapply(dgpMAM_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["MdAsE",1]-x$sample002$result.bm["MdAsE",1])/x$sample002$result.bm["MdAsE",1]))

name.model <- c("ANN)", "AAN)", "ANA)", "AAA)", "AAdN)", "AAdA)", "MAM)")
matresult.mdase <- array(NA, c(7,7), dimnames = list(paste0("DGP: ETS(", name.model), paste0("MDL: ETS(", name.model)))
matresult.mdase["DGP: ETS(ANN)",] <- colMeans(mdase.dgpANN_UR_l2)*100
matresult.mdase["DGP: ETS(AAN)",] <- colMeans(mdase.dgpAAN_UR_l2)*100
matresult.mdase["DGP: ETS(ANA)",] <- colMeans(mdase.dgpANA_UR_l2)*100
matresult.mdase["DGP: ETS(AAdA)",] <- colMeans(mdase.dgpAAdA_UR_l2)*100
matresult.mdase["DGP: ETS(MAM)",] <- colMeans(mdase.dgpMAM_UR_l2, na.rm = TRUE)*100
matresult.mdase

# MsE
mse.dgpANN_UR_l2 <- sapply(dgpANN_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["MsE",1]-x$sample002$result.bm["MsE",1])/x$sample002$result.bm["MsE",1]))
mse.dgpAAN_UR_l2 <- sapply(dgpAAN_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["MsE",1]-x$sample002$result.bm["MsE",1])/x$sample002$result.bm["MsE",1]))
mse.dgpANA_UR_l2 <- sapply(dgpANA_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["MsE",1]-x$sample002$result.bm["MsE",1])/x$sample002$result.bm["MsE",1]))
mse.dgpAAdA_UR_l2 <- sapply(dgpAAdA_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["MsE",1]-x$sample002$result.bm["MsE",1])/x$sample002$result.bm["MsE",1]))
mse.dgpMAM_UR_l2 <- sapply(dgpMAM_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["MsE",1]-x$sample002$result.bm["MsE",1])/x$sample002$result.bm["MsE",1]))

name.model <- c("ANN)", "AAN)", "ANA)", "AAA)", "AAdN)", "AAdA)", "MAM)")
matresult.mse <- array(NA, c(7,7), dimnames = list(paste0("DGP: ETS(", name.model), paste0("MDL: ETS(", name.model)))
matresult.mse["DGP: ETS(ANN)",] <- colMeans(mse.dgpANN_UR_l2)*100
matresult.mse["DGP: ETS(AAN)",] <- colMeans(mse.dgpAAN_UR_l2)*100
matresult.mse["DGP: ETS(ANA)",] <- colMeans(mse.dgpANA_UR_l2)*100
matresult.mse["DGP: ETS(AAdA)",] <- colMeans(mse.dgpAAdA_UR_l2)*100
matresult.mse["DGP: ETS(MAM)",] <- colMeans(mse.dgpMAM_UR_l2, na.rm = TRUE)*100
matresult.mse

# MIS_95
mis95.dgpANN_UR_l2 <- sapply(dgpANN_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["MIS_95",1]-x$sample002$result.bm["MIS_95",1])/x$sample002$result.bm["MIS_95",1]))
mis95.dgpAAN_UR_l2 <- sapply(dgpAAN_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["MIS_95",1]-x$sample002$result.bm["MIS_95",1])/x$sample002$result.bm["MIS_95",1]))
mis95.dgpANA_UR_l2 <- sapply(dgpANA_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["MIS_95",1]-x$sample002$result.bm["MIS_95",1])/x$sample002$result.bm["MIS_95",1]))
mis95.dgpAAdA_UR_l2 <- sapply(dgpAAdA_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["MIS_95",1]-x$sample002$result.bm["MIS_95",1])/x$sample002$result.bm["MIS_95",1]))
mis95.dgpMAM_UR_l2 <- sapply(dgpMAM_UR_l2, function(xx) 
  sapply(xx, function(x) (x$sample002$result.reg["MIS_95",1]-x$sample002$result.bm["MIS_95",1])/x$sample002$result.bm["MIS_95",1]))

name.model <- c("ANN)", "AAN)", "ANA)", "AAA)", "AAdN)", "AAdA)", "MAM)")
matresult.mis95 <- array(NA, c(7,7), dimnames = list(paste0("DGP: ETS(", name.model), paste0("MDL: ETS(", name.model)))
matresult.mis95["DGP: ETS(ANN)",] <- colMeans(mis95.dgpANN_UR_l2)*100
matresult.mis95["DGP: ETS(AAN)",] <- colMeans(mis95.dgpAAN_UR_l2)*100
matresult.mis95["DGP: ETS(ANA)",] <- colMeans(mis95.dgpANA_UR_l2)*100
matresult.mis95["DGP: ETS(AAdA)",] <- colMeans(mis95.dgpAAdA_UR_l2)*100
matresult.mis95["DGP: ETS(MAM)",] <- colMeans(mis95.dgpMAM_UR_l2)*100
matresult.mis95

# Parameter: DGP: ETS(AAdA)

param.dgpAAdA.modelAAdA.reg <- t(sapply(dgpAAdA_modelAAdA_UR_l2, function(x) x$sample002$param.reg))
param.dgpAAdA.modelAAdA.bm <- t(sapply(dgpAAdA_modelAAdA_UR_l2, function(x) x$sample002$param.bm))

boxplot(cbind(param.dgpAAdA.modelAAdA.reg[,1], param.dgpAAdA.modelAAdA.bm[,1],
              param.dgpAAdA.modelAAdA.reg[,2], param.dgpAAdA.modelAAdA.bm[,2]))
boxplot(cbind(param.dgpAAdA.modelAAdA.reg[,3], param.dgpAAdA.modelAAdA.bm[,3]), ylim = c(0, 0.02))
boxplot(cbind(param.dgpAAdA.modelAAdA.reg[,4], param.dgpAAdA.modelAAdA.bm[,4]), ylim = c(0.90, 1))

boxplot(cbind(param.dgpAAdA.modelAAdA.reg[,5], param.dgpAAdA.modelAAdA.bm[,5],
              param.dgpAAdA.modelAAdA.reg[,6], param.dgpAAdA.modelAAdA.bm[,6]))

boxplot(cbind(rowSums(param.dgpAAdA.modelAAdA.reg[,7:12]), rowSums(param.dgpAAdA.modelAAdA.bm[,7:12])))
boxplot(cbind(param.dgpAAdA.modelAAdA.reg[,1], param.dgpAAdA.modelAAdA.bm[,1],
              param.dgpAAdA.modelAAdA.reg[,2], param.dgpAAdA.modelAAdA.bm[,2],
              param.dgpAAdA.modelAAdA.reg[,3], param.dgpAAdA.modelAAdA.bm[,3],
              param.dgpAAdA.modelAAdA.reg[,4], param.dgpAAdA.modelAAdA.bm[,4],
              param.dgpAAdA.modelAAdA.reg[,5], param.dgpAAdA.modelAAdA.bm[,5],
              param.dgpAAdA.modelAAdA.reg[,6], param.dgpAAdA.modelAAdA.bm[,6],
              param.dgpAAdA.modelAAdA.reg[,7], param.dgpAAdA.modelAAdA.bm[,7]))

# Tuning parameter
# AAdA v AAdA
tunpar.dgpAAdA.modelAAdA <- sapply(dgpAAdA_modelAAdA_UR_l2, function(x) x$sample002$tunpar)
drmse.dgpAAdA.modelAAdA <- rmse.dgpAAdA_UR_l2[,1]
drmse.dgpAAdA.modelAAdA <- sapply(dgpAAdA_modelAAdA_UR_l2, function(x) x$sample002$result.reg[1,1])

plot(tunpar.dgpAAdA.modelAAdA, drmse.dgpAAdA.modelAAdA,
     pch = 20, col = "grey", bg = "grey",
     ylab = "RMSE from ETS wih Regularisation", xlab = "Tuning Parameter (Lambda)",
     main = "Relationship between Accuracy and Tuning Parameter")
lo <- lowess(drmse.dgpAAdA.modelAAdA~tunpar.dgpAAdA.modelAAdA)
lines(lo, col = "red", lwd = 2)

# AAN v AAN
tunpar.dgpAAN.modelAAN <- sapply(dgpAAN_modelAAN_UR_l2, function(x) x$sample002$tunpar)
drmse.dgpAAN.modelAAN <- rmse.dgpAAN_UR_l2[,1]
drmse.dgpAAN.modelAAN <- sapply(dgpAAN_modelAAN_UR_l2, function(x) x$sample002$result.reg[1,1])

plot(tunpar.dgpAAN.modelAAN, drmse.dgpAAN.modelAAN,
     pch = 20, col = "grey", bg = "grey",
     ylab = "RMSE from ETS wih Regularisation", xlab = "Tuning Parameter (Lambda)",
     main = "Relationship between Accuracy and Tuning Parameter")
lo <- lowess(drmse.dgpAAN.modelAAN~tunpar.dgpAAN.modelAAN)
lines(lo, col = "red", lwd = 2)

# PRESENTATION, in percentage already
round(matresult.rmsse, 2)
round(matresult.mdase, 2)
round(matresult.mse, 2)
round(matresult.mis95, 2)

# DGP: ETS(ANN) --> "overmodelled"
modelname <- c("ANN", "AAN", "ANA", "AAA", "AAdN", "AAdA", "MAM")
boxplot(rmsse.dgpANN_UR_l2*100, ylim = c(-100, 100), names = modelname, 
        col = c("pink", rep("white",5), "brown"), outcex = 0.5, outcol = "grey", pch = 20,
        main = "DGP: ETS(ANN), overmodelled", ylab = "difference in RMSSE in %", xlab = "Model")
lines(apply(mdase.dgpANN_UR_l2*100, 2, mean), col = "red", pch = 19, bg = "red", type = "o")
abline(h = 0, col = "blue")
legend("bottomleft", col = c("pink", "white", "grey", "brown"), fill = c("pink", "white", "grey", "brown"), 
       c("correct model", "overmodelled", "undermodelled", "wrong model"), cex = 0.8, horiz = TRUE)

colnames(rmsse.dgpANN_UR_l2) <- c("ANN", "AAN", "ANA", "AAA", "AAdN", "AAdA", "MAM")
nemenyi(rmsse.dgpANN_UR_l2, plottype = "vmcb", main = "DGP: ETS(ANN)")

# DGP: ETS(AAN) --> "some undermodelled, some overmodelled"
modelname <- c("ANN", "AAN", "ANA", "AAA", "AAdN", "AAdA", "MAM")
boxplot(rmsse.dgpAAN_UR_l2*100, ylim = c(-100, 100), names = modelname, 
        col = c("grey", "pink", rep("white", 4), "brown"), outcex = 0.5, outcol = "grey", pch = 20,
        main = "DGP: ETS(AAN), some undermodelled, some overmodelled", ylab = "difference in RMSSE in %", xlab = "Model")
lines(apply(mdase.dgpAAN_UR_l2*100, 2, mean), col = "red", pch = 19, bg = "red", type = "o")
abline(h = 0, col = "blue")
legend("bottomleft", col = c("pink", "white", "grey", "brown"), fill = c("pink", "white", "grey", "brown"), 
       c("correct model", "overmodelled", "undermodelled", "wrong model"), cex = 0.8, horiz = TRUE)

colnames(rmsse.dgpAAN_UR_l2) <- c("ANN", "AAN", "ANA", "AAA", "AAdN", "AAdA", "MAM")
nemenyi(rmsse.dgpAAN_UR_l2, plottype = "vmcb", main = "DGP: ETS(AAN)")

# DGP: ETS(ANA) --> "some undermodelled, some overmodelled"
modelname <- c("ANN", "AAN", "ANA", "AAA", "AAdN", "AAdA", "MAM")
boxplot(rmsse.dgpANA_UR_l2*100, ylim = c(-100, 100), names = modelname, 
        col = c("grey", "grey", "pink", rep("white", 3), "brown"), outcex = 0.5, outcol = "grey", pch = 20,
        main = "DGP: ETS(ANA), undermodelled, some overmodelled", ylab = "difference in RMSSE in %", xlab = "Model")
lines(apply(mdase.dgpANA_UR_l2*100, 2, mean), col = "red", pch = 19, bg = "red", type = "o")
abline(h = 0, col = "blue")
legend("bottomleft", col = c("pink", "white", "grey", "brown"), fill = c("pink", "white", "grey", "brown"), 
       c("correct model", "overmodelled", "undermodelled", "wrong model"), cex = 0.8, horiz = TRUE)

colnames(rmsse.dgpANA_UR_l2) <- c("ANN", "AAN", "ANA", "AAA", "AAdN", "AAdA", "MAM")
nemenyi(rmsse.dgpANA_UR_l2, plottype = "vmcb", main = "DGP: ETS(ANA)")

# DGP: ETS(AAdA) --> "mostly undermodelled"
modelname <- c("ANN", "AAN", "ANA", "AAA", "AAdN", "AAdA", "MAM")
boxplot(rmsse.dgpAAdA_UR_l2*100, ylim = c(-100, 100), names = modelname, 
        col = c(rep("grey",5), "pink",  "brown"), outcex = 0.5, outcol = "grey", pch = 20,
        main = "DGP: ETS(AAdA), mostly undermodelled", ylab = "difference in RMSSE in %", xlab = "Model")
lines(apply(mdase.dgpAAdA_UR_l2*100, 2, mean), col = "red", pch = 19, bg = "red", type = "o")
abline(h = 0, col = "blue")
legend("bottomleft", col = c("pink", "white", "grey", "brown"), fill = c("pink", "white", "grey", "brown"), 
       c("correct model", "overmodelled", "undermodelled", "wrong model"), cex = 0.8, horiz = TRUE)

colnames(rmsse.dgpAAdA_UR_l2) <- c("ANN", "AAN", "ANA", "AAA", "AAdN", "AAdA", "MAM")
nemenyi(rmsse.dgpAAdA_UR_l2, plottype = "vmcb", main = "DGP: ETS(AAdA), all models")
nemenyi(rmsse.dgpAAdA_UR_l2[,1:6], plottype = "vmcb", main = "DGP: ETS(AAdA), except MAM")


# parameter model: MAM, dgp = ANA
param.dgpANA.modelMAM.URl2.reg <- t(sapply(dgpANA_modelMAM_UR_l2, function(x) c(x$sample002$param.reg[1:5], sum(x$sample002$param.reg[6:11]))))
colnames(param.dgpANA.modelMAM.URl2.reg)[6] <- "sumSeason"

param.dgpANA.modelMAM.URl2.bm <- t(sapply(dgpANA_modelMAM_UR_l2, function(x) c(x$sample002$param.bm[1:5], sum(x$sample002$param.bm[6:11]))))
colnames(param.dgpANA.modelMAM.URl2.bm)[6] <- "sumSeason"

boxplot(
cbind(param.dgpANA.modelMAM.URl2.reg[,1], param.dgpANA.modelMAM.URl2.bm[,1],
      param.dgpANA.modelMAM.URl2.reg[,2], param.dgpANA.modelMAM.URl2.bm[,2],
      param.dgpANA.modelMAM.URl2.reg[,3], param.dgpANA.modelMAM.URl2.bm[,3]), 
names = c("alpha.reg", "alpha.bm", "beta.reg", "beta.bm", "gamma.reg", "gamma.bm"), main = "Smoothing parameters \n DGP: ETS(ANA), Model: ETS(MAM)")

boxplot(
  cbind(param.dgpANA.modelMAM.URl2.reg[,4], param.dgpANA.modelMAM.URl2.bm[,4],
        param.dgpANA.modelMAM.URl2.reg[,5], param.dgpANA.modelMAM.URl2.bm[,5],
        param.dgpANA.modelMAM.URl2.reg[,6], param.dgpANA.modelMAM.URl2.bm[,6]), 
  names = c("lvl.reg", "lvl.bm", "trnd.reg", "trnd.bm", "seas.reg", "seas.bm"), main = "Initial values \n DGP: ETS(ANA), Model: ETS(MAM)")

# parameter model: ANA, dgp = ANA

param.dgpANA.modelANA.URl2.reg <- t(sapply(dgpANA_modelANA_UR_l2, function(x) c(x$sample002$param.reg[1:3], sum(x$sample002$param.reg[4:9]))))
colnames(param.dgpANA.modelANA.URl2.reg)[4] <- "sumSeason"

param.dgpANA.modelANA.URl2.bm <- t(sapply(dgpANA_modelANA_UR_l2, function(x) c(x$sample002$param.bm[1:3], sum(x$sample002$param.bm[4:9]))))
colnames(param.dgpANA.modelANA.URl2.bm)[4] <- "sumSeason"

boxplot(
  cbind(param.dgpANA.modelANA.URl2.reg[,1], param.dgpANA.modelANA.URl2.bm[,1],
        param.dgpANA.modelANA.URl2.reg[,2], param.dgpANA.modelANA.URl2.bm[,2]), 
  names = c("alpha.reg", "alpha.bm", "gamma.reg", "gamma.bm"), main = "Smoothing parameters \n DGP: ETS(ANA), Model: ETS(ANA)")

boxplot(
  cbind(param.dgpANA.modelANA.URl2.reg[,3], param.dgpANA.modelANA.URl2.bm[,3],
        param.dgpANA.modelANA.URl2.reg[,4], param.dgpANA.modelANA.URl2.bm[,4]), 
  names = c("lvl.reg", "lvl.bm", "seas.reg", "seas.bm"), main = "Initial values \n DGP: ETS(ANA), Model: ETS(ANA)")
