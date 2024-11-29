################################# get results #################################
nsim = 1e+5
lag_max = 150





load("logitGaussIS.RData")
estLogitGaussIS = colSums(
  t(logitGaussIS$res$sample) * logitGaussIS$res$normWeights
)
essLogitGaussIS = 1 / sum(logitGaussIS$res$normWeights^2)
esspsLogitGaussIS = essLogitGaussIS / as.numeric(logitGaussIS$time)
timeLogitGaussIS = as.numeric(logitGaussIS$time)
essLogitGaussIS = summary(essLogitGaussIS)



load("logitGauss_pSUN.RData")
estLogitGauss_pSUN = rowMeans(logitGauss_pSUN$res)
essLogitGauss_pSUN = double(length(estLogitGauss_pSUN))
for (iess in 1:length(essLogitGauss_pSUN)) {
  essLogitGauss_pSUN[iess] = sns::ess(logitGauss_pSUN$res[iess, ])
}
esspsLogitGauss_pSUN = summary(
  essLogitGauss_pSUN / as.numeric(logitGauss_pSUN$time)
)
timeLogitGauss_pSUN = as.numeric(logitGauss_pSUN$time)
essLogitGauss_pSUN[essLogitGauss_pSUN > nsim] = nsim
essLogitGauss_pSUN = summary(essLogitGauss_pSUN)



load("logitGaussPG.RData")
estLogitGaussPG = rowMeans(logitGaussPG$res)
essLogitGaussPG = double(length(estLogitGaussPG))
for (iess in 1:length(essLogitGaussPG)) {
  essLogitGaussPG[iess] = sns::ess(logitGaussPG$res[iess, ])
}
esspsLogitGaussPG = summary(
  essLogitGaussPG / as.numeric(logitGaussPG$time)
)
timeLogitGaussPG = as.numeric(logitGaussPG$time)
essLogitGaussPG = summary(essLogitGaussPG)



load("logitGaussUPG.RData")
estLogitGaussUPG = rowMeans(logitGaussUPG$res)
essLogitGaussUPG = double(length(estLogitGaussUPG))
for (iess in 1:length(essLogitGaussUPG)) {
  essLogitGaussUPG[iess] = sns::ess(logitGaussUPG$res[iess, ])
}
esspsLogitGaussUPG = summary(
  essLogitGaussUPG / as.numeric(logitGaussUPG$time)
)
timeLogitGaussUPG = as.numeric(logitGaussUPG$time)
essLogitGaussUPG = summary(essLogitGaussUPG)





load("logitLaplacit_pSUN.RData")
estLogitLaplacit_pSUN = rowMeans(logitLaplacit_pSUN$res)
essLogitLaplacit_pSUN = double(length(estLogitLaplacit_pSUN))
for (iess in 1:length(essLogitLaplacit_pSUN)) {
  essLogitLaplacit_pSUN[iess] = sns::ess(logitLaplacit_pSUN$res[iess, ])
}
esspsLogitLaplacit_pSUN = summary(
  essLogitLaplacit_pSUN / as.numeric(logitLaplacit_pSUN$time)
)
timeLogitLaplacit_pSUN = as.numeric(logitLaplacit_pSUN$time)
essLogitLaplacit_pSUN[essLogitLaplacit_pSUN > nsim] = nsim
essLogitLaplacit_pSUN = summary(essLogitLaplacit_pSUN)



load("logitLaplacitPG.RData")
estLogitLaplacitPG = rowMeans(logitLaplacitPG$res)
essLogitLaplacitPG = double(length(estLogitLaplacitPG))
for (iess in 1:length(essLogitLaplacitPG)) {
  essLogitLaplacitPG[iess] = sns::ess(logitLaplacitPG$res[iess, ])
}
esspsLogitLaplacitPG = summary(
  essLogitLaplacitPG / as.numeric(logitLaplacitPG$time)
)
timeLogitLaplacitPG = as.numeric(logitLaplacitPG$time)
essLogitLaplacitPG = summary(essLogitLaplacitPG)



load("logitLaplacitUPG.RData")
estLogitLaplacitUPG = rowMeans(logitLaplacitUPG$res)
essLogitLaplacitUPG = double(length(estLogitLaplacitUPG))
for (iess in 1:length(essLogitLaplacitUPG)) {
  essLogitLaplacitUPG[iess] = sns::ess(logitLaplacitUPG$res[iess, ])
}
esspsLogitLaplacitUPG = summary(
  essLogitLaplacitUPG / as.numeric(logitLaplacitUPG$time)
)
timeLogitLaplacitUPG = as.numeric(logitLaplacitUPG$time)
essLogitLaplacitUPG = summary(essLogitLaplacitUPG)




load("logitDirLapl_pSUN.RData")
estLogitDirLapl_pSUN = rowMeans(logitDirLapl_pSUN$res)
essLogitDirLapl_pSUN = double(length(estLogitDirLapl_pSUN))
for (iess in 1:length(essLogitDirLapl_pSUN)) {
  essLogitDirLapl_pSUN[iess] = sns::ess(logitDirLapl_pSUN$res[iess, ])
}
esspsLogitDirLapl_pSUN = summary(
  essLogitDirLapl_pSUN / as.numeric(logitDirLapl_pSUN$time)
)
timeLogitDirLapl_pSUN = as.numeric(logitDirLapl_pSUN$time)
essLogitDirLapl_pSUN[essLogitDirLapl_pSUN > nsim] = nsim
essLogitDirLapl_pSUN = summary(essLogitDirLapl_pSUN)



load("logitDirLaplPG.RData")
estLogitDirLaplPG = rowMeans(logitDirLaplPG$res)
essLogitDirLaplPG = double(length(estLogitDirLaplPG))
for (iess in 1:length(essLogitDirLaplPG)) {
  essLogitDirLaplPG[iess] = sns::ess(logitDirLaplPG$res[iess, ])
}
esspsLogitDirLaplPG = summary(
  essLogitDirLaplPG / as.numeric(logitDirLaplPG$time)
)
timeLogitDirLaplPG = as.numeric(logitDirLaplPG$time)
essLogitDirLaplPG = summary(essLogitDirLaplPG)



load("logitDirLaplUPG.RData")
estLogitDirLaplUPG = rowMeans(logitDirLaplUPG$res)
essLogitDirLaplUPG = double(length(estLogitDirLaplUPG))
for (iess in 1:length(essLogitDirLaplUPG)) {
  essLogitDirLaplUPG[iess] = sns::ess(logitDirLaplUPG$res[iess, ])
}
esspsLogitDirLaplUPG = summary(
  essLogitDirLaplUPG / as.numeric(logitDirLaplUPG$time)
)
timeLogitDirLaplUPG = as.numeric(logitDirLaplUPG$time)
essLogitDirLaplUPG = summary(essLogitDirLaplUPG)



load("probitGauss.RData")
estProbitGauss = rowMeans(probitGauss$res$sample)
essProbitGauss = nsim
esspsProbitGauss = essProbitGauss / as.numeric(probitGauss$time)
timeProbitGauss = as.numeric(probitGauss$time)
essProbitGauss = summary(essProbitGauss)
esspsProbitGauss = summary(esspsProbitGauss)



load("probitLaplacit.RData")
estProbitLaplacit = rowMeans(probitLaplacit$res)
essProbitLaplacit = double(length(estProbitLaplacit))
for (iess in 1:length(essProbitLaplacit)) {
  essProbitLaplacit[iess] = sns::ess(probitLaplacit$res[iess, ])
}
esspsProbitLaplacit = essProbitLaplacit / as.numeric(probitLaplacit$time)
timeProbitLaplacit = as.numeric(probitLaplacit$time)
essProbitLaplacit = summary(essProbitLaplacit)
esspsProbitLaplacit = summary(esspsProbitLaplacit)



load("probitDirLapl.RData")
estProbitDirLapl = rowMeans(probitDirLapl$res)
essProbitDirLapl = double(length(estProbitDirLapl))
for (iess in 1:length(essProbitDirLapl)) {
  essProbitDirLapl[iess] = sns::ess(probitDirLapl$res[iess, ])
}
esspsProbitDirLapl = essProbitDirLapl / as.numeric(probitDirLapl$time)
timeProbitDirLapl = as.numeric(probitDirLapl$time)
essProbitDirLapl = summary(essProbitDirLapl)
esspsProbitDirLapl = summary(esspsProbitDirLapl)



load("logitCauchy_pSUN.RData")
load("logitCauchyPG.RData")
load("logitCauchyUPG.RData")
load("probitCauchy.RData")
medCauchyLogit_pSUN = matrixStats::rowMedians(logitCauchy_pSUN$res)
medCauchyLogitPG = matrixStats::rowMedians(logitCauchyPG$res)
medCauchyLogitUPG = matrixStats::rowMedians(logitCauchyUPG$res)
medCauchyProbit = matrixStats::rowMedians(probitCauchy$res)

logitEmpCDF = ecdf(c(
  logitCauchy_pSUN$res[1 ,], logitCauchyPG$res[1, ], logitCauchyUPG$res[1, ]
))
probitEmpCDF = ecdf(probitCauchy$res[1, ])

quantLogitCauchy_pSUN = logitEmpCDF(logitCauchy_pSUN$res[1, ])
quantLogitCauchyPG = logitEmpCDF(logitCauchyPG$res[1, ])
quantLogitCauchyUPG = logitEmpCDF(logitCauchyUPG$res[1, ])
quantProbitCauchy = probitEmpCDF(probitCauchy$res[1, ])

quantLogitCauchy_pSUN = acf(
  quantLogitCauchy_pSUN, lag.max = lag_max, plot = FALSE
)$acf[2:(lag_max + 1)]
quantLogitCauchyPG = acf(
  quantLogitCauchyPG, lag.max = lag_max, plot = FALSE
)$acf[2:(lag_max + 1)]
quantLogitCauchyUPG = acf(
  quantLogitCauchyUPG, lag.max = lag_max, plot = FALSE
)$acf[2:(lag_max + 1)]
quantProbitCauchy = acf(
  quantProbitCauchy, lag.max = lag_max, plot = FALSE
)$acf[2:(lag_max + 1)]


plot(estLogitGaussIS, pch = 4)
points(estLogitGauss_pSUN, col = "red", pch = 4)
points(estLogitGaussPG, col = "blue", pch = 4)
points(estLogitGaussUPG, col = "green4")

plot(estLogitLaplacit_pSUN, col = "red", pch = 4)
points(estLogitLaplacitPG, col = "blue", pch = 4)
points(estLogitLaplacitUPG, col = "green4")

plot(estLogitDirLapl_pSUN, col = "red", pch = 4)
points(estLogitDirLaplPG, col = "blue", pch = 4)
points(estLogitDirLaplUPG, col = "green4")





#################################### plots ####################################

########## Estimates using logit and probit links and different priors
png("CancerSAGE_est.png", width = 750, height = 950, res = 125)
par(mfrow = c(2, 3))
par(mar = c(4.5, 4, 1, 1))

seq_x = c(1, 100, 200, 300, 400, 500)
seq_y = c(-10, -5, 0, 5, 10, 15, 20)

plot(
  estLogitGauss_pSUN, pch = 4, cex = 0.5, ylim = c(-11.5, 21.5),
  xlab = expression(beta), ylab = "Values", xaxt = "n", yaxt = "n"
)
axis(1, seq_x)
axis(2, seq_y)
abline(h = seq_y, lwd = 0.075)
abline(v = seq_x, lwd = 0.075)

plot(
  estLogitLaplacit_pSUN, pch = 4, cex = 0.5, ylim = c(-11.5, 21.5),
  xlab = expression(beta), ylab = "Values", xaxt = "n", yaxt = "n"
)
axis(1, seq_x)
axis(2, seq_y)
abline(h = seq_y, lwd = 0.075)
abline(v = seq_x, lwd = 0.075)

plot(
  estLogitDirLapl_pSUN, pch = 4, cex = 0.5, ylim = c(-11.5, 21.5),
  xlab = expression(beta), ylab = "Values", xaxt = "n", yaxt = "n"
)
axis(1, seq_x)
axis(2, seq_y)
abline(h = seq_y, lwd = 0.075)
abline(v = seq_x, lwd = 0.075)

seq_x = c(1, 100, 200, 300, 400, 500)
seq_y = c(-5, -2.5, 0, 2.5, 5, 7.5, 10)

plot(
  estProbitGauss, pch = 4, cex = 0.5, ylim = c(-6.5, 12),
  xlab = expression(beta), ylab = "Values", xaxt = "n", yaxt = "n"
)
axis(1, seq_x)
axis(2, seq_y, labels = paste(seq_y))
abline(h = seq_y, lwd = 0.075)
abline(v = seq_x, lwd = 0.075)

plot(
  estProbitLaplacit, pch = 4, cex = 0.5, ylim = c(-6.5, 12),
  xlab = expression(beta), ylab = "Values", xaxt = "n", yaxt = "n"
)
axis(1, seq_x)
axis(2, seq_y, labels = paste(seq_y))
abline(h = seq_y, lwd = 0.075)
abline(v = seq_x, lwd = 0.075)

plot(
  estProbitDirLapl, pch = 4, cex = 0.5, ylim = c(-6.5, 12),
  xlab = expression(beta), ylab = "Values", xaxt = "n", yaxt = "n"
)
axis(1, seq_x)
axis(2, seq_y, labels = paste(seq_y))
abline(h = seq_y, lwd = 0.075)
abline(v = seq_x, lwd = 0.075)

dev.off()




########## Estimates using Dirichlet-Laplace prior using different algorithms
png("CancerSAGE_logitEst_DirLapl.png", width = 750, height = 500, res = 125)
par(mfrow = c(1, 3))
par(mar = c(4.5, 4, 1, 1))
seq_x = c(1, 100, 200, 300, 400, 500)
seq_y = c(-10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10)

plot(
  estLogitDirLaplUPG, pch = 20, cex = 0.5, ylim = c(-11.5, 10), col = "red",
  xlab = expression(beta), ylab = "Values", xaxt = "n", yaxt = "n"
)
axis(1, seq_x)
axis(2, seq_y, labels = paste(seq_y))
abline(h = seq_y, lwd = 0.075)
abline(v = seq_x, lwd = 0.075)

plot(
  estLogitDirLaplPG, pch = 20, cex = 0.5, ylim = c(-11.5, 10), col = "blue",
  xlab = expression(beta), ylab = "Values", xaxt = "n", yaxt = "n"
)
axis(1, seq_x)
axis(2, seq_y, labels = paste(seq_y))
abline(h = seq_y, lwd = 0.075)
abline(v = seq_x, lwd = 0.075)

plot(
  estLogitDirLapl_pSUN, pch = 4, cex = 0.5, ylim = c(-11.5, 10),
  xlab = expression(beta), ylab = "Values", xaxt = "n", yaxt = "n"
)
axis(1, seq_x)
axis(2, seq_y, labels = paste(seq_y))
abline(h = seq_y, lwd = 0.075)
abline(v = seq_x, lwd = 0.075)

dev.off()




########## Median using Cauchy prior with logit and probit links

##### Logit
png("CancerSAGE_medLogitCauchy.png", width = 750, height = 500, res = 125)
par(mfrow = c(1, 3))
par(mar = c(4.5, 4, 1, 1))

seq_x = c(1, 100, 200, 300, 400, 500)
seq_y = c(-1, 0, 1, 2, 3, 4)

plot(
  medCauchyLogitUPG, pch = 20, cex = 0.5, ylim = c(-2, 5), col = "red",
  xlab = expression(beta), ylab = "Values", xaxt = "n", yaxt = "n"
)
axis(1, seq_x)
axis(2, seq_y)
abline(h = seq_y, lwd = 0.075)
abline(v = seq_x, lwd = 0.075)

plot(
  medCauchyLogitPG, pch = 20, cex = 0.5, ylim = c(-2, 5), col = "blue",
  xlab = expression(beta), ylab = "Values", xaxt = "n", yaxt = "n"
)
axis(1, seq_x)
axis(2, seq_y)
abline(h = seq_y, lwd = 0.075)
abline(v = seq_x, lwd = 0.075)

plot(
  medCauchyLogit_pSUN, pch = 4, cex = 0.5, ylim = c(-2, 5),
  xlab = expression(beta), ylab = "Values", xaxt = "n", yaxt = "n"
)
axis(1, seq_x)
axis(2, seq_y)
abline(h = seq_y, lwd = 0.075)
abline(v = seq_x, lwd = 0.075)

dev.off()



##### Probit
png("CancerSAGE_medProbitCauchy.png", res = 125)
par(mfrow = c(1, 1))
par(mar = c(4.5, 4, 1, 1))

seq_x = c(1, 100, 200, 300, 400, 500)
seq_y = c(-0.5, 0, 0.5, 1, 1.5, 2)

plot(
  medCauchyProbit, pch = 4, cex = 0.5, ylim = c(-1, 2.5), col = "black",
  xlab = expression(beta), ylab = "Values", xaxt = "n", yaxt = "n"
)
axis(1, seq_x)
axis(2, seq_y, labels = paste(seq_y))
abline(h = seq_y, lwd = 0.075)
abline(v = seq_x, lwd = 0.075)

dev.off()




########## Cauchy case: a sort of ranks autocorrelation function
png("CancerSAGE_RACF_cauchy.png", width = 750, res = 125)
par(mfrow = c(1, 2))
par(mar = c(4.5, 4, 1, 1))

seq_lag = c(1, 50, 100, 150)

plot(
  x = c(1:lag_max), y = quantLogitCauchy_pSUN,
  type = "l", lty = 1, lwd = 2,
  xlab = "Lag", ylab = "ACF", ylim = c(0, 1), xaxt = "n"
)
points(
  x = c(1:lag_max), y = quantLogitCauchyPG,
  type = "l", lty = 2, lwd = 2, col = "blue"
)
points(
  x = c(1:lag_max), y = quantLogitCauchyUPG,
  type = "l", lty = 3, lwd = 2, col = "red"
)
axis(1, at = seq_lag)
#abline(h = 1.96 / sqrt(nsim), lty = 5, lwd = 2, col = "darkgray")
#abline(h = -1.96 / sqrt(nsim), lty = 5, lwd = 2, col = "darkgray")
abline(h = c(0:10) / 10, lwd = 0.075)
abline(v = c(1, 25, 50, 75, 100, 125, 150), lwd = 0.075)

plot(
  x = c(1:lag_max), y = quantProbitCauchy,
  type = "l", lty = 1, lwd = 2,
  xlab = "Lag", ylab = "ACF", ylim = c(0, 1), xaxt = "n"
)
axis(1, at = seq_lag)
#abline(h = 1.96 / sqrt(nsim), lty = 5, lwd = 2, col = "darkgray")
#abline(h = -1.96 / sqrt(nsim), lty = 5, lwd = 2, col = "darkgray")
abline(h = c(0:10) / 10, lwd = 0.075)
abline(v = c(1, 25, 50, 75, 100, 125, 150), lwd = 0.075)

dev.off()





#################################### tables ###################################

########## effective sample size per second logit case
esspsLogit = matrix(nrow = 10, ncol = 7)

esspsLogit[1, -1] = esspsLogitGaussUPG
esspsLogit[2, -1] = esspsLogitGaussPG
esspsLogit[3, -1] = esspsLogitGauss_pSUN
esspsLogit[4, -1] = summary(esspsLogitGaussIS)

esspsLogit[5, -1] = esspsLogitLaplacitUPG
esspsLogit[6, -1] = esspsLogitLaplacitPG
esspsLogit[7, -1] = esspsLogitLaplacit_pSUN

esspsLogit[8, -1] = esspsLogitDirLaplUPG
esspsLogit[9, -1] = esspsLogitDirLaplPG
esspsLogit[10, -1] = esspsLogitDirLapl_pSUN

xtable::xtable(esspsLogit, digits = 2)




########## effective sample size per second probit case
esspsProbit = matrix(nrow = 3, ncol = 6)

esspsProbit[1, ] = esspsProbitGauss
esspsProbit[2, ] = esspsProbitLaplacit
esspsProbit[3, ] = esspsProbitDirLapl

xtable::xtable(esspsProbit, digits = 2)




########## effective sample size logit case
essLogit = matrix(nrow = 10, ncol = 7)

essLogit[1, -1] = essLogitGaussUPG
essLogit[2, -1] = essLogitGaussPG
essLogit[3, -1] = essLogitGauss_pSUN
essLogit[4, -1] = essLogitGaussIS

essLogit[5, -1] = essLogitLaplacitUPG
essLogit[6, -1] = essLogitLaplacitPG
essLogit[7, -1] = essLogitLaplacit_pSUN

essLogit[8, -1] = essLogitDirLaplUPG
essLogit[9, -1] = essLogitDirLaplPG
essLogit[10, -1] = essLogitDirLapl_pSUN

xtable::xtable(essLogit, digits = 2)




########## effective sample size probit case
essProbit = matrix(nrow = 3, ncol = 6)

essProbit[1, ] = essProbitGauss
essProbit[2, ] = essProbitLaplacit
essProbit[3, ] = essProbitDirLapl

xtable::xtable(essProbit, digits = 2)




########## computational times in seconds
times = matrix(nrow = 13, ncol = 3)

times[1, 3] = logitGaussUPG$time
times[2, 3] = logitGaussPG$time
times[3, 3] = logitGauss_pSUN$time
times[4, 3] = logitGaussIS$time
times[5, 3] = logitLaplacitUPG$time
times[6, 3] = logitLaplacitPG$time
times[7, 3] = logitLaplacit_pSUN$time
times[8, 3] = logitDirLaplUPG$time
times[9, 3] = logitDirLaplPG$time
times[10, 3] = logitDirLapl_pSUN$time
times[11, 3] = probitGauss$time
times[12, 3] = probitLaplacit$time
times[13, 3] = probitDirLapl$time

xtable::xtable(times, digits = 2)
