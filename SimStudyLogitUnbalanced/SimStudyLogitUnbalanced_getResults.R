################################# get results #################################

########## ACFs for the different methods and sample sizes
png("logitUnb_ACF.png", width = 750, height = 500, res = 125)
par(mfrow = c(1, 3))
par(mar = c(4.5, 4, 1, 1))
lag_max = 150
seq_lag = c(1, 50, 100, 150)

load(paste(
  "./SimStudyLogitUnbalanced_50n1pGauss",
  "/SimStudyLogitUnbalanced_50n1pGauss.RData",
  sep = ""
))
plot(
  x = c(1:lag_max), y = res[[1]]$acf[2:(lag_max + 1), 1],
  type = "l", lty = 1, lwd = 2,
  xlab = "Lag", ylab = "ACF", ylim = c(0, 1), xaxt = "n"
)
points(
  x = c(1:lag_max), y = res[[2]]$acf[2:(lag_max + 1), 1],
  type = "l", lty = 2, lwd = 2, col = "blue"
)
points(
  x = c(1:lag_max), y = res[[3]]$acf[2:(lag_max + 1), 1],
  type = "l", lty = 3, lwd = 2, col = "red"
)
axis(1, at = seq_lag)
#abline(h = 1.96 / sqrt(nsim), lty = 5, lwd = 2, col = "darkgray")
#abline(h = -1.96 / sqrt(nsim), lty = 5, lwd = 2, col = "darkgray")
abline(h = c(0:10) / 10, lwd = 0.075)
abline(v = c(1, 25, 50, 75, 100, 125, 150), lwd = 0.075)

load(paste(
  "./SimStudyLogitUnbalanced_200n1pGauss",
  "/SimStudyLogitUnbalanced_200n1pGauss.RData",
  sep = ""
))
plot(
  x = c(1:lag_max), y = res[[1]]$acf[2:(lag_max + 1), 1],
  type = "l", lty = 1, lwd = 2,
  xlab = "Lag", ylab = "ACF", ylim = c(0, 1), xaxt = "n"
)
points(
  x = c(1:lag_max), y = res[[2]]$acf[2:(lag_max + 1), 1],
  type = "l", lty = 2, lwd = 2, col = "blue"
)
points(
  x = c(1:lag_max), y = res[[3]]$acf[2:(lag_max + 1), 1],
  type = "l", lty = 3, lwd = 2, col = "red"
)
axis(1, at = seq_lag)
#abline(h = 1.96 / sqrt(nsim), lty = 5, lwd = 2, col = "darkgray")
#abline(h = -1.96 / sqrt(nsim), lty = 5, lwd = 2, col = "darkgray")
abline(h = c(0:10) / 10, lwd = 0.075)
abline(v = c(1, 25, 50, 75, 100, 125, 150), lwd = 0.075)

load(paste(
  "./SimStudyLogitUnbalanced_1000n1pGauss",
  "/SimStudyLogitUnbalanced_1000n1pGauss.RData",
  sep = ""
))
plot(
  x = c(1:lag_max), y = res[[1]]$acf[2:(lag_max + 1), 1],
  type = "l", lty = 1, lwd = 2,
  xlab = "Lag", ylab = "ACF", ylim = c(0, 1), xaxt = "n"
)
points(
  x = c(1:lag_max), y = res[[2]]$acf[2:(lag_max + 1), 1],
  type = "l", lty = 2, lwd = 2, col = "blue"
)
points(
  x = c(1:lag_max), y = res[[3]]$acf[2:(lag_max + 1), 1],
  type = "l", lty = 3, lwd = 2, col = "red"
)
axis(1, at = seq_lag)
#abline(h = 1.96 / sqrt(nsim), lty = 5, lwd = 2, col = "darkgray")
#abline(h = -1.96 / sqrt(nsim), lty = 5, lwd = 2, col = "darkgray")
abline(h = c(0:10) / 10, lwd = 0.075)
abline(v = c(1, 25, 50, 75, 100, 125, 150), lwd = 0.075)

dev.off()




########## table ESSps, ESS, and times
logitUnb_tab = matrix(nrow = 9, ncol = 5)

load(paste(
  "./SimStudyLogitUnbalanced_50n1pGauss",
  "/SimStudyLogitUnbalanced_50n1pGauss.RData",
  sep = ""
))
logitUnb_tab[1, -c(1:2)] = c(
  res[[3]]$ess / as.numeric(res[[3]]$time),
  res[[2]]$ess / as.numeric(res[[2]]$time),
  res[[1]]$ess / as.numeric(res[[1]]$time)
)
logitUnb_tab[2, -c(1:2)] = c(
  res[[3]]$ess,
  res[[2]]$ess,
  res[[1]]$ess
)
logitUnb_tab[3, -c(1:2)] = c(
  res[[3]]$time,
  res[[2]]$time,
  res[[1]]$time
)

load(paste(
  "./SimStudyLogitUnbalanced_200n1pGauss",
  "/SimStudyLogitUnbalanced_200n1pGauss.RData",
  sep = ""
))
logitUnb_tab[4, -c(1:2)] = c(
  res[[3]]$ess / as.numeric(res[[3]]$time),
  res[[2]]$ess / as.numeric(res[[2]]$time),
  res[[1]]$ess / as.numeric(res[[1]]$time)
)
logitUnb_tab[5, -c(1:2)] = c(
  res[[3]]$ess,
  res[[2]]$ess,
  res[[1]]$ess
)
logitUnb_tab[6, -c(1:2)] = c(
  res[[3]]$time,
  res[[2]]$time,
  res[[1]]$time
)

load(paste(
  "./SimStudyLogitUnbalanced_1000n1pGauss",
  "/SimStudyLogitUnbalanced_1000n1pGauss.RData",
  sep = ""
))
logitUnb_tab[7, -c(1:2)] = c(
  res[[3]]$ess / as.numeric(res[[3]]$time),
  res[[2]]$ess / as.numeric(res[[2]]$time),
  res[[1]]$ess / as.numeric(res[[1]]$time)
)
logitUnb_tab[8, -c(1:2)] = c(
  res[[3]]$ess,
  res[[2]]$ess,
  res[[1]]$ess
)
logitUnb_tab[9, -c(1:2)] = c(
  res[[3]]$time,
  res[[2]]$time,
  res[[1]]$time
)

xtable::xtable(logitUnb_tab, digits = 2)
