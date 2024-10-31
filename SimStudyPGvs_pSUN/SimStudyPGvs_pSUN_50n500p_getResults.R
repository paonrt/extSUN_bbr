################################# get results #################################
PGvs_pSUN_getRes = function() {
  groupsBeta = list()
  breaksBeta = c(1, c(2:length(trueBeta))[diff(trueBeta) != 0], p + 1)

  groupsBeta[["intercept"]] = 1
  groupsBeta[["0.05"]] = c(breaksBeta[2]:(breaksBeta[3] - 1))
  groupsBeta[["0.10"]] = c(breaksBeta[3]:(breaksBeta[4] - 1))
  groupsBeta[["0.20"]] = c(breaksBeta[4]:(breaksBeta[5] - 1))
  groupsBeta[["0.30"]] = c(breaksBeta[5]:(breaksBeta[6] - 1))
  groupsBeta[["0.40"]] = c(breaksBeta[6]:(breaksBeta[7] - 1))
  groupsBeta[["0.50"]] = c(breaksBeta[12]:(breaksBeta[13] - 1))
  groupsBeta[["0.60"]] = c(breaksBeta[7]:(breaksBeta[8] - 1))
  groupsBeta[["0.70"]] = c(breaksBeta[8]:(breaksBeta[9] - 1))
  groupsBeta[["0.80"]] = c(breaksBeta[9]:(breaksBeta[10] - 1))
  groupsBeta[["0.90"]] = c(breaksBeta[10]:(breaksBeta[11] - 1))
  groupsBeta[["0.95"]] = c(breaksBeta[11]:(breaksBeta[12] - 1))


  pSUN_acf = matrix(
    0, nrow = length(groupsBeta), ncol = nrow(res[[1]]$pSUN_acf)
  )
  PG_acf = matrix(0, nrow = length(groupsBeta), ncol = nrow(res[[1]]$PG_acf))
  UPG_acf = matrix(0, nrow = length(groupsBeta), ncol = nrow(res[[1]]$UPG_acf))

  pSUN_ess = double(length = length(groupsBeta))
  PG_ess = double(length = length(groupsBeta))
  UPG_ess = double(length = length(groupsBeta))

  pSUN_essps = double(length = length(groupsBeta))
  PG_essps = double(length = length(groupsBeta))
  UPG_essps = double(length = length(groupsBeta))

  pSUN_mse = double(length = length(groupsBeta))
  PG_mse = double(length = length(groupsBeta))
  UPG_mse = double(length = length(groupsBeta))

  pSUN_mpm = double(length = length(groupsBeta))
  PG_mpm = double(length = length(groupsBeta))
  UPG_mpm = double(length = length(groupsBeta))

  pSUN_time = double(length = niter)
  PG_time = double(length = niter)
  UPG_time = double(length = niter)

  for (iiter in 1:niter) {

    pSUN_time[iiter] = as.numeric(res[[iiter]]$pSUN_time)
    PG_time[iiter] = as.numeric(res[[iiter]]$PG_time)
    UPG_time[iiter] = as.numeric(res[[iiter]]$UPG_time)

    for (igroup in 1:length(groupsBeta)) {

      pSUN_acf[igroup, ] = pSUN_acf[igroup, ] +
        tryCatch(
          expr = rowMeans(res[[iiter]]$pSUN_acf[, groupsBeta[[igroup]]]),
          error = function(err) res[[iiter]]$pSUN_acf[, groupsBeta[[igroup]]]
        )

      PG_acf[igroup, ] = PG_acf[igroup, ] +
        tryCatch(
          expr = rowMeans(res[[iiter]]$PG_acf[, groupsBeta[[igroup]]]),
          error = function(err) res[[iiter]]$PG_acf[, groupsBeta[[igroup]]]
        )

      UPG_acf[igroup, ] = UPG_acf[igroup, ] +
        tryCatch(
          expr = rowMeans(res[[iiter]]$UPG_acf[, groupsBeta[[igroup]]]),
          error = function(err) res[[iiter]]$UPG_acf[, groupsBeta[[igroup]]]
        )


      pSUN_ess[igroup] = pSUN_ess[igroup] + mean(
        res[[iiter]]$pSUN_ess[groupsBeta[[igroup]]]
      )

      PG_ess[igroup] = PG_ess[igroup] + mean(
        res[[iiter]]$PG_ess[groupsBeta[[igroup]]]
      )

      UPG_ess[igroup] = UPG_ess[igroup] + mean(
        res[[iiter]]$UPG_ess[groupsBeta[[igroup]]]
      )


      pSUN_essps[igroup] = pSUN_essps[igroup] + mean(
        res[[iiter]]$pSUN_ess[groupsBeta[[igroup]]] / as.numeric(
          res[[iiter]]$pSUN_time
        )
      )

      PG_essps[igroup] = PG_essps[igroup] + mean(
        res[[iiter]]$PG_ess[groupsBeta[[igroup]]] / as.numeric(
          res[[iiter]]$PG_time
        )
      )

      UPG_essps[igroup] = UPG_essps[igroup] + mean(
        res[[iiter]]$UPG_ess[groupsBeta[[igroup]]] / as.numeric(
          res[[iiter]]$UPG_time
        )
      )# iiter = 172, igroup = 4 and 49-th component


      pSUN_mse[igroup] = pSUN_mse[igroup] + mean(
        res[[iiter]]$pSUN_err[groupsBeta[[igroup]]]^2
      )

      PG_mse[igroup] = PG_mse[igroup] + mean(
        res[[iiter]]$PG_err[groupsBeta[[igroup]]]^2
      )

      UPG_mse[igroup] = UPG_mse[igroup] + mean(
        res[[iiter]]$UPG_err[groupsBeta[[igroup]]]^2
      )

      pSUN_mpm[igroup] = pSUN_mpm[igroup] + mean(
        res[[iiter]]$pSUN_err[groupsBeta[[igroup]]] +
          trueBeta[groupsBeta[[igroup]]]
      )

      PG_mpm[igroup] = PG_mpm[igroup] + mean(
        res[[iiter]]$PG_err[groupsBeta[[igroup]]] +
          trueBeta[groupsBeta[[igroup]]]
      )

      UPG_mpm[igroup] = UPG_mpm[igroup] + mean(
        res[[iiter]]$UPG_err[groupsBeta[[igroup]]] +
          trueBeta[groupsBeta[[igroup]]]
      )

    }

  }

  return(list(
    pSUN_acf = pSUN_acf / niter,
    PG_acf = PG_acf / niter,
    UPG_acf = UPG_acf / niter,
    pSUN_ess = pSUN_ess / niter,
    PG_ess = PG_ess / niter,
    UPG_ess = UPG_ess / niter,
    pSUN_essps = pSUN_essps / niter,
    PG_essps = PG_essps / niter,
    UPG_essps = UPG_essps / niter,
    pSUN_mse = pSUN_mse / niter,
    PG_mse = PG_mse / niter,
    UPG_mse = UPG_mse / niter,
    pSUN_time = summary(pSUN_time),
    PG_time = summary(PG_time),
    UPG_time = summary(UPG_time),
    pSUN_mpm = pSUN_mpm / niter,
    PG_mpm = PG_mpm / niter,
    UPG_mpm = UPG_mpm / niter
  ))

}

IS_getResults = function() {
  groupsBeta = list()
  breaksBeta = c(1, c(2:length(trueBeta))[diff(trueBeta) != 0], p + 1)

  groupsBeta[["intercept"]] = 1
  groupsBeta[["0.05"]] = c(breaksBeta[2]:(breaksBeta[3] - 1))
  groupsBeta[["0.10"]] = c(breaksBeta[3]:(breaksBeta[4] - 1))
  groupsBeta[["0.20"]] = c(breaksBeta[4]:(breaksBeta[5] - 1))
  groupsBeta[["0.30"]] = c(breaksBeta[5]:(breaksBeta[6] - 1))
  groupsBeta[["0.40"]] = c(breaksBeta[6]:(breaksBeta[7] - 1))
  groupsBeta[["0.50"]] = c(breaksBeta[12]:(breaksBeta[13] - 1))
  groupsBeta[["0.60"]] = c(breaksBeta[7]:(breaksBeta[8] - 1))
  groupsBeta[["0.70"]] = c(breaksBeta[8]:(breaksBeta[9] - 1))
  groupsBeta[["0.80"]] = c(breaksBeta[9]:(breaksBeta[10] - 1))
  groupsBeta[["0.90"]] = c(breaksBeta[10]:(breaksBeta[11] - 1))
  groupsBeta[["0.95"]] = c(breaksBeta[11]:(breaksBeta[12] - 1))

  IS_ess = double(length = length(groupsBeta))
  IS_essps = double(length = length(groupsBeta))
  IS_mse = double(length = length(groupsBeta))
  IS_time = double(length = niter)

  for (iiter in 1:niter) {

    IS_time[iiter] = as.numeric(res[[iiter]]$IS_time)

    for (igroup in 1:length(groupsBeta)) {

      IS_ess[igroup] = IS_ess[igroup] + mean(
        res[[iiter]]$IS_ess[groupsBeta[[igroup]]]
      )


      IS_essps[igroup] = IS_essps[igroup] + mean(
        res[[iiter]]$IS_ess[groupsBeta[[igroup]]] / as.numeric(
          res[[iiter]]$IS_time
        )
      )

      IS_mse[igroup] = IS_mse[igroup] + mean(
        res[[iiter]]$IS_err[groupsBeta[[igroup]]]^2 #+
          #trueBeta[groupsBeta[[igroup]]]
      )

    }

  }

  return(list(
    IS_ess = IS_ess / niter,
    IS_essps = IS_essps / niter,
    IS_mse = IS_mse / niter,
    IS_time = summary(IS_time)
  ))
}



load("./SimStudyPGvs_pSUN_50n500pGauss/SimStudyPGvs_pSUN_50n500pGauss.RData")
Gauss = PGvs_pSUN_getRes()
IS = IS_getResults()

load(
  "./SimStudyPGvs_pSUN_50n500pLaplacit/SimStudyPGvs_pSUN_50n500pLaplacit.RData"
)
Laplacit = PGvs_pSUN_getRes()

load(
  "./SimStudyPGvs_pSUN_50n500pDirLapl/SimStudyPGvs_pSUN_50n500pDirLapl.RData"
)
DirLapl = PGvs_pSUN_getRes()





#################################### plots ####################################

########## Mean of ACF for intercept for different priors and algorithms
png("PGvspSUN_interceptACF.png", width = 750, height = 500, res = 125)
par(mfrow = c(1, 3))
par(mar = c(4.5, 4, 1, 1))
lag_max = 150
seq_lag = c(1, 50, 100, 150)

plot(
  x = c(1:lag_max), y = Gauss$pSUN_acf[1, 2:(lag_max + 1)],
  type = "l", lty = 1, lwd = 2,
  xlab = "Lag", ylab = "ACF", ylim = c(0, 1), xaxt = "n"
)
points(
  x = c(1:lag_max), y = Gauss$PG_acf[1, 2:(lag_max + 1)],
  type = "l", lty = 2, lwd = 2, col = "blue"
)
points(
  x = c(1:lag_max), y = Gauss$UPG_acf[1, 2:(lag_max + 1)],
  type = "l", lty = 3, lwd = 2, col = "red"
)
axis(1, at = seq_lag)
#abline(h = 1.96 / sqrt(nsim), lty = 5, lwd = 2, col = "darkgray")
#abline(h = -1.96 / sqrt(nsim), lty = 5, lwd = 2, col = "darkgray")
abline(h = c(0:10) / 10, lwd = 0.075)
abline(v = c(1, 25, 50, 75, 100, 125, 150), lwd = 0.075)

plot(
  x = c(1:lag_max), y = Laplacit$pSUN_acf[1, 2:(lag_max + 1)],
  type = "l", lty = 1, lwd = 2,
  xlab = "Lag", ylab = "ACF", ylim = c(0, 1), xaxt = "n"
)
points(
  x = c(1:lag_max), y = Laplacit$PG_acf[1, 2:(lag_max + 1)],
  type = "l", lty = 2, lwd = 2, col = "blue"
)
points(
  x = c(1:lag_max), y = Laplacit$UPG_acf[1, 2:(lag_max + 1)],
  type = "l", lty = 3, lwd = 2, col = "red"
)
axis(1, at = seq_lag)
#abline(h = 1.96 / sqrt(nsim), lty = 5, lwd = 2, col = "darkgray")
#abline(h = -1.96 / sqrt(nsim), lty = 5, lwd = 2, col = "darkgray")
abline(h = c(0:10) / 10, lwd = 0.075)
abline(v = c(1, 25, 50, 75, 100, 125, 150), lwd = 0.075)

plot(
  x = c(1:lag_max), y = DirLapl$pSUN_acf[1, 2:(lag_max + 1)],
  type = "l", lty = 1, lwd = 2,
  xlab = "Lag", ylab = "ACF", ylim = c(0, 1), xaxt = "n"
)
points(
  x = c(1:lag_max), y = DirLapl$PG_acf[1, 2:(lag_max + 1)],
  type = "l", lty = 2, lwd = 2, col = "blue"
)
points(
  x = c(1:lag_max), y = DirLapl$UPG_acf[1, 2:(lag_max + 1)],
  type = "l", lty = 3, lwd = 2, col = "red"
)
axis(1, at = seq_lag)
#abline(h = 1.96 / sqrt(nsim), lty = 5, lwd = 2, col = "darkgray")
#abline(h = -1.96 / sqrt(nsim), lty = 5, lwd = 2, col = "darkgray")
abline(h = c(0:10) / 10, lwd = 0.075)
abline(v = c(1, 25, 50, 75, 100, 125, 150), lwd = 0.075)
dev.off()




########## Mean of Posterior Means by different beta's groups and algorithms
#           for the Dirichlet-Laplace prior
MPM = matrix(nrow = 12, ncol = 9)

MPM[, 1] = Gauss$UPG_mpm
MPM[, 2] = Gauss$PG_mpm
MPM[, 3] = Gauss$pSUN_mpm
MPM[, 4] = Laplacit$UPG_mpm
MPM[, 5] = Laplacit$PG_mpm
MPM[, 6] = Laplacit$pSUN_mpm
MPM[, 7] = DirLapl$UPG_mpm
MPM[, 8] = DirLapl$PG_mpm
MPM[, 9] = DirLapl$pSUN_mpm

png("PGvspSUN_DirLaplMPM.png", width = 900, height = 500, res = 125)
par(mfrow = c(1, 2))
par(mar = c(4.5, 4, 1, 1))
plot(
  x = c(1:length(DirLapl$pSUN_mpm)), y = DirLapl$pSUN_mpm,
  type = "p", pch = 4, cex = 2,
  xlab = "Group", ylab = "Value", ylim = c(-0.4, 0.4), col = "black",
  xaxt = "n", yaxt = "n"
)
points(DirLapl$UPG_mpm, pch = 20, col = "red")
axis(1, at = c(2, 4, 6, 8, 10, 12))
axis(2, at = c(-0.35, -0.25, -0.15, 0, 0.15, 0.25, 0.35))
abline(h = c(-0.35, -0.25, -0.15, 0, 0.15, 0.25, 0.35), lwd = 0.075)
abline(v = c(2, 4, 6, 8, 10, 12), lwd = 0.075)
plot(
  x = c(1:length(DirLapl$pSUN_mpm)), y = DirLapl$pSUN_mpm,
  type = "p", pch = 4, cex = 2,
  xlab = "Group", ylab = "Value", ylim = c(-0.4, 0.4), col = "black",
  xaxt = "n", yaxt = "n"
)
points(DirLapl$PG_mpm, pch = 20, col = "blue")
axis(1, at = c(2, 4, 6, 8, 10, 12))
axis(2, at = c(-0.35, -0.25, -0.15, 0, 0.15, 0.25, 0.35))
abline(h = c(-0.35, -0.25, -0.15, 0, 0.15, 0.25, 0.35), lwd = 0.075)
abline(v = c(2, 4, 6, 8, 10, 12), lwd = 0.075)
dev.off()





#################################### tables ###################################

########## mean of ESSps by beta's groups
# ESSps = matrix(nrow = 10, ncol = 13)
# ESSps[1, -1] = Gauss$UPG_essps
# ESSps[2, -1] = Gauss$PG_essps
# ESSps[3, -1] = Gauss$pSUN_essps
# ESSps[4, -1] = IS$IS_essps
# ESSps[5, -1] = Laplacit$UPG_essps
# ESSps[6, -1] = Laplacit$PG_essps
# ESSps[7, -1] = Laplacit$pSUN_essps
# ESSps[8, -1] = DirLapl$UPG_essps
# ESSps[9, -1] = DirLapl$PG_essps
# ESSps[10, -1] = DirLapl$pSUN_essps
# xtable::xtable(ESSps, digits = 2)

### Gaussian prior
Gauss_ESSps = matrix(nrow = 13, ncol = 4)

Gauss_ESSps[-1, 1] = Gauss$UPG_essps
Gauss_ESSps[-1, 2] = Gauss$PG_essps
Gauss_ESSps[-1, 3] = Gauss$pSUN_essps
Gauss_ESSps[-1, 4] = IS$IS_essps

xtable::xtable(Gauss_ESSps, digits = 2)


### Laplacit prior
Laplacit_ESSps = matrix(nrow = 13, ncol = 3)

Laplacit_ESSps[-1, 1] = Laplacit$UPG_essps
Laplacit_ESSps[-1, 2] = Laplacit$PG_essps
Laplacit_ESSps[-1, 3] = Laplacit$pSUN_essps

xtable::xtable(Laplacit_ESSps, digits = 2)


### Dirichlet-Laplace prior
DirLapl_ESSps = matrix(nrow = 13, ncol = 3)

DirLapl_ESSps[-1, 1] = DirLapl$UPG_essps
DirLapl_ESSps[-1, 2] = DirLapl$PG_essps
DirLapl_ESSps[-1, 3] = DirLapl$pSUN_essps

xtable::xtable(DirLapl_ESSps, digits = 2)




########## mean of ESS by beta's groups

### Gaussian prior
Gauss_ESS = matrix(nrow = 13, ncol = 4)

Gauss_ESS[-1, 1] = Gauss$UPG_ess
Gauss_ESS[-1, 2] = Gauss$PG_ess
Gauss_ESS[-1, 3] = Gauss$pSUN_ess
Gauss_ESS[-1, 4] = IS$IS_ess

xtable::xtable(Gauss_ESS, digits = 2)


### Laplacit prior
Laplacit_ESS = matrix(nrow = 13, ncol = 3)

Laplacit_ESS[-1, 1] = Laplacit$UPG_ess
Laplacit_ESS[-1, 2] = Laplacit$PG_ess
Laplacit_ESS[-1, 3] = Laplacit$pSUN_ess

xtable::xtable(Laplacit_ESS, digits = 2)


### Dirichlet-Laplace prior
DirLapl_ESS = matrix(nrow = 13, ncol = 3)

DirLapl_ESS[-1, 1] = DirLapl$UPG_ess
DirLapl_ESS[-1, 2] = DirLapl$PG_ess
DirLapl_ESS[-1, 3] = DirLapl$pSUN_ess

xtable::xtable(DirLapl_ESS / 5, digits = 2)




########## summaries of computational times
times = matrix(nrow = 10, ncol = 8)
times[1, -c(1, 2)] = Gauss$UPG_time
times[2, -c(1, 2)] = Gauss$PG_time
times[3, -c(1, 2)] = Gauss$pSUN_time
times[4, -c(1, 2)] = IS$IS_time
times[5, -c(1, 2)] = Laplacit$UPG_time
times[6, -c(1, 2)] = Laplacit$PG_time
times[7, -c(1, 2)] = Laplacit$pSUN_time
times[8, -c(1, 2)] = DirLapl$UPG_time / 5
times[9, -c(1, 2)] = DirLapl$PG_time / 5
times[10, -c(1, 2)] = DirLapl$pSUN_time / 5

xtable::xtable(times, digits = 2)




########## MSE using pSUN-Gibbs and pSUN-IS by different beta's groups and
#           priors
MSE = matrix(nrow = 12, ncol = 10)

MSE[, 1] = Gauss$UPG_mse
MSE[, 2] = Gauss$PG_mse
MSE[, 3] = Gauss$pSUN_mse
MSE[, 4] = IS$IS_mse
MSE[, 5] = Laplacit$UPG_mse
MSE[, 6] = Laplacit$PG_mse
MSE[, 7] = Laplacit$pSUN_mse
MSE[, 8] = DirLapl$UPG_mse
MSE[, 9] = DirLapl$PG_mse
MSE[, 10] = DirLapl$pSUN_mse

xtable::xtable(MSE[, c(3, 7, 10, 4)])





### check why UPG and PG underestimate the MSE using Dirichlet-Laplace
groupsBeta = list()
breaksBeta = c(1, c(2:length(trueBeta))[diff(trueBeta) != 0], p + 1)

groupsBeta[["intercept"]] = 1
groupsBeta[["0.05"]] = c(breaksBeta[2]:(breaksBeta[3] - 1))
groupsBeta[["0.10"]] = c(breaksBeta[3]:(breaksBeta[4] - 1))
groupsBeta[["0.20"]] = c(breaksBeta[4]:(breaksBeta[5] - 1))
groupsBeta[["0.30"]] = c(breaksBeta[5]:(breaksBeta[6] - 1))
groupsBeta[["0.40"]] = c(breaksBeta[6]:(breaksBeta[7] - 1))
groupsBeta[["0.50"]] = c(breaksBeta[12]:(breaksBeta[13] - 1))
groupsBeta[["0.60"]] = c(breaksBeta[7]:(breaksBeta[8] - 1))
groupsBeta[["0.70"]] = c(breaksBeta[8]:(breaksBeta[9] - 1))
groupsBeta[["0.80"]] = c(breaksBeta[9]:(breaksBeta[10] - 1))
groupsBeta[["0.90"]] = c(breaksBeta[10]:(breaksBeta[11] - 1))
groupsBeta[["0.95"]] = c(breaksBeta[11]:(breaksBeta[12] - 1))

pSUN_pm = matrix(nrow = length(groupsBeta), ncol = niter)
PG_pm = matrix(nrow = length(groupsBeta), ncol = niter)
UPG_pm = matrix(nrow = length(groupsBeta), ncol = niter)

for (iiter in 1:niter) {

  for (igroup in 1:length(groupsBeta)) {

    pSUN_pm[igroup, iiter] = mean(
      res[[iiter]]$pSUN_err[groupsBeta[[igroup]]] +
        trueBeta[groupsBeta[[igroup]]]
    )

    PG_pm[igroup, iiter] = mean(
      res[[iiter]]$PG_err[groupsBeta[[igroup]]] +
        trueBeta[groupsBeta[[igroup]]]
    )

    UPG_pm[igroup, iiter] = mean(
      res[[iiter]]$UPG_err[groupsBeta[[igroup]]] +
        trueBeta[groupsBeta[[igroup]]]
    )

  }

}

plot(pSUN_pm[2, ], ylim = c(-3, 0.6))
points(PG_pm[2, ], col = "blue", pch = 20, cex = 1)
points(UPG_pm[2, ], col = "red", pch = 20, cex = 0.5)

plot(pSUN_pm[2, ], ylim = c(-3, 0.6))
points(PG_pm[2, ], col = "blue", pch = 20, cex = 1)

plot(pSUN_pm[2, ], ylim = c(-3, 0.6))
points(UPG_pm[2, ], col = "red", pch = 20, cex = 0.5)

plot(pSUN_pm[2, 1:500], ylim = c(-3, 0.6))
points(PG_pm[2, 1:500], col = "blue", pch = 20, cex = 1)
points(UPG_pm[2, 1:500], col = "red", pch = 20, cex = 0.5)
abline(h = qlogis(0.05), col = "green3") # true value

mean(abs(pSUN_pm[2, ] - qlogis(0.05)))
mean(abs(PG_pm[2, ] - qlogis(0.05)))
mean(abs(UPG_pm[2, ] - qlogis(0.05)))