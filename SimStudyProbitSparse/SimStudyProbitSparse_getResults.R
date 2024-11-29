################################# get results #################################
ProbitSparse_getRes = function() {

  groupsBeta = list()
  breaksBeta = c(1, 2, c(2:length(trueBeta))[diff(trueBeta) != 0], p + 1)

  groupsBeta[["intercept"]] = 1
  groupsBeta[["0.05"]] = c(breaksBeta[3]:(breaksBeta[4] - 1))
  groupsBeta[["0.10"]] = c(breaksBeta[4]:(breaksBeta[5] - 1))
  groupsBeta[["0.20"]] = c(breaksBeta[5]:(breaksBeta[6] - 1))
  groupsBeta[["0.30"]] = c(breaksBeta[6]:(breaksBeta[7] - 1))
  groupsBeta[["0.40"]] = c(breaksBeta[7]:(breaksBeta[8] - 1))
  groupsBeta[["0.50"]] = c(breaksBeta[2]:(breaksBeta[3] - 1))
  groupsBeta[["0.60"]] = c(breaksBeta[8]:(breaksBeta[9] - 1))
  groupsBeta[["0.70"]] = c(breaksBeta[9]:(breaksBeta[10] - 1))
  groupsBeta[["0.80"]] = c(breaksBeta[10]:(breaksBeta[11] - 1))
  groupsBeta[["0.90"]] = c(breaksBeta[11]:(breaksBeta[12] - 1))
  groupsBeta[["0.95"]] = c(breaksBeta[12]:(breaksBeta[13] - 1))

  acf = matrix(0, nrow = length(groupsBeta), ncol = nrow(res[[1]]$acf))
  ess = double(length = length(groupsBeta))
  essps = double(length = length(groupsBeta))
  mse = double(length = length(groupsBeta))
  mpm = double(length = length(groupsBeta))
  time = double(length = niter)

  for (iiter in 1:niter) {

    time[iiter] = as.numeric(res[[iiter]]$time)

    for (igroup in 1:length(groupsBeta)) {

      acf[igroup, ] = acf[igroup, ] +
        tryCatch(
          expr = rowMeans(res[[iiter]]$acf[, groupsBeta[[igroup]]]),
          error = function(err) res[[iiter]]$acf[, groupsBeta[[igroup]]]
        )

      ess[igroup] = ess[igroup] + mean(
        res[[iiter]]$ess[groupsBeta[[igroup]]]
      )

      essps[igroup] = essps[igroup] + mean(
        res[[iiter]]$ess[groupsBeta[[igroup]]] / as.numeric(
          res[[iiter]]$time
        )
      )

      mse[igroup] = mse[igroup] + mean(
        res[[iiter]]$err[groupsBeta[[igroup]]]^2
      )

      mpm[igroup] = mpm[igroup] + mean(
        res[[iiter]]$err[groupsBeta[[igroup]]] +
          trueBeta[groupsBeta[[igroup]]]
      )

    }

  }

  return(list(
    acf = acf / niter,
    ess = ess / niter,
    essps = essps / niter,
    mse = mse / niter,
    time = summary(time),
    mpm = mpm / niter
  ))

}



load("./SimStudyProbitSparse_50n500pGauss/SimStudyProbitSparse_50n500pGauss.RData")
Gauss = ProbitSparse_getRes()

time = double(length = niter)
for (iiter in 1:niter) {
  time[iiter] = as.numeric(res[[iiter]]$time)
}
Gauss$essps = rep(mean(nsim / time), length(Gauss$essps))

load("./SimStudyProbitSparse_50n500pLaplacit/SimStudyProbitSparse_50n500pLaplacit.RData")
Laplacit = ProbitSparse_getRes()

load("./SimStudyProbitSparse_50n500pDirLapl/SimStudyProbitSparse_50n500pDirLapl.RData")
DirLapl = ProbitSparse_getRes()





#################################### tables ####################################

mse = matrix(nrow = 13, ncol = 3)

mse[-1, 1] = Gauss$mse
mse[-1, 2] = Laplacit$mse
mse[-1, 3] = DirLapl$mse

xtable::xtable(mse, digits = 2)


essps = matrix(nrow = 13, ncol = 3)

essps[-1, 1] = Gauss$essps
essps[-1, 2] = Laplacit$essps
essps[-1, 3] = DirLapl$essps

xtable::xtable(essps, digits = 2)


ess = matrix(nrow = 13, ncol = 3)

ess[-1, 1] = rep(1e+4, length(Gauss$ess))
ess[-1, 2] = Laplacit$ess / 5
ess[-1, 3] = DirLapl$ess

xtable::xtable(ess, digits = 2)


time = matrix(nrow = 3, ncol = 6)

time[1, ] = Gauss$time
time[2, ] = Laplacit$time / 5
time[3, ] = DirLapl$time

xtable::xtable(time, digits = 2)