#################################### plots ####################################
plot_int = function() {
  intercept = double(length = length(seq_quant))
  for (iquant in 1:length(seq_quant)) {
    intercept[iquant] =  mean(
      trueBeta["intercept", ] <= empQuant["intercept", iquant, ]
    )
  }
  plot(
    ylim = c(0.075, 0.925), xlim = c(0.075, 0.925), x = 0, col = "white",
    xlab = "Quantile", ylab = "Coverage", xaxt = "n", yaxt = "n"
  )
  abline(h = c(0, 0.2, 0.4, 0.6, 0.8), lwd = 0.1)
  abline(v = c(0, 0.2, 0.4, 0.6, 0.8), lwd = 0.1)
  abline(a = 0, b = 1, col = "grey50", lty = 2, lwd = 2.5)
  points(seq_quant, intercept, type = "l", col = "black", lwd = 0.75)
  axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8))
  axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8))
  return(intercept)
}

plot_not = function() {
  notIntercept = double(length = length(seq_quant))
  for (iquant in 1:length(seq_quant)) {
    notIntercept[iquant] = mean(
      trueBeta["notIntercept", ] <= empQuant["notIntercept", iquant, ]
    )
  }
  plot(
    ylim = c(0.075, 0.925), xlim = c(0.075, 0.925), x = 0, col = "white",
    xlab = "Quantile", ylab = "Coverage", xaxt = "n", yaxt = "n"
  )
  abline(h = c(0, 0.2, 0.4, 0.6, 0.8), lwd = 0.1)
  abline(v = c(0, 0.2, 0.4, 0.6, 0.8), lwd = 0.1)
  abline(a = 0, b = 1, col = "grey50", lty = 2, lwd = 2.5)
  points(seq_quant, notIntercept, type = "l", col = "black", lwd = 0.75)
  axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8))
  axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8))
  return(notIntercept)
}

##### intercept
png("IScoverage_intercept.png", width = 750, height = 750, res = 125)
par(mar = c(4, 4, 1, 1))
par(mfrow = c(3, 2))
load("./SimStudyIS_n50p500/SimStudyIS_n50p500.RData")
plot_int()
load("./SimStudyIS_n50p1000/SimStudyIS_n50p1000.RData")
plot_int()
load("./SimStudyIS_n100p500/SimStudyIS_n100p500.RData")
plot_int()
load("./SimStudyIS_n100p1000/SimStudyIS_n100p1000.RData")
plot_int()
load("./SimStudyIS_n200p500/SimStudyIS_n200p500.RData")
plot_int()
load("./SimStudyIS_n200p1000/SimStudyIS_n200p1000.RData")
plot_int()
dev.off()

##### not intercept
png("IScoverage_notIntercept.png", width = 750, height = 750, res = 125)
par(mar = c(4, 4, 1, 1))
par(mfrow = c(3, 2))
load("./SimStudyIS_n50p500/SimStudyIS_n50p500.RData")
plot_not()
load("./SimStudyIS_n50p1000/SimStudyIS_n50p1000.RData")
plot_not()
load("./SimStudyIS_n100p500/SimStudyIS_n100p500.RData")
plot_not()
load("./SimStudyIS_n100p1000/SimStudyIS_n100p1000.RData")
plot_not()
load("./SimStudyIS_n200p500/SimStudyIS_n200p500.RData")
plot_not()
load("./SimStudyIS_n200p1000/SimStudyIS_n200p1000.RData")
plot_not()
dev.off()





################################### tables ####################################

########## coverage
cov_int = function() {
  intercept = double(length = length(seq_quant))
  for (iquant in 1:length(seq_quant)) {
    intercept[iquant] =  mean(
      trueBeta["intercept", ] <= empQuant["intercept", iquant, ]
    )
  }
  return(intercept)
}

cov_not = function() {
  notIntercept = double(length = length(seq_quant))
  for (iquant in 1:length(seq_quant)) {
    notIntercept[iquant] = mean(
      trueBeta["notIntercept", ] <= empQuant["notIntercept", iquant, ]
    )
  }
  return(notIntercept)
}


tmp_int = matrix(nrow = length(seq_quant), ncol = 7)
tmp_int[, 1] = seq_quant

load("./SimStudyIS_n50p500/SimStudyIS_n50p500.RData")
tmp_int[, 2] = cov_int()
load("./SimStudyIS_n50p1000/SimStudyIS_n50p1000.RData")
tmp_int[, 3] = cov_int()
load("./SimStudyIS_n100p500/SimStudyIS_n100p500.RData")
tmp_int[, 4] = cov_int()
load("./SimStudyIS_n100p1000/SimStudyIS_n100p1000.RData")
tmp_int[, 5] = cov_int()
load("./SimStudyIS_n200p500/SimStudyIS_n200p500.RData")
tmp_int[, 6] = cov_int()
load("./SimStudyIS_n200p1000/SimStudyIS_n200p1000.RData")
tmp_int[, 7] = cov_int()

xtable::xtable(tmp_int, digits = 4)


tmp_not = matrix(nrow = length(seq_quant), ncol = 7)
tmp_not[, 1] = seq_quant

load("./SimStudyIS_n50p500/SimStudyIS_n50p500.RData")
tmp_not[, 2] = cov_not()
load("./SimStudyIS_n50p1000/SimStudyIS_n50p1000.RData")
tmp_not[, 3] = cov_not()
load("./SimStudyIS_n100p500/SimStudyIS_n100p500.RData")
tmp_not[, 4] = cov_not()
load("./SimStudyIS_n100p1000/SimStudyIS_n100p1000.RData")
tmp_not[, 5] = cov_not()
load("./SimStudyIS_n200p500/SimStudyIS_n200p500.RData")
tmp_not[, 6] = cov_not()
load("./SimStudyIS_n200p1000/SimStudyIS_n200p1000.RData")
tmp_not[, 7] = cov_not()

xtable::xtable(tmp_not, digits = 4)




########## effective sample size per second
tmp_ESSps = matrix(nrow = 6, ncol = 8)

tmp_ESSps[, 1] = c(rep(50, 2), rep(100, 2), rep(200, 2))
tmp_ESSps[, 2] = rep(c(500, 1000), 3)

load("./SimStudyIS_n50p500/SimStudyIS_n50p500.RData")
tmp_ESSps[1, -c(1, 2)] = summary(effSamplSize / timesSecs)
load("./SimStudyIS_n50p1000/SimStudyIS_n50p1000.RData")
tmp_ESSps[2, -c(1, 2)] = summary(effSamplSize / timesSecs)
load("./SimStudyIS_n100p500/SimStudyIS_n100p500.RData")
tmp_ESSps[3, -c(1, 2)] = summary(effSamplSize / timesSecs)
load("./SimStudyIS_n100p1000/SimStudyIS_n100p1000.RData")
tmp_ESSps[4, -c(1, 2)] = summary(effSamplSize / timesSecs)
load("./SimStudyIS_n200p500/SimStudyIS_n200p500.RData")
tmp_ESSps[5, -c(1, 2)] = summary(effSamplSize / timesSecs)
load("./SimStudyIS_n200p1000/SimStudyIS_n200p1000.RData")
tmp_ESSps[6, -c(1, 2)] = summary(effSamplSize / timesSecs)

xtable::xtable(tmp_ESSps, digits = 2)




########## computational time (second)
tmp_time = matrix(nrow = 6, ncol = 8)

tmp_time[, 1] = c(rep(50, 2), rep(100, 2), rep(200, 2))
tmp_time[, 2] = rep(c(500, 1000), 3)

load("./SimStudyIS_n50p500/SimStudyIS_n50p500.RData")
tmp_time[1, -c(1, 2)] = summary(timesSecs)
load("./SimStudyIS_n50p1000/SimStudyIS_n50p1000.RData")
tmp_time[2, -c(1, 2)] = summary(timesSecs)
load("./SimStudyIS_n100p500/SimStudyIS_n100p500.RData")
tmp_time[3, -c(1, 2)] = summary(timesSecs)
load("./SimStudyIS_n100p1000/SimStudyIS_n100p1000.RData")
tmp_time[4, -c(1, 2)] = summary(timesSecs)
load("./SimStudyIS_n200p500/SimStudyIS_n200p500.RData")
tmp_time[5, -c(1, 2)] = summary(timesSecs)
load("./SimStudyIS_n200p1000/SimStudyIS_n200p1000.RData")
tmp_time[6, -c(1, 2)] = summary(timesSecs)

xtable::xtable(tmp_time, digits = 2)




########## effective sample size
tmp_ESS = matrix(nrow = 6, ncol = 8)

tmp_ESS[, 1] = c(rep(50, 2), rep(100, 2), rep(200, 2))
tmp_ESS[, 2] = rep(c(500, 1000), 3)

load("./SimStudyIS_n50p500/SimStudyIS_n50p500.RData")
tmp_ESS[1, -c(1, 2)] = summary(effSamplSize)
load("./SimStudyIS_n50p1000/SimStudyIS_n50p1000.RData")
tmp_ESS[2, -c(1, 2)] = summary(effSamplSize)
load("./SimStudyIS_n100p500/SimStudyIS_n100p500.RData")
tmp_ESS[3, -c(1, 2)] = summary(effSamplSize)
load("./SimStudyIS_n100p1000/SimStudyIS_n100p1000.RData")
tmp_ESS[4, -c(1, 2)] = summary(effSamplSize)
load("./SimStudyIS_n200p500/SimStudyIS_n200p500.RData")
tmp_ESS[5, -c(1, 2)] = summary(effSamplSize)
load("./SimStudyIS_n200p1000/SimStudyIS_n200p1000.RData")
tmp_ESS[6, -c(1, 2)] = summary(effSamplSize)

xtable::xtable(tmp_ESS, digits = 2)
