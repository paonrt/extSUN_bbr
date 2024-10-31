################################# get results #################################
load("LungCancer_ModSel.RData")

logitSelected = logitVarSel$postInclProbs >= 0.5
probitSelected = probitVarSel$postInclProbs >= 0.5

X_medianLogit = X[, logitSelected]
X_medianProbit = X[, probitSelected >= 0.5]

medianLogitSigma = logitSigma[logitSelected, logitSelected]
medianProbitSigma = probitSigma[probitSelected, probitSelected]

medianLogitMu = mu[logitSelected]
medianProbitMu = mu[probitSelected]

set.seed(1994)

medianLogit = pSUN::pSUN_LogitGauss_IS(
  nsim = 2.5e+4, nuSUT = 100,
  A = B %*% X_medianLogit %*% diag(
    sqrt(diag(medianLogitSigma)), nrow = ncol(X_medianLogit)
  ),
  b = B %*% X_medianLogit %*% medianLogitMu, 
  xi = medianLogitMu, Omega = medianLogitSigma, plim = 1e-4
)

medianProbit = pSUN::pSUN_ProbitGauss_RNG(
  nsim = 2.5e+4, A = B %*% X_medianProbit %*% diag(
    sqrt(diag(medianProbitSigma)), nrow = ncol(X_medianProbit)
  ), b = B %*% X_medianProbit %*% medianProbitMu, xi = medianProbitMu,
  Omega = medianProbitSigma, plim = 1e-4
)

1 / (1 + exp(medianProbit$logNormConst - medianLogit$logNormConst)) # 0.7474743
1 / (1 + exp(medianLogit$logNormConst - medianProbit$logNormConst)) # 0.2525257
sum(logitSelected) # 49
sum(probitSelected) # 45

sum((logitSelected == TRUE) & (probitSelected == TRUE)) # 40

medianLogit$logNormConst # -18.31429
medianProbit$logNormConst # -19.39948
