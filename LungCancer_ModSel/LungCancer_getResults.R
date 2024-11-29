################################# get results #################################
load("logitVarSel.RData")
load("probitVarSel.RData")
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

logit_covmedProbit = pSUN::pSUN_LogitGauss_IS(
  nsim = 2.5e+4, nuSUT = 100,
  A = B %*% X_medianProbit %*% diag(
    sqrt(diag(logitSigma[probitSelected, probitSelected])),
    nrow = ncol(X_medianProbit)
  ),
  b = B %*% X_medianProbit %*% mu[probitSelected],
  xi = mu[probitSelected], Omega = logitSigma[probitSelected, probitSelected],
  plim = 1e-4
)

probit_covmedProbit = pSUN::pSUN_ProbitGauss_RNG(
  nsim = 2.5e+4,
  A = B %*% X_medianLogit %*% diag(
    sqrt(diag(probitSigma[logitSelected, logitSelected])),
    nrow = ncol(X_medianLogit)
  ),
  b = B %*% X_medianLogit %*% mu[logitSelected],
  xi = mu[logitSelected], Omega = probitSigma[logitSelected, logitSelected],
  plim = 1e-4
)

1 / (1 + exp(medianProbit$logNormConst - medianLogit$logNormConst)) # 0.9192447
1 / (1 + exp(medianLogit$logNormConst - medianProbit$logNormConst)) # 0.0807553
sum(logitSelected) # 38
sum(probitSelected) # 29

sum((logitSelected == TRUE) & (probitSelected == TRUE)) # 18

medianLogit$logNormConst #  -17.10568
medianProbit$logNormConst # -19.5378

### BayesFactor Logit over Probit MAD
exp(medianLogit$logNormConst - medianProbit$logNormConst) # 11.3831


### logit / probit with covariates from MAD of probit
1 / (1 + exp(medianProbit$logNormConst - logit_covmedProbit$logNormConst))
# logit  0.7084
# probit 0.2916
### Bayes factor
exp(logit_covmedProbit$logNormConst - medianProbit$logNormConst)
# 2.4293

### logit / probit with covariates from MAD of logit
1 / (1 + exp(probit_covmedProbit$logNormConst - medianLogit$logNormConst))
# logit  0.7048
# probit 0.2952
### Bayes factor
exp(medianLogit$logNormConst - probit_covmedProbit$logNormConst)
# 2.3871