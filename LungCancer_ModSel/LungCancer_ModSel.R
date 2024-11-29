###############################################################################
#
#
# Lung Cancer: Models Evaluation Via Bayes Factor
#
#
###############################################################################





############################### Initialization ################################

### clear IDE
rm(list = ls()); gc(); cat("\14")

### set sequential computation for all workers
# Sys.setenv(OPENBLAS_NUM_THREADS = 1)
# Sys.setenv(MKL_NUM_THREADS = 1)
# Sys.setenv(OMP_NUM_THREADS = 1)
# Sys.setenv(BLIS_NUM_THREADS = 1)
# Sys.setenv(FFTW_NUM_THREADS = 1)

### set the seed
set.seed(1994)



##### load and pre-process data

### get function rescale from arm package
arm_rescale = function (x, binary.inputs = "center") {
  if (!is.numeric(x)) {
    x <- as.numeric(factor(x))
    x.obs <- x[!is.na(x)]
  }
  x.obs <- x[!is.na(x)]
  if (length(unique(x.obs)) == 2) {
    if (binary.inputs == "0/1") {
      x <- (x - min(x.obs))/(max(x.obs) - min(x.obs))
      return(x)
    } else if (binary.inputs == "-0.5,0.5") {
      return(x - 0.5)
    } else if (binary.inputs == "center") {
      return(x - mean(x.obs))
    } else if (binary.inputs == "full") {
      return((x - mean(x.obs))/(2 * sd(x.obs)))
    }
  } else {
    return((x - mean(x.obs)) / (2 * sd(x.obs)))
  }
}

### load raw-data from pSUN package
lung_cancer = pSUN::lung_cancer

### get response variables
y = lung_cancer[, 1]

### get standardized design matrix (with intercept)
X = cbind(
  rep(1, length(y)),
  apply(lung_cancer[, -1], 2, arm_rescale)
)



### number of scans in MCMC for posterior draws from model space
nsim = 1e+4

### number of simulations for normalizing constant computation
ntim = 1e+4

### likelihood matrix
B = 2 * diag(y) - diag(1, nrow = length(y))

### prior mean vector
mu = rep(0, ncol(X))

### prior variance-covariance matrix
# logit
logitSigma = diag(52.6379, nrow = ncol(X))
# probit
probitSigma = diag(16, nrow = ncol(X))

### verbose parameter
verbose = 1





################## Real Data Analysis of Lung Cancer Dataset ##################

########## variable selection

##### logit case

### run posterior simulation algorithm
logitVarSel = pSUN::LSBR_LogitGauss_VarSel(
  nsim, X = X, y = y, mu = mu, Sigma = logitSigma,
  ntim = ntim, verbose = verbose
)

### save results
save(list = c("logitVarSel"), file = "logitVarSel.RData")



##### probit case

### run posterior simulation algorithm
probitVarSel = pSUN::LSBR_ProbitGauss_VarSel(
  nsim, X = X, y = y, mu = mu, Sigma = probitSigma,
  ntim = ntim, verbose = verbose
)

### save results
save(list = c("probitVarSel"), file = "probitVarSel.RData")




########## median probability models comparison

##### logit case

### get variables selected variables into the median probability model
logitVarSelMedProbMod = logitVarSel$postInclProbs >= 0.5

### evaluate the median probability model
logitMedProbMod = pSUN::pSUN_LogitGauss_IS(
  ntim,
  A = B %*% X[, logitVarSelMedProbMod] %*%
    diag(
      sqrt(diag(logitSigma)[logitVarSelMedProbMod]),
      nrow = sum(logitVarSelMedProbMod)
    ),
  b = B %*% X[, logitVarSelMedProbMod] %*% mu[logitVarSelMedProbMod],
  xi = mu[logitVarSelMedProbMod],
  Omega = matrix(
    logitSigma[logitVarSelMedProbMod, logitVarSelMedProbMod],
    nrow = sum(logitVarSelMedProbMod)
  )
)



##### probit case

### get variables selected variables into the median probability model
probitVarSelMedProbMod = probitVarSel$postInclProbs >= 0.5

### evaluate the median probability model
probitMedProbMod = pSUN::pSUN_ProbitGauss_RNG(
  ntim,
  A = B %*% X[, probitVarSelMedProbMod] %*%
    diag(
      sqrt(diag(probitSigma)[probitVarSelMedProbMod]),
      nrow = sum(probitVarSelMedProbMod)
    ),
  b = B %*% X[, probitVarSelMedProbMod] %*% mu[probitVarSelMedProbMod],
  xi = mu[probitVarSelMedProbMod],
  Omega = matrix(
    probitSigma[probitVarSelMedProbMod, probitVarSelMedProbMod],
    nrow = sum(probitVarSelMedProbMod)
  )
)



### Bayes factor logit / probit
BayesFactor_logitVSprobit = exp(
  logitMedProbMod$logNormConst - probitMedProbMod$logNormConst
)

### posterior probabilities
# logit
postProb_logitMedProbMod = 1 / (
  1 + exp(probitMedProbMod$logNormConst - logitMedProbMod$logNormConst)
)
# probit
postProb_probitMedProbMod = 1 / (
  1 + exp(logitMedProbMod$logNormConst - probitMedProbMod$logNormConst)
)





################################### Closure ###################################
save(list = ls(), file = "LungCancer_ModSel.RData")