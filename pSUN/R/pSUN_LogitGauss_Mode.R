################################################################################
#
# Perturbed Unified Skew Normal (pSUN), Logit Gauss: Mode Finder
#
################################################################################


############### description
# This function computes mode from a pSUN distribution with:
#
#   1. m selective variables composed by:
#      a. joint normal variable N_m(0,I),
#      b. independent scale variables distributed as a logistic Kolmogorov
#         distribution;
#
#   2. d observed variables composed by:
#      a. joint normal variable N_p(xi, Omega);
#      b. constant scale variables equal to 1.
#
# Notice that the logistic Kolmogorov is 4 times the square of Kolmogorov
# distribution.

############### function
pSUN_LogitGauss_Mode = function(
  A, b, xi, Omega, beta_start = rep(0, ncol(A)), maxit_mode = +Inf,
  gradtol_mode = sqrt(.Machine$double.eps), details = FALSE
) {

  # A is the matrix in the constraint;
  # b is the vector in the constraint;
  # xi is the location vector;
  # Omega is the scale matrix;
  # beta_start is the starting point for the hyper-parameters;
  # maxit_mode is the maximum number of iterations for mode finding for the
  #  Laplace approximation;
  # gradtol_mode is the convergence tolerance w.r.t. gradient for mode finding;
  # details is a logical variable, if TRUE then also the negative inverse of
  #  the Hessian matrix and the Laplace approximation of the log-normalizing
  #  constant are returned


  ###############--------------- function script

  ### get number of latent variables and dimension of the pSUN distribution
  # number of latent variables
  m = nrow(A)
  # dimension of the pSUN distribution
  d = ncol(A)

  ### compute correlation matrix of Omega and get square root of its diagonal
  if (is.null(dim(Omega))) {
    sqrtDiagOmega = sqrt(Omega)
    OmegaBar = matrix(1, nrow = 1, ncol = 1)
  } else {
    sqrtDiagOmega = sqrt(diag(Omega))
    OmegaBar = cov2cor(Omega)
  }


  ### compute transpose of Cholesky factor of OmegaBar (L^T)
  tCholOmegaBar = chol(OmegaBar)

  ### get starting point
  betaTilde = (beta_start - xi) / sqrtDiagOmega

  ### counter for mode finding
  counter_mode = 0

  ### compute logit CDF
  logitCDF = c(plogis(A %*% betaTilde + b))

  ### compute gradient of logit CDF
  gradLogitCDF = 1 - logitCDF

  ### compute gradient of Lambda
  gradLambda = t(A) %*% gradLogitCDF




  ########## run mode finder
  while (TRUE) {

    ### compute negative Hessian of logit CDF
    negHessLogitCDF = logitCDF * gradLogitCDF

    ### compute decomposition of negative Hessian matrix (M^T)
    decompNegHes = t(sqrt(negHessLogitCDF) * A)

    ### compute Cholesky factor of I + M %*% OmegaBar %*% M^T
    cholL_ras = t(chol(
      diag(1, nrow = m) + crossprod(tCholOmegaBar %*% decompNegHes)
    ))

    ### compute inverse of the above Cholesky factor
    invCholL_ras = forwardsolve(cholL_ras, diag(1, nrow = m))

    ### compute (M^T %*% M) %*% beta_tile + gradLambda
    b_ras = tcrossprod(decompNegHes) %*% betaTilde + gradLambda

    ### compute solve(OmegaBar) %*% betaTilde
    a_ras = b_ras - crossprod(invCholL_ras %*% t(decompNegHes)) %*%
      OmegaBar %*% b_ras

    ### update betaTilde
    betaTilde = OmegaBar %*% a_ras

    ### update counter
    counter_mode = counter_mode + 1

    ### update logit CDF
    logitCDF = c(plogis(A %*% betaTilde + b))

    ### update gradient of logit CDF
    gradLogitCDF = 1 - logitCDF

    ### update gradient of Lambda
    gradLambda = t(A) %*% gradLogitCDF

    ### check gradient norm
    if (sqrt(sum((gradLambda - a_ras)^2)) <= gradtol_mode) break

    ### check number of iterations
    if (counter_mode >= maxit_mode) break

  }




  ##### return results

  ### compute mode
  beta = xi + sqrtDiagOmega * betaTilde

  ### compute details or return just the mode
  if (details) {

    ### update Cholesky factor of I + M %*% OmegaBar %*% M^T
    cholL_ras = t(chol(
      diag(1, nrow = m) + crossprod(tCholOmegaBar %*% decompNegHes)
    ))

    ### update inverse of the above Cholesky factor
    invCholL_ras = forwardsolve(cholL_ras, diag(1, nrow = m))

    ### compute Laplace approximation of normalizing constant
    approxLogNormConst = - 0.5 * crossprod(a_ras, betaTilde) +
      sum(log(logitCDF)) - sum(log(diag(cholL_ras))) +
      sum(log(sqrtDiagOmega))

    ### compute negative inverse of the Hessian
    negInvHess = (
      OmegaBar - crossprod(invCholL_ras %*% t(decompNegHes) %*% OmegaBar)
    ) * (
      tcrossprod(sqrtDiagOmega)
    )

    ### return mode and details
    return(list(
      mode = beta,
      approxLogNormConst = c(approxLogNormConst),
      negInvHess = negInvHess
    ))

  } else {

    ### return only the mode
    return(list(mode = beta))
  }

}





############### appendix

if (FALSE) { # appendix script is not run

  # debugging
  {
    rm(list = ls()); gc(); cat("\14")
    source(paste(.libPaths()[1], "/myFUN/00000.R", sep = ""))
    source(paste(.libPaths()[1], "/myFUN/pSUN_LogitGauss_Mode.R", sep = ""))
    set.seed(1984)
    d = 500
    m = 30
    A = matrix(rnorm(m*d), nrow = m, ncol = d)
    b = rep(0, m)
    xi = rep(0, d)
    Omega = rnorm(d)
    Omega = Omega %*% t(Omega) + diag(1, nrow = d)
    beta_start = rep(0, ncol(A))
    maxit_mode = +Inf
    gradtol_mode = sqrt(.Machine$double.eps)
    details = TRUE
    #stop("it's ok!")

    t0 = Sys.time()
    res = pSUN_LogitGauss_Mode(
      A, b, xi, Omega, beta_start, maxit_mode, gradtol_mode, details = details
    )
    tn = Sys.time()
    print(difftime(tn, t0, units = "sec"))
    #plot(res$mode)
  }


  ###############


  # test CancerSAGE
  {
    rm(list = ls()); gc(); cat("\14")
    source(paste(.libPaths()[1], "/myFUN/00000.R", sep = ""))
    source(paste(.libPaths()[1], "/myFUN/pSUN_LogitGauss_Mode.R", sep = ""))
    load(paste(
      .libPaths()[1], "/myFUN/datasets/processedDataset_CancerSAGE.RData",
      sep = ""
    ))
    set.seed(1984)
    X = X_CancerSAGE
    y = y_CancerSAGE
    B = 2 * diag(y) - diag(1, nrow = length(y))
    Sigma = diag((pi * 4 / sqrt(3))^2, nrow = ncol(X) )
    mu = rep(0, ncol(X))

    A = B %*% X %*% diag(sqrt(diag(Sigma)), nrow = ncol(X))
    b = B %*% X %*% mu
    xi = mu
    Omega = Sigma

    t0 = Sys.time()
    res_pSUN_mode = pSUN_LogitGauss_Mode(A, b, xi, Omega)
    tn = Sys.time()
    print(difftime(tn, t0, units = "sec"))
    #plot(res_pSUN_mode$mode)
  }


  ###############


  # test PimaIndians
  {
    rm(list = ls()); gc(); cat("\14")
    source(paste(.libPaths()[1], "/myFUN/00000.R", sep = ""))
    source(paste(.libPaths()[1], "/myFUN/pSUN_LogitGauss_Mode.R", sep = ""))
    load(paste(
      .libPaths()[1], "/myFUN/datasets/processedDataset_PimaIndians.RData",
      sep = ""
    ))
    set.seed(1984)
    X = X_PimaIndians
    y = y_PimaIndians
    B = 2 * diag(y) - diag(1, nrow = length(y))
    Sigma = diag((pi * 4 / sqrt(3))^2, nrow = ncol(X) )
    mu = rep(0, ncol(X))

    A = B %*% X %*% diag(sqrt(diag(Sigma)), nrow = ncol(X))
    b = B %*% X %*% mu
    xi = mu
    Omega = Sigma

    t0 = Sys.time()
    res_pSUN_mode = pSUN_LogitGauss_Mode(A, b, xi, Omega)
    tn = Sys.time()
    print(difftime(tn, t0, units = "sec"))
    #plot(res_pSUN_mode$mode)
  }


  ###############


  # test Imbalanced Data
  {
    rm(list = ls()); gc(); cat("\14")
    source(paste(.libPaths()[1], "/myFUN/00000.R", sep = ""))
    source(paste(.libPaths()[1], "/myFUN/pSUN_LogitGauss_Mode.R", sep = ""))
    set.seed(1984)
    n = 2.5 * 1e+2
    X = matrix(1, nrow = n, ncol = 1)
    y = c(1, rep(0, n - 1))
    B = 2 * diag(y) - diag(1, nrow = length(y))
    Sigma = diag((pi * 4 / sqrt(3))^2, nrow = ncol(X) )
    mu = rep(0, ncol(X))

    A = B %*% X %*% diag(sqrt(diag(Sigma)), nrow = ncol(X))
    b = B %*% X %*% mu
    xi = mu
    Omega = Sigma

    t0 = Sys.time()
    res_pSUN_mode = pSUN_LogitGauss_Mode(A, b, xi, Omega)
    tn = Sys.time()
    print(difftime(tn, t0, units = "sec"))
    #plot(res_pSUN_mode$mode)
  }


  ###############


  # test ProstateCancer
  {
    rm(list = ls()); gc(); cat("\14")
    source(paste(.libPaths()[1], "/myFUN/00000.R", sep = ""))
    source(paste(.libPaths()[1], "/myFUN/pSUN_LogitGauss_Mode.R", sep = ""))
    load(paste(
      .libPaths()[1], "/myFUN/datasets/processedDataset_ProstateCancer.RData",
      sep = ""
    ))
    set.seed(1984)
    X = X_ProstateCancer
    y = y_ProstateCancer
    B = 2 * diag(y) - diag(1, nrow = length(y))
    Sigma = diag((pi * 4 / sqrt(3))^2, nrow = ncol(X) )
    mu = rep(0, ncol(X))

    A = B %*% X %*% diag(sqrt(diag(Sigma)), nrow = ncol(X))
    b = B %*% X %*% mu
    xi = mu
    Omega = Sigma

    t0 = Sys.time()
    res_pSUN_mode = pSUN_LogitGauss_Mode(A, b, xi, Omega)
    tn = Sys.time()
    print(difftime(tn, t0, units = "sec"))
    #plot(res_pSUN_mode$mode)
  }

}
