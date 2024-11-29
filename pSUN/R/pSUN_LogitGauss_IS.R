#' @title Perturbed Unified Skew Normal Logit Gaussian: Importance Sampling

#' @description
#' Importance sampling for a perturbed unified skew normal (pSUN) logit Gaussian
#'  distribution.

#' @usage
#' pSUN_LogitGauss_IS = function(nsim, A, b, xi, Omega, nuSUT = NULL,
#'                               V = NULL, plim = 1 / nsim)
#'
#' @details
#' It is said \eqn{X} follows a pSUN logit Gaussian distribution if
#'  \deqn{
#'   X = \xi + \mathrm{diag}^\frac{1}{2} (\Omega) \, Z \vert T \le A Z + b \\
#'   Z \vert W \sim N(0, \bar{\Omega}_W) \perp\!\!\!\perp T \vert V \sim N(0,
#'    \Theta_V) \\
#'   W_i = 1 \, , \, V_i \overset{i.i.d.}{\sim} \mathrm{LK}(\cdot)
#'  }
#'  where \eqn{\bar{\Omega}_W = \mathrm{diag}^\frac{1}{2}(W) \bar{\Omega}
#'  \mathrm{diag}^\frac{1}{2}(W)}, \eqn{\Theta_V = \mathrm{diag}^\frac{1}{2}(V)
#'  \Theta \mathrm{diag}^\frac{1}{2}(V)}, \eqn{\bar{\Omega}} is the correlation
#'  matrix of \eqn{\Omega} and \eqn{\mathrm{LK}(\cdot)} is the CDF of a logistic
#'  Kolmogorov distribution.
#'
#' Notice that a logistic Kolmogorov distribution is four times the square of
#'  a Kolmogorov random variable.
#'
#' The importance density is a Unified Skew-t distribution (SUT). It is said
#'  \eqn{X} follows a SUT distribution if
#'  \deqn{
#'   X = \xi + \sqrt{S} \mathrm{diag}^\frac{1}{2} (\Omega) \, Z  \\
#'   Z \sim \mathrm{SUN}(\Theta, A, b, 0, \bar{\Omega}) \perp\!\!\!\perp S \sim
#'   \mathrm{inv.gamma}(\nu / 2, \nu / 2)
#'  }
#'  where \eqn{\bar{\Omega}} is the correlation matrix of \eqn{\Omega}. In this
#'  function, \eqn{\Theta} is restricted to be the identity matrix.
#'  See Onorati and Liseo (2022) for more details.

#' @references
#' Onorati, P., & Liseo, B. (2022). An extension of the unified skew-normal
#'  family of distributions and application to Bayesian binary regression.
#'  arXiv preprint arXiv:2209.03474.

#' @param nsim The number of simulations.
#' @param A The slant matrix.
#' @param b The slant vector.
#' @param xi The location vector.
#' @param Omega The scale matrix.
#' @param nuSUT The shape parameter of SUT distribution.
#' @param V The vector of constant values for the scales of selective variable.
#' @param plim Minimum empirical acceptance probability such that the
#'  function does not stop with an error.
#' @return Sample, normalized weights and logarithm of normalizing constant of
#'  the pSUN distribution.

#' @import matrixStats
#' @export
pSUN_LogitGauss_IS = function(
  nsim, A, b, xi, Omega, nuSUT = NULL, V = NULL, plim = 1 / nsim
) {

  ### check xi
  xi = matrix(xi, ncol = 1)

  ### check b
  b = matrix(b, ncol = 1)

  ### check Omega
  Omega = matrix(Omega, nrow = ncol(A))

  ### dimension of the distribution
  d = ncol(A)

  ### dimension of selective latent variable
  m = nrow(A)



  ##### parameters of normal of observed variable (i.e. z)

  ### vector of square root of diagonal of Omega
  diagOmegaOh = sqrt(diag(Omega))

  ### correlation matrix associated to Omega
  OmegaBar = cov2cor(Omega)

  ### Cholesky factor of OmegaBar
  invOmegaBar = chol(OmegaBar)

  ### logarithm of determinant of Omega
  logDetOmega =  2 * sum(
    log(diagOmegaOh) + log(diag(invOmegaBar))
  )

  ### inverse of OmegaBar
  invOmegaBar = chol2inv(invOmegaBar)

  ### inverse of diagOmegaOh
  invDiagOmegaOh = 1 / diagOmegaOh

  ### inverse of Omega
  invOmega = tcrossprod(invDiagOmegaOh) * invOmegaBar

  ### compute logarithm of inverse normalizing constant of Gaussian
  logInvNormConst_gau_PDF = -0.5 * d * log(2 * pi) - 0.5 * logDetOmega




  ########## parameters of SUT

  ### get nuSUT if not provided
  if (is.null(nuSUT)) nuSUT = max(m - d, 100)

  ### get the mode of the pSUN
  logitMode = c(
    pSUN_LogitGauss_Mode(A = A, b = b, xi = xi, Omega = Omega)$mode
  )

  ### Quasi Monte Carlo sequence
  dx = c(1:1e+2) / (1e+2 + 1)



  ##### if V is not provided then compute expected value of V given X equal to
  #      its mode
  if (is.null(V)) {

    ### get upper bounds of t given z equal to its mode
    upper = A %*% (invDiagOmegaOh * (logitMode - xi)) + c(b)

    ### initialize V
    V = double(length = m)

    ### get truncated logistic random variables via Quasi Monte Carlo
    t = qlogis(dx %o% plogis(upper))

    ### take absolute value of t
    t = abs(t)

    ### pre-calculus for mean of V given t
    t = 1 + exp(-t)

    ### compute mean of V given t
    t = t * t + t * (t + exp(-log(t - 1) + log(log(t))))

    ### compute mean of V given x equal to its mode
    V = colMeans(t)

  }

  ### get inverse square root of V
  invSqrt_V = sqrt(1 / c(V))

  ### get slant matrix of SUT
  ASUT = invSqrt_V * A

  ### get slant vector of SUT
  bSUT = invSqrt_V * b

  ### get the mode of the approximating SUT before affine transformation
  probitMode = c(pSUN_ProbitGauss_Mode(
    A = ASUT, b = bSUT, xi = xi, Omega = Omega, beta_start = logitMode
  )$mode)

  ### get location of the approximating SUT
  xiSUT = logitMode - probitMode




  ########## simulation and log-weights computation

  ##### parameters of truncated normal variable (i.e. eps = t - Az)

  ### minus cross-covariance matrix between z and eps
  cov_zEps = OmegaBar %*% t(ASUT)

  ### variance-covariance matrix of eps
  Sigma_eps = ASUT %*% cov_zEps + diag(1, nrow = m)

  ### get cross-covariance matrix between z and eps
  cov_zEps = -cov_zEps

  ### get parameters of the exponentially tilted measure
  expTilt_par = TruncatedNormal_mvrandn_getExpTilt(
    l = rep(-Inf, m), u = bSUT, Sig = Sigma_eps, mu = NULL,
    x0 = NULL
  )

  ### inverse of variance-covariance matrix
  # permuted version of inverse of Sigma_eps
  invSigma_eps = chol2inv(t(expTilt_par$output$Lfull))
  # get reverse permutation order
  rev_perm = sort(
    expTilt_par$output$perm,
    decreasing = FALSE, index.return = TRUE
  )$ix
  # reverse permutation
  invSigma_eps = invSigma_eps[rev_perm, rev_perm]



  ##### parameters of z given eps

  ### regression coefficients matrix of the mean
  regrMean_zGNeps = cov_zEps %*% invSigma_eps 

  ### decomposition of variance-covariance matrix with nugget
  decompVar_zGNeps = tryCatch(
    t(chol(
      OmegaBar - regrMean_zGNeps %*% t(cov_zEps) +
        diag(sqrt(.Machine$double.eps), d)
    )),
    error = function(err) {
      return(backsolve(
        r = chol(
          invOmegaBar + crossprod(ASUT) +
            diag(sqrt(.Machine$double.eps), d)
        ),
        x = diag(1, nrow = d)
      ))
    }
  )



  ##### sample from SUT and compute

  ### sample from truncated normal
  eps = TruncatedNormal_mvrandn_runExpTilt(
    mvrandnObj = expTilt_par, n = nsim, plim = plim
  )

  ### sample from standardized SUN
  z = regrMean_zGNeps %*% eps$sample + decompVar_zGNeps %*% matrix(
    rnorm(d * nsim), nrow = d, ncol = nsim
  )

  ### sample from standardized SUT
  z = t(sqrt(1 / rgamma(nsim, shape = 0.5 * nuSUT, rate = 0.5 * nuSUT)) * t(z))

  ### compute standardized argument for pSUN
  z_pSUN =  invDiagOmegaOh * (xiSUT - c(xi)) + z

  ### compute un-normalized log-density of Gaussian
  log_gau_PDF = - 0.5 * colSums(crossprod(invOmegaBar, z_pSUN) * z_pSUN)

  ### compute logarithm of logistic CDF
  log_logis_CDF = colSums(plogis(A %*% z_pSUN + c(b), log.p = TRUE))

  ### compute quadratic form of SUT
  quadForm_SUT = colSums(crossprod(invOmegaBar, z) * z)

  ### compute un-normalized log-density of Student t
  log_stu_PDF = - 0.5 * (nuSUT + d) * log(1 + quadForm_SUT / nuSUT)

  ### compute logarithm of non-central t CDF
  log_nct_CDF = pmvnctID(
    x = sqrt((nuSUT + d) / (nuSUT + quadForm_SUT)) * ASUT %*% z,
    nu = nuSUT + d,
    delta = -bSUT, dx = dx, n = nsim,
    maxComp = max((1 * 1e+9) %/% (d * length(dx) * 8), 1), # max GB * 1e+9
    log.p = TRUE
  )

  ### compute logarithm of inverse normalizing constant of Student t
  logInvNormConst_stu_PDF = lgamma(0.5 * (nuSUT + d)) - lgamma(0.5 * nuSUT) -
    0.5 * d * log(nuSUT * pi) - 0.5 * logDetOmega #logDetOmegaSUT

  ### compute log-ratio
  logWeights = log_gau_PDF + log_logis_CDF - log_stu_PDF - log_nct_CDF +
    eps$logNormConst + logInvNormConst_gau_PDF - logInvNormConst_stu_PDF




  ########## get results

  ### get maximum of log-weights
  max_logWeights = max(logWeights)

  ### compute exponential of scaled log-weights
  exp_SlogWeights = exp(logWeights - max_logWeights)

  ### compute logarithm of normalizing constant of pSUN
  logNormConst_pSUN = max_logWeights + log(sum(exp_SlogWeights)) - log(nsim)

  ### initialize
  res = list()

  ### compute final values of the sample
  res[["sample"]] = c(xiSUT) + c(diagOmegaOh) * z # diagOmegaOhSUT

  ### compute normalized weights
  res[["normWeights"]] = exp_SlogWeights / sum(exp_SlogWeights)

  ### get normalizing constant
  res[["logNormConst"]] = logNormConst_pSUN

  ### return result
  return(res)

}