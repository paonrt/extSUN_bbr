#' @title Perturbed Unified Skew Normal Probit Gaussian: Random Number Generator

#' @description
#' Exact random number generator for a perturbed unified skew normal (pSUN)
#'  probit Gaussian distribution. It is actually a SUN distribution
#'  (Arellano-Valle and Azzalini, 2006) with a different parametrization.

#' @usage
#' pSUN_ProbitGauss_RNG = function(nsim, A, b, xi, Omega, plim = 1 / nsim)

#' @details
#' It is said \eqn{X} follows a pSUN probit Gaussian distribution if
#'  \deqn{
#'   X = \xi + \mathrm{diag}^\frac{1}{2} (\Omega) \, Z \vert T \le A Z + b \\
#'   Z \vert W \sim N(0, \bar{\Omega}_W) \perp\!\!\!\perp T \vert V \sim N(0,
#'    \Theta_V) \\
#'   W_i = 1 \, , \, V_i = 1
#'  }
#'  where \eqn{\bar{\Omega}_W = \mathrm{diag}^\frac{1}{2}(W) \bar{\Omega}
#'  \mathrm{diag}^\frac{1}{2}(W)}, \eqn{\Theta_V = \mathrm{diag}^\frac{1}{2}(V)
#'  \Theta \mathrm{diag}^\frac{1}{2}(V)}, \eqn{\bar{\Omega}} is the correlation
#'  matrix of \eqn{\Omega}.
#'
#' Draws from truncated Gaussian are obtained using Botev (2017) algorithm.

#' @references
#' Arellano‐Valle, R. B., & Azzalini, A. (2006). On the unification of families
#'  of skew‐normal distributions. Scandinavian Journal of Statistics, 33(3),
#'  561-574.
#'
#' Botev, Z. I. (2017). The normal law under linear restrictions: simulation and
#'  estimation via minimax tilting. Journal of the Royal Statistical Society
#'  Series B: Statistical Methodology, 79(1), 125-148.

#' @param nsim The number of simulations.
#' @param A The slant matrix.
#' @param b The slant vector.
#' @param xi The location vector.
#' @param Omega The scale matrix.
#' @param plim Minimum empirical acceptance probability such that the
#'  function does not stop with an error.
#' @return Sample and logarithm of normalizing constant of
#'  the pSUN distribution.

#' @export
pSUN_ProbitGauss_RNG = function(nsim, A, b, xi, Omega, plim = 1 / nsim) {

  ### check xi
  xi = matrix(xi, ncol = 1)

  ### check b
  b = matrix(b, ncol = 1)

  ### check Omega
  Omega = matrix(Omega, nrow = ncol(A))




  ########## compute useful parameters for sampling

  ### dimension of the distribution
  d = ncol(A)

  ### dimension of selective latent variable
  m = nrow(A)



  ##### parameters of normal of observed variable (i.e. z)

  ### vector of square root of diagonal of Omega
  diagOmegaOh = sqrt(diag(Omega))

  ### correlation matrix associated to Omega
  OmegaBar = cov2cor(Omega)



  ##### parameters of truncated normal variable (i.e. eps = t - Az)

  ### minus cross-covariance matrix between z and eps
  cov_zEps = OmegaBar %*% t(A)

  ### variance-covariance matrix of eps
  Sigma_eps = A %*% cov_zEps + diag(1, nrow = m)

  ### inverse of variance-covariance matrix
  invSigma_eps = chol2inv(chol(Sigma_eps))

  ### get cross-covariance matrix between z and eps
  cov_zEps = - cov_zEps



  ##### parameters of z given eps

  ### regression coefficients matrix of the mean
  regrMean_zGNeps = cov_zEps %*% invSigma_eps

  ### decomposition of variance-covariance matrix with nugget
  decompVar_zGNeps = tryCatch(
    t(chol(
      OmegaBar - regrMean_zGNeps %*% t(cov_zEps) +
        diag(sqrt(.Machine$double.eps), d) 
    )),
    error = function(err){ # if chol fails then use svd
      #svdObj = svd(
      # OmegaBar - regrMean_zGNeps %*% t(cov_zEps) +
      #  diag(sqrt(.Machine$double.eps), d),
      # nv = 0
      #)
      #return(t(sqrt(svdObj$d) * t(svdObj$u)))
      return(backsolve(
        r = chol(
          chol2inv(chol(OmegaBar)) + crossprod(A) +
            diag(sqrt(.Machine$double.eps), d)
        ),
        x = diag(1, nrow = d)
      ))
    }
  )




  ########## sampling and get log-normalizing constant

  ##### sample eps

  ### get parameters of the exponentially tilted measure
  expTilt_par = TruncatedNormal_mvrandn_getExpTilt(
    l = rep(-Inf, m), u = b, Sig = Sigma_eps, mu = NULL
  )
  ### run simulations of the exponentially tilted measure
  eps = TruncatedNormal_mvrandn_runExpTilt(
    mvrandnObj = expTilt_par, n = nsim, plim = plim
  )



  ### sample Z
  Z = regrMean_zGNeps %*% eps$sample + decompVar_zGNeps %*% matrix(
    rnorm(d * nsim), nrow = d, ncol = nsim
  )




  ########## get results

  ### initialize
  res = list()

  ### compute final values of the sample
  res[["sample"]] = c(xi) + diagOmegaOh * Z

  ### get normalizing constant
  res[["logNormConst"]] = eps$logNormConst

  ### return result
  return(res)

}