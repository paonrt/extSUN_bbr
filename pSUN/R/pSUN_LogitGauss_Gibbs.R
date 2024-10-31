#' @title Perturbed Unified Skew Normal Logit Gaussian: Gibbs Sampler

#' @description
#' Gibbs Sampler for a perturbed unified skew normal (pSUN) logit Gaussian
#'  distribution.

#' @usage
#' pSUN_LogitGauss_Gibbs = function(nsim, A, b, xi, Omega, x0 = NULL,
#'                                  burn = 1, thin = 1, plim = 1 / nsim,
#'                                  verbose = +Inf)

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
#' Draws from truncated Gaussian are obtained using Botev (2017) algorithm. See
#'  Onorati and Liseo (2022) for more details.
#'
#' \code{x0}: the starting point is generated in the following way:
#' \itemize{
#'  \item if \code{x0 == NULL} then the starting point for \code{X} is the
#'   null vector,
#'  \item \code{else} set the starting point for \code{X} equal to \code{x0}.
#' }
#'
#' Only one value every \code{thin} values is kept in the chain, so the true
#'  number of complete scans will be \code{nsim * thin + burn}. By default
#'  \code{thin = 1}, that is no thinning.
#'
#' The current time and the current number of iteration are printed one every
#' \code{verbose} iterations. Furthermore:
#' \itemize{
#'  \item if \code{verbose == +-Inf} then there is no printing,
#'  \item if \code{verbose != +-Inf} then at least start and end of simulation
#'   are reported.
#' }

#' @references
#' Botev, Z. I. (2017). The normal law under linear restrictions: simulation and
#'  estimation via minimax tilting. Journal of the Royal Statistical Society
#'  Series B: Statistical Methodology, 79(1), 125-148.
#'
#' Onorati, P., & Liseo, B. (2022). An extension of the unified skew-normal
#'  family of distributions and application to Bayesian binary regression.
#'  arXiv preprint arXiv:2209.03474.

#' @param nsim The number of simulations.
#' @param A The slant matrix.
#' @param b The slant vector.
#' @param xi The location vector.
#' @param Omega The scale matrix.
#' @param x0 The starting point.
#' @param burn The number of draws to be discarded as burn-in.
#' @param thin The thinning parameter.
#' @param plim Minimum empirical acceptance probability such that the
#'  function does not stop with an error.
#' @param verbose The period for printing status of the chain.
#' @return Sample of the pSUN distribution.

#' @export
pSUN_LogitGauss_Gibbs = function(
  nsim, A, b, xi, Omega,
  x0 = NULL, burn = 1, thin = 1, plim = 1 / nsim, verbose = +Inf
) {

  ### print start time if required
  if (!is.infinite(verbose)) {
    print(paste("pSUN_LogitGauss_Gibbs: start time at ", Sys.time(), sep = ""))
  }

  ### check verbose
  if (verbose == 0) verbose = nsim + 1

  ### check burn
  if (burn <= 0) burn = 1

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



  ##### parameters of normal of observed variable (i.e. Z)

  ### vector of square root of diagonal of Omega
  diagOmegaOh = sqrt(diag(Omega))

  ### correlation matrix associated to Omega
  OmegaBar = cov2cor(Omega)

  ### inverse of OmegaBar
  invOmegaBar = chol2inv(chol(OmegaBar))

  ### minus cross-covariance matrix between Z and eps
  cov_zEps = OmegaBar %*% t(A)

  ### variance-covariance matrix of linear transformation Az
  AOmegaBartA = A %*% cov_zEps

  ### get cross-covariance matrix between Z and eps
  cov_zEps = -cov_zEps




  ########## initialize MCMC

  ### set Z
  if (is.null(x0)) Z = rep(0, d) else Z = (x0 - xi) / diagOmegaOh

  ### sample T
  T = qlogis(runif(m) * plogis(A %*% Z + b))

  ### compute eps
  eps = T - A %*% Z

  ### store matrix
  sample = matrix(nrow = d, ncol = nsim)

  ### starting point for exponential tilting optimization
  expTilt_par = list(output = list(x = NULL, muV = NULL))




  ########## run burn-in (if required)
  for (iburn in 1:burn) {

    ### update V
    V = c(rpostlogiskolmo(n = m, x = T))



    ##### update T and Z

    ### sample eps
    # variance-covariance matrix
    Sigma_eps = diag(V) + AOmegaBartA
    # get parameters of the exponentially tilted measure
    expTilt_par = TruncatedNormal_mvrandn_getExpTilt(
      l = rep(-Inf, m), u = b, Sig = Sigma_eps, mu = NULL,
      x0 = c(expTilt_par$output$x, expTilt_par$output$muV)
    )
    # run simulation of the exponentially tilted measure
    eps = TruncatedNormal_mvrandn_runExpTilt(
      mvrandnObj = expTilt_par, n = 1, plim = plim
    )$sample

    ### update Z
    # permuted version of inverse of Sigma_eps
    invSigma_eps = chol2inv(t(expTilt_par$output$Lfull))
    # get reverse permutation order
    rev_perm = sort(
      expTilt_par$output$perm,
      decreasing = FALSE, index.return = TRUE
    )$ix
    # reverse permutation
    invSigma_eps = invSigma_eps[rev_perm, rev_perm]
    # regression coefficients matrix of the mean
    regrMean_zGNeps = cov_zEps %*% invSigma_eps
    # decomposition of variance-covariance matrix with nugget
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
            invOmegaBar + crossprod(A, A / sqrt(V)) +
              diag(sqrt(.Machine$double.eps), d)
          ),
          x = diag(1, nrow = d)
        ))
      }
    )
    # get Z
    Z = regrMean_zGNeps %*% eps + decompVar_zGNeps %*% rnorm(d)

    ### update T 
    T = eps + A %*% Z



  }




  ########## run scans to be collected
  for (insim in 1:nsim) {

    ### thinning
    for (ithin in 1:thin) {

      ### update V
      V = c(rpostlogiskolmo(n = m, x = T))



      ##### update T and Z

      ### sample eps
      # variance-covariance matrix
      Sigma_eps = diag(V) + AOmegaBartA
      # get parameters of the exponentially tilted measure
      expTilt_par = TruncatedNormal_mvrandn_getExpTilt(
        l = rep(-Inf, m), u = b, Sig = Sigma_eps, mu = NULL,
        x0 = c(expTilt_par$output$x, expTilt_par$output$muV)
      )
      # run simulation of the exponentially tilted measure
      eps = TruncatedNormal_mvrandn_runExpTilt(
        mvrandnObj = expTilt_par, n = 1, plim = plim
      )$sample

      ### update Z
      # permuted version of inverse of Sigma_eps
      invSigma_eps = chol2inv(t(expTilt_par$output$Lfull))
      # get reverse permutation order
      rev_perm = sort(
        expTilt_par$output$perm,
        decreasing = FALSE, index.return = TRUE
      )$ix
      # reverse permutation
      invSigma_eps = invSigma_eps[rev_perm, rev_perm]
      # regression coefficients matrix of the mean
      regrMean_zGNeps = cov_zEps %*% invSigma_eps
      # decomposition of variance-covariance matrix with nugget
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
              invOmegaBar + crossprod(A, A / sqrt(V)) +
                diag(sqrt(.Machine$double.eps), d)
            ),
            x = diag(1, nrow = d)
          ))
        }
      )
      # get Z
      Z = regrMean_zGNeps %*% eps + decompVar_zGNeps %*% rnorm(d)

      ### update T 
      T = eps + A %*% Z



    }

    ### after thinning collect last values
    sample[, insim] = Z

    ### print status of the chain
    if (insim %% verbose == 0) {
      print(paste(
        "iteration ", insim, " of ", nsim, " completed at time ", Sys.time(),
        sep = ""
      ))
    }

  }




  ########## get results

  ### affine transformation of Z
  sample = c(xi) + diagOmegaOh * sample

  ### print end time if required
  if (!is.infinite(verbose)) {
    print(paste("pSUN_LogitGauss_Gibbs: end time at ", Sys.time(), sep = ""))
  }

  ### return result
  return(sample)

}