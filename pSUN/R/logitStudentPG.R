#' @title Logit Model with Student-t Prior: Polya-Gamma Gibbs Sampler

#' @description
#' Polya-Gamma Gibbs sampler for the logit model with Student-t prior. It is
#'  used the code from package \pkg{BayesLogit} with slight modifications.

#' @usage
#' logitStudentPG = function(nsim, X, y, nu, mu, Sigma, x0 = NULL, burn = 1,
#'                           thin = 1, verbose = +Inf)

#' @details
#' \code{x0}: the starting point is generated in the following way:
#' \itemize{
#'  \item if \code{x0 == NULL} then the starting point for \code{beta} is the
#'   null vector,
#'  \item \code{else} set the starting point for \code{beta} equal to \code{x0}.
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
#'
#' The algorithm described in Polson et al. (2013) is used. The mixture density
#'  is a \eqn{\mathrm{inv.gamma} \big( \nu/2, \nu/2 \big)}.
#'
#' Pure R implementation of package \pkg{BayesLogit} is used for sampling from
#'  Polya-Gamma distribution.

#' @references
#' Polson, N. G., Scott, J. G., & Windle, J. (2013). Bayesian
#'  inference for logistic models using Pólya–Gamma latent variables. Journal
#'  of the American statistical Association, 108(504), 1339-1349.

#' @param nsim The number of simulations.
#' @param X The design matrix.
#' @param y The response vector.
#' @param nu The shape parameter.
#' @param mu The prior location vector.
#' @param Sigma The prior scale matrix.
#' @param x0 The starting point.
#' @param burn The number of draws to be discarded as burn-in.
#' @param thin The thinning parameter.
#' @param verbose The period for printing status of the chain.
#' @return Posterior sample of the parameters.

#' @export
logitStudentPG = function(
  nsim, X, y, nu, mu, Sigma, x0 = NULL, burn = 1, thin = 1, verbose = +Inf
) {

  ### print start time if required
  if (!is.infinite(verbose)) {
    print(paste("logitStudentPG: start time at ", Sys.time(), sep = ""))
  }

  ### check verbose
  if (verbose == 0) verbose = nsim + 1

  ### check burn
  if (burn <= 0) burn = 1

  ### check Sigma
  Sigma = matrix(Sigma, nrow = ncol(X))




  ########## compute useful parameters for sampling and initialize MCMC

  ### number of parameters
  p = ncol(X)

  ### number of observations
  n = nrow(X)

  ### inverse of the Cholesky factor of prior scale matrix
  invCholSigma = forwardsolve(
    l = t(chol(Sigma)),
    x = diag(1, nrow = p)
  )

  ### inverse of the prior scale matrix
  invSigma = crossprod(invCholSigma)

  ### compute t(X) %*% (y - 1/2)
  tXyOh = t(X) %*% (y-0.5)

  ### compute invSigma %*% mu
  invSigmaMu = invSigma %*% mu

  ### set beta
  if (is.null(x0)) beta = rep(0, p) else beta = x0

  ### store matrix
  sample = matrix(nrow = p, ncol = nsim)




  ########## run burn-in (if required)
  for (iburn in 1:burn) {

    ### update omega
    omega = rpg.devroye.R(n, 1, drop(X %*% beta))

    ### update W
    W = 1 / rgamma(
      n = 1, shape = 0.5 * (nu + p),
      rate = 0.5 * (nu + crossprod(invCholSigma %*% (beta - mu)))
    )



    ##### update beta

    ### decomposition of variance-covariance matrix of beta given omega and y
    #    (V_omega)
    # upper Cholesky factor of the inverse
    decompV_omega = chol(t(X) %*% (X * omega) + invSigma / W)
    # inverse of the upper Cholesky factor of the inverse
    decompV_omega = backsolve(
      r = decompV_omega,
      x = diag(1, nrow = p)
    )

    ### get V_omega
    V_omega = tcrossprod(decompV_omega)

    ### sample beta
    beta = V_omega %*% (tXyOh + invSigmaMu / W) + decompV_omega %*% rnorm(p)

  }




  ########## run scans to be collected
  for (insim in 1:nsim) {

    ### thinning
    for (ithin in 1:thin) {

      ### update omega
      omega = rpg.devroye.R(n, 1, drop(X %*% beta))

      ### update W
      W = 1 / rgamma(
        n = 1, shape = 0.5 * (nu + p),
        rate = 0.5 * (nu + crossprod(invCholSigma %*% (beta - mu)))
      )



      ##### update beta

      ### decomposition of variance-covariance matrix of beta given omega and y
      #   (V_omega)
      # upper Cholesky factor of the inverse
      decompV_omega = chol(t(X) %*% (X * omega) + invSigma / W)
      # inverse of the upper Cholesky factor of the inverse
      decompV_omega = backsolve(
        r = decompV_omega,
        x = diag(1, nrow = p)
      )

      ### get V_omega
      V_omega = tcrossprod(decompV_omega)

      ### sample beta
      beta = V_omega %*% (tXyOh + invSigmaMu / W) + decompV_omega %*% rnorm(p)



    }

    ### after thinning collect last values
    sample[, insim] = beta

    ### print status of the chain
    if (insim %% verbose == 0) {
      print(paste(
        "iteration ", insim, " of ", nsim, " completed at time ", Sys.time(),
        sep = ""
      ))
    }

  }




  ### print end time if required
  if (!is.infinite(verbose)) {
    print(paste("logitStudentPG: end time at ", Sys.time(), sep = ""))
  }

  ### return result
  return(sample)

}