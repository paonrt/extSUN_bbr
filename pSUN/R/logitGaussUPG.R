#' @title Logit Model with Gaussian Prior: Ultimate Polya-Gamma Gibbs Sampler

#' @description
#' Ultimate Polya-Gamma Gibbs sampler for the logit model with Gaussian prior.
#' It is used the code from package \pkg{UPG} with slight modifications.

#' @usage
#' logitGaussUPG = function(nsim, X, y, Sigma, x0 = NULL, burn = 1,
#'                          thin = 1, verbose = +Inf)

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
#' The algorithm described in Zens et al. (2023) is used.
#'
#' Pure R implementation of package \pkg{BayesLogit} is used for sampling from
#'  Polya-Gamma distribution.

#' @references
#' Zens, G., Frühwirth-Schnatter, S., & Wagner, H. (2023). Ultimate Pólya Gamma
#'  Samplers–Efficient MCMC for possibly imbalanced binary and categorical
#'  data. Journal of the American Statistical Association, 1-12.

#' @param nsim The number of simulations.
#' @param X The design matrix.
#' @param y The response vector.
#' @param Sigma The prior scale matrix, it must be diagonal.
#' @param x0 The starting point.
#' @param burn The number of draws to be discarded as burn-in.
#' @param thin The thinning parameter.
#' @param verbose The period for printing status of the chain.
#' @return Posterior sample of the parameters.

#' @export
logitGaussUPG = function(
  nsim, X, y, Sigma,
  x0 = NULL, burn = 1, thin = 1, verbose = +Inf
) {


  ### print start time if required
  if (!is.infinite(verbose)) {
    print(paste("logitGaussUPG: start time at ", Sys.time(), sep = ""))
  }

  ### check verbose
  if (verbose == 0) verbose = nsim + 1

  ### check burn
  if (burn <= 0) burn = 1

  ### check Sigma
  # it must be a matrix
  Sigma = matrix(Sigma, nrow = ncol(X))
  # it must be diagonal
  if (!all(diag(diag(Sigma), nrow = nrow(Sigma)) == Sigma)) {
    stop("Sigma must be a diagonal matrix")
  }

  ### prior shape for working parameter delta, 2.5 is the default.
  d0 = 2.5

  ###	prior rate for working parameter delta, 1.5 is the default.
  D0 = 1.5

  ### prior variance for working parameter gamma, 100 is the default.
  G0 = 100




  ########## compute useful parameters for sampling and initialize MCMC

  ### number of parameters
  p = ncol(X)

  ### number of observations
  n = nrow(X)

  ### diagonal of Sigma
  diagSigma = diag(Sigma)

  ### inverse of the prior scale matrix
  invSigma = diag(1 / diagSigma, nrow = p)

  ### functions for sampling latent utilities
  Fe_ = function(p){log(p)-log(1-p)}
  Fe  = function(x) plogis(x)#{exp(x) / (1+exp(x))}

  ### some precomputations
  Sy0   = sum(y) == 0
  SyN   = sum(y) == n
  y0    = y == 0
  y1    = y == 1

  ### set beta
  if (is.null(x0)) beta = rep(0, p) else beta = x0

  ### store matrix
  sample = matrix(nrow = p, ncol = nsim)




  ########## run burn-in (if required)
  for (iburn in 1:burn) {

    ### update latent utilities
    U.i = runif(n)
    lambda = X %*% beta
    z = lambda + Fe_(y + U.i * (1 - y - (1 - Fe(-lambda))))

    ### update omega
    #omega = rpg.devroye.R(n, 1, drop(X %*% beta))
    #omega  = pgdraw::pgdraw(2, abs(z - X %*% beta.draw))
    omega = rpg.devroye.R(n, 2, z - drop(X %*% beta))

    ### draw from working prior
    gamma.star = rnorm(1, 0, sqrt(G0))

    ### shift utilities
    ztilde = z + gamma.star

    ### posterior quantities for gamma conditional
    XW         = X * omega
    tXW        = t(XW)
    Bn         = chol2inv(chol(invSigma + tXW %*% X))
    mb         = colSums(XW)
    mg         = sum(ztilde * omega)
    mn         = tXW %*% ztilde
    tmbBn      = t(mb) %*% Bn
    Gn         = 1 / ( (1 / G0) + sum(omega) - tmbBn %*% mb)
    gn         = Gn * (mg - tmbBn %*% mn)

    ### truncation points for gamma conditional
    U          = ifelse(Sy0,  Inf,  min(ztilde[y1]))
    L          = ifelse(SyN, -Inf,  max(ztilde[y0]))

    ### simulate from gamma posterior
    gamma = truncnorm::rtruncnorm(1, a = L, b = U, mean = gn, sd = sqrt(Gn))

    ### reverse shift
    z.shift = ztilde - gamma

    ### draw from working prior
    delta.star = 1 / rgamma(1, d0, D0)

    ### posterior quantities for delta conditional
    mn = tXW %*% z.shift
    bn = Bn %*% mn

    ### simulate from delta posterior
    delta = 1 / rgamma(
      1, d0 + 0.5 * n,
      D0 + 0.5 * delta.star * (t(bn) %*% invSigma %*% bn + sum(
        omega * (z.shift - X %*% bn)^2
      ))
    )

    ### sample coefficient
    sqrtBn    = t(chol(Bn))
    beta = sqrt(delta.star / delta) * bn + sqrtBn %*% rnorm(p)

  }




  ########## run scans to be collected
  for (insim in 1:nsim) {

    ### thinning
    for (ithin in 1:thin) {

      ### update latent utilities
      U.i = runif(n)
      lambda = X %*% beta
      z = lambda + Fe_(y + U.i * (1 - y - (1 - Fe(-lambda))))

      ### update omega
      #omega = rpg.devroye.R(n, 1, drop(X %*% beta))
      #omega  = pgdraw::pgdraw(2, abs(z - X %*% beta.draw))
      omega = rpg.devroye.R(n, 2, z - drop(X %*% beta))

      ### draw from working prior
      gamma.star = rnorm(1, 0, sqrt(G0))

      ### shift utilities
      ztilde = z + gamma.star

      ### posterior quantities for gamma conditional
      XW         = X * omega
      tXW        = t(XW)
      Bn         = chol2inv(chol(invSigma + tXW %*% X))
      mb         = colSums(XW)
      mg         = sum(ztilde * omega)
      mn         = tXW %*% ztilde
      tmbBn      = t(mb) %*% Bn
      Gn         = 1 / ( (1 / G0) + sum(omega) - tmbBn %*% mb)
      gn         = Gn * (mg - tmbBn %*% mn)

      ### truncation points for gamma conditional
      U          = ifelse(Sy0,  Inf,  min(ztilde[y1]))
      L          = ifelse(SyN, -Inf,  max(ztilde[y0]))

      ### simulate from gamma posterior
      gamma = truncnorm::rtruncnorm(1, a = L, b = U, mean = gn, sd = sqrt(Gn))

      ### reverse shift
      z.shift = ztilde - gamma

      ### draw from working prior
      delta.star = 1 / rgamma(1, d0, D0)

      ### posterior quantities for delta conditional
      mn = tXW %*% z.shift
      bn = Bn %*% mn

      ### simulate from delta posterior
      delta = 1 / rgamma(
        1, d0 + 0.5 * n,
        D0 + 0.5 * delta.star * (t(bn) %*% invSigma %*% bn + sum(
          omega * (z.shift - X %*% bn)^2
        ))
      )

      ### sample coefficient
      sqrtBn    = t(chol(Bn))
      beta = sqrt(delta.star / delta) * bn + sqrtBn %*% rnorm(p)



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
    print(paste("logitGaussUPG: end time at ", Sys.time(), sep = ""))
  }

  ### return result
  return(sample)

}