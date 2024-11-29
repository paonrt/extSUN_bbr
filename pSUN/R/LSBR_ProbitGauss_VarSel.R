#' @title Probit Model with Gaussian Prior: Variables Selection

#' @description
#' Variable Selection using Bayes Factor and Metropolis steps in a Linear
#'  Symmetric Binary Regression (LSBR) model space.

#' @usage
#' LSBR_ProbitGauss_VarSel = function(nsim, X, y, mu, Sigma,
#'                                    ntim = 1e+4, x0 = NULL, verbose = +Inf)

#' @details
#' A posterior sample from the parameter space is drawn and a
#'  Metropolis-within-Gibbs algorithm is used. Given a model, the marginal data
#'  density is computed by using the procedure described in Durante (2019).
#'  The update from the model space is the same as in Onorati and Liseo (2022).
#'
#' The prior on model space is the discrete uniform on models including the
#'  intercept. Models without the intercept are not included in the space model.

#' The prior on parameters is a Gaussian distribution with mean vector \code{mu}
#'  and variance-covariance matrix \code{Sigma}.
#'
#' If \code{x0 == NULL} then the starting point is get using the posterior
#'  credible interval at 0.95 level. The i-th variable is included in the
#'  starting point if zero is not included
#'
#' The current time and the current number of iteration are printed one every
#' \code{verbose} iterations. Furthermore:
#' \itemize{
#'  \item if \code{verbose == +-Inf} then there is no printing,
#'  \item if \code{verbose != +-Inf} then at least start and end of simulation
#'   are reported.
#' }

#' @references
#' Durante, D. (2019). Conjugate Bayes for probit regression via unified
#'  skew-normal distributions. Biometrika, 106(4), 765-779.
#'
#' Onorati, P., & Liseo, B. (2022). An extension of the unified skew-normal
#'  family of distributions and application to Bayesian binary regression.
#'  arXiv preprint arXiv:2209.03474.

#' @param nsim The number of simulations.
#' @param X The design matrix.
#' @param y The response vector.
#' @param mu The prior location vector.
#' @param Sigma The prior scale matrix.
#' @param ntim The number of draws for the normalizing constant computation.
#' @param x0 The starting point.
#' @param verbose The period for printing status of the chain.
#' @return Posterior inclusion probabilities and a posterior sample from model
#'  space.

#' @export
LSBR_ProbitGauss_VarSel = function(
  nsim, X, y, mu, Sigma, ntim = 1e+4, x0 = NULL, verbose = +Inf
) {

  ### print start time if required
  if (!is.infinite(verbose)) {
    print(paste(
      "LSBR_ProbitGauss_VarSel: start time at ", Sys.time(), sep = ""
    ))
  }

  ### check
  if (!all(X[, 1] == 1)) {
    stop("The first column of X must be composed only by 1")
  }

  ### check verbose
  if (verbose == 0) verbose = nsim + 1

  ### check mu
  mu = matrix(mu, ncol = 1)

  ### check Sigma
  Sigma = matrix(Sigma, nrow = ncol(X))




  ########## compute useful parameters for sampling

  ### number of the parameters
  p = ncol(X)

  ### number of observations
  n = nrow(X)

  ### likelihood matrix
  B = 2 * diag(y) - diag(1, nrow = length(y))

  ### parameters of the SUN distribution
  # slant matrix
  A = B %*% X %*% diag(sqrt(diag(Sigma)), nrow = ncol(X))
  # slant vector
  b = B %*% X %*% mu
  # location vector
  xi = mu
  # scale matrix
  Omega = Sigma




  ########## initialize MCMC

  ### posterior sample matrix
  sample = matrix(nrow = p, ncol = nsim)



  ##### starting point

  ### if it is null then get from posterior credible interval at 0.95
  if (is.null(x0)) {
    currMod = logical(p)
    currMod[1] = TRUE
    x0 = pSUN_ProbitGauss_RNG(
      ntim, A = A, b = b, xi = xi, Omega = Omega
    )
    for (i in 2:p) {
      credInt = quantile(
        x0$sample[i, ], probs = c(0.025, 0.975)
      )
      if (credInt[1] <= 0 && credInt[2] >= 0) {
        currMod[i] = FALSE
      } else {
        currMod[i] = TRUE
      }
    }
  }

  ### get model size
  currModSize = sum(currMod)

  ### compute logarithm of posterior density
  logPostCurrMod = pSUN_ProbitGauss_RNG(
    ntim,
    A = matrix(A[, currMod], ncol = currModSize),
    b = b, xi = xi[currMod],
    Omega = matrix(Omega[currMod, currMod], ncol = currModSize)
  )$logNormConst



  #### initialize proposal
  propMod = currMod




  ########## run MCMC
  for (insim in 1:nsim) {

    ### Metropolis steps
    for (imetro in 2:p) {



      ##### proposal

      ### get proposal
      propMod[imetro] = !currMod[imetro]

      ### get model size
      propModSize = sum(propMod)

      ### compute logarithm of posterior density
      logPostPropMod = pSUN_ProbitGauss_RNG(
        ntim, A = matrix(A[, propMod], ncol = propModSize),
        b = b, xi = xi[propMod],
        Omega = matrix(Omega[propMod,propMod], ncol = propModSize)
      )$logNormConst



      ##### acceptance test

      ### compute acceptance log-probability
      accLogprob = logPostPropMod - logPostCurrMod

      ### sample from exponential
      expo = rgamma(n = 1, shape = 1, rate = 1)

      ### check
      if (expo >= -accLogprob) {
        # update
        currMod = propMod
        currModSize = propModSize
        logPostCurrMod = logPostPropMod 
      } else {
        # reverse proposal
        propMod[imetro] = !currMod[imetro]
      }



    }

    ### after Metropolis steps collect last value
    sample[, insim] = currMod

    ### print status of the chain
    if (insim %% verbose == 0) {
      print(paste(
        "iteration ", insim, " of ", nsim, " completed at time ", Sys.time(),
        sep = ""
      ))
    }

  }




  ########## get results

  ### compute posterior inclusion probabilities
  postInclProbs = rowMeans(sample)

  ### print end time if required
  if (!is.infinite(verbose)) {
    print(paste("LSBR_ProbitGauss_VarSel: end time at ", Sys.time(), sep = ""))
  }

  ### return results
  return(list(
    sample = sample,
    postInclProbs = postInclProbs
  ))

}