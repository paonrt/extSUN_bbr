#' @title Logit Model with Dirichlet-Laplace Prior: Polya-Gamma Gibbs Sampler

#' @description
#' Polya-Gamma Gibbs sampler for the logit model with Dirichlet-Laplace prior as
#'  described in Bhattacharya at al. (2015). It is used the code from package
#'  \pkg{BayesLogit} with slight modifications.

#' @usage
#' logitDirLaplPG = function(nsim, X, y, mu, Sigma, x0 = NULL, burn = 1,
#'                           thin = 1, verbose = +Inf)

#' @details
#' \code{x0}: the starting point is generated in the following way:
#' \itemize{
#'  \item if \code{x0 == NULL} then sample \code{x0} from unit circle. Thus set
#'   the starting value of \code{beta} equal to \code{x0} times the square root
#'   of the expected value of the scaling random variable,
#'  \item \code{else} set the starting point for \code{beta} equal to \code{x0}.
#' }

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
#' The algorithm described in Polson et al. (2013) is used. A uniform discrete
#'  prior on the support \eqn{ \{ \frac{1}{300}, \frac{2}{300}, \dots,
#'  \frac{299}{300}, 1 \} } is used for the parameter of the Dirichlet-Laplace
#'  prior and, in the Gibbs sampler, its starting points is always \eqn{0.5}.
#'
#' Pure R implementation of package \pkg{BayesLogit} is used for sampling from
#'  Polya-Gamma distribution.

#' @references
#' Bhattacharya, A., Pati, D., Pillai, N. S., & Dunson, D. B. (2015).
#'  Dirichlet–Laplace priors for optimal shrinkage. Journal of the American
#'  Statistical Association, 110(512), 1479-1490.
#'
#' Polson, N. G., Scott, J. G., & Windle, J. (2013). Bayesian
#'  inference for logistic models using Pólya–Gamma latent variables. Journal
#'  of the American statistical Association, 108(504), 1339-1349.

#' @param nsim The number of simulations.
#' @param X The design matrix.
#' @param y The response vector.
#' @param mu The prior location vector.
#' @param Sigma The prior scale matrix, it must be diagonal.
#' @param x0 The starting point.
#' @param burn The number of draws to be discarded as burn-in.
#' @param thin The thinning parameter.
#' @param verbose The period for printing status of the chain.
#' @return Posterior sample of the parameters.

#' @import GIGrvg
#' @import statmod
#' @export
logitDirLaplPG = function(
  nsim, X, y, mu, Sigma,
  x0 = NULL, burn = 1, thin = 1, verbose = +Inf
) {

  ### print start time if required
  if (!is.infinite(verbose)) {
    print(paste("logitDirLaplPG: start time at ", Sys.time(), sep = ""))
  }

  ### check verbose
  if (verbose == 0) verbose = nsim+1

  ### check burn
  if (burn <= 0) burn = 1

  ### check Sigma
  # it must be a matrix
  Sigma = matrix(Sigma, nrow = ncol(X))
  # it must be diagonal
  if (!all(diag(diag(Sigma), nrow = nrow(Sigma)) == Sigma)) {
    stop("Sigma must be a diagonal matrix")
  }




  ########## compute useful parameters for sampling and initialize MCMC

  ### number of parameters
  p = ncol(X)

  ### number of observations
  n = nrow(X)

  ### diagonal of Sigma
  diagSigma = diag(Sigma)

  ### square root of the diagonal of Sigma
  diagSigmaOh = sqrt(diagSigma)

  ### inverse of the prior scale matrix
  invSigma = diag(1 / diagSigma, nrow = p)

  ### get discrete support of the prior of psi
  seq_psi = c(1:300) / 300

  ### compute constant values of log-density of full conditional of psi
  logConstFullPsi = -p * (seq_psi * log(2) + lgamma(seq_psi))

  ### compute t(X) %*% (y - 1/2)
  tXyOh = t(X) %*% (y - 0.5)

  ### compute invSigma %*% mu
  invSigmaMu = mu / diagSigma

  ### set beta (or sample it) and Z
  if (is.null(x0)) {
    x0 = rnorm(p)
    Z = sqrt(2 * mean(4 * seq_psi * (1 + seq_psi))) * x0 / sqrt(sum(x0^2))
    beta = mu + diagSigmaOh * Z
  } else {
    beta = x0
    Z = (beta - mu) / diagSigmaOh
  }

  ### set psi
  psi = seq_psi[ceiling(length(seq_psi) / 2)]

  ### sample zeta from its full conditional distribution
  # unnormalized zeta
  zeta = apply(matrix(2 * abs(Z)), 1, FUN = function(xappl) {
    tryCatch(
      expr = GIGrvg::rgig(n = 1, lambda = psi - 1, chi = xappl, psi = 1),
      error = function(err) {
        GIGrvg::rgig(
          n = 1, lambda = psi - 1,
          chi = xappl + sqrt(.Machine$double.eps),
          psi = 1
        )
      }
    )
  })
  # normalized zeta
  zeta = zeta / sum(zeta)

  ### sample tau from its full conditional distribution
  tau = tryCatch(
    GIGrvg::rgig(
      n = 1, lambda = (psi-1) * p, chi = 2 * sum(abs(Z) / zeta), psi = 1
    ),
    error = function(err){
      GIGrvg::rgig(
        n = 1, lambda = (psi - 1) * p,
        chi = 2 * sum(abs(Z) / zeta) + sqrt(.Machine$double.eps),
        psi = 1
      )
    }
  )

  ### store matrix
  sample = matrix(nrow = p, ncol = nsim)




  ########## run burn-in (if required)
  for (iburn in 1:burn) {

    ### update omega
    omega = rpg.devroye.R(n, 1, drop(X %*% beta))



    ##### update W

    ### update psi
    # compute weights of discrete full conditional of psi
    weiPsi = logConstFullPsi + (1 - p * seq_psi) * log(tau) + (seq_psi - 1) *
      sum(log(zeta))
    weiPsi = exp(weiPsi - max(weiPsi))
    weiPsi = weiPsi / sum(weiPsi)
    # sample psi
    psi = sample(seq_psi, 1, prob = weiPsi)

    ### update zeta
    # un-normalized zeta
    zeta = apply(matrix(2 * abs(Z)), 1, FUN = function(xappl) {
      tryCatch(
        expr = GIGrvg::rgig(n = 1, lambda = psi - 1, chi = xappl, psi = 1),
        error = function(err) {
          GIGrvg::rgig(
            n = 1, lambda = psi - 1,
            chi = xappl + sqrt(.Machine$double.eps),
            psi = 1
          )
        }
      )
    })
    # normalized zeta
    zeta = zeta / sum(zeta)

    ### update tau
    tau = tryCatch(
      GIGrvg::rgig(
        n = 1, lambda = (psi - 1) * p, chi = 2 * sum(abs(Z) / zeta), psi = 1
      ),
      error = function(err) {
        GIGrvg::rgig(
          n = 1, lambda = (psi - 1) * p,
          chi = 2 * sum(abs(Z) / zeta) + sqrt(.Machine$double.eps),
          psi = 1
        )
      }
    )

    ### update eta
    # inverse of eta
    eta = statmod::rinvgauss(n = p, mean = tau * zeta / abs(Z), shape = 1) +
      sqrt(.Machine$double.eps)
    # get eta
    eta = 1 / eta

    ### get W
    W = eta * (zeta * tau)^2



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

    ### get Z
    Z = (beta-mu) / diagSigmaOh

  }




  ########## run scans to be collected
  for (insim in 1:nsim) {

    ### thinning
    for (ithin in 1:thin) {

      ### update omega
      omega = rpg.devroye.R(n, 1, drop(X %*% beta))

      ##### update W

      ### update psi
      # compute weights of discrete full conditional of psi
      weiPsi = logConstFullPsi + (1 - p * seq_psi) * log(tau) + (seq_psi - 1) *
        sum(log(zeta))
      weiPsi = exp(weiPsi - max(weiPsi))
      weiPsi = weiPsi / sum(weiPsi)
      # sample psi
      psi = sample(seq_psi, 1, prob = weiPsi)

      ### update zeta
      # un-normalized zeta
      zeta = apply(matrix(2 * abs(Z)), 1, FUN = function(xappl) {
        tryCatch(
          expr = GIGrvg::rgig(n = 1, lambda = psi - 1, chi = xappl, psi = 1),
          error = function(err) {
            GIGrvg::rgig(
              n = 1, lambda = psi - 1,
              chi = xappl + sqrt(.Machine$double.eps),
              psi = 1
            )
          }
        )
      })
      # normalized zeta
      zeta = zeta / sum(zeta)

      ### update tau
      tau = tryCatch(
        GIGrvg::rgig(
          n = 1, lambda = (psi - 1) * p, chi = 2 * sum(abs(Z) / zeta), psi = 1
        ),
        error = function(err) {
          GIGrvg::rgig(
            n = 1, lambda = (psi - 1) * p,
            chi = 2 * sum(abs(Z) / zeta) + sqrt(.Machine$double.eps),
            psi = 1
          )
        }
      )

      ### update eta
      # inverse of eta
      eta = statmod::rinvgauss(n = p, mean = tau * zeta / abs(Z), shape = 1) +
        sqrt(.Machine$double.eps)
      # get eta
      eta = 1 / eta

      ### get W
      W = eta * (zeta * tau)^2



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

      ### get Z
      Z = (beta - mu) / diagSigmaOh



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
    print(paste("logitDirLaplPG: end time at ", Sys.time(), sep = ""))
  }

  ### return result
  return(sample)

}