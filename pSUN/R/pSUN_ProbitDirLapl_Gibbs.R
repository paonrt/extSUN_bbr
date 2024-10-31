#' @title Perturbed Unified Skew Normal Probit Dirichlet-Laplace: Gibbs Sampler

#' @description
#' Gibbs Sampler for a perturbed unified skew normal (pSUN) probit
#'  Dirichlet-Laplace distribution.

#' @usage
#' pSUN_ProbitDirLapl_Gibbs = function(nsim, A, b, xi, Omega, x0 = NULL,
#'                                    burn = 1, thin = 1, plim = 1 / nsim,
#'                                    verbose = +Inf)

#' @details
#' It is said \eqn{X} follows a pSUN probit Dirichlet-Laplace distribution if
#'  \deqn{
#'   X = \xi + \mathrm{diag}^\frac{1}{2} (\Omega) \, Z \vert T \le A Z + b \\
#'   Z \vert W \sim N(0, \bar{\Omega}_W) \perp\!\!\!\perp T \vert V \sim N(0,
#'    \Theta_V) \\
#'   W_i \overset{i.i.d.}{\sim} \mathrm{SDL}(\cdot; \psi) \, , \,
#'    V_i = 1
#'  }
#'  where \eqn{\bar{\Omega}_W = \mathrm{diag}^\frac{1}{2}(W) \bar{\Omega}
#'  \mathrm{diag}^\frac{1}{2}(W)}, \eqn{\Theta_V = \mathrm{diag}^\frac{1}{2}(V)
#'  \Theta \mathrm{diag}^\frac{1}{2}(V)}, \eqn{\bar{\Omega}} is the correlation
#'  matrix of \eqn{\Omega}. \eqn{\mathrm{SDL}(\cdot; \psi)} is the CDF of the
#'  scaling mixture distribution of the Dirichlet-Laplace prior with parameter
#'  \eqn{\psi} (Bhattacharya at al., 2015).
#'
#' In this function, \eqn{\Theta} is restricted to be the identity matrix.
#'  Furthermore, \eqn{W_i \overset{i.i.d.}{\sim}
#'  \mathrm{SDL}(\cdot; \psi)} implies
#'  \deqn{
#'   W_i =  G^2_i E_i \, , \, G_i \overset{i.i.d.}{\sim}
#'   \mathrm{gamma}(\psi, 1 / 2) \perp\!\!\!\perp E_i \overset{i.i.d.}{\sim}
#'   \exp(1 / 2)
#'  }
#'  or, equivalently
#'  \deqn{
#'   W = \widetilde{G}_i U E \\
#'   \widetilde{G}_i \overset{i.i.d.}{\sim} \mathrm{gamma}(d \psi, 1 / 2) \\
#'   U \sim \mathrm{Dirichlet}_d(\psi, \psi, \dots, \psi) \\
#'   E_i \overset{i.i.d.}{\sim} \exp(1 / 2)
#'  }
#'  a uniform discrete prior on the support \eqn{ \{ \frac{1}{300},
#'  \frac{2}{300}, \dots, \frac{299}{300}, 1 \} } is used for \eqn{\psi}.
#'
#' Draws from truncated Gaussian are obtained using Botev (2017) algorithm. See
#'  Onorati and Liseo (2022) for more details.
#'
#' \code{x0}: the starting point is generated in the following way:
#' \itemize{
#'  \item if \code{x0 == NULL} then sample \code{x0} from unit circle. Thus set
#'   the starting value of \code{beta} equal to \code{x0} times the square root
#'   of the expected value of the scaling random variable,
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

#' @references
#' Bhattacharya, A., Pati, D., Pillai, N. S., & Dunson, D. B. (2015).
#'  Dirichlet–Laplace priors for optimal shrinkage. Journal of the American
#'  Statistical Association, 110(512), 1479-1490.
#'
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
#' @param Omega The scale matrix, it must be diagonal.
#' @param x0 The starting point.
#' @param burn The number of draws to be discarded as burn-in.
#' @param thin The thinning parameter.
#' @param plim Minimum empirical acceptance probability such that the
#'  function does not stop with an error.
#' @param verbose The period for printing status of the chain.
#' @return Sample of the pSUN distribution.

#' @import GIGrvg
#' @import statmod
#' @import truncnorm
#' @export
pSUN_ProbitDirLapl_Gibbs = function(
  nsim, A, b, xi, Omega,
  x0 = NULL, burn = 1, thin = 1, plim = 1 / nsim, verbose = +Inf
) {

  ### print start time if required
  if (!is.infinite(verbose)) {
    print(paste(
      "pSUN_ProbitDirLapl_Gibbs: start time at ", Sys.time(), sep = ""
    ))
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
  # it must be a matrix
  Omega = matrix(Omega, nrow = ncol(A))
  # it must be diagonal
  if (!all(diag(diag(Omega), nrow = nrow(Omega)) == Omega)) {
    stop("Omega must be a diagonal matrix")
  }




  ########## compute useful parameters for sampling

  ### dimension of the distribution
  d = ncol(A)

  ### dimension of selective latent variable
  m = nrow(A)

  ### vector of square root of diagonal of Omega
  diagOmegaOh = sqrt(diag(Omega))

  ### inner product of matrix A
  tAA = t(A) %*% A

  ### get discrete support of the prior of psi
  seq_psi = c(1:300) / 300

  ### compute constant values of log-density of full conditional of psi
  logConstFullPsi = -d * (seq_psi * log(2) + lgamma(seq_psi))




  ########## initialize MCMC

  ### set Z (or sample it)
  if (is.null(x0)) {
    x0 = rnorm(d)
    Z = sqrt(2 * mean(4 * seq_psi * (1 + seq_psi))) * x0 / sqrt(sum(x0^2))
  } else {
    Z = (x0 - xi) / diagOmegaOh
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
      n = 1, lambda = (psi - 1) * d, chi = 2 * sum(abs(Z) / zeta), psi = 1
    ),
    error = function(err) {
      GIGrvg::rgig(
        n = 1, lambda = (psi - 1) * d,
        chi = 2 * sum(abs(Z) / zeta) + sqrt(.Machine$double.eps),
        psi = 1
      )
    }
  )

  ### sample T
  T = truncnorm::rtruncnorm(n = m, b = A %*% Z + b)

  ### compute eps
  eps = T - A %*% Z

  ### store matrix
  sample = matrix(nrow = d, ncol = nsim)

  ### starting point for exponential tilting optimization
  expTilt_par = list(output = list(x = NULL, muV = NULL ))




  ########## run burn-in (if required)
  for (iburn in 1:burn) {

    ##### update W

    ### update psi
    # compute weights of discrete full conditional of psi
    weiPsi = logConstFullPsi + (1 - d * seq_psi) * log(tau) +
      (seq_psi - 1)*sum(log(zeta))
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
        n = 1, lambda = (psi - 1) * d, chi = 2 * sum(abs(Z) / zeta), psi = 1
      ),
      error = function(err) {
        GIGrvg::rgig(
          n = 1, lambda = (psi - 1) * d,
          chi = 2 * sum(abs(Z) / zeta) + sqrt(.Machine$double.eps),
          psi = 1
        )
      }
    )

    ### update eta
    # inverse of eta
    eta = statmod::rinvgauss(n = d, mean = tau * zeta / abs(Z), shape = 1) +
      sqrt(.Machine$double.eps)
    # get eta
    eta = 1 / eta

    ### get W
    W = eta * (zeta * tau)^2



    ##### update T and Z

    ### sample eps
    # minus cross-covariance matrix between Z and eps
    cov_zEps = W * t(A)
    # variance-covariance matrix of linear transformation Az
    AOmegaBarWtA = A %*% cov_zEps
    # get cross-covariance matrix between Z and eps
    cov_zEps = -cov_zEps
    # variance-covariance matrix
    Sigma_eps = diag(1, nrow = m) + AOmegaBarWtA
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
        diag(W, nrow = d) - regrMean_zGNeps %*% t(cov_zEps) +
          diag(sqrt(.Machine$double.eps), nrow = d)
      )),
      error = function(err) { # if chol fails then use svd
        #svdObj = svd(
        # OmegaBar - regrMean_zGNeps %*% t(cov_zEps) +
        #  diag(sqrt(.Machine$double.eps), d),
        # nv = 0
        #)
        #return(t(sqrt(svdObj$d) * t(svdObj$u)))
        return(backsolve(
          r = chol(diag(1 / W, nrow = d) + tAA +
                     diag(sqrt(.Machine$double.eps), d)),
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

      ##### update W

      ### update psi
      # compute weights of discrete full conditional of psi
      weiPsi = logConstFullPsi + (1 - d * seq_psi) * log(tau) +
        (seq_psi - 1) * sum(log(zeta))
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
          n = 1, lambda = (psi - 1) * d, chi = 2 * sum(abs(Z) / zeta), psi = 1
        ),
        error = function(err) {
          GIGrvg::rgig(
            n = 1, lambda = (psi - 1) * d,
            chi = 2 * sum(abs(Z) / zeta) + sqrt(.Machine$double.eps),
            psi = 1
          )
        }
      )

      ### update eta
      # inverse of eta
      eta = statmod::rinvgauss(n = d, mean = tau * zeta / abs(Z), shape = 1) +
        sqrt(.Machine$double.eps)
      # get eta
      eta = 1 / eta

      ### get W
      W = eta * (zeta * tau)^2



      ##### update T and Z

      ### sample eps
      # minus cross-covariance matrix between Z and eps
      cov_zEps = W * t(A)
      # variance-covariance matrix of linear transformation Az
      AOmegaBarWtA = A %*% cov_zEps
      # get cross-covariance matrix between Z and eps
      cov_zEps = -cov_zEps
      # variance-covariance matrix
      Sigma_eps = diag(1, nrow = m) + AOmegaBarWtA
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
          diag(W, nrow = d) - regrMean_zGNeps %*% t(cov_zEps) +
            diag(sqrt(.Machine$double.eps), nrow = d)
        )),
        error = function(err) { # if chol fails then use svd
          #svdObj = svd(
          # OmegaBar - regrMean_zGNeps %*% t(cov_zEps) +
          #  diag(sqrt(.Machine$double.eps), d),
          # nv = 0
          #)
          #return(t(sqrt(svdObj$d) * t(svdObj$u)))
          return(backsolve(
            r = chol(diag(1 / W, nrow = d) + tAA +
                       diag(sqrt(.Machine$double.eps), d)),
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
    if(insim %% verbose == 0) {
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
    print(paste("pSUN_ProbitDirLapl_Gibbs: end time at ", Sys.time(), sep = ""))
  }

  ### return result
  return(sample)

}