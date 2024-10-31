#' @title The Multivariate Non-Central Student-t Distribution with Identity
#'  Matrix

#' @description
#' Distribution function for the multivariate non-central Student-t distribution
#'  with identity matrix as scale matrix.

#' @usage
#' pmvnctID = function(x, nu = NULL, delta, n = NULL, dx = NULL, log.p = FALSE,
#'                     maxComp = 1e+9, stringLength = 100)

#' @details
#' It is said \eqn{X} follows a multivariate non-central Student-t distribution
#'  with shape \eqn{\nu}, scale matrix \eqn{\Omega} and non centrality vector
#'  \eqn{\delta} if
#'  \deqn{
#'    X = \sqrt{W} (Z + \delta) \, ,
#'  }
#'  where \eqn{Z \sim N(0, \Omega) \perp\!\!\!\perp W \sim \mathrm{inv.gamma}
#'  (\nu / 2, \nu / 2)}. Notice that this function requires \eqn{\Omega = I}.
#'
#' The distribution function is computed via a quasi Monte Carlo algorithm as
#' described in Onorati and Liseo (2022).

#' @references
#' Onorati, P., & Liseo, B. (2022). An extension of the unified skew-normal
#'  family of distributions and application to Bayesian binary regression.
#'  arXiv preprint arXiv:2209.03474.

#' @param x Matrix of arguments.
#' @param nu Shape parameter.
#' @param delta Non centrality vector.
#' @param n Dimension of multivariate non-central Student-t.
#' @param dx Vector of pseudo-sampled points.
#' @param log.p If it is true then the logarythm is returned.
#' @param maxComp Maximum number of rows of x to be computed for each iteration.
#' @param stringLength Number of pseudo-sapled points to be computed. It is used
#'  only if the points are not provided.
#' @return Vector of distribution function values.

#' @import matrixStats
#' @export
pmvnctID = function(
  x, nu = NULL, delta, n = ncol(x), dx = NULL, log.p = FALSE,
  maxComp = 1e+9, stringLength = 100
) {

  ### get points for quasi Monte Carlo ( they are 1 / sqrt(quantileInvGamma) )
  if (is.null(dx)) {
    if (is.null(nu)) stop("dx is NULL, nu cannot be NULL as well")
    dx = sqrt(qgamma(
      1 - c(1:stringLength) / (stringLength + 1),
      shape = 0.5 * nu, rate = 0.5 * nu
    ))
  }

  ### compute number of iterations
  nIter = ceiling(n / maxComp)



  ##### do iterations

  ### initialize results
  log_nct_CDF = double(length = n)

  for (inIter in 1:nIter) {

    ### bounds
    low = (inIter - 1) * maxComp + 1
    upp = min(inIter * maxComp, n)

    ### compute logarithm of Gaussian CDF
    log_gauss_CDF = pnorm(x[, low:upp] %o% dx + c(delta), log.p = TRUE)

    ### compute sums of Gaussian log-CDF
    sumsLog_gauss_CDF = matrix(colSums(log_gauss_CDF), nrow = upp - low + 1)

    ### get max of sumsLog_gauss_CDF
    maxSumsLog_gauss_CDF = matrixStats::rowMaxs(sumsLog_gauss_CDF)

    ### compute logarithm of non-central t CDF
    log_nct_CDF[low:upp] = maxSumsLog_gauss_CDF -
      log(length(dx)) + log(rowSums(
        exp(sumsLog_gauss_CDF-maxSumsLog_gauss_CDF)
      ))

    ### clear memory
    invisible(gc())

  }



  ### get exponential of log-CDF if required
  if (!log.p) log_nct_CDF = exp(log_nct_CDF)

  ### return result
  return(log_nct_CDF)

}