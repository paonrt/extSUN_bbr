# Multivariate Truncated Normal: Get Parameters for Exponential Tilting
#  (from TruncatedNormal)

# This function comes from mvrandn of TruncatedNormal package. It searches
# the saddlepoint, and the output should be used only as input of
# function TruncatedNormal_mvrandn_runExpTilt.
#
# It is NOT exported in the namespace.

############### function
# no help and no export in the namespace
#' @import TruncatedNormal
#' @import nleqslv
TruncatedNormal_mvrandn_getExpTilt <- function(
  l, u, Sig, mu = NULL, x0 = NULL
) {

  # parameters are the same of mvrandn of TruncatedNormal package, see it for
  #  help

  ##### reload namespace
  #evalCpp = Rcpp::evalCpp
  #auglag = alabama::auglag
  #nleqslv = nleqslv::nleqslv
  #randtoolbox::sobol
  #dnorm = stats::dnorm
  #pnorm = stats::pnorm
  #qnorm = stats::qnorm
  #rnorm = stats::rnorm
  #pt = stats::pt
  #qt = stats::qt
  #rexp = stats::rexp
  #gradpsi
  #jacpsi
  #psy
  #cholperm
  #trandn
  #mvnrnd

  ex_l = l
  ex_u = u
  ex_Sig = Sig
  ex_mu = mu



  d <- length(l); # basic input check
  if (length(u)!= d | d != sqrt(length(Sig)) | any(l > u)) {
    stop("l, u, and Sig have to match in dimension with u>l")
  }
  if (!is.null(mu)) {
    l <- l - mu
    u <- u - mu
  }
  #Univariate case
  if (d == 1) {
    stop("Use another function for univariate case")
  }
  # Cholesky decomposition of matrix
  out <- TruncatedNormal::cholperm(
    as.matrix(Sig), as.numeric(l), as.numeric(u)
  );
  Lfull <- out$L;
  l <- out$l;
  u <- out$u;
  D <- diag(Lfull);
  perm <- out$perm;
  #if (any(D < 1e-10)) {
  #  warning("Method may fail as covariance matrix is singular!")
  #}
  L <- Lfull / D;
  u <- u / D;
  l <- l / D; # rescale
  L <- L - diag(d); # remove diagonal
  # find optimal tilting parameter via non-linear equation solver
  if (is.null(x0)) x0 <- rep(0, 2 * length(l) - 2)
  solvneq <- nleqslv::nleqslv(
    x0, fn = TruncatedNormal:::gradpsi, jac = TruncatedNormal:::jacpsi,
    L = L, l = l, u = u, global = "pwldog", method = "Broyden",
    control = list(maxit = 500L, allowSingular = TRUE)
  )
  xmu <- solvneq$x
  exitflag <- solvneq$termcd
  #if (
  #   !(exitflag %in% c(1,2)) ||
  #     !isTRUE(all.equal(solvneq$fvec, rep(0, length(x0)), tolerance = 1e-6))
  #   ) {
  #     warning('Did not find a solution to the nonlinear system in `mvrandn`!')
  #   }
  #xmu <- nleq(l,u,L) # nonlinear equation solver (not used, 40 times slower!)
  x <- xmu[1:(d - 1)];
  muV <- xmu[d:(2 * d - 2)]; # assign saddlepoint x* and mu*
  # compute psi star
  psistar <- TruncatedNormal:::psy(x, L, l, u, muV);

  ### initialize result list
  res = list()

  ### get input
  res[["input"]] = list(
    l = ex_l, u = ex_u, Sig = ex_Sig, mu = ex_mu
  )

  ### get output
  res[["output"]] = list(
    d = d, L = L, l = l, u = u, x = x, muV = muV,
    psistar = psistar, perm = perm, Lfull = Lfull
  )

  return(res)
}