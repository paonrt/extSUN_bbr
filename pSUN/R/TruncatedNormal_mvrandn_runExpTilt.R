# Multivariate Truncated Normal: RNG and log-normalizing constant computation
#  via Exponential Tilting (from TruncatedNormal)

# This function comes from mvrandn of TruncatedNormal package. It samples from a
# multivariate truncated Gaussian random variables via acceptance-rejection
# algorithm and it computes the logarithm of normalizing constant via
# importance sampling. In both cases the instrumental density is given by
# exponential tilting obtained by function TruncatedNormal_mvrandn_getExpTilt.
#
# It is NOT exported in the namespace.

############### function
# no help and no export in the namespace
TruncatedNormal_mvrandn_runExpTilt <- function(mvrandnObj, n, plim = 1 / n) {

  # mvrandnObj is the result object from TruncatedNormal_mvrandn_getPar;
  # n is the sample size to be obtain from AR algorithm;
  # plim is the minimum empirical acceptance probability such that the function
  #  does not stop.

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

  mu = mvrandnObj$input$mu
  d = mvrandnObj$output$d
  L = mvrandnObj$output$L
  l = mvrandnObj$output$l
  u = mvrandnObj$output$u
  muV = mvrandnObj$output$muV
  psistar = mvrandnObj$output$psistar
  perm = mvrandnObj$output$perm
  Lfull = mvrandnObj$output$Lfull
  # minimum initial trials
  min_nsim = max(1e+2, n)




  # start acceptance rejection sampling
  Z <- matrix(0, nrow = d, ncol = n)
  accept <- 0L; iter <- 0L; nsim <- min(c(min_nsim, 1e+4)); ntotsim <- 0L;
  logNormConst = -Inf; reachMax = FALSE; reachAcc = FALSE; smallProbFlag = TRUE
  while (accept < n) { # while # of accepted is less than n
    call <- TruncatedNormal:::mvnrnd(n = nsim, L = L, l = l, u = u, mu = muV);
    lLR = call$logpr - psistar
    if (is.infinite(logNormConst)) {
      logNormConst = psistar + log(sum(exp(lLR)))
    } else {
      logNormConst = logNormConst + log1p(exp(
        psistar + log(sum(exp(lLR))) - logNormConst
      ))
    }
    ntotsim <- ntotsim + nsim
    idx <- stats::rexp(nsim) > -lLR; # acceptance tests
    m <- sum(idx)
    if (m > n - accept) {
      m <- n - accept
      idx <- which(idx)[1:m]
    }
    if(m > 0){
      Z[,(accept+1):(accept+m)] <- call$Z[,idx];  # accumulate accepted
    }
    accept <- accept + m; # keep track of accepted
    iter <- iter + 1L;  # keep track of while loop iterations
    nsim <- min(c(
      1e+4, ceiling((n - accept) * ntotsim / accept)
    ))
    # if iterations are getting large, give warning
    if (
      (ntotsim > max(2e+4, 1 / plim)) &
        (accept / ntotsim < plim) &
        smallProbFlag
    ) {
      #warning("Acceptance probability smaller than 0.001")
      #print("Acceptance probability smaller than 0.001")
      stop(paste( "Acceptance probability smaller than ", plim, sep = "" ))
      smallProbFlag = FALSE
    } else if (iter > 1e+5) { # if iterations too large, seek approximation
      reachAcc = accept > 0
      accept = n + 1
      reachMax = TRUE
      if (accept == 0) {
        warning("Could not sample from truncated Normal - check input")
      } else if(accept > 1) {
        Z <- Z[,1:accept]
        warning("Sample of size smaller than n returned.")
      }
    }
  }

  logNormConst = logNormConst - log(ntotsim)

  if (!reachMax | reachAcc) {
    ## finish sampling; postprocessing
    out = sort(perm, decreasing = FALSE, index.return = TRUE)
    order = out$ix
    Z = Lfull %*% Z # reverse scaling of L
    Z = Z[order, ] # reverse the Cholesky permutation

    if (!is.null(mu)) {
      return(list(
        sample = Z + mu,
        logNormConst = logNormConst,
        ntotsim = ntotsim
      ))
    } else {
      return(list(
        sample = Z,
        logNormConst = logNormConst,
        ntotsim = ntotsim
      ))
    }

  } else {
    return(list(
      sample = NA,
      logNormConst = logNormConst,
      ntotsim = ntotsim
    ))
  }
}
