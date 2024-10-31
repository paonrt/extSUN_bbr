###############################################################################
#
#
# Simulation Study: Probit Model and Sparse Parameters
#
#
###############################################################################





############################### Initialization ################################

### clear IDE
rm(list = ls()); gc(); cat("\14")

### set sequential computation for all workers
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(BLIS_NUM_THREADS = 1)
Sys.setenv(FFTW_NUM_THREADS = 1)

### set the seed
set.seed(1994)

### set number of iterations
niter = 2400

### set number of simulations
nsim = 5e+4

### set number of observations
n = 50

### set number of parameters
p = 500

### set prior parameters
mu = rep(0, p)
Sigma = diag(16, nrow = p)

### set true beta vector
# fraction of non zero parameters
not0 = 0.2
# set non zero parameters
trueBeta = c(
  rep(qnorm(0.05), p * not0 * 0.1),
  rep(qnorm(0.1), p * not0 * 0.1),
  rep(qnorm(0.2), p * not0 * 0.1),
  rep(qnorm(0.3), p * not0 * 0.1),
  rep(qnorm(0.4), p * not0 * 0.1),
  rep(qnorm(0.6), p * not0 * 0.1),
  rep(qnorm(0.7), p * not0 * 0.1),
  rep(qnorm(0.8), p * not0 * 0.1),
  rep(qnorm(0.9), p * not0 * 0.1),
  rep(qnorm(0.95), p * not0 * 0.1)
)
# set the others to zero
trueBeta = c(
  rep(0, p - length(trueBeta)), trueBeta
)





############### Sim Study: sample size of 50 vs 500 parameters ################

##### parallel environment setup

### get cluster
cluster = parallel::makeCluster(24)

### get true number of iterations according to cluster size
niter = ceiling(niter / length(cluster)) * length(cluster)

### export variables
parallel::clusterExport(cl = cluster, varlist = ls(all.names = TRUE))

### identify cores
for (id_cluster in 1:length(cluster)) {
  parallel::clusterExport(cl = cluster[id_cluster], "id_cluster")
}

### open report file
parallel::clusterEvalQ(
  cl = cluster, sink(paste("report", id_cluster, ".txt", sep = ""))
)

### set seed
parallel::clusterSetRNGStream(cl = cluster)




########## do iterations in parallel
res = parallel::clusterApply(
  cl = cluster, x = matrix(1:niter), fun = function(xapp) {

    ### print start of iteration
    print(paste(
      "iteration ", xapp, " start at time ", Sys.time()
    ))

    ### sample the design matrix (except the intercept)
    # draw from independent standard Gaussian distributions
    X = matrix(
      c(rep(1, n), rnorm(n * (p - 1), sd = 0.5)),
      nrow = n, ncol = p
    )
    # scale
    X[, -1] = 0.5 * scale(X[, -1])

    ### sample response variable
    y = rbinom(n, 1, plogis(X %*% trueBeta))

    ### compute BX
    BX = (2 * y - 1) * X



    ##### pSUN Gibbs sampler

    ### start time
    time = Sys.time()

    ### run
    output = pSUN::pSUN_ProbitBessit_Gibbs(
      nsim, A = t(sqrt(diag(Sigma)) * t(BX)), b = BX %*% mu, nu = 1,
      xi = mu, Omega = Sigma
    )

    ### end time
    time = difftime(Sys.time(), time, units = "secs")

    ### autocorrelation function
    acf = apply(matrix(1:p), 1, FUN = function(yapp) {
      acf(output[yapp, ], lag.max = 250, plot = FALSE)$acf
    })

    ### effective sample size
    ess = apply(matrix(1:p), 1, FUN = function(yapp) {
      sns::ess(output[yapp, ])
    })

    ### compute difference between estimated posterior means and true values
    err = rowMeans(output) - trueBeta



    ### get results
    return(list(
      time = time,
      acf = acf,
      ess = ess,
      err = err
    ))

  }
)




##### parallel environment strip

### close report file
parallel::clusterEvalQ(cl = cluster, sink())

### leave cluster
parallel::stopCluster(cluster)





################################### Closure ###################################
save(list = ls(), file = "SimStudyProbitSparse_50n500pLaplacit.RData")