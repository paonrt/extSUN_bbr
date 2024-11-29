###############################################################################
#
#
# Simulation Study: Logit Model and Unbalanced Data
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

### set number of simulations
nsim = 1e+5

### set sample size
n = 200

### set number of parameters
p = 1

### set prior parameters
mu = rep(0, p)
Sigma = diag(52.6379, nrow = p)

### set response variables
y = c(1, rep(0, n - 1))

### set design matrix
X = matrix(rep(1, n), ncol = 1)

### likelihood matrix
BX = (2 * y - 1) * X

### posterior parameters
A = t(sqrt(diag(Sigma)) * t(BX))
b = BX %*% mu

### verbose parameter
verbose = ceiling(nsim / 10)

### prior and algorithm cases
cases = c(
  "pSUN", "PG", "UPG"
)





################ Sim Study: sample size of 200 vs 1 parameter #################

##### set parallel environment

### get cluster
cluster = parallel::makeCluster(3)

### export variables
parallel::clusterExport(cl = cluster, varlist = ls())

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
  cl = cluster, x = matrix(cases), fun = function(xapp) {

    ### Polya-Gamma Gibbs sampler
    if (xapp == "PG") {
      # start time
      time = Sys.time()
      # run
      resWorker = pSUN::logitGaussPG(
        nsim, X = X, y = y, mu = mu, Sigma = Sigma, verbose = verbose
      )
      # end time
      time = difftime(Sys.time(), time, units = "secs")
    }

    ### pSUN Gibbs sampler
    if (xapp == "pSUN") {
      # start time
      time = Sys.time()
      # run
      resWorker = pSUN::pSUN_LogitGauss_Gibbs(
        nsim, A = A, b = b, xi = mu, Omega = Sigma, verbose = verbose
      )
      # end time
      time = difftime(Sys.time(), time, units = "secs")
    }

    ### Ultimate Polya-Gamma Gibbs sampler
    if (xapp == "UPG") {
      # start time
      time = Sys.time()
      # run
      resWorker = pSUN::logitGaussUPG(
        nsim, X = X, y = y, Sigma = Sigma, verbose = verbose
      )
      # end time
      time = difftime(Sys.time(), time, units = "secs")
    }

    ### get results
    # acf
    acf = apply(matrix(1:p), 1, FUN = function(yapp) {
      acf(resWorker[yapp, ], lag.max = 150, plot = FALSE)$acf
    })
    # effective sample size
    ess = apply(
      matrix(1:p), 1, FUN = function(yapp) sns::ess(resWorker[yapp, ])
    )
    # take sample
    sam = resWorker
    # return
    return(list(
      time = time,
      acf = acf,
      ess = ess,
      sam = sam
    ))

  }
)




##### parallel environment strip

### close report file
parallel::clusterEvalQ(cl = cluster, sink())

### leave cluster
parallel::stopCluster(cluster)





################################### Closure ###################################
save(list = ls(), file = "SimStudyLogitUnbalanced_200n1pGauss.RData")