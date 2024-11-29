###############################################################################
#
#
# Simulation Study: Polya-Gamma vs pSUN Methods
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

### set minimum number of iterations
niter = 2400

### set number of simulations
nsim = 1e+4

### set number of observations
n = 50

### set number of parameters
p = 500

### set prior parameters
mu = rep(0, p)
Sigma = diag(52.6379, nrow = p)

### set true beta vector
trueBeta = c(
  0,
  rep(qlogis(0.05), p * 0.05),
  rep(qlogis(0.1), p * 0.1),
  rep(qlogis(0.2), p * 0.1),
  rep(qlogis(0.3), p * 0.1),
  rep(qlogis(0.4), p * 0.1),
  rep(qlogis(0.6), p * 0.1),
  rep(qlogis(0.7), p * 0.1),
  rep(qlogis(0.8), p * 0.1),
  rep(qlogis(0.9), p * 0.1),
  rep(qlogis(0.95), p * 0.05),
  rep(0, p * 0.1 - 1)
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

    ### stopping variable for while loop
    toDo = TRUE

    ### run until get all successes
    while (toDo) {

      ### reset stopping variable
      toDo = FALSE

      ### sample the design matrix (except the intercept)
      # draw from independent standard Gaussian distributions
      X = matrix(
        c(rep(1, n), rnorm(n * (p - 1), sd = 0.5)), nrow = n, ncol = p
      )
      # scale
      X[, -1] = 0.5 * scale(X[, -1])

      ### sample response variable
      y = rbinom(n, 1, plogis(X %*% trueBeta))

      ### compute BX
      BX = (2 * y - 1) * X



      ##### importance sampling (pSUN)

      ### start time
      IS_time = Sys.time()

      ### run
      IS_output = tryCatch(
        expr = pSUN::pSUN_LogitGauss_IS(
          nsim, A = t(sqrt(diag(Sigma)) * t(BX)), b = BX %*% mu,
          xi = mu, Omega = Sigma
        ),
        error = function(err) return(FALSE)
      )

      ### end time
      IS_time = difftime(Sys.time(), IS_time, units = "secs")

      ### check success
      if (isFALSE(IS_output[1])) {

        ### redo iteration with new data if there is a fail
        toDo = TRUE

        ### print the error
        print("IS failed, redo the iteration with new data")

      } else {

        ### effective sample size
        IS_ess = rep(1 / sum(IS_output$normWeights^2), p)

        ### compute difference between estimated posterior means and true values
        IS_err = colSums(t(IS_output$sample) * IS_output$normWeights) - trueBeta

      }



      ##### pSUN

      ### start time
      pSUN_time = Sys.time()

      ### run
      pSUN_output = tryCatch(
        expr = pSUN::pSUN_LogitGauss_Gibbs(
          nsim, A = t(sqrt(diag(Sigma)) * t(BX)), b = BX %*% mu,
          xi = mu, Omega = Sigma
        ),
        error = function(err) return(FALSE)
      )

      ### end time
      pSUN_time = difftime(Sys.time(), pSUN_time, units = "secs")

      ### check success
      if (isFALSE(pSUN_output[1])) {

        ### redo iteration with new data if there is a fail
        toDo = TRUE

        ### print the error
        print("pSUN failed, redo the iteration with new data")

      } else {

        ### autocorrelation function
        pSUN_acf = apply(matrix(1:p), 1, FUN = function(yapp) {
          acf(pSUN_output[yapp, ], lag.max = 250, plot = FALSE)$acf
        })

        ### effective sample size
        pSUN_ess = apply(matrix(1:p), 1, FUN = function(yapp) {
          sns::ess(pSUN_output[yapp, ])
        })

        ### compute difference between estimated posterior means and true values
        pSUN_err = rowMeans(pSUN_output) - trueBeta

      }



      ##### Polya-Gamma

      ### start time
      PG_time = Sys.time()

      ### run
      PG_output = tryCatch(
        expr = pSUN::logitGaussPG(
          nsim, X = X, y = y, mu = mu, Sigma = Sigma
        ),
        error = function(err) return(FALSE)
      )

      ### end time
      PG_time = difftime(Sys.time(), PG_time, units = "secs")

      ### check success
      if (isFALSE(PG_output[1])) {

        ### redo iteration with new data if there is a fail
        toDo = TRUE

        ### print the error
        print("PG failed, redo the iteration with new data")

      } else {

        ### autocorrelation function
        PG_acf = apply(matrix(1:p), 1, FUN = function(yapp) {
           acf(PG_output[yapp, ], lag.max = 250, plot = FALSE)$acf
        })

        ### effective sample size
        PG_ess = apply(matrix(1:p), 1, FUN = function(yapp) {
          sns::ess(PG_output[yapp, ])
        })

        ### compute difference between estimated posterior means and true values
        PG_err = rowMeans(PG_output) - trueBeta

      }



      ##### Ultimate Polya-Gamma

      ### start time
      UPG_time = Sys.time()

      ### run
      UPG_output = tryCatch(
        expr = pSUN::logitGaussUPG(
          nsim, X = X, y = y, Sigma = Sigma
        ),
        error = function(err) return(FALSE)
      )

      ### end time
      UPG_time = difftime(Sys.time(), UPG_time, units = "secs")

      ### check success
      if (isFALSE(UPG_output[1])) {

        ### redo iteration with new data if there is a fail
        toDo = TRUE

        ### print the error
        print("UPG failed, redo the iteration with new data")

      } else {

        ### autocorrelation function
        UPG_acf = apply(matrix(1:p), 1, FUN = function(yapp) {
          acf(UPG_output[yapp, ], lag.max = 250, plot = FALSE)$acf
        })

        ### effective sample size
        UPG_ess = apply(matrix(1:p), 1, FUN = function(yapp) {
          sns::ess(UPG_output[yapp, ])
        })

        ### compute difference between estimated posterior means and true values
        UPG_err = rowMeans(UPG_output) - trueBeta

      }

    }



    ### get results
    return(list(
      # importance sampling (pSUN)
      IS_time = IS_time,
      IS_ess = IS_ess,
      IS_err = IS_err,
      # pSUN
      pSUN_time = pSUN_time,
      pSUN_acf = pSUN_acf,
      pSUN_ess = pSUN_ess,
      pSUN_err = pSUN_err,
      # Polya-Gamma
      PG_time = PG_time,
      PG_acf = PG_acf,
      PG_ess = PG_ess,
      PG_err = PG_err,
      # Ultimate Polya-Gamma
      UPG_time = UPG_time,
      UPG_acf = UPG_acf,
      UPG_ess = UPG_ess,
      UPG_err = UPG_err
    ))

  }
)




##### parallel environment strip

### close report file
parallel::clusterEvalQ(cl = cluster, sink())

### leave cluster
parallel::stopCluster(cluster)





################################### Closure ###################################
save(list = ls(), file = "SimStudyPGvs_pSUN_50n500pGauss.RData")