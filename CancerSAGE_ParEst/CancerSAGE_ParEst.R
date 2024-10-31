###############################################################################
#
#
# Cancer SAGE Dataset: Parameters Estimation
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



##### load and pre-process data

### load raw-data from pSUN package
cancer_sage = pSUN::cancer_sage

### get response variables
y = cancer_sage[, 1]

### get standardized design matrix (with intercept)
X = cbind(
  rep(1, length(y)),
  apply(cancer_sage[, -1], 2, arm::rescale)
)



### number of scans in MCMC and IS algorithms
nsim = 1e+5

### likelihood matrix
B = 2 * diag(y) - diag(1, nrow = length(y))




################## Real Data Analysis of Cancer SAGE Dataset ##################

##### set parallel environment

### get cluster
cluster = parallel::makeCluster(17)

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




########## do computation with different links and priors in parallel
parallel::clusterEvalQ(cl = cluster, expr = {

  ##### probit model with Gaussian prior
  if (id_cluster == 1) {

    ### get prior parameters
    # scale matrix
    Sigma = diag(16, nrow = ncol(X))
    # location vector
    mu = rep(0, ncol(X))

    ### compute posterior parameters
    # slant matrix
    A = B %*% X %*% diag(sqrt(diag(Sigma)), nrow = ncol(X))
    # slant vector
    b = B %*% X %*% mu
    # location vector
    xi = mu
    # scale matrix
    Omega = Sigma

    ### run posterior simulation algorithm
    t0 = Sys.time()
    res = pSUN::pSUN_ProbitGauss_RNG(nsim, A = A, b = b, xi = xi, Omega = Omega)
    tn = Sys.time()
    print(difftime(tn, t0, units = "sec"))

    ### save results
    # list
    probitGauss = list(res = res, time = difftime(tn, t0, units = "sec"))
    # save
    save(list = c("probitGauss"), file = "probitGauss.RData")

  }



  ##### probit model with Cauchy prior
  if (id_cluster == 2) {

    ### get prior parameters
    # scale matrix
    Sigma = diag(0.3807, nrow = ncol(X))
    # location vector
    mu = rep(0, ncol(X))

    ### compute posterior parameters
    # slant matrix
    A = B %*% X %*% diag(sqrt(diag(Sigma)), nrow = ncol(X))
    # slant vector
    b = B %*% X %*% mu
    # location vector
    xi = mu
    # scale matrix
    Omega = Sigma

    ### run posterior simulation algorithm
    t0 = Sys.time()
    res = pSUN::pSUN_ProbitStudent_Gibbs(
      nsim, A = A, b = b, nu = 1, xi = xi, Omega = Omega, verbose = 1000
    )
    tn = Sys.time()
    print( difftime(tn, t0, units = "sec") )

    ### save results
    # list
    probitCauchy = list(res = res, time = difftime(tn, t0, units = "sec"))
    # save
    save(list = c("probitCauchy"), file = "probitCauchy.RData")

  }



  ##### probit model with independent Laplace priors
  if (id_cluster == 3) {

    ### get prior parameters
    # scale matrix
    Sigma = diag(6.8487, nrow = ncol(X))
    # location vector
    mu = rep(0, ncol(X))

    ### compute posterior parameters
    # slant matrix
    A = B %*% X %*% diag(sqrt(diag(Sigma)), nrow = ncol(X))
    # slant vector
    b = B %*% X %*% mu
    # location vector
    xi = mu
    # scale matrix
    Omega = Sigma

    ### run posterior simulation algorithm
    t0 = Sys.time()
    res = pSUN::pSUN_ProbitBessit_Gibbs(
      nsim, A = A, b = b, nu = 1, xi = xi, Omega = Omega, verbose = 1000
    )
    tn = Sys.time()
    print(difftime(tn, t0, units = "sec"))

    ### save results
    # list
    probitLaplacit = list(res = res, time = difftime(tn, t0, units = "sec"))
    # save
    save(list = c("probitLaplacit"), file = "probitLaplacit.RData")

  }



  ##### probit model with Dirichlet-Laplace prior
  if (id_cluster == 4) {

    ### get prior parameters
    # scale matrix
    Sigma = diag(2.7655, nrow = ncol(X))
    # location vector
    mu = rep(0, ncol(X))

    ### compute posterior parameters
    # slant matrix
    A = B %*% X %*% diag(sqrt(diag(Sigma)), nrow = ncol(X))
    # slant vector
    b = B %*% X %*% mu
    # location vector
    xi = mu
    # scale matrix
    Omega = Sigma

    ### run posterior simulation algorithm
    t0 = Sys.time()
    res = pSUN::pSUN_ProbitDirLapl_Gibbs(
      nsim, A = A, b = b, xi = xi, Omega = Omega, verbose = 1000
    )
    tn = Sys.time()
    print(difftime(tn, t0, units = "sec"))

    ### save results
    # list
    probitDirLapl = list(res = res, time = difftime(tn, t0, units = "sec"))
    # save
    save(list = c("probitDirLapl"), file = "probitDirLapl.RData")

  }



  ##### logit model with Gaussian prior and importance sampling
  if (id_cluster == 5) {

    ### get prior parameters
    # scale matrix
    Sigma = diag(52.6379, nrow = ncol(X))
    # location vector
    mu = rep(0, ncol(X))

    ### compute posterior parameters
    # slant matrix
    A = B %*% X %*% diag(sqrt(diag(Sigma)), nrow = ncol(X))
    # slant vector
    b = B %*% X %*% mu
    # location vector
    xi = mu
    # scale matrix
    Omega = Sigma

    ### run posterior simulation algorithm
    t0 = Sys.time()
    res = pSUN::pSUN_LogitGauss_IS(nsim, A = A, b = b, xi = xi, Omega = Omega)
    tn = Sys.time()
    print(difftime(tn, t0, units = "sec"))

    ### save results
    # list
    logitGaussIS = list(res = res, time = difftime(tn, t0, units = "sec"))
    # save
    save(list = c("logitGaussIS"), file = "logitGaussIS.RData")

  }



  ##### logit model with Gaussian prior and pSUN Gibbs sampler
  if (id_cluster == 6) {

    ### get prior parameters
    # scale matrix
    Sigma = diag(52.6379, nrow = ncol(X))
    # location vector
    mu = rep(0, ncol(X))

    ### compute posterior parameters
    # slant matrix
    A = B %*% X %*% diag(sqrt(diag(Sigma)), nrow = ncol(X))
    # slant vector
    b = B %*% X %*% mu
    # location vector
    xi = mu
    # scale matrix
    Omega = Sigma

    ### run posterior simulation algorithm
    t0 = Sys.time()
    res = pSUN::pSUN_LogitGauss_Gibbs(
      nsim, A = A, b = b, xi = xi, Omega = Omega, verbose = 1000
    )
    tn = Sys.time()
    print(difftime(tn, t0, units = "sec"))

    ### save results
    # list
    logitGauss_pSUN = list(res = res, time = difftime(tn, t0, units = "sec"))
    # save
    save(list = c("logitGauss_pSUN"), file = "logitGauss_pSUN.RData")

  }



  ##### logit model with Cauchy prior and pSUN Gibbs sampler
  if (id_cluster == 7) {

    ### get prior parameters
    # scale matrix
    Sigma = diag(1.2525, nrow = ncol(X))
    # location vector
    mu = rep(0, ncol(X))

    ### compute posterior parameters
    # slant matrix
    A = B %*% X %*% diag(sqrt(diag(Sigma)), nrow = ncol(X))
    # slant vector
    b = B %*% X %*% mu
    # location vector
    xi = mu
    # scale matrix
    Omega = Sigma

    ### run posterior simulation algorithm
    t0 = Sys.time()
    res = pSUN::pSUN_LogitStudent_Gibbs(
      nsim, A = A, b = b, nu = 1, xi = xi, Omega = Omega, verbose = 1000
    )
    tn = Sys.time()
    print(difftime(tn, t0, units = "sec"))

    ### save results
    # list
    logitCauchy_pSUN = list(res = res, time = difftime(tn, t0, units = "sec"))
    # save
    save(list = c("logitCauchy_pSUN"), file = "logitCauchy_pSUN.RData")

  }



  ##### logit model with independent Laplace priors and pSUN Gibbs sampler
  if (id_cluster == 8) {

    ### get prior parameters
    # scale matrix
    Sigma = diag(22.5314, nrow = ncol(X))
    # location vector
    mu = rep(0, ncol(X))

    ### compute posterior parameters
    # slant matrix
    A = B %*% X %*% diag(sqrt(diag(Sigma)), nrow = ncol(X))
    # slant vector
    b = B %*% X %*% mu
    # location vector
    xi = mu
    # scale matrix
    Omega = Sigma

    ### run posterior simulation algorithm
    t0 = Sys.time()
    res = pSUN::pSUN_LogitBessit_Gibbs(
      nsim, A = A, b = b, nu = 1, xi = xi, Omega = Omega, verbose = 1000
    )
    tn = Sys.time()
    print(difftime(tn, t0, units = "sec"))

    ### save results
    # list
    logitLaplacit_pSUN = list(res = res, time = difftime(tn, t0, units = "sec"))
    # save
    save(list = c("logitLaplacit_pSUN"), file = "logitLaplacit_pSUN.RData")

  }



  ##### logit model with Dirichlet-Laplace prior and pSUN Gibbs sampler
  if (id_cluster == 9) {

    ### get prior parameters
    # scale matrix
    Sigma = diag(9.0984, nrow = ncol(X))
    # location vector
    mu = rep(0, ncol(X))

    ### compute posterior parameters
    # slant matrix
    A = B %*% X %*% diag(sqrt(diag(Sigma)), nrow = ncol(X))
    # slant vector
    b = B %*% X %*% mu
    # location vector
    xi = mu
    # scale matrix
    Omega = Sigma

    ### run posterior simulation algorithm
    t0 = Sys.time()
    res = pSUN::pSUN_LogitDirLapl_Gibbs(
      nsim, A = A, b = b, xi = xi, Omega = Omega, verbose = 1000
    )
    tn = Sys.time()
    print(difftime(tn, t0, units = "sec"))

    ### save results
    # list
    logitDirLapl_pSUN = list(res = res, time = difftime(tn, t0, units = "sec"))
    # save
    save(list = c("logitDirLapl_pSUN"), file = "logitDirLapl_pSUN.RData")

  }



  ##### logit model with Gaussian prior and Polya-Gamma Gibbs sampler
  if (id_cluster == 10) {

    ### get prior parameters
    # scale matrix
    Sigma = diag(52.6379, nrow = ncol(X))
    # location vector
    mu = rep(0, ncol(X))

    ### run posterior simulation algorithm
    t0 = Sys.time()
    res = pSUN::logitGaussPG(
      nsim, X = X, y = y, mu = mu, Sigma = Sigma, verbose = 1000
    )
    tn = Sys.time()
    print(difftime(tn, t0, units = "sec"))

    ### save results
    # list
    logitGaussPG = list(res = res, time = difftime(tn, t0, units = "sec"))
    # save
    save(list = c("logitGaussPG"), file = "logitGaussPG.RData")

  }



  ##### logit model with Cauchy prior and Polya-Gamma Gibbs sampler
  if (id_cluster == 11) {

    ### get prior parameters
    # scale matrix
    Sigma = diag(1.2525, nrow = ncol(X))
    # location vector
    mu = rep(0, ncol(X))

    ### run posterior simulation algorithm
    t0 = Sys.time()
    res = pSUN::logitStudentPG(
      nsim, X = X, y = y, nu = 1, mu = mu, Sigma = Sigma, verbose = 1000
    )
    tn = Sys.time()
    print(difftime(tn, t0, units = "sec"))

    ### save results
    # list
    logitCauchyPG = list(res = res, time = difftime(tn, t0, units = "sec"))
    # save
    save(list = c("logitCauchyPG"), file = "logitCauchyPG.RData")
  }



  ##### logit model with independent Laplace priors and Polya-Gamma Gibbs
  #      sampler
  if (id_cluster == 12) {

    ### get prior parameters
    # scale matrix
    Sigma = diag(22.5314, nrow = ncol(X))
    # location vector
    mu = rep(0, ncol(X))

    ### run posterior simulation algorithm
    t0 = Sys.time()
    res = pSUN::logitBessitPG(
      nsim, X = X, y = y, nu = 1, mu = mu, Sigma = Sigma, verbose = 1000
    )
    tn = Sys.time()
    print(difftime(tn, t0, units = "sec"))

    ### save results
    # list
    logitLaplacitPG = list(res = res, time = difftime(tn, t0, units = "sec"))
    # save
    save(list = c("logitLaplacitPG"), file = "logitLaplacitPG.RData")

  }




  ##### logit model with Dirichlet-Laplace prior and Polya-Gamma Gibbs sampler
  if (id_cluster == 13) {

    ### get prior parameters
    # scale matrix
    Sigma = diag(9.0984, nrow = ncol(X))
    # location vector
    mu = rep(0, ncol(X))

    ### run posterior simulation algorithm
    t0 = Sys.time()
    res = pSUN::logitDirLaplPG(
      nsim, X = X, y = y, mu = mu, Sigma = Sigma, verbose = 1000
    )
    tn = Sys.time()
    print(difftime(tn, t0, units = "sec"))

    ### save results
    # list
    logitDirLaplPG = list(res = res, time = difftime(tn, t0, units = "sec"))
    # save
    save(list = c("logitDirLaplPG"), file = "logitDirLaplPG.RData")

  }



  ##### logit model with Gaussian prior and Ultimate Polya-Gamma Gibbs sampler
  if (id_cluster == 14) {

    ### get prior parameters
    # scale matrix
    Sigma = diag(52.6379, nrow = ncol(X))

    ### run posterior simulation algorithm
    t0 = Sys.time()
    res = pSUN::logitGaussUPG(
      nsim, X = X, y = y, Sigma = Sigma, verbose = 1000
    )
    tn = Sys.time()
    print(difftime(tn, t0, units = "sec"))

    ### save results
    # list
    logitGaussUPG = list(res = res, time = difftime(tn, t0, units = "sec"))
    # save
    save(list = c("logitGaussUPG"), file = "logitGaussUPG.RData")

  }



  ##### logit model with Cauchy prior and Ultimate Polya-Gamma Gibbs sampler
  if (id_cluster == 15) {

    ### get prior parameters
    # scale matrix
    Sigma = diag(1.2525, nrow = ncol(X))

    ### run posterior simulation algorithm
    t0 = Sys.time()
    res = pSUN::logitStudentUPG(
      nsim, X = X, y = y, nu = 1, Sigma = Sigma, verbose = 1000
    )
    tn = Sys.time()
    print(difftime(tn, t0, units = "sec"))

    ### save results
    # list
    logitCauchyUPG = list(res = res, time = difftime(tn, t0, units = "sec"))
    # save
    save(list = c("logitCauchyUPG"), file = "logitCauchyUPG.RData")
  }



  ##### logit model with independent Laplace priors and Ultimate Polya-Gamma
  #      Gibbs sampler
  if (id_cluster == 16) {

    ### get prior parameters
    # scale matrix
    Sigma = diag(22.5314, nrow = ncol(X))

    ### run posterior simulation algorithm
    t0 = Sys.time()
    res = pSUN::logitBessitUPG(
      nsim, X = X, y = y, nu = 1, Sigma = Sigma, verbose = 1000
    )
    tn = Sys.time()
    print(difftime(tn, t0, units = "sec"))

    ### save results
    # list
    logitLaplacitUPG = list(res = res, time = difftime(tn, t0, units = "sec"))
    # save
    save(list = c("logitLaplacitUPG"), file = "logitLaplacitUPG.RData")

  }




  ##### logit model with Dirichlet-Laplace prior and Ultimate Polya-Gamma Gibbs
  #      sampler
  if (id_cluster == 17) {

    ### get prior parameters
    # scale matrix
    Sigma = diag(9.0984, nrow = ncol(X))

    ### run posterior simulation algorithm
    t0 = Sys.time()
    res = pSUN::logitDirLaplUPG(
      nsim, X = X, y = y, Sigma = Sigma, verbose = 1000
    )
    tn = Sys.time()
    print(difftime(tn, t0, units = "sec"))

    ### save results
    # list
    logitDirLaplUPG = list(res = res, time = difftime(tn, t0, units = "sec"))
    # save
    save(list = c("logitDirLaplUPG"), file = "logitDirLaplUPG.RData")

  }

})




##### parallel environment strip

### close report file
parallel::clusterEvalQ(cl = cluster, sink())

### leave cluster
parallel::stopCluster(cluster)





################################### Closure ###################################
save(list = ls(), file = "CancerSAGE_ParEst.RData")