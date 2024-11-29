###############################################################################
#
#
# Simulation Study: Importance Sampling Method
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




########## simulation study framework

### sample size
n = 100

### number of parameters
p = 500

### minimum total number of iterations (it will be ceiled according to number
#    of workers)
niter = 1e+4

### prior variance-covariance matrix
Sigma = diag(52.6379, nrow = p)

### prior mean vector
mu = rep(0, p)

### number of simulations
nsim = 1e+4

### sequence of empirical quantiles to be stored for each iteration
seq_quant = c(1:19) * 0.05




######################## Sim Study: coverage analysis #########################

########## parallel environment setup

### get cluster
cluster = parallel::makeCluster(24)

### number of iterations per worker
niterWorker = ceiling(niter / length(cluster))

### true number of iterations after ceiling
niter = niterWorker * length(cluster)

### export variables
parallel::clusterExport(cl = cluster, varlist = ls(all.names = TRUE))

### identify cluster
for (worker in 1:length(cluster)) {
  parallel::clusterExport(cl = cluster[worker], "worker")
}; rm("worker")

### open report files
invisible(
  parallel::clusterEvalQ(cl = cluster, {
    sink(paste("reportWorker", worker, ".txt", sep = ""))
  })
)

### set seed
parallel::clusterSetRNGStream(cl = cluster)




### run in parallel the external script in the previous upper level
parallel::clusterEvalQ(
  cl = cluster, expr = source("../SimStudyIS_ExeScript.R")
)




########## parallel environment strip

##### initialize merging results objects

### time in secs vector
timesSecs = double(length = niter)

### effective sample size vector
effSamplSize = double(length = niter)

### true beta (intercept and not intercept) matrix
trueBeta = matrix(nrow = 2, ncol = niter, dimnames = list(
  beta = c("intercept", "notIntercept"), iter = 1:niter
))

### empirical quantiles array
empQuant = array(
  data = NA,
  dim = c(2, length(seq_quant), niter),
  dimnames = list(
    beta = c("intercept", "notIntercept"),
    quantiles = seq_quant, iter = 1:niter
  )
)

### temporary environment for workers variables
tmpEnvWorker = new.env()



##### merge results
for (iworker in 0:(length(cluster)-1)) {

  ### load worker environment
  load(
    paste("tmpWorker", iworker + 1, ".RData", sep = ""),
    envir = tmpEnvWorker
  )

  ### get times
  timesSecs[
    (iworker * niterWorker + 1) : (iworker * niterWorker + niterWorker)
  ] = tmpEnvWorker$timesSecs

  ### get effective sample size
  effSamplSize[
    (iworker * niterWorker + 1) : (iworker * niterWorker + niterWorker)
  ] = tmpEnvWorker$effSamplSize

  ### get true beta
  trueBeta[
    , (iworker * niterWorker + 1) : (iworker * niterWorker + niterWorker)
  ] = tmpEnvWorker$trueBeta

  ### get empirical quantiles
  empQuant[
    , , (iworker * niterWorker + 1) : (iworker * niterWorker + niterWorker)
  ] = tmpEnvWorker$empQuant

}; rm(list = c("iworker", "tmpEnvWorker"))



### close report file
parallel::clusterEvalQ(cl = cluster, sink())

### leave cluster
parallel::stopCluster(cl = cluster)





################################### Closure ###################################
save(list = ls(), file = paste("SimStudyIS_n", n, "p", p, ".RData", sep = ""))