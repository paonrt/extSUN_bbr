###############################################################################
#
#
# Simulation Study: Importance Sampling Execution Script
#
#
###############################################################################

##### initialize variables

### get name for temporary .RData file of the worker
tmpName = paste("tmpWorker", worker, sep = "")

### get number of iterations per worker
nIter = niterWorker

### time in secs vector
timesSecs = double(length = nIter)

### effective sample size vector
effSamplSize = double(length = nIter)

### true beta (intercept and not intercept) matrix
trueBeta = matrix(nrow = 2, ncol = nIter, dimnames = list(
  beta = c("intercept", "notIntercept"), iter = 1:nIter
))

### empirical quantiles array
empQuant = array(
  data = NA,
  dim = c(2, length(seq_quant), nIter),
  dimnames = list(
    beta = c("intercept", "notIntercept"),
    quantiles = seq_quant,
    iter = 1:nIter
  )
)

### print start time
cat(paste("start time at ", Sys.time(), "\n", sep = ""))




########## main body execution
for (inIter in 1:nIter) {

  ### stopping variable for while loop
  toDo = TRUE

  while (toDo) {

    ### reset stopping variable
    toDo = FALSE

    ### sample true beta
    betaTrue = pi * 4 / sqrt(3) * rnorm(p)

    ### sample and scale design matrix
    # sampling
    X = matrix(rnorm(n * p), nrow = n, ncol = p)
    # scaling
    X = 0.5 * scale(X)

    ### sample binary outputs
    y = rbinom(n , 1, plogis(X %*% betaTrue))

    ### likelihood matrix
    B = 2 * diag(y) - diag(1, nrow = length(y))

    ### posterior parameters
    A = B %*% X %*% diag(sqrt(diag(Sigma)), nrow = ncol(X))
    b = B %*% X %*% mu
    xi = mu
    Omega = Sigma

    ### try posterior simulation
    t0 = Sys.time()
    res_pSUN_IS = tryCatch(
      pSUN::pSUN_LogitGauss_IS(nsim, A = A, b = b, xi = xi, Omega = Omega),
      error = function(err) return(FALSE)
    )
    tn = Sys.time()

    ### check if it is ok
    if (isFALSE(res_pSUN_IS[1])) toDo = TRUE

  }

  ### get time in secs
  timesSecs[inIter] = difftime(tn, t0, units = "secs")

  ### get effective sample size
  effSamplSize[inIter] = 1 / sum( res_pSUN_IS$normWeights^2 )




  ########## quantiles

  ### sort for the intercept
  sortObjInt = sort.int( res_pSUN_IS$sample[1,], index.return = TRUE)
  sortWeiInt = res_pSUN_IS$normWeights[sortObjInt$ix]
  sortDraInt = res_pSUN_IS$sample[1, sortObjInt$ix]

  ### sort for the not intercept
  sortObjNot = sort.int( res_pSUN_IS$sample[2,], index.return = TRUE)
  sortWeiNot = res_pSUN_IS$normWeights[sortObjNot$ix]
  sortDraNot = res_pSUN_IS$sample[2, sortObjNot$ix]

  ### empirical CDF of the intercept
  empCDF_Int = cumsum(sortWeiInt)

  ### empirical CDF of the not intercept
  empCDF_Not = cumsum(sortWeiNot)



  ##### get quantile
  for (iquant in 1:length(seq_quant)) {

    ### intercept
    # get brackets
    id = max(c(1:nsim)[empCDF_Int <= seq_quant[iquant]])
    # get quantile via linear interpolation
    empQuant[1, iquant, inIter] = sortDraInt[id] +
      (seq_quant[iquant] - empCDF_Int[id]) /
        (empCDF_Int[id + 1] - empCDF_Int[id]) *
        (sortDraInt[id + 1] - sortDraInt[id])

    ### not intercept
    # get brackets
    id = max(c(1:nsim)[empCDF_Not <= seq_quant[iquant]])
    # get quantile via linear interpolation
    empQuant[2, iquant, inIter] = sortDraNot[id] +
      (seq_quant[iquant] - empCDF_Not[id]) /
        (empCDF_Not[id + 1] - empCDF_Not[id]) *
        (sortDraNot[id + 1] - sortDraNot[id])

  }

  ### get true beta
  trueBeta[, inIter] = betaTrue[1:2]

  ### print time of completation
  cat(paste(
    "iteration ", inIter, " of ", nIter, " complete at time ", Sys.time(), "\n",
    sep = ""
  ))
}




### after for loop save data
save(list = ls(all.names = TRUE), file = paste(tmpName, ".RData", sep = ""))

### print end time
cat(paste(
  "end time at ", Sys.time(), "\n", sep = ""
))