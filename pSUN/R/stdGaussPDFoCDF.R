# Standard Gaussian: Computation of the Ratio of PDF Over CDF

# This function compute the value of dnorm(x) / pnorm(x). If x < -5e+6 then,
# in order to avoid numerical issue, the asymptotic approximation is used.
#
# It is NOT exported in the namespace

############### function
# no help and no export in the namespace
stdGaussPDFoCDF = function(x) {

  # x is the vector of arguments whose the ratio must be computed.


  ###############--------------- function script

  ### initialize the results vector
  res = double(length = length(x))

  ### chek if there are small values for arguments
  ctrl = x <= -5e+6

  ### compute the ratio using logarithms if x > -5e+6
  res[!ctrl] = exp(
    dnorm(x[!ctrl], log = TRUE) - pnorm(x[!ctrl], log.p = TRUE)
  )

  ### compute the ratio using asymptotic approximation if x <= -5+6
  res[ctrl] = -x[ctrl]

  ### return the results vector
  return(res)

}