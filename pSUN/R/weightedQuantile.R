# Quantile Computation for a Weighted Sample

# This function compute the quantiles of a weighted sample.
#
# It is NOT exported in the namespace

############### function
# no help and no export in the namespace
weightedQuantile = function(x, w, probs = seq(0, 1, 0.25)) {

  # x numeric vector of values of sample,
  # w numeric vector of normalized weights,
  # probs numeric vector of probabilities with values in [0, 1].

  res = double(length(probs))
  names(res) = paste(probs * 100, "%", sep = "")

  sorted = sort.int(x, index.return = TRUE)
  sorted[["w"]] = w[sorted$ix]
  emp_cdf = cumsum(sorted$w)

  for (i in 1:length(probs)) {

    if (probs[i] == 0) {
      res[i] = sorted$x[1]
    } else if (probs[i] == 1) {
      res[i] = tail(sorted$x, 1)
    } else {
      low = tail(c(1:length(x))[emp_cdf <= probs[i]], 1)
      res[i] = sorted$x[low] +
        (probs[i] - emp_cdf[low]) / (emp_cdf[low + 1] - emp_cdf[low]) *
          (sorted$x[low + 1] - sorted$x[low])
    }

  }

  return(res)

}
