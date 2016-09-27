library("mvtnorm")

# Returns a list with two components:
#   $operon.data: a matrix of simulated expression data for the operon
#   $operon.calls: a matrix of 'gold standard' calls for the simulated expression data
generate.sim.data.for.operon <- function(numExps, mixing, mu0, mu1, var) {
  # numExps: the number of experiments' worth of data to generate for the operon (the number of samples from the distribution)
  # mixing: the mixing parameter to use
  # mu0: Vector with the low-component mean expression value for each gene; length must be the number of genes in the operon
  # mu1: Vector (the same length as mu0) with the high-component mean expression value for each gene
  # var: Within-state covariance matrix for the operon (N x N, where N is the number of genes in the operon; must be a matrix, even if N==1 (in which case we have a 1x1 matrix))

  colOrder <- sample(numExps)  # randomly choose a column order
  num.from.high <- as.integer(mixing*numExps)  # how many samples will come from the 'active' component

  # Generate the raw data
  if(num.from.high==0) {
    samples <- t(rmvnorm(numExps, mu0, var))
    operon.data <- rbind(samples[,colOrder,drop=FALSE])
  } else if(num.from.high==numExps) {
    samples <- t(rmvnorm(numExps, mu1, var))
    operon.data <- rbind(samples[,colOrder,drop=FALSE])
  } else {
    high.samples <- t(rmvnorm(num.from.high, mu1, var))
    low.samples <- t(rmvnorm(numExps-num.from.high, mu0, var))
    operon.data <- cbind(high.samples, low.samples)[,colOrder,drop=FALSE]
  }
  
  # Generate the true calls for the data
  operon.calls <- matrix((colOrder <= num.from.high), nr=nrow(operon.data), nc=ncol(operon.data), byrow=TRUE)

  # Return the data and calls
  return(list(operon.data=operon.data, operon.calls=operon.calls))
}
