# IMPORTANT: IN THIS CODE, PI IS THE PROBABILITY OF BEING INACTIVE, or (1 - the probability of being active)

library("MCMCpack")
library("mvtnorm")
library(bayesm)

find.posterior <- function(exp.data, n, mu0.mean, mu1.mean, mu0.var=3*(mu1.mean-mu0.mean)^2, mu1.var=mu0.var, wishart_scale=(mu1.mean-mu0.mean)^2)
{
# mu0 refers to the low component mean
# mu1 refers to the high component mean
# For each, .mean is the mean of their prior distribution
# and .var is the variance of their prior distribution
# wishart_scale is a single number, which will be repeated down the main diagonal of a matrix to get the "scale" parameter for the inverse Wishart distirbution describing the prior for the variance

cat("Running Gibbs sampler with", n, "iterations...\n")
    
if (is.vector(exp.data)) {
  x <- t(matrix(exp.data))  # Convert a vector argument into a 1-row matrix, not a column vector
} else {
  x <- exp.data  # leave a matrix argument alone 
}

m <- dim(x)[1]

mu0.a <- rep(mu0.mean, m)
mu1.a <- rep(mu1.mean, m)
mu0.b <- mu0.var*diag(m)
mu1.b <- mu1.var*diag(m)
df_wishart <- m+2
scale_wishart <- wishart_scale*diag(m)
pi.a <- 5
pi.b <- 5

mu0 <- list()
mu1 <- list()
var <- list()
pi <- list()
z <- matrix(nrow = n, ncol = ncol(x))
colnames(z) <- names(x)

mu0[[1]] <- mu0.a
mu1[[1]] <- mu1.a
var[[1]] <- scale_wishart
pi[[1]] <- .5

z[1,] <- generate.ind(x, mu0[[1]], mu1[[1]], var[[1]], pi[[1]])

for(i in 2:n) {
  index <- which(z[i-1,] == 0)
  if(length(index)==0) {
    cat("All assigned to cluster 2 on iteration",i,"\n")
    index <- sample(1:ncol(z), 3)  # randomly assign 3 to the other cluster
  }
  if(length(index)==ncol(z)) {
    cat("All assigned to cluster 1 on iteration",i,"\n")
    index <- sample(index, ncol(z)-3)  # randomly assign 3 to the other cluster
  }
  cluster1 <- x[,index]   # a matrix; just the columns of x which were assigned to the lower mixture component
  cluster2 <- x[,-index]  # likewise a matrix for the upper mixture component
  if(is.vector(cluster1)) {   # Two reasons this may have happened.  
                              # 1) There is only 1 row in x (one gene in the operon), or 
                              # 2) there was only one column of x (one experiment) assigned to this cluster. 
                              # This code handles either case and converts the cluster back to the appropriate-dimensioned matrix (a 1-row matrix in the first case, or a 1-column matrix in the second case.  Or even a 1x1 matrix if both cases apply.)
    cluster1 <- matrix(cluster1, nrow=m)
  }
  if(is.vector(cluster2)) {
    cluster2 <- matrix(cluster2, nrow=m)
  }
  mu0[[i]] <- update.mu(cluster1, var[[i-1]], mu0.a, mu0.b)
  mu1[[i]] <- update.mu(cluster2, var[[i-1]], mu1.a, mu1.b)
  
  if(is.na(mean(mu0[[i]]) > mean(mu1[[i]]))) {cat("i =", i, "mu0[[i]] =", mu0[[i]], "mu1[[i]] =", mu1[[i]], "\n")}  # debugging
  
  if(mean(mu0[[i]]) > mean(mu1[[i]])) {
    temp <- mu1[i]
    mu1[i] <- mu0[i]
    mu0[i] <- temp
    temp <- cluster1
    cluster1 <- cluster2
    cluster2 <- temp
  }
  
  var[[i]] <- update.var(cluster1, cluster2, mu0[[i]], mu1[[i]], df_wishart, scale_wishart)
  
  pi[[i]] <- rbeta(1, pi.a + length(cluster1), pi.b + length(cluster2))
  z[i,] <- generate.ind(x, mu0[[i]], mu1[[i]], var[[i]], pi[[i]])
}

return(list(ind=z,mu0=mu0,mu1=mu1,var=var,pi=pi))

}

generate.ind <- function(data, mu0, mu1, var, pi){
  # we calculate p, a vector of probabilities
  # where p[i] is the probability that datapoint i is inactive
  
  rooti <- backsolve(chol(var),diag(nrow(var)))
  num <- pi * exp(apply(data, 2, lndMvn, mu0, rooti))
  denom <- num + (1-pi) * exp(apply(data, 2, lndMvn, mu1, rooti))

  #  the above three lines are equivalent to these two lines:
  #   num <- pi * apply(data, 2, dmvnorm, mu0, var)
  #   denom <- num + (1-pi) * apply(data, 2, dmvnorm, mu1, var)
  #  but the three-line version is more computationally efficient (by a large factor, perhaps 50x)
  #   due to only computing the Cholesky decomposition and backsolve once, 
  #   and/or the fact that it uses a different R package which may be coded more efficiently.   
  
  p <- num/denom
  u <- runif(length(p), 0, 1)
  indicators <- ifelse(u <= p, 0, 1)
  
  return(indicators)
}

update.mu <- function(data, var, a, b){
  n <- dim(data)[2]
  inv.b <- solve(b)
  inv.var <- solve(var)
  temp <- solve(inv.b + n*inv.var)
  mean <- temp %*% (inv.b%*%a + n*inv.var%*%rowMeans(data))
  var <- temp
  mu <- rmvnorm(n =1, mean = mean, sigma = var)
  return(mu)
}

update.var <- function(cluster1, cluster2, mu0, mu1, a, b){
  n1 <- dim(cluster1)[2]
  n2 <- dim(cluster2)[2]
  m <- dim(cluster1)[1]
  df <- a + n1 + n2
  mult.by.t <- function(mat) mat %*% t(mat)
  scale <- b + mult.by.t(as.matrix(cluster1 - matrix(rep(mu0, n1), nrow=m))) + mult.by.t(as.matrix(cluster2 - matrix(rep(mu1, n2), nrow=m)))
  var <- riwish(df, scale)
  return(var)
}
