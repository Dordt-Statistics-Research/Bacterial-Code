
m <- as.matrix(temp)

XGene = "fig|83333.1.peg.1018"
YGene = "fig|83333.1.peg.1020"
  
plot.genes(XGene, YGene)
    
#for (i in 100:50) {
#  plot.genes(i, i+2)
#}

plot.genes <- function(genex, geney, hardlimits=TRUE) {
  xlab <- paste(c("Gene ", as.character(genex)), collapse=" ")
  ylab <- paste(c("Gene ", as.character(geney)), collapse=" ")
  if (hardlimits) {
    plot(m[genex,], m[geney,], xlim=c(0,16), ylim=c(0,16), ylab=ylab, xlab=xlab)
  } else {
    plot(m[genex,], m[geney,], ylab=ylab, xlab=xlab)
  }
  reg <- lm(m[geney,]~m[genex,])
  title(main=as.character(summary(reg)$r.squared))
  abline(reg)
}

#TODO
# Run all possible graphs and find r^2 values and list them from highest to lowest excluding 1's
# Get list of all genes that only take on one expression level for all of the experiments

# For each gene i
  # For each gene j
    # get slope and r^2 for lm(j~i)
    # in a separate matrices of size m x m place these values

get.r_squared_matrix <- function(data) {
  mat <- as.matrix(data)
  
  ###...TEMPORARY...###
  mat <- mat[1:100,]
  
  head(mat)
  num_genes = dim(mat)[1]
  print(num_genes)
  return_mat <- matrix(nrow=num_genes, ncol=num_genes)
  for (i in 2:num_genes) {
    for (j in 1:(i-1)) {
      model <- lm(mat[j,]~mat[i,])
      return_mat[i,j] <- summary(model)$r.squared
    }
  }
  return(return_mat)
}

r_squared_matrix <- get.r_squared_matrix(temp)


sum <- 0
thresh <- 0.95
for (i in 2:100) {
  for (j in 1:(i-1)) {
    if (r_squared_matrix[i,j] > thresh) {
      sum <- sum + 1
    }
  }
}
print(sum)
