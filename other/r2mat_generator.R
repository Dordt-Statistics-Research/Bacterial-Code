get.r_squared_matrix <- function(data) {
  mat <- as.matrix(data)
  
  
  
  num_genes = dim(mat)[1]
  return_mat <- matrix(nrow=num_genes, ncol=num_genes)
  for (i in 2:num_genes) {
    for (j in 1:(i-1)) {
      model <- lm(mat[j,]~mat[i,])
      return_mat[i,j] <- summary(model)$r.squared
    }
  }
  return(return_mat)
}

temp = readRDS(file="83333.1.Rds")
r_squared_matrix <- get.r_squared_matrix(temp)
saveRDS(r_squared_matrix, file="r2mat.Rds")