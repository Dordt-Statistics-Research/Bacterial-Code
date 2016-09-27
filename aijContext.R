# Methods for aijContexts

# Gives a vector of gene names; the same as the rownames of the expression data, and guaranteed to be in the same order
get.genenames <- function(context) UseMethod("get.genenames")

# Gives a vector of experiment names; the same as the colnames of the expression data, and guaranteed to be in the same order
get.expnames <- function(context) UseMethod("get.expnames")

# Two methods giving the expression data matrix. For both, each row is a gene and each column is an experiment, named appropriately.
# The difference is in how outliers are reported. 
# Which data points are outliers, is defined when the context is first formed, by the outlierHandling argument passed to get.aij.context()
# The first method reports these outliers as NAs in the expression data matrix
get.expression.data.outliersNA <- function(context) UseMethod("get.expression.data.outliersNA")
# The second method reports these outliers as Inf (if they were high outliers) or -Inf (if they were low outliers)
get.expression.data.outliersInf <- function(context) UseMethod("get.expression.data.outliersInf")

# Gives a list of character vectors, where each character vector represents one operon, and 
#   each element of each character vector is a gene name. All genes in the context will be present in exactly one operon.
get.operons <- function(context) UseMethod("get.operons")

# Gives a filename containing R code required for the proper use of the context
#   (for aijContext, that's this file; it will presumably be a different file for subclasses)
get.sourcefile <- function(context) UseMethod("get.sourcefile")

# Gives the univariate or multivariate aijs assuming for the moment that all genes/operons are 2-component
get.uni.raw.aijs <- function(context) UseMethod("get.uni.raw.aijs")
get.multi.raw.aijs <- function(context) UseMethod("get.multi.raw.aijs")

# Gives a named vector with the fitted mu0 parameter for each gene/operon respectively
get.uni.mu0 <- function(context) UseMethod("get.uni.mu0")
get.multi.mu0 <- function(context) UseMethod("get.multi.mu0")

# Gives a named vector with the fitted mu1 parameter for each gene/operon respectively
get.uni.mu1 <- function(context) UseMethod("get.uni.mu1")
get.multi.mu1 <- function(context) UseMethod("get.multi.mu1")

# Gives a named vector with the fitted variance for each gene
get.uni.var <- function(context) UseMethod("get.uni.var")

# Gives a list containing the fitted covariance matrix for each operon
get.multi.var <- function(context) UseMethod("get.multi.var")

# Gives a named vector with the fitted mixing parameter for each gene/operon respectively
get.uni.pi <- function(context) UseMethod("get.uni.pi")
get.multi.pi <- function(context) UseMethod("get.multi.pi")

# Gives a named vector with a value for each gene/operon respectively
# Positive BICdiff indicates the 1-component model was a better fit
get.BICdiffs.by.gene <- function(context) UseMethod("get.BICdiffs.by.gene")
get.BICdiffs.by.operon <- function(context) UseMethod("get.BICdiffs.by.operon")

# Univariate. Gives gene names from the expression data in whatever context it is passed.
# BIC.bonus: A bonus applied to the 2-component BIC before comparing it to the 1-component BIC
# Higher values of BIC.bonus will give fewer 1-component genes and more 2-component genes.
get.1comp.genenames <- function(context, BIC.bonus) UseMethod("get.1comp.genenames")

# Multivariate. Gives indices in the operons list from whatever context it is passed.
# BIC.bonus: A bonus applied to the 2-component BIC before comparing it to the 1-component BIC
# Higher values of BIC.bonus will give fewer 1-component operons and more 2-component operons.
get.1comp.operons <- function(context, BIC.bonus) UseMethod("get.1comp.operons")

# Univariate. Gives the results of Mclust on each gene's expression data, assuming the gene is 1-component
get.1comp.fits <- function(context) UseMethod("get.1comp.fits")

# Gives the natural log of the Bayes Factor (calculated by gene/operon respectively)
# Higher means more likely to be 1-component
get.bayesfactors.by.gene <- function(context) UseMethod("get.bayesfactors.by.gene")
get.bayesfactors.by.operon <- function(context) UseMethod("get.bayesfactors.by.operon")

# Gives C-values calculated by gene/operon respectively (% confidences each gene/operon is 2-component)
# B is the NBF value at which C reaches 0
# The interpretation of A differs dramatically based on whether exponential=TRUE/FALSE
#   when exponential=TRUE, A is a horizontal scaling coefficient
#   when exponential=FALSE, A is the NBF value at which C reaches 1
# For get.Cvalues.by.gene or get.Cvalues.by.operon.oldmethod, A and B must be single values
# For get.Cvalues.by.operon, A and B can be either single values, 
#   or expressions, quoted using quote() or bquote(), involving 'p' (operon size)
get.Cvalues.by.gene <- function(context, A, B, exponential) UseMethod("get.Cvalues.by.gene")
get.Cvalues.by.operon <- function(context, A, B, exponential) UseMethod("get.Cvalues.by.operon")
get.Cvalues.by.operon.oldmethod <- function(context, A, B, exponential) UseMethod("get.Cvalues.by.operon.oldmethod")
# oldmethod uses the simple Ln(K)/N normalization which does not account for differences in operon size

# Gives a named vector with the C-value for each gene, but using the calculated-by-operon C-values
# Arguments are as to get.Cvalues.by.operon() or get.Cvalues.by.operon.oldmethod() respectively
get.operon.Cvalues.by.gene <- function(context, A, B, exponential) UseMethod("get.operon.Cvalues.by.gene")
get.operon.Cvalues.by.gene.oldmethod <- function(context, A, B, exponential) UseMethod("get.operon.Cvalues.by.gene.oldmethod")

###############################################################

# Wraps up expression data, operons, and the results of some initial calculations
#   into an "aijContext" which can be passed to other methods to generate aijs
get.aij.context <- function(expression.data, operons=list(), outlierHandling="None", numGibbsSamplerRuns=500, numBurnedRuns=50, mu0.mean=8, mu1.mean=9, ...) {
  # For explanations of the format of expression.data and operons, see get.aijs.all.methods() in Aij_generation.R
  # outlierHandling: One of
  #   "None": Don't detect outliers or do anything special with them
  #   "Pct x" (where x is a number between 0 and 100 inclusive): Remove the top and bottom x% of expression values from each gene
  #   "Num x" (where x is a number between 0 and ncol(expression.data) inclusive): Remove the top and bottom x expression values from each gene
  #   "SD x" (where x is a nonnegative number): For each gene, remove all expression values outside x standard deviations from the mean
  #   "Wins x" (where x is a nonnegative number): For each gene, winsorize all expression values outside x standard deviations from the mean
  #   Importantly, outliers are always computed by gene.  However, methods which operate on an operon level (notably MultiMM
  #     and its friends) require complete operons; so for those methods, any experiment which is dropped from any gene will be
  #     dropped from all genes in that operon, not just the gene(s) it was an outlier for.
  #     This may lead to much more data being dropped than it would seem at first glance. 
  
  # remove genes referenced in operons that are not in rownames(expression.data)
  operons <- lapply(operons, function(operon) operon[operon %in% rownames(expression.data)])
  # remove any resulting zero-length operons
  operons <- operons[sapply(operons, length)!=0]
  # add single-gene operons for any genes not already appearing in operons
  genesInOperons <- unlist(operons)
  operons <- c(operons, as.list(rownames(expression.data)[!(rownames(expression.data) %in% genesInOperons)]))
  # check for duplicates
  if(anyDuplicated(genesInOperons)) {
    dups <- genesInOperons[duplicated(genesInOperons)]
    stop(paste("The following gene(s) appear more than once in operons:", paste(dups, collapse=", ")))
  }
  
  expression.data <- handleOutliers(expression.data, outlierHandling)
  
  multiAijResults <- get.aijResults(expression.data, operons, numGibbsSamplerRuns, numBurnedRuns, mu0.mean, mu1.mean, ..., multivariate=TRUE)
  uniAijResults <- get.aijResults(expression.data, operons, numGibbsSamplerRuns, numBurnedRuns, mu0.mean, mu1.mean, ..., multivariate=FALSE)
  expression.data <- expression.data[names(multiAijResults$mu0),]  # make sure rows are in the same order
  returnme <- list(multiAijResults=multiAijResults, uniAijResults=uniAijResults, 
                   expression.data=expression.data, operons=operons)
  class(returnme) <- c("aijContext", "list")  # inherit from list
  return(returnme)
}

### 'Public' S3 methods for class 'aijContext'

get.genenames.aijContext <- function(context) rownames(context$expression.data)
get.expnames.aijContext <- function(context) colnames(context$expression.data)
get.expression.data.outliersNA.aijContext <- function(context) {x <- context$expression.data; x[is.infinite(x)] <- NA; x}
get.expression.data.outliersInf.aijContext <- function(context) context$expression.data
get.operons.aijContext <- function(context) context$operons
get.sourcefile.aijContext <- function(context) "aijContext.R"
get.uni.raw.aijs.aijContext <- function(context) context$uniAijResults$aijs
get.multi.raw.aijs.aijContext <- function(context) context$multiAijResults$aijs
get.uni.mu0.aijContext <- function(context) context$uniAijResults$mu0
get.multi.mu0.aijContext <- function(context) context$multiAijResults$mu0
get.uni.mu1.aijContext <- function(context) context$uniAijResults$mu1
get.multi.mu1.aijContext <- function(context) context$multiAijResults$mu1
get.uni.var.aijContext <- function(context) context$uniAijResults$var
get.multi.var.aijContext <- function(context) context$multiAijResults$var
get.uni.pi.aijContext <- function(context) context$uniAijResults$pi
get.multi.pi.aijContext <- function(context) context$multiAijResults$pi

get.BICdiffs.by.gene.aijContext <- function(context) {
  bics <- apply(get.expression.data.outliersNA(context), 1, function(genedata) mclustBIC(genedata[!is.na(genedata)], G=c(1:2), modelNames="E"))
  return(bics[1,]-bics[2,])
}

get.BICdiffs.by.operon.aijContext <- function(context) {
  bics <- sapply(get.operons(context), function(operon) {
    operon.data <- get.expression.data.outliersNA(context)[unlist(operon),,drop=FALSE]
    operon.data <- operon.data[,!apply(operon.data, 2, function(col) any(is.na(col))),drop=FALSE]  # remove any experiment where any of this operon's genes were dropped due to outlier
    return(mclustBIC(t(operon.data), G=c(1:2), modelNames=ifelse(length(operon)==1,"E","EEE")))
  })
  return(bics[1,]-bics[2,])
}

get.1comp.genenames.aijContext <- function(context, BIC.bonus) {
  return(get.genenames(context)[get.BICdiffs.by.gene(context) > BIC.bonus])
}

get.1comp.operons.aijContext <- function(context, BIC.bonus) {
  return(which(get.BICdiffs.by.operon(context) > BIC.bonus))
}

get.1comp.fits.aijContext <- function(context) {
  apply(get.expression.data.outliersNA(context), 1, function(genedata) Mclust(genedata[!is.na(genedata)], G=1, modelNames="X"))
}

get.bayesfactors.by.gene.aijContext <- function(context) {
  returnme <- apply(get.expression.data.outliersNA(context), 1, function(genedata) {
    genedata <- genedata[!is.na(genedata)]
    log.likelihood.1 <- Mclust(data=genedata, G=1, modelNames="E")$loglik
    log.likelihood.2 <- Mclust(data=genedata, G=2, modelNames="E")$loglik
    return(log.likelihood.1-log.likelihood.2)
  })
  names(returnme) <- get.genenames(context)
  return(returnme)
}

get.bayesfactors.by.operon.aijContext <- function(context) {
  sapply(get.operons(context), function(operon) {
    operon.data <- get.expression.data.outliersNA(context)[unlist(operon),,drop=FALSE]
    operon.data <- operon.data[,!apply(operon.data, 2, function(col) any(is.na(col))), drop=FALSE]  # remove any experiment where any of this operon's genes were dropped due to outlier
    log.likelihood.1 <- Mclust(data=t(operon.data), G=1, modelNames=ifelse(length(operon)==1,"E","EEE"))$loglik
    log.likelihood.2 <- Mclust(data=t(operon.data), G=2, modelNames=ifelse(length(operon)==1,"E","EEE"))$loglik
    return(log.likelihood.1-log.likelihood.2)
  })
}

get.Cvalues.by.gene.aijContext <- function(context, A, B, exponential) {
  if(exponential && A <= 0) stop("Exponential method requires A > 0")
  if(!exponential && A >= B) stop("Piecewise-linear method requires A < B")
  nbf <- get.bayesfactors.by.gene(context)/length(get.expnames(context))
  if(exponential) C <- -exp((nbf-B)*A)+1
  else C <- (nbf-B)/(A-B)
  C <- pmax(0, pmin(1, C))  # bound C between 0 and 1 inclusive
  names(C) <- names(nbf)
  return(C)
}

get.Cvalues.by.operon.aijContext <- function(context, A, B, exponential) {
  bf <- get.bayesfactors.by.operon(context)
  operons <- get.operons(context)
  N <- length(get.expnames(context))
  C <- sapply(1:length(operons), function(opNum) {
    p <- length(operons[[opNum]])
    nbf <- -exp(-2*bf[opNum]/N)/p
    actual.A <- eval(A, envir=list(p=length(operons[[opNum]])))  # Get the actual value for A, given this value of p
    actual.B <- eval(B, envir=list(p=length(operons[[opNum]])))  # Get the actual value for B, given this value of p
    if(exponential && actual.A <= 0) stop("Exponential method requires A > 0 for all positive integer values of p")
    if(!exponential && actual.A >= actual.B) stop("Piecewise-linear method requires A < B for all positive integer values of p")
    if(exponential) return(-exp((nbf-actual.B)*actual.A)+1)
    else return((nbf-actual.B)/(actual.A-actual.B))
  })
  C <- pmax(0, pmin(1, C))  # bound C between 0 and 1 inclusive
  names(C) <- names(bf)
  return(C)
}

get.Cvalues.by.operon.oldmethod.aijContext <- function(context, A, B, exponential) {
  if(exponential && A <= 0) stop("Exponential method requires A > 0")
  if(!exponential && A >= B) stop("Piecewise-linear method requires A < B")
  nbf <- get.bayesfactors.by.operon(context)/length(get.expnames(context))
  if(exponential) C <- -exp((nbf-B)*A)+1
  else C <- (nbf-B)/(A-B)
  C <- pmax(0, pmin(1, C))  # bound C between 0 and 1 inclusive
  names(C) <- names(nbf)
  return(C)
}

get.operon.Cvalues.by.gene.aijContext <- function(context, A, B, exponential) {
  get.values.by.gene(context, get.Cvalues.by.operon(context, A, B, exponential))
}

get.operon.Cvalues.by.gene.oldmethod.aijContext <- function(context, A, B, exponential) {
  get.values.by.gene(context, get.Cvalues.by.operon.oldmethod(context, A, B, exponential))
}

### Private methods to be called only in this file ###

# Given a vector of values corresponding to operons, return a named vector of the same values where names are genes 
#   and each gene has the value that was assigned to its operon
get.values.by.gene <- function(context, values.by.operon) {
  values.by.gene <- rep(NA, length(get.genenames(context)))
  names(values.by.gene) <- get.genenames(context)
  operons <- get.operons(context)
  for(opNum in 1:length(operons)) values.by.gene[operons[[opNum]]] <- values.by.operon[opNum]
  if(anyNA(values.by.gene, recursive=TRUE)) stop("NAs found in values.by.gene")
  return(values.by.gene)
}

handleOutliers <- function(expression.data, method) {
  # method: See comments on "outlierHandling" argument to get.aij.context; its options and behavior are identical
  
  #   "None": Don't detect outliers or do anything special with them
  #   "Pct x" (where x is a number between 0 and 100): Remove the top and bottom x% of expression values from each gene
  #   "Num x" (where x is a number): Remove the top and bottom x expression values from each gene
  #   "SD x" (where x is a number): For each gene, remove all expression values outside x standard deviations from the mean
  #   "Wins x" (where x is a number): For each gene, winsorize all expression values outside x standard deviations from the mean
  #   Importantly, outliers are always computed by gene.  However, methods which operate on an operon level (notably MultiMM
  #     and its friends) require complete operons; so for those methods, any experiment which is dropped from any gene will be
  #     dropped from all genes in that operon, not just the gene(s) it was an outlier for.
  #     This may lead to much more data being dropped than it would seem at first glance. 
  
  if(method=="None") return(expression.data)
  if(any(is.na(expression.data))) stop("expression.data cannot contain NAs prior to outlier handling")
  
  first.space <- regexpr(" ", method)[1]
  method.name <- substr(method, 1, first.space-1)
  method.parameter <- as.numeric(substr(method, first.space+1, nchar(method)))
  if(is.na(method.parameter)) stop(paste("Invalid method:", method))
  
  handled <- t(apply(expression.data, 1, function(gene) {
    if(method.name=="Pct") {
      cutoffRank.low <- round((method.parameter/100)*length(gene))
      cutoffRank.high <- round(((100-method.parameter)/100)*length(gene))+1  # e.g. if method.parameter is 100*(1/length(gene)), then we want to remove only the highest rank, so we want cutoffRank to be length(gene); but (100-method.parameter)/100 is 1-1/length(gene); multiplying by length(gene) gives length(gene)-1; but we actually desire length(gene), so we see why the +1 is necessary
      ranks <- rank(gene)
      gene[ranks<=cutoffRank.low] <- -Inf
      gene[ranks>=cutoffRank.high] <- Inf
      return(gene)
    } else if(method.name=="Num") {
      ranks <- rank(gene)
      gene[ranks<=method.parameter] <- -Inf
      gene[ranks>=length(gene)-method.parameter+1] <- Inf  # e.g. if method.parameter is 1, this is gene[ranks>=length(gene)] <- Inf which turns only the top rank into Inf
      return(gene)
    } else if(method.name=="SD") {
      gene.mean <- mean(gene)
      gene.sd <- sd(gene)
      gene[gene > gene.mean+method.parameter*gene.sd] <- Inf
      gene[gene < gene.mean-method.parameter*gene.sd] <- -Inf
      return(gene)
    } else if(method.name=="Wins") {
      gene.mean <- mean(gene)
      gene.sd <- sd(gene)
      upperBound <- gene.mean+method.parameter*gene.sd
      lowerBound <- gene.mean-method.parameter*gene.sd
      gene[gene>upperBound] <- upperBound
      gene[gene<lowerBound] <- lowerBound
      return(gene)
    } else {
      stop(paste("Unrecognized method:", method))
    }
  }))
  
  rownames(handled) <- rownames(expression.data)
  colnames(handled) <- colnames(expression.data)
  return(handled)
  
}

# Gets aijs and associated fitting parameters, assuming for the moment that all genes are 2-component
get.aijResults <- function(expression.data, operons, numGibbsSamplerRuns=500, numBurnedRuns=50, mu0.mean=8, mu1.mean=9, ..., multivariate=TRUE) {
  # expression.data: should be in its final form (post-outlier-handling) as would be returned by get.expression.data.outliersInf() later
  
  cluster <- getCluster()
  # Spread all needed data to the other processes
  clusterExport(cluster, c("expression.data", "operons", "numGibbsSamplerRuns", "numBurnedRuns", "multivariate"), envir = environment())
  clusterCall(cluster, source, "gibbs_sampler_multivariate.R")
  
  resultsList <- parLapplyLB(cluster, operons, function(operon) {
    if(multivariate) {
      operon.data <- expression.data[operon,,drop=FALSE]
      operon.data <- operon.data[,!apply(operon.data, 2, function(col) any(is.infinite(col))), drop=FALSE]  # remove any experiment where any of this operon's genes were dropped due to outlier
      if(anyNA(operon.data, recursive=T)) {
        cat("operon.data was:\n"); print(operon.data)
        stop("NA found in operon.data")
      }
      MGSresults <- find.posterior(operon.data, numGibbsSamplerRuns, mu0.mean, mu1.mean, ...)
      operon.indicators <- rbind(MGSresults$ind, 0, 1)  # The 'Greco correction': average in a 0 and a 1 in addition to the indicators. Suppose we're doing 100 iterations and they all come up 0, we really shouldn't report aij=0, we should report aij<0.01 instead. Same logic applies for other aijs; let's not overstate our confidence.
      operon.aijs <- apply(operon.indicators[-(1:numBurnedRuns),],2,mean)
      names(operon.aijs) <- colnames(operon.data)
      aijs <- matrix(NA, nr=nrow(operon.data), nc=ncol(expression.data), dimnames=list(rownames(operon.data),colnames(expression.data)))
      aijs[,colnames(operon.data)] <- t(as.matrix(sapply(operon, function(gene) operon.aijs)))  # leaves the excluded experiments as NA for now
      for(gene in operon) {
        aijs[,expression.data[gene,]==-Inf] <- min(operon.aijs)  # low outlier in any gene -> all genes in the operon get the lowest generated aij for the operon
        aijs[,expression.data[gene,]==Inf] <- max(operon.aijs)  # high outlier in any gene -> all genes in the operon get the highest generated aij for the operon
        # Note: the above assumes that two different genes in the same operon won't have different-direction outliers in the same experiment.  This seems a reasonable assumption to me. 
      }
      if(anyNA(aijs, recursive=T)) {
        cat("aijs was:\n"); print(aijs)
        cat("operon.data was:\n"); print(operon.data)
        cat("operon.aijs was:\n"); print(operon.aijs)
        stop("NA found in aijs")
      }
      mu0mat <- do.call(rbind, MGSresults$mu0[-(1:numBurnedRuns)])
      mu1mat <- do.call(rbind, MGSresults$mu1[-(1:numBurnedRuns)])
      mu0 <- colMeans(mu0mat)
      mu1 <- colMeans(mu1mat)
      varsums <- matrix(0, nr=length(operon), nc=length(operon))
      nonBurnedIndices <- (numBurnedRuns+1):(length(MGSresults$var))
      lapply(MGSresults$var[-(1:numBurnedRuns)], function(covmat) {
        varsums <<- varsums + covmat
      })
      var <- varsums / length(nonBurnedIndices)
      rownames(var) <- operon
      colnames(var) <- operon
      pi <- rep(1-mean(as.vector(unlist(MGSresults$pi[-(1:numBurnedRuns)]))), length(operon))  # The 1- is because the Gibbs sampler's pi is the probability of being inactive
    } else {
      rows <- lapply(operon, function(gene) {
        gene.data <- expression.data[gene,]
        gene.data.nonOutliers <- gene.data[!is.infinite(gene.data)]  # remove any experiments which were dropped for this gene due to outlier
        if(anyNA(gene.data.nonOutliers, recursive=T)) {
          cat("gene.data.nonOutliers was:\n"); print(gene.data.nonOutliers)
          stop("NA found in gene.data.nonOutliers")
        }
        MGSresults <- find.posterior(gene.data.nonOutliers, numGibbsSamplerRuns, mu0.mean, mu1.mean, ...)
        indicators <- rbind(MGSresults$ind, 0, 1)  # The 'Greco correction': average in a 0 and a 1 in addition to the indicators. Suppose we're doing 100 iterations and they all come up 0, we really shouldn't report aij=0, we should report aij<0.01 instead. Same logic applies for other aijs; let's not overstate our confidence.
        aijs <- rep(NA, length(gene.data))
        names(aijs) <- names(gene.data)
        generated.aijs <- apply(indicators[-(1:numBurnedRuns),],2,mean)
        aijs[names(gene.data.nonOutliers)] <- generated.aijs  # leaves the excluded experiments as NA for now
        aijs[gene.data==-Inf] <- min(generated.aijs)  # low outliers get the lowest generated aij for that gene
        aijs[gene.data==Inf] <- max(generated.aijs)  # high outliers get the highest generated aij for that gene
        if(anyNA(aijs, recursive=T)) {
          cat("aijs was:\n"); print(aijs)
          cat("gene.data was:\n"); print(gene.data)
          cat("generated.aijs was:\n"); print(generated.aijs)
          stop("NA found in aijs")
        }
        mu0 <- mean(as.vector(unlist(MGSresults$mu0[-(1:numBurnedRuns)])))
        mu1 <- mean(as.vector(unlist(MGSresults$mu1[-(1:numBurnedRuns)])))
        var <- mean(as.vector(unlist(MGSresults$var[-(1:numBurnedRuns)])))
        names(mu0) <- gene
        names(mu1) <- gene
        names(var) <- gene
        pi <- 1-mean(as.vector(unlist(MGSresults$pi[-(1:numBurnedRuns)])))  # The 1- is because the Gibbs sampler's pi is the probability of being inactive
        names(pi) <- gene
        return(list(aijs=aijs, mu0=mu0, mu1=mu1, var=var, pi=pi))
      })
      aijs <- do.call(rbind, lapply(rows, function(row) row$aijs))
      rownames(aijs) <- unlist(operon)
      mu0 <- do.call(c, lapply(rows, function(row) row$mu0))
      mu1 <- do.call(c, lapply(rows, function(row) row$mu1))
      var <- do.call(c, lapply(rows, function(row) row$var))
      pi <- do.call(c, lapply(rows, function(row) row$pi))
      names(var) <- operon
    }
    names(mu0) <- operon
    names(mu1) <- operon
    names(pi) <- operon
    return(list(aijs=aijs, mu0=mu0, mu1=mu1, var=var, pi=pi))
  })
  
  aijs <- do.call(rbind, lapply(resultsList, function(operon) operon$aijs))
  mu0 <- do.call(c, lapply(resultsList, function(operon) operon$mu0))
  mu1 <- do.call(c, lapply(resultsList, function(operon) operon$mu1))
  if(multivariate) {
    var <- do.call(list, lapply(resultsList, function(operon) operon$var)) 
  } else {
    var <- do.call(c, lapply(resultsList, function(operon) operon$var)) 
  }
  pi <- do.call(c, lapply(resultsList, function(operon) operon$pi))
  
  # Clean up parallel resources  
  stopCluster(cluster)
  
  if(anyNA(aijs, recursive=T)) stop("Ended up with NAs in the final aijs for get.aijResults")
  
  return(list(aijs=aijs, mu0=mu0, mu1=mu1, var=var, pi=pi))
  
}

# This function will utilize the Rmpi package (if available) for multi-node cluster computation.
# Else no worries, still works with built-in R packages (e.g. "parallel") on a normal computer, just slower
getCluster <- function(useParallel=TRUE) {
  # useParallel: if FALSE, don't try to do anything in parallel at all, even multicore on a normal computer
  #   (this is even slower than with useParallel=TRUE on a normal computer)
  
  if(useParallel) {
    library(parallel)
    haveRMPI <- require("Rmpi", quietly=T)
    if(haveRMPI) {
      # if Rmpi is available, we will run on as many cluster nodes as we've been assigned, as long as that is more than one.
      universeSize <- mpi.universe.size()
      if(universeSize>1) {useMPI <- TRUE; numSlaves=universeSize-1}  # one process is the master, so we can spawn (n-1) slaves
      else if(universeSize==1) useMPI <- FALSE  # We have only one node, so we won't use MPI for this run
      else if(universeSize==0) {
        # Using a version of MPI which doesn't support the mpi.universe.size() call
        universeSize <- as.integer(system2("echo", args="$NUMPROCS", stdout=TRUE))  # user must set this environment variable to the total number of processes being used for this run
        if(is.na(universeSize)) useMPI <- FALSE  # user didn't set the environment variable
        else if(universeSize>1) {useMPI <- TRUE; numSlaves=universeSize-1}
        else useMPI <- FALSE  # We have only one node, so we won't use MPI for this run
      }
    } else {
      useMPI <- FALSE  # Rmpi not available
    }
    
    if(useMPI) {
      # Initialize parallel resources
      cluster <- makeCluster(numSlaves, type="MPI", outfile="")
      
      # Test that parallel resources are up and running - i.e. we
      # actually are running on all the nodes we think we are.
      getInfo <- function() { return(paste("Process", mpi.comm.rank(), "running on node", Sys.info()["nodename"])) }
      cat(sep="", paste(clusterCall(cluster, getInfo), "\n"))
    } else {
      # Else, we will run using multicore computation on a single machine.
      cat("Rmpi not found and/or only one node detected. Running on single node / single machine.\n")
      cluster <- makePSOCKcluster(getOption("mc.cores", 2L))
    }
    
    return(cluster)
    
  } else {
    # Define stubs for the functions in library(parallel) that we'd actually use
    clusterExport <<- function(...) NULL  # do nothing
    clusterCall <<- function(cluster, func, ...) func(...)  # ignore the 'cluster' argument and execute locally
    stopCluster <<- function(...) NULL  # do nothing
    parLapply <<- function(cluster, ..) lapply(...)  # ignore the 'cluster' argument and execute locally
    parLapplyLB <<- function(cluster, ...) lapply(...)  # ignore the 'cluster' argument and execute locally
    parSapply <<- function(cluster, ..) sapply(...)  # ignore the 'cluster' argument and execute locally
    parSapplyLB <<- function(cluster, ...) sapply(...)  # ignore the 'cluster' argument and execute locally
    return(NULL)  # no cluster to return
  }
}
