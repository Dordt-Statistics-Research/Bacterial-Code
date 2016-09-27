source("aijContext.R")
library(mclust)

# The code in this script can be very compute-intensive and therefore take a long time to run.
# Best run on a computing cluster. 

get.aijs.all.methods <- function(expression.data, operons=list(), outlierHandling="None") {
  # expression.data: A matrix of expression data, each row is a gene and each column is an experiment
  #   Both rows and columns should be named appropriately
  # operons: A list of character vectors
  #   Each character vector represents one operon, and each element of a character vector
  #   should be a gene name which matches one of the gene names in expression.data.
  #   Any gene names in operons which are not also rownames of expression.data will be ignored.
  #   No gene should appear in two different operons (or more than once in the same operon).
  #   Any genes present in expression.data but not in operons will be assumed to be in single-gene operons.
  #   Therefore the default behavior (operons=list()) assumes every gene is in its own operon. 
  #   With operons=list(), there is no difference between the UniMM and MultiMM methods, as the 
  #   distinction between those methods is only that MultiMM incorporates operon structure. 
  # outlierHandling: Default "None", for other options see comments in aijContext.R on get.aij.context()
  
  aijContext <- get.aij.context(expression.data, operons, outlierHandling)
  return(list(
    MT = MT(aijContext),
    TT = TT(aijContext),
    RB = RB(aijContext),
    UniMM_12 = UniMM(aijContext, 12),
    MultiMM_12 = MultiMM(aijContext, 12),
    UniMM_0 = UniMM(aijContext, 0),
    MultiMM_0 = MultiMM(aijContext, 0),
    UniMM_NBF = UniMM(aijContext),
    MultiMM_NBF = MultiMM(aijContext)
  ))
}

coinflip <- function(context) UseMethod("coinflip")  # An additional aij-generating method that literally coinflips to determine calls of either 1 or 0
MT <- function(context) UseMethod("MT")
TT <- function(context, bufferSize=0.1) UseMethod("TT")
RB <- function(context) UseMethod("RB")
UniMM <- function(context, BIC.bonus) UseMethod("UniMM")
MultiMM <- function(context, BIC.bonus) UseMethod("MultiMM")
ConfMM <- function(context, BIC.bonus, confidence) UseMethod("ConfMM")  # 'Confidence' should be a vector the same length as get.operons(context), giving confidences between 0 and 1 inclusive
UniMM_NBF <- function(context) UseMethod("UniMM_NBF")
MultiMM_NBF <- function(context) UseMethod("MultiMM_NBF")
ConfMM_NBF <- function(context, confidence) UseMethod("ConfMM_NBF")  # 'Confidence' should be a vector the same length as get.operons(context), giving confidences between 0 and 1 inclusive

# S3 method for class 'aijContext'
coinflip.aijContext <- function(context) {
  expdata <- get.expression.data.outliersNA(context)
  returnme <- matrix(rbinom(length(expdata), 1, 0.5), nr=nrow(expdata), nc=ncol(expdata), dimnames=dimnames(expdata))
  return(returnme)
}

# S3 method for class 'aijContext'
MT.aijContext <- function(context) {
  aijs <- 1*(RB(context) >= 0.5)  # simple dichotomization of RB. Note that any outliers that were removed, become 0.5's for RB and were automatically rounded up to 1 here
  aijs[get.expression.data.outliersInf(context) == -Inf] <- 0  # revise the above so that low outliers get 0 not 1
  return(aijs)
}

# S3 method for class 'aijContext'
TT.aijContext <- function(context, bufferSize=0.1) {
  # Expression values with percentiles within bufferSize of 0.5 are given aij of 0.5. Others are called 0 or 1 a la MT.
  aijs <- RB(context)
  aijs[aijs > 0.5+bufferSize] <- 1
  aijs[aijs < 0.5-bufferSize] <- 0
  aijs[aijs != 0 & aijs != 1] <- 0.5  # everything remaining gets 0.5
  return(aijs)
}

# S3 method for class 'aijContext'
RB.aijContext <- function(context) {
  ranks <- apply(get.expression.data.outliersNA(context), 2, rank, na.last="keep")  # na.last="keep" means the rank of any NA is NA
  aijs <- ranks/nrow(ranks)
  aijs[is.na(aijs)] <- 0.5  # expression data which was excluded for outlier reasons just gets a 0.5 "don't know"
  return(aijs)
}

# S3 method for class 'aijContext'
UniMM.aijContext <- function(context, BIC.bonus) {
  aij.matrix <- get.uni.raw.aijs(context)
  one.comp.genenames <- get.1comp.genenames(context, BIC.bonus)
  two.comp.genenames <- setdiff(get.genenames(context), one.comp.genenames)
  
  get.1comp.fits(context)  # Get and throw away the 1-component fits for the context.  This ensures
                            # it's cached and ready to go, before we start parallel processing (so we 
                            # avoid the scenario where each slave is independently trying to generate
                            # and then cache the 1-component fits).  Getting the fits from cache is
                            # cheap, so if it's already cached then this check doesn't take much time.
  
  # Compute imputed aijs in parallel
  cluster <- getCluster()
  clusterExport(cluster, c("context", "two.comp.genenames"), envir=environment())
  clusterCall(cluster, source, get.sourcefile(context))
  imputed <- parLapplyLB(cluster, one.comp.genenames, function(one.comp.genename) {
    get.imputed.aijs(context, two.comp.genenames, one.comp.genename)
  })
  if(anyNA(imputed, recursive=T)) stop("NAs found in imputed")
  stopCluster(cluster)
  names(imputed) <- one.comp.genenames

  # Now assign all of the pre-computed imputed aijs into the appropriate rows
  for(one.comp.genename in one.comp.genenames) aij.matrix[one.comp.genename,] <- imputed[[one.comp.genename]]
  if(anyNA(aij.matrix, recursive=T)) stop("Ended up with NAs in the result for UniMM")
  return(aij.matrix)
}

# S3 method for class 'aijContext'
MultiMM.aijContext <- function(context, BIC.bonus) {
  aij.matrix <- get.multi.raw.aijs(context)
  one.comp.operons <- get.1comp.operons(context, BIC.bonus)
  print("Got 1-comp operons")
  genes.in.1comp.operons <- unlist(get.operons(context)[one.comp.operons])
  genes.in.2comp.operons <- setdiff(get.genenames(context), genes.in.1comp.operons)
  
  get.1comp.fits(context)  # Get and throw away the 1-component fits for the context.  This ensures
                            # it's cached and ready to go, before we start parallel processing (so we 
                            # avoid the scenario where each slave is independently trying to generate
                            # and then cache the 1-component fits).  Getting the fits from cache is
                            # cheap, so if it's already cached then this check doesn't take much time.
  
  # Compute imputed aijs in parallel
  cluster <- getCluster()
  clusterExport(cluster, c("context", "genes.in.2comp.operons"), envir=environment())
  clusterCall(cluster, source, get.sourcefile(context))
  imputed <- parLapplyLB(cluster, one.comp.operons, function(operonnum) {
    genes.in.operon <- get.operons(context)[[operonnum]]
    get.imputed.aijs(context, genes.in.2comp.operons, genes.in.operon)
  })
  stopCluster(cluster)
  # Now assign all of the pre-computed imputed aijs into the appropriate rows
  if(length(one.comp.operons) > 0) {
    for(index in 1:length(one.comp.operons)) {
      operonnum <- one.comp.operons[[index]]
      genes.in.operon <- get.operons(context)[[operonnum]]
      aij.matrix[genes.in.operon,] <- imputed[[index]]
    }
  }
  if(anyNA(aij.matrix, recursive=T)) stop("Ended up with NAs in the result for MultiMM")
  return(aij.matrix)
}

# S3 method for class 'aijContext'
ConfMM.aijContext <- function(context, BIC.bonus, confidence) {
  if(length(get.operons(context))!=length(confidence)) {
    stop(paste0("Confidence vector length (",length(confidence),") does not equal operon list length (",length(get.operons(context)),")"))
  }
  if(any(confidence<0, confidence>1)) stop("Confidence values must be between 0 and 1 inclusive.")
  
  UniMM <- UniMM(context, BIC.bonus)
  MultiMM <- MultiMM(context, BIC.bonus)
  
  aij.matrix <- UniMM
  for(i in 1:length(get.operons(context)))
    aij.matrix[get.operons(context)[[i]],] <- confidence[[i]] * MultiMM[get.operons(context)[[i]],] + (1 - confidence[[i]]) * UniMM[get.operons(context)[[i]],]
  
  return(aij.matrix)
}

# S3 method for class 'aijContext'
UniMM_NBF.aijContext <- function(context) {
  N <- length(get.expnames(context))
  confs.by.gene <- get.Cvalues.by.gene(context, A=-5/N, B=0, exponential=FALSE)
  
  get.1comp.fits(context)  # Get and throw away the 1-component fits for the context.  This ensures
                           # it's cached and ready to go, before we start parallel processing (so we 
                           # avoid the scenario where each slave is independently trying to generate
                           # and then cache the 1-component fits).  Getting the fits from cache is
                           # cheap, so if it's already cached then this check doesn't take much time.
  
  # Compute imputed aijs in parallel
  cluster <- getCluster()
  clusterExport(cluster, c("context", "confs.by.gene"), envir=environment())
  clusterCall(cluster, source, get.sourcefile(context))
  imputed <- parLapplyLB(cluster, get.genenames(context), function(genename) get.imputed.aijs.bayesian2(context, genename, confs.by.gene))
  names(imputed) <- get.genenames(context)
  stopCluster(cluster)
  
  returnme <- do.call(rbind, lapply(get.genenames(context), function(genename) confs.by.gene[genename]*get.uni.raw.aijs(context)[genename,] + (1-confs.by.gene[genename])*imputed[[genename]]))
  dimnames(returnme) <- dimnames(get.uni.raw.aijs(context))
  if(anyNA(returnme, recursive=T)) stop("Ended up with NAs in the result for UniMM_NBF")
  return(returnme)
}

# S3 method for class 'aijContext'
MultiMM_NBF.aijContext <- function(context) {
  N <- length(get.expnames(context))
  confs.by.gene <- get.operon.Cvalues.by.gene(context, A=bquote(-exp(2*(3+2*p)/.(N))/p), B=quote(-1/p), exponential=FALSE)

  get.1comp.fits(context)  # Get and throw away the 1-component fits for the context.  This ensures
                            # it's cached and ready to go, before we start parallel processing (so we 
                            # avoid the scenario where each slave is independently trying to generate
                            # and then cache the 1-component fits).  Getting the fits from cache is
                            # cheap, so if it's already cached then this check doesn't take much time.
  
  # Compute imputed aijs in parallel
  cluster <- getCluster()
  clusterExport(cluster, c("context", "confs.by.gene"), envir=environment())
  clusterCall(cluster, source, get.sourcefile(context))
  imputed <- parLapplyLB(cluster, 1:length(get.operons(context)), function(opNum) {
    genes.in.operon <- get.operons(context)[[opNum]]
    get.imputed.aijs.bayesian2(context, genes.in.operon, confs.by.gene)
  })
  stopCluster(cluster)
  
  aij.matrix <- get.multi.raw.aijs(context)
  final.aij.matrix <- matrix(NA, nc=ncol(aij.matrix), nr=nrow(aij.matrix), dimnames=dimnames(aij.matrix))
  operons <- get.operons(context)
  for(opNum in 1:length(operons)) {
    final.aij.matrix[operons[[opNum]],] <- confs.by.gene[opNum]*aij.matrix[operons[[opNum]],] + (1-confs.by.gene[opNum])*imputed[[opNum]]
  }
  if(anyNA(final.aij.matrix, recursive=T)) stop("Ended up with NAs in the result for MultiMM_NBF")
  return(final.aij.matrix)
}

# S3 method for class 'aijContext'
ConfMM_NBF.aijContext <- function(context, confidence) {
  
  if(length(get.operons(context))!=length(confidence)) {
    stop(paste0("Confidence vector length (",length(confidence),") does not equal operon list length (",length(get.operons(context)),")"))
  }
  if(any(confidence<0, confidence>1)) stop("Confidence values must be between 0 and 1 inclusive.")
  
  UniMM <- UniMM_NBF(context)
  MultiMM <- MultiMM_NBF(context)
  
  aij.matrix <- UniMM
  for(i in 1:length(get.operons(context)))
    aij.matrix[get.operons(context)[[i]],] <- confidence[[i]] * MultiMM[get.operons(context)[[i]],] + (1 - confidence[[i]]) * UniMM[get.operons(context)[[i]],]
  
  return(aij.matrix)
}

### Private methods ###

# is value in the closed interval [target-tol, target+tol]?
isWithinTol <- function(target, tol, value) {
  if(value < target-tol) return(FALSE)
  if(value > target+tol) return(FALSE)
  return(TRUE)
}

get.imputed.aijs <- function(context, two.comp.genenames, one.comp.operon) {
  # two.comp.genenames: Gene names of all genes identified as 2-component
  # one.comp.operon: Gene names of a set of genes which are
  #   identified as 1-component and all belong to the same operon.  (Call this function once per operon.)
  #   For a single-gene operon, pass only the single gene name.
  
  one.comp.fits <- get.1comp.fits(context)
  
  get.raw.imputed.aijs <- function(genename) {
    onecompmean <- one.comp.fits[[genename]]$parameters$mean
    onecompsd <- sqrt(one.comp.fits[[genename]]$parameters$variance$sigmasq)
    
    components <- lapply(two.comp.genenames, function(twocompgenename) {
      # get a list of possible other components
      # each 'component' in the list, is itself a list (genename, 0|1) where 0|1 is 0 if the lower component matched or 1 if the upper component matched
      # therefore this function will return either a 0-length list (if neither component of this gene matched),
      # a 1-length list (if one component of this gene matched), or a 2-length list (if both components of this gene matched)
      twocompsd <- sqrt(get.uni.var(context)[twocompgenename])
      if(!isWithinTol(twocompsd, 0.1, onecompsd)) return(list())  # not within sd tolerance, so neither component matches.
      twocompmeans <- c(get.uni.mu0(context)[twocompgenename], get.uni.mu1(context)[twocompgenename])
      return(list(list(twocompgenename,0),list(twocompgenename,1))[c(isWithinTol(twocompmeans[1],0.1,onecompmean), isWithinTol(twocompmeans[2],0.1,onecompmean))])
    })
    components <- unlist(components, recursive=F)  # collect all matching components into one list of components
    
    if(length(components)==0) return(rep(0.5, length(get.expnames(context))))  # no matching components: we can't determine aijs, so guess 0.5

    if(anyNA(components, recursive = T)) stop(paste("NAs found in components for", paste(one.comp.operon,collapse=" ")))
    
    # In this matrix, each column is the aijs for this gene according to the specified component match
    compaijs <- sapply(components, function(comp) {
      mu0 <- get.uni.mu0(context)[comp[[1]]]
      mu1 <- get.uni.mu1(context)[comp[[1]]]
      sd <- sqrt(get.uni.var(context)[comp[[1]]])
      pi <- get.uni.pi(context)[comp[[1]]]
      names(pi) <- NULL  # otherwise all the aijs end up named with this gene name; we want them named only with the experiment name
      if(anyNA(c(mu0,mu1,sd,pi),rec=T)) {
        cat("mu0 =", mu0, "; mu1 =", mu1, "; sd =", sd, "; pi =", pi, "\n")
        stop("NAs found in one of the above")
      }
      aijs <- sapply(get.expression.data.outliersInf(context)[genename,], function(eij) {
        if(is.infinite(eij)) return(eij)  # we'll leave +/-Inf marked as such and process at the end
        on.dnorm <- pi*dnorm(eij, mean=mu1, sd=sd)
        off.dnorm <- (1-pi)*dnorm(eij, mean=mu0, sd=sd)
        return(on.dnorm/(on.dnorm+off.dnorm))
      })
      aijs[aijs==Inf] <- max(aijs[!is.infinite(aijs)])  # high outliers get the highest generated aij for the gene
      aijs[aijs==-Inf] <- min(aijs[!is.infinite(aijs)])  # low outliers get the lowest generated aij for the gene
      return(aijs)
    })

    return(apply(compaijs, 1, mean))
  }
  
  raw.imputed.aij.matrix <- t(sapply(one.comp.operon, get.raw.imputed.aijs))
  if(anyNA(raw.imputed.aij.matrix, recursive = T)) {
    cat("one.comp.operon was", one.comp.operon, "; raw.imputed.aij.matrix was"); str(raw.imputed.aij.matrix)
    stop("NAs found in raw.imputed.aij.matrix")
  }
  return(matrix(colMeans(raw.imputed.aij.matrix), nr=length(one.comp.operon), nc=length(get.expnames(context)), byrow=TRUE))  # average the aijs for each gene to get the aijs for the operon
}

get.imputed.aijs.bayesian2 <- function(context, operon, confs.by.gene) {
  # operon: Names of genes which all belong to the same operon.
  #   (Call this function once per operon.)
  #   For a single-gene operon, pass only the single gene name. 
  # confs.by.gene: Vector of confidences (between 0 and 1 inclusive) that each gene is 2-component.
  #   This should be the same length as get.genenames(context). 
  
  one.comp.fits <- get.1comp.fits(context)
  
  get.raw.imputed.aijs <- function(genename) {
    onecompmean <- one.comp.fits[[genename]]$parameters$mean
    onecompsd <- sqrt(one.comp.fits[[genename]]$parameters$variance$sigmasq)
    
    components <- lapply(get.genenames(context), function(othergene) {
      # get a list of possible other components
      # each 'component' in the list, is itself a list (genename, 0|1) where 0|1 is 0 if the lower component matched or 1 if the upper component matched
      # therefore this function will return either a 0-length list (if neither component of this gene matched),
      # a 1-length list (if one component of this gene matched), or a 2-length list (if both components of this gene matched)
      twocompsd <- sqrt(get.uni.var(context)[othergene])
      if(!isWithinTol(twocompsd, 0.1, onecompsd)) return(list())  # not within sd tolerance, so neither component matches.
      twocompmeans <- c(get.uni.mu0(context)[othergene], get.uni.mu1(context)[othergene])
      return(list(list(othergene,0),list(othergene,1))[c(isWithinTol(twocompmeans[1],0.1,onecompmean), isWithinTol(twocompmeans[2],0.1,onecompmean))])
    })
    components <- unlist(components, recursive=F)  # collect all matching components into one list of components

    if(length(components)==0) return(rep(0.5, length(get.expnames(context))))  # no matching components: we can't determine aijs, so guess 0.5
    
    weights <- confs.by.gene[sapply(components, `[[`, 1)]
    if(all(weights==0)) return(rep(0.5, length(get.expnames(context))))  # the only matching components are genes we know aren't 2-component, so again we can't determine aijs

    # In this matrix, each column is the aijs for this gene according to the specified component match
    compaijs <- sapply(components, function(comp) {
      mu0 <- get.uni.mu0(context)[comp[[1]]]
      mu1 <- get.uni.mu1(context)[comp[[1]]]
      sd <- sqrt(get.uni.var(context)[comp[[1]]])
      pi <- get.uni.pi(context)[comp[[1]]]
      names(pi) <- NULL  # otherwise all the aijs end up named with this gene name; we want them named only with the experiment name
      aijs <- sapply(get.expression.data.outliersInf(context)[genename,], function(eij) {
        if(is.infinite(eij)) return(eij)  # we'll leave +/-Inf marked as such and process at the end
        on.dnorm <- pi*dnorm(eij, mean=mu1, sd=sd)
        off.dnorm <- (1-pi)*dnorm(eij, mean=mu0, sd=sd)
        return(on.dnorm/(on.dnorm+off.dnorm))
      })
      aijs[aijs==Inf] <- max(aijs[!is.infinite(aijs)])  # high outliers get the highest generated aij for the gene
      aijs[aijs==-Inf] <- min(aijs[!is.infinite(aijs)])  # low outliers get the lowest generated aij for the gene
      return(aijs)
    })

    return(apply(compaijs, 1, function(row) sum(weights*row))/sum(weights))
  }
  
  raw.imputed.aij.matrix <- t(sapply(operon, get.raw.imputed.aijs))
  return(matrix(colMeans(raw.imputed.aij.matrix), nr=length(operon), nc=length(get.expnames(context)), byrow=TRUE))  # average the aijs for each gene to get the aijs for the operon
}
