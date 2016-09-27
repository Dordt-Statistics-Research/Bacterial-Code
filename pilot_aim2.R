source("wrAijContext.R")
source("Aij_validation.R")
source("analyze_scenarios_complete.R")

getZScores <- function(vec) {
  (vec-mean(vec))/sd(vec)
}

gene <- function(pegnums) {
  # pegnums: A vector (possibly length 1) of pegnums
  sapply(pegnums, function(pegnum) {
    if(substr(pegnum, 1, 3)=="fig")  return(pegnum)
    else return(paste0("fig|83333.1.peg.", pegnum))
  })
}

do.analysis <- function(wrAijContext, BIC.bonus, Bayesian2=FALSE, A, B, exponential) {
  # If Bayesian2 is TRUE, BIC.bonus is ignored
  # If Bayesian2 is FALSE, A, B, and exponential are ignored
  
  expression.data <- get.expression.data.outliersNA(wrAijContext)
  genenames <- get.genenames(wrAijContext)
  if(Bayesian2) aijs <- MultiMM_Bayesian2(wrAijContext, A, B, exponential)
  else aijs <- MultiMM(wrAijContext, BIC.bonus)
  cat("Got aijs\n")

  meanEij <- rowMeans(expression.data, na.rm=T)
  sdEij <- apply(expression.data, 1, sd, na.rm=T)
  #multiBICdiff.by.operon <- get.BICdiffs.by.operon(wrAijContext)
  #multiBICdiff <- rep(NA, length(genenames))
  #names(multiBICdiff) <- genenames
  #for(index in 1:length(get.operons(wrAijContext))) {
  #  operon <- get.operons(wrAijContext)[[index]]
  #  multiBICdiff[operon] <- multiBICdiff.by.operon[index]
  #}
  #uniBICdiff <- get.BICdiffs.by.gene(wrAijContext)
  logMultiBayesFactor.by.operon <- get.bayesfactors.by.operon(wrAijContext)
  logMultiBayesFactor <- rep(NA, length(genenames))
  names(logMultiBayesFactor) <- genenames
  for(index in 1:length(get.operons(wrAijContext))) {
    operon <- get.operons(wrAijContext)[[index]]
    logMultiBayesFactor[operon] <- logMultiBayesFactor.by.operon[index]
  }
  meanAij <- rowMeans(aijs)
  pctAijsOff <- rowMeans(aijs<0.5)
  
  # Find the logarithm of the product of the numbers in vec, attempting to avoid numerical problems
  logprod <- function(vec) {
    sum(log(vec))  # as opposed to log(prod(vec)), which first computes the product and risks a 0 or Inf in doing so
  }
  
  logSStat <- apply(aijs, 1, function(row) logprod(1-row))
  
  if(Bayesian2) {
    C.values.by.gene <- get.operon.Cvalues.by.gene(wrAijContext, A, B, exponential)
  } else {
    C.values.by.gene <- as.integer(!(genenames %in% unlist(get.operons(wrAijContext)[get.1comp.operons(wrAijContext, BIC.bonus)])))
    names(C.values.by.gene) <- genenames
  }
  
  if(file.exists("overallmap.Rdata")) {
    load("overallmap.Rdata")
    genenameMap <- overallmap  # a more thorough mapping that combines E_coli_common_genenames_edited.txt and Carrera_genenames.txt
    remove(overallmap)
  } else {
    source("genenameMaps.R")
    genenameMap <- getGenenameMap("inputs/E_coli_common_genenames_edited.txt")
  }
  common.names <- genenameMap[match(genenames,genenameMap$fig),]$common
  common.names[is.na(common.names)] <- genenameMap[match(genenames[is.na(common.names)],genenameMap$b),]$common
  
  if(grepl("Hope Ecoli", get.data.source(wrAijContext))) {
    scenarios <- get.scenarios("inputs/83333.1.scenarios")
    scenariosByGene <- lapply(genenames, function(genename) scenarios[sapply(scenarios, function(scenario) genename %in% scenario)])
    remove(scenarios)
    numScenarios <- sapply(scenariosByGene, length)
    scenarioNames <- sapply(scenariosByGene, function(scenariosList) paste(names(scenariosList), collapse=", "))
  } else {
    numScenarios <- rep(NA, length(genenames))
    scenarioNames <- rep(NA, length(genenames))
  }
  
  otherGenesInOperon <- sapply(genenames, function(genename) setdiff(unlist(get.operons(wrAijContext)[sapply(get.operons(wrAijContext), function(operon) genename %in% operon)]), genename))
  numGenesInOperon <- sapply(otherGenesInOperon, function(otherGenes) length(otherGenes)+1)
  
  corMatrix <- cor(t(expression.data), t(expression.data), use="pairwise.complete.obs")
  diag(corMatrix) <- NA  # put NAs on the main diagonal because we want to get maxes/mins etc of the other numbers
  
  maxCor <- apply(corMatrix, 1, max, na.rm=T)
  minCor <- apply(corMatrix, 1, min, na.rm=T)
  #meanCor <- rowMeans(corMatrix, na.rm=T)
  #sdCor <- apply(corMatrix, 1, sd, na.rm=T)
  #totalMeanCor <- mean(meanCor)
  #totalSDCor <- sd(corMatrix, na.rm=T)
  #maxCorStdGene <- (maxCor-meanCor)/sdCor
  #minCorStdGene <- (minCor-meanCor)/sdCor
  #maxCorStdOverall <- (maxCor-totalMeanCor)/totalSDCor
  #minCorStdOverall <- (minCor-totalMeanCor)/totalSDCor
  #maxCSGp <- rank(maxCorStdGene)/length(maxCorStdGene)
  #maxCSGp[maxCSGp>0.5] <- 1-maxCSGp[maxCSGp>0.5]
  #minCSGp <- rank(minCorStdGene)/length(minCorStdGene)
  #minCSGp[minCSGp>0.5] <- 1-minCSGp[minCSGp>0.5]
  #maxCSOp <- rank(maxCorStdOverall)/length(maxCorStdOverall)
  #maxCSOp[maxCSOp>0.5] <- 1-maxCSOp[maxCSOp>0.5]
  #minCSOp <- rank(minCorStdOverall)/length(minCorStdOverall)
  #minCSOp[minCSOp>0.5] <- 1-minCSOp[minCSOp>0.5]
  maxCorWithinOperon <- sapply(1:length(genenames), function(geneIndex) { 
    operonmates <- otherGenesInOperon[[geneIndex]]
    if(length(operonmates)==0) return(NA) 
    else return(max(corMatrix[geneIndex, otherGenesInOperon[[geneIndex]]]))
  })

  # Data frame columns are named according to predictor names in "Aim 2 predictors.docx"
  results <- data.frame(
    PegID = get_pegnum(genenames),
    ComName = common.names,
    #WScore = WeightedScore,
    #MeanEij = meanEij,
    #SDEij = sdEij,
    #MultiBICDiff = multiBICdiff,
    #UniBICDiff = uniBICdiff,
    LogMBF = logMultiBayesFactor,
    #MeanAij = meanAij,
    LogSStat = logSStat,
    PctOff = pctAijsOff,
    Cvalue = C.values.by.gene,
    NumSc = numScenarios,
    ScName = scenarioNames,
    NumOp = numGenesInOperon,
    OtherPegsInOperon = sapply(otherGenesInOperon, function(othergenes) ifelse(length(othergenes)==0, "(none)", paste(get_pegnum(othergenes), collapse=", ")))
    #MaxCor = maxCor,
    #PegMaxCor = get_pegnum(colnames(corMatrix)[apply(corMatrix, 1, which.max)]),
    #MaxCorStdOverall = maxCorStdOverall,
    #MaxCorOverallPseudoP = maxCSOp,
    #MaxCorStdGene = maxCorStdGene,
    #MaxCorGenePseudoP = maxCSGp,
    #MinCor = minCor,
    #PegMinCor = get_pegnum(colnames(corMatrix)[apply(corMatrix, 1, which.min)]),
    #MinCorStdOverall = minCorStdOverall,
    #MinCorOverallPseudoP = minCSOp,
    #MinCorStdGene = minCorStdGene,
    #MinCorGenePseudoP = minCSGp,
    #MaxCorOp = maxCorWithinOperon
  )
  
  # Remove fig|83333.1.rna.x genes, leaving only the fig|83333.1.peg.x genes
  results <- results[!is.na(suppressWarnings(as.numeric(results$PegID))),]
  
  #return(results[order(results$WScore, decreasing = TRUE),])
  return(results)

}

makeSubsetDropTable <- function(data.source.1, data.source.2) {
  # Doesn't actually work with arbitrary data sources (at the moment), only Hope Ecoli / Sim Hope Ecoli and subsets of those
  results1 <- do.analysis(get.wrapped.aij.context.from.ds(data.source.1))
  results2 <- do.analysis(get.wrapped.aij.context.from.ds(data.source.2))
  
  commonCols <- c("PegID", "ComName", "NumSc", "ScName", "OtherPegsInOperon")  # the columns that don't change regardless of subsetting
  commonColIndices <- match(commonCols, colnames(results1))
  diffColIndices <- (1:ncol(results1))[-commonColIndices]  # the other columns (the ones that do change with subsetting)
  colnames(results1)[diffColIndices] <- paste0(colnames(results1)[diffColIndices],"_ds1")
  colnames(results2)[diffColIndices] <- paste0(colnames(results2)[diffColIndices],"_ds2")
  
  allResults <- merge(results1[,c(commonColIndices,diffColIndices)], results2[,c(commonColIndices,diffColIndices)])  # put the common columns first
  allResults$Diff_WScore <- allResults$WScore_ds1-allResults$WScore_ds2
  allResults <- allResults[order(abs(allResults$Diff_WScore), decreasing = TRUE),]
  
  # Now interleave the columns to juxtapose ds1 and ds2
  ds1ColsStart <- length(commonCols)+1  # the first of the specifically-ds1 columns
  ds2ColsStart <- ncol(results1)+1  # the first of the specifically-ds2 columns
  colSequence <- c(1:(ds1ColsStart-1), ds2ColsStart, ds1ColsStart, ncol(allResults))  # Start with all the shared columns, then both weighted scores, then the weighted score difference
  for(i in 1:(ds2ColsStart-ds1ColsStart-1)) colSequence <- c(colSequence, ds2ColsStart+i, ds1ColsStart+i)  # add all other columns in pairs
  allResults <- allResults[,colSequence]
  write.table(allResults, file="Subset_analysis.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
}

cleanResults <- function(results, valContext) {
  actuallyInactive <- apply(!get.gold.calls(valContext), 1, all)
  resultsgenes <- gene(results$PegID)
  genesWithCalls <- resultsgenes %in% names(actuallyInactive)
  resultsgenes <- resultsgenes[genesWithCalls]  # only count the genes we have calls for
  results <- results[genesWithCalls,]  # remove the same genes from results
  actuallyInactive <- actuallyInactive[resultsgenes]  # likewise, remove the calls we don't have genes for
  
  # Now randomly permute all the genes.  This prevents the coming order() calls with many ties from accidentally falling back on a useful ranking
  permutation <- sample(length(actuallyInactive))  # generate a random permutation
  results <- results[permutation,]
  resultsgenes <- resultsgenes[permutation]
  actuallyInactive <- actuallyInactive[permutation]
  
  results$PegID <- as.integer(results$PegID)
  results$ComName <- order(results$ComName)  # convert to integer ranking
  results$ScName <- order(results$ScName)  # convert to integer ranking
  results$OtherPegsInOperon <- NULL  # remove this for this analysis
  #results$PegMaxCor <- as.integer(results$PegMaxCor)
  #results$PegMinCor <- as.integer(results$PegMinCor)
  #results$MaxCorOp[is.na(results$MaxCorOp)] <- 1  # replace NAs with 1
  return(list(results=results, resultsgenes=resultsgenes, actuallyInactive=actuallyInactive))
}

getScores <- function(results, valContext) {
  cleaned <- cleanResults(results, valContext)
  results <- cleaned$results
  resultsgenes <- cleaned$resultsgenes
  actuallyInactive <- cleaned$actuallyInactive
  remove(cleaned)
  
  computeScores <- function(ordering) {
    resultsgenes <- resultsgenes[ordering]
    actuallyInactive <- actuallyInactive[ordering]
    top50score <- sum(actuallyInactive[resultsgenes[1:50]])  # how many of the top 50 genes in results were actually always-inactive?
    avgIndex <- mean(match(names(actuallyInactive)[actuallyInactive], resultsgenes), na.rm=T)  # average index of the always-inactive genes in results
    return(invisible(c(top50=top50score, avgIndex=avgIndex)))
  }
  
  return(do.call(rbind, lapply(colnames(results), function(colname) {
    incScores <- computeScores(order(results[[colname]]))
    decScores <- computeScores(order(results[[colname]], decreasing=TRUE))
    return(matrix(c(incScores, decScores), nr=2, nc=2, byrow=TRUE, dimnames=list(paste0(colname, "_", c("inc", "dec")), c("top50", "avgIndex"))))
  })))
    
}

getGLMresults <- function(cleanedResults) {
  results <- cleanedResults$results
  actuallyInactive <- cleanedResults$actuallyInactive
  remove(cleanedResults)
  results$WScore <- NULL  # Remove WScore
  
  baseequation <- paste("actuallyInactive ~", paste(colnames(results), collapse=" + "))
  baseglm <- glm(as.formula(baseequation), data=results, family=binomial(logit))
  
  return(c(list(base=baseglm), combn(colnames(results), 2, FUN=function(pair) {
    interaction <- paste(pair, collapse=":")
    glm <- list(glm(as.formula(paste(baseequation, "+", interaction)), data=results, family=binomial(logit)))
    names(glm) <- interaction
    return(glm)
  })))
}

getPvals <- function(GLMresults) {
  GLMsummaries <- lapply(GLMresults, summary)
  pvals <- lapply(GLMsummaries, function(summary) summary$coefficients[,"Pr(>|z|)"])
  pvals[[1]]["base"] <- NA  # first pval is missing an interaction term, so we give it an NA. Furthermore, that interaction term is used later as the name of the row, so we give it the name "base"
  pvals <- lapply(1:length(pvals), function(index) {
    vec <- pvals[[index]]
    if(index==1) return(vec)  # don't change the first vector
    if(length(vec) < 19) vec <- c(vec, Interaction=NA)  # I believe adding the NA at the end (assuming last term missing) is always correct
    names(vec) <- names(GLMresults[[index]]$coefficients)
    return(vec)
  })
  names(pvals) <- sapply(pvals, function(vec) names(vec)[19])  # 19 being the column number for the interaction term
  pvals <- lapply(pvals, function(vec) {names(vec)[19] <- "Interaction"; return(vec)})  # Rename the 19th column to Interaction for all GLMs
  return(do.call(rbind, pvals))
}

getIndivGLMresults <- function(cleanedResults) {
  results <- cleanedResults$results
  actuallyInactive <- cleanedResults$actuallyInactive
  remove(cleanedResults)
  
  glms <- lapply(results, function(predictor) {
    if(is.logical(predictor)) {
      # add 1 to each cell count in the contingency table, to avoid numerical problems
      actuallyInactive <- c(actuallyInactive, FALSE, FALSE, TRUE, TRUE)
      predictor <- c(predictor, FALSE, TRUE, FALSE, TRUE)
    }
    glm(actuallyInactive ~ predictor, family=binomial(logit))
  })
  names(glms) <- colnames(results)
  return(glms)
}

getPairwiseGLMresults <- function(cleanedResults, withInteraction=FALSE) {
  results <- cleanedResults$results
  actuallyInactive <- cleanedResults$actuallyInactive
  remove(cleanedResults)
  results$WScore <- NULL  # Remove WScore
  
  glms <- combn(colnames(results), 2, simplify=FALSE, FUN=function(pair) {
    if(withInteraction) rhs <- paste(pair, collapse=" * ")
    else rhs <- paste(pair, collapse=" + ")
    glm(as.formula(paste("actuallyInactive ~ ", rhs)), data=results, family=binomial(logit))
  })
  names(glms) <- combn(colnames(results), 2, FUN=paste, collapse=",")
  return(glms)
}

library(lmtest)
getLRTresults <- function(anyGLMresults) {
  lrts <- lapply(anyGLMresults, function(glmresults) {
    lrtest(glmresults, 1:(length(glmresults$coefficients)-1))  # remove all predictors except intercept to get the null model
  })
  names(lrts) <- names(anyGLMresults)
  return(lrts)
}

#wrAijContext <- get.wrapped.aij.context.from.ds("Hope Ecoli")
#valContext <- get.aij.validation.context(wrAijContext, list(MultiMM(wrAijContext)), "mattE")

#results <- do.analysis(wrAijContext)
#write.table(results, file="aim2_3-15-16.txt", quote=FALSE, sep="\t", row.names=FALSE)

#scenarios <- get.scenarios("inputs/83333.1.scenarios")
#scenario.means <- cbind(
#  t(sapply(scenarios, function(scenario) colMeans(results[scenario, c("LogMBF", "LogSStat", "PctOff", "NumSc", "NumOp")]))),
#  matrix(sapply(scenarios, function(scenario) mean(exp(results[scenario, "LogMBF"]))), ncol = 1, dimnames=list(names(scenarios), "MBF"))
#)

#cleanedResults <- cleanResults(results, valContext)

#print(getScores(cleanedResults))

#GLMresults <- getGLMresults(cleanedResults)
#View(getPvals(GLMresults))

#indivGLMresults <- getIndivGLMresults(cleanedResults)
#indivLRTresults <- getLRTresults(indivGLMresults)
#indivpvals <- sort(sapply(indivLRTresults, function(x) x[[5]][2]))

#pairwiseGLMresults <- getPairwiseGLMresults(cleanedResults)
#pairwiseLRTresults <- getLRTresults(pairwiseGLMresults)
#pairwisepvals <- sort(sapply(pairwiseLRTresults, function(x) x[[5]][2]))

#fruitfulness <- -log10(pairwisepvals/sapply(strsplit(names(pairwisepvals), ","), function(preds) min(indivpvals[preds])))

#interactionGLMresults <- getPairwiseGLMresults(cleanedResults, withInteraction = TRUE)
#interactionLRTresults <- getLRTresults(interactionGLMresults)
#pvals <- sort(sapply(interactionLRTresults, function(x) x[[5]][2]))
