source("simulate_data.R")
library(mclust)

lbf.gene <- function(expdata) {
  log.likelihood.1 <- Mclust(data=expdata, G=1, modelNames="E")$loglik
  log.likelihood.2 <- Mclust(data=expdata, G=2, modelNames="E")$loglik
  return(log.likelihood.1-log.likelihood.2)
}
lbf.operon <- function(operon.expdata) {
  if(nrow(operon.expdata)==1) return(lbf.gene(as.vector(operon.expdata)))
  log.likelihood.1 <- Mclust(data=t(operon.expdata), G=1, modelNames="EEE")$loglik
  log.likelihood.2 <- Mclust(data=t(operon.expdata), G=2, modelNames="EEE")$loglik
  return(log.likelihood.1-log.likelihood.2)
}

nbf.gene <- function(expdata) lbf.gene(expdata)/length(expdata)
nbf.operon <- function(operon.expdata) lbf.operon(operon.expdata)/ncol(operon.expdata)

plot.NBF.by.meangap <- function(operonSize=1, numPoints=300, meangapMin=0, meangapMax=6, numExps=1000, mixing=0.5) {
  meangaps <- runif(n = numPoints, min = meangapMin, max = meangapMax)
  plot(meangaps, sapply(meangaps, function(meangap) {
    nbf.operon(generate.sim.data.for.operon(numExps, mixing, rep(0,operonSize), rep(meangap,operonSize), matrix(1))$operon.data)
  }), xlab="Mean gap", ylab="NBF", main=paste0("NBF vs. mean gap\n(",operonSize,"-gene operons with ",numExps," experiments and mixing ",mixing,")"))
}

plot.NBF.by.mixing <- function(operonSize=1, numPoints=300, mixingMin=0.2, mixingMax=0.8, meangap=4, numExps=1000) {
  mixings <- runif(n = numPoints, min = mixingMin, max = mixingMax)
  plot(mixings, sapply(mixings, function(mixing) {
    nbf.operon(generate.sim.data.for.operon(numExps, mixing, rep(0,operonSize), rep(meangap,operonSize), matrix(1))$operon.data)
  }), cex=0.5, xlab="Mixing", ylab="NBF", main=paste0("NBF vs. mixing\n(",operonSize,"-gene operons with mean gap ",meangap," and ",numExps," experiments)"))
}

plot.LBF.by.experiments <- function(operonSize=1, numPoints=300, numExperimentsMax=5000, meangap=4, mixing=0.5) {
  numexps <- round(runif(n = numPoints, min = 1, max = numExperimentsMax))
  plot(numexps, sapply(numexps, function(numexp) {
    lbf.operon(generate.sim.data.for.operon(numexp, mixing, rep(0,operonSize), rep(meangap,operonSize), diag(operonSize))$operon.data)
  }), xlim=c(0,numExperimentsMax), xlab="Number of experiments", ylab="Ln(Bayes Factor)", main=paste0("Ln(Bayes Factor) vs. number of experiments\n(",operonSize,"-gene operons, mean gap ",meangap,")"))
}

plot.LBF.by.experiments.colorByOperonSize <- function(numPoints=300, operonSizeMin=1, operonSizeMax=8, numExperimentsMax=5000, meangap=4, mixing=0.5) {
  numexps <- round(runif(n = numPoints, min = 1, max = numExperimentsMax))
  operonSize <- round(runif(n = numPoints, min=operonSizeMin-0.5, max=operonSizeMax+0.5))
  lbfs <- sapply(1:numPoints, function(index) {
    cat("Generating point", index, "; operonSize =", operonSize[index], "and numexps =", numexps[index],"\n")
    lbf.operon(generate.sim.data.for.operon(numexps[index], mixing, rep(0,operonSize[index]), rep(meangap,operonSize[index]), diag(operonSize[index]))$operon.data)
  })
  colors <- rainbow(operonSizeMax-operonSizeMin+1, end=0.8)
  plot(numexps, lbfs, col=colors[operonSize], xlim=c(0,numExperimentsMax), 
       xlab="Number of experiments", ylab="Ln(Bayes Factor)", main=paste0("Ln(Bayes Factor) vs. number of experiments\n(mean gap ", meangap,")"))
  legend("bottomleft", legend=paste0(operonSizeMin:operonSizeMax,"-gene operons"), fill=colors)
  return(invisible(data.frame(lbf=lbfs,numexps=numexps,operonSize=operonSize)))
}

plot.NBF.by.operonSize <- function(numPointsPerSize=50, operonSizeMin=1, operonSizeMax=10, numExps=1000, meangap=4, mixing=0.5) {
  returnme <- lapply(operonSizeMin:operonSizeMax, function(operonSize) {
    cat("Starting points for", operonSize, "\n")
    sapply(1:numPointsPerSize, function(x) nbf.operon(generate.sim.data.for.operon(numExps, mixing, rep(0,operonSize), rep(meangap,operonSize), diag(operonSize))$operon.data))
  })
  boxplot(returnme, names=operonSizeMin:operonSizeMax, ylab="NBF", main=paste0("NBF vs. operon size\n(mean gap ",meangap," and ",numExps," experiments)"))
  return(invisible(returnme))
}
