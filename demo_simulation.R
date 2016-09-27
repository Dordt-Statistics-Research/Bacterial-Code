source("simulate_data.R")
source("wrAijContext.R")

run.tests <- function() {
  operon.data <- lapply(c(2, 3, 4), function(x) test.simulated.data(x, 907, .7, 100))
  operon.deviations <- do.call(rbind, lapply(operon.data, function(x) x$mean.deviations["all",]))
  operon.squared.deviations <- do.call(rbind, lapply(operon.data, function(x) x$mean.squared.deviations["all",]))
  rownames(operon.deviations) <- c("2gene", "3gene", "4gene")
  rownames(operon.squared.deviations) <- c("2gene", "3gene", "4gene")
  write.csv(1-operon.deviations, file="operon_deviations.csv")
  write.csv(1-operon.squared.deviations, file="operon_squared_deviations.csv")
  write.csv(1-sqrt(operon.squared.deviations), file="operon_sqrt_squared_deviations.csv")
  
  conf.data <- lapply(c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1), function(x) test.simulated.data(3, 907, x, 100))
  conf.deviations <- do.call(rbind, lapply(conf.data, function(x) x$mean.deviations["all",]))
  conf.squared.deviations <- do.call(rbind, lapply(conf.data, function(x) x$mean.squared.deviations["all",]))
  rownames(conf.deviations) <- c(".0conf", ".1conf", ".2conf", ".3conf", ".4conf", ".5conf", ".6conf", ".7conf", ".8conf", ".9conf", "1.0conf")
  rownames(conf.squared.deviations) <- c(".0conf", ".1conf", ".2conf", ".3conf", ".4conf", ".5conf", ".6conf", ".7conf", ".8conf", ".9conf", "1.0conf")
  write.csv(1-conf.deviations, file="conf_deviations.csv")
  write.csv(1-conf.squared.deviations, file="conf_squared_deviations.csv")
  write.csv(1-sqrt(conf.squared.deviations), file="conf_sqrt_squared_deviations.csv")
}


# Returns a list with two components, "mean.deviations" and "mean.squared.deviations"
# Each component is a 3x3 data.frame with
# columns "UniMM", "MultiMM", "ComboMM" and
# rows "all", "operon", and "independent"
test.simulated.data <- function(operonSize, numExps, percentConf, numTrials, dataFromSaved=TRUE) {
  # operonSize: Number of genes in the operon
  # percentConf: Percentage of the trials which should be actual operons.  The others will be independent genes.
  # numTrials: Number of trials to run
  
  savefile <- paste0("Bayesian2-0_Sim_Temp/",operonSize,"gene_",numExps,"exps_",
                     as.integer(percentConf*100),"conf_",numTrials,"trial",".Rdata")
  
  if(dataFromSaved && file.exists(savefile)) {
    load(savefile)
  } else {
    data <- get.simulated.data(operonSize, numExps, percentConf, numTrials)
    aijContexts <- lapply(data$operon.data, get.aij.context, list(paste0("gene",1:operonSize)))
    all.aijs <- lapply(aijContexts, CompareMM, percentConf)
    save(data, all.aijs, file=savefile)
  }
  
  num.actual.operon.trials <- as.integer(percentConf*numTrials)
  actual.operon.aijs <- head(all.aijs, num.actual.operon.trials)
  independent.aijs <- tail(all.aijs, numTrials-num.actual.operon.trials)
  
  all.calls <- do.call(rbind, data$operon.calls)
  all.aijs <- list(
    do.call(rbind, lapply(all.aijs, function(x) x$UniMM)),
    do.call(rbind, lapply(all.aijs, function(x) x$MultiMM)),
    do.call(rbind, lapply(all.aijs, function(x) x$ComboMM))
  )
  
  actual.operon.calls <- do.call(rbind, head(data$operon.calls, num.actual.operon.trials))
  actual.operon.aijs <- list(
    do.call(rbind, lapply(actual.operon.aijs, function(x) x$UniMM)),
    do.call(rbind, lapply(actual.operon.aijs, function(x) x$MultiMM)),
    do.call(rbind, lapply(actual.operon.aijs, function(x) x$ComboMM))
  )
  
  independent.calls <- do.call(rbind, tail(data$operon.calls, numTrials-num.actual.operon.trials))
  independent.aijs <- list(
    do.call(rbind, lapply(independent.aijs, function(x) x$UniMM)),
    do.call(rbind, lapply(independent.aijs, function(x) x$MultiMM)),
    do.call(rbind, lapply(independent.aijs, function(x) x$ComboMM))
  )
  
  deviations <- c(
    lapply(all.aijs, get.deviations, all.calls),
    lapply(actual.operon.aijs, get.deviations, actual.operon.calls),
    lapply(independent.aijs, get.deviations, independent.calls)
  )
  
  methodnames <- c("UniMM", "MultiMM", "ComboMM")
  setnames <- c("all", "operon", "independent")
  
  mean.deviations <- matrix(unlist(lapply(deviations, mean)), nrow=3, ncol=3, byrow=TRUE)
  colnames(mean.deviations) <- methodnames
  rownames(mean.deviations) <- setnames
  
  mean.squared.deviations <- matrix(unlist(lapply(deviations, function(x) mean(x*x))), nrow=3, ncol=3, byrow=TRUE)
  colnames(mean.squared.deviations) <- methodnames
  rownames(mean.squared.deviations) <- setnames
  
  return(list(
    mean.deviations = mean.deviations,
    mean.squared.deviations = mean.squared.deviations
  ))
}


# Returns a vector of the deviations between the two lists
get.deviations <- function(aijs, gold.calls) {
  # aijs: an NxM list of predictions
  # calls: an NxM list with actual values
  # where the rows (R) are genes
  # where the columns (C) are expiriments
  if(!(is.null(nrow(aijs)) && is.null(nrow(gold.calls)))) {
    if(nrow(aijs)!=nrow(gold.calls))
      stop(paste0("The number of aij genes (",nrow(aijs),
                  ") does not equal the number of gold_call genes (",nrow(gold.calls),")"))
    if(ncol(aijs)!=ncol(gold.calls))
      stop(paste0("The number of aij expiriments (",ncol(aijs),
                  ") does not equal the number of gold_call experiments (",ncol(gold.calls),")"))
  }
  return(abs(aijs - gold.calls))
}


# Returns a list with two components, "operon.data" and "operon.calls"
# "operon.data" is a list numTrials long
#   Each component of operon.data is a matrix of simulated expression data operonSize x numExps
# "operon.calls" is a list numTrials long
#   Each component of operon.calls is a matrix of 'gold standard' calls for the corresponding matrix in operon.data
get.simulated.data <- function(operonSize, numExps, percentConf, numTrials) {
  # operonSize: Number of genes in the operon
  # percentConf: Percentage of the trials which should be actual operons.  The others will be independent genes.
  # numTrials: Number of trials to run
  
  num.actual.operon.trials <- as.integer(percentConf*numTrials)
  
  if(num.actual.operon.trials > 0) {
    actual.operon.data <- lapply(1:num.actual.operon.trials, function(i) {
      generate.sim.data.for.operon(
        numExps = numExps,
        mixing = runif(1, min=0.2, max=0.8),
        mu0 = rep(8, operonSize),
        mu1 = rep(10, operonSize),
        var = diag(operonSize)  # an identity matrix operonSize x operonSize
      )
    })
  } else {
    actual.operon.data <- list()
  }
  
  if(num.actual.operon.trials < numTrials) {
    independent.data <- lapply(1:(numTrials-num.actual.operon.trials), function(i) {
      data <- lapply(1:operonSize, function(j) {
        generate.sim.data.for.operon(
          numExps = numExps,
          mixing = runif(1, min=0.2, max=0.8),
          mu0 = 8,
          mu1 = 10,
          var = matrix(1)  # a 1x1 matrix
        )
      })
      
      return(list(
        operon.data = do.call(rbind, lapply(data, function(x) x$operon.data)),
        operon.calls = do.call(rbind, lapply(data, function(x) x$operon.calls))
      ))
      
    })
  } else {
    independent.data <- list()
  }
  
  all.data <- c(actual.operon.data, independent.data)
  
  return(list(
    operon.data = lapply(all.data, function(x) {
      rownames(x$operon.data) <- paste0("gene",1:nrow(x$operon.data))
      return(x$operon.data)
    }),
    operon.calls = lapply(all.data, function(x) {
      rownames(x$operon.calls) <- paste0("gene",1:nrow(x$operon.calls))
      return(x$operon.calls)
    })
  ))
}
