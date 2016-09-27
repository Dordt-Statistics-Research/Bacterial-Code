source("validation_plots.R")
source("analyze_scenarios_complete.R")
source("craig_utility_funcs.R")

hope.context <- get.wrapped.aij.context.from.ds("Hope Ecoli")
carrera.context <- get.wrapped.aij.context.from.ds("Carrera Ecoli")

scenarios <- get.scenarios("inputs/83333.1.scenarios")

get.merged.results <- function(BIC.bonus) {
  nonconstantColnames <- c("PegID", "ComName", "LogMBF", "LogSStat", "PctOff", "onecomp", "OtherPegsInOperon")  # columns whose information differs based on whether we're using Hope or Carrera
  load(paste0("hope_results_",BIC.bonus,".Rdata"))
  hope.results$PegID <- as.character(hope.results$PegID)
  hope.results <- renameByName(hope.results, nonconstantColnames, paste0("Hope.", nonconstantColnames))
  hope.results$Hope.ComName <- NULL
  load(paste0("carrera_results_",BIC.bonus,".Rdata"))
  carrera.results$PegID <- as.character(carrera.results$PegID)
  carrera.results <- renameByName(carrera.results, nonconstantColnames, paste0("Carrera.", nonconstantColnames))
  carrera.results$Carrera.ComName <- NULL
  carrera.results$NumSc <- NULL
  carrera.results$ScName <- NULL
  carrera.results$NumOp <- NULL
  
  carrera.results$Hope.genename <- b_to_fig(carrera.results$Carrera.PegID)
  carrera.results$Hope.genename[is.na(carrera.results$Hope.genename)] <- paste(carrera.results$Carrera.PegID[is.na(carrera.results$Hope.genename)], "(no fig name found)")
  hope.results$Hope.genename <- hope.gene(hope.results$Hope.PegID)
  merged.results <- merge(hope.results, carrera.results, by="Hope.genename", all=TRUE)
  rownames(merged.results) <- merged.results$Hope.genename
  return(merged.results)
}

get.all.results <- function() {
  all.results.12 <- get.merged.results(12)
  all.results.0 <- get.merged.results(0)
  
  nonconstantColnames.bonus <- c("Hope.LogMBF", "Carrera.LogMBF", "Hope.LogSStat", "Carrera.LogSStat", "Hope.PctOff", "Carrera.PctOff", "Hope.onecomp", "Carrera.onecomp")  # columns whose information differs based on the BIC.bonus parameter
  all.results.12 <- renameByName(all.results.12, nonconstantColnames.bonus, paste0(nonconstantColnames.bonus, ".12"))
  all.results.0 <- renameByName(all.results.0, nonconstantColnames.bonus, paste0(nonconstantColnames.bonus, ".0"))
  all.results <- merge(all.results.12, all.results.0, all=TRUE)
  rownames(all.results) <- all.results$Hope.genename
  return(all.results)
}

get.scenario.means <- function(all.results) {
  scenario.means <- as.data.frame(do.call(rbind, lapply(scenarios, function(scenario) {
    colMeans(all.results[all.results$Hope.genename %in% scenario,  c("Hope.LogMBF.12", "Hope.LogSStat.12", "Hope.PctOff.12", "Carrera.LogMBF.12", "Carrera.LogSStat.12", "Carrera.PctOff.12", "NumSc", "NumOp")])
  })))
  rownames(scenario.means) <- names(scenarios)
  return(scenario.means)
}

# whichfig: E.g. 1 if the 'top' experiments are to be determined based on the first fig
hope_exp_without_top <- function(figs, howmanytop, whichfig=1) {
  expdata <- get.expression.data.outliersNA(hope.context)[figs,,drop=FALSE]
  expdata[,-order(expdata[whichfig,], decreasing = TRUE, na.last=NA)[1:howmanytop]]
}
hope_exp_without_bot <- function(figs, howmanybot, whichfig=1) {
  expdata <- get.expression.data.outliersNA(hope.context)[figs,,drop=FALSE]
  expdata[,-order(expdata[whichfig,], decreasing = FALSE, na.last=NA)[1:howmanybot]]
}
carrera_exp_without_top <- function(figs, howmanytop, whichfig=1) {
  expdata <- get.expression.data.outliersNA(carrera.context)[fig_to_b(figs),,drop=FALSE]
  expdata[,-order(expdata[whichfig,], decreasing = TRUE, na.last=NA)[1:howmanytop]]
}
carrera_exp_without_bot <- function(figs, howmanybot, whichfig=1) {
  expdata <- get.expression.data.outliersNA(carrera.context)[fig_to_b(figs),,drop=FALSE]
  expdata[,-order(expdata[whichfig,], decreasing = FALSE, na.last=NA)[1:howmanybot]]
}

myBICdiff_gene <- function(expdata) {
  bic <- mclustBIC(expdata, G=1:2, modelNames="E")
  bic[1] - bic[2]
}
myBICdiff_operon <- function(op_expdata) {
  bic <- mclustBIC(t(op_expdata), G=1:2, modelNames="EEE")
  bic[1] - bic[2]
}

myBayesFactor_gene <- function(expdata) {
  log.likelihood.1 <- Mclust(data=t(expdata), G=1, modelNames="E")$loglik
  log.likelihood.2 <- Mclust(data=t(expdata), G=2, modelNames="E")$loglik
  log.likelihood.1-log.likelihood.2
}
myBayesFactor_operon <- function(op_expdata) {
  log.likelihood.1 <- Mclust(data=t(op_expdata), G=1, modelNames="EEE")$loglik
  log.likelihood.2 <- Mclust(data=t(op_expdata), G=2, modelNames="EEE")$loglik
  log.likelihood.1-log.likelihood.2
}

stats_for_scenario <- function(genevec) {
  # genevec: Vector of gene names (fig format)
  print(all.results[all.results$Hope.genename %in% genevec,-1])  # leave out "Hope.genename" as it's the same as the rownames
  readline()
  for(gene in sort(genevec)) {
    hist.exp.Bayesian1(hope.context, gene, 0)
    readline()
    hist.exp.Bayesian1(carrear.context, gene, 0)
    readline()
  }
}
