load("nonCompleteExpts_hope_ecoli.Rdata")
source("craig_utility_funcs.R")

source("expnameMaps.R")
expnameMap <- get.expname.map("inputs/Ecoli_Cel_to_chip_fromClaire_edited.tab")

get_pegnum <- function(genenames) {
  # genenames: Single genename, or vector of genenames
  return(sapply(genenames, function(genename) {
    if(substr(genename, 1, 16)=="fig|83333.1.peg.") return(substr(genename, 17, nchar(genename)))
    else return(genename)
  }))
}

# Accepts any matrix, e.g. aijs, aij calls, gold calls, whatever
remove_CEL_gz <- function(mat) {
  ncharColnames <- nchar(colnames(mat))  # just so we don't have to compute it over and over
  if(all(substr(colnames(mat), ncharColnames-6, ncharColnames) == ".CEL.gz")) {
    colnames(mat) <- substr(colnames(mat), 1, ncharColnames-7)  # Remove the .CEL.gz from each colname
  }
  return(mat)
}

analyze.scenario <- function(scenario, aijs) {
  # scenario: A vector of gene names
  # aijs: A named matrix
  
  avgAijs <- colMeans(aijs[scenario,])
  isComplete <- !(sapply(colnames(aijs), get_chipname, expnameMap=expnameMap) %in% nonCompleteExpts)

  return(list(
    pval = t.test(avgAijs ~ isComplete)$p.value, 
    meanAijNonComplete = mean(avgAijs[!isComplete]), 
    sdAijNonComplete = sd(avgAijs[!isComplete]),
    meanAijComplete = mean(avgAijs[isComplete]),
    sdAijComplete = sd(avgAijs[isComplete])
  ))
}

get.scenarios <- function(scenariofile) {
  scenariodata <- read.table(scenariofile, header=FALSE, sep="\t", row.names=NULL, col.names=c("Pathway", "Genes"), as.is=TRUE)
  scenarios <- lapply(scenariodata[["Genes"]], function(genedata) strsplit(genedata, ",")[[1]])
  names(scenarios) <- scenariodata[["Pathway"]]
  return(scenarios)
}

runAnalysis <- function(aijs, resultsFilename) {
  scenarios <- get.scenarios("inputs/83333.1.scenarios")
  scenarioresults <- lapplyWithNames(scenarios, analyze.scenario, remove_CEL_gz(aijs))
  meanAijNonComp <- sapply(scenarioresults, function(results) results$meanAijNonComplete)
  meanAijComp <- sapply(scenarioresults, function(results) results$meanAijComplete)
  all.results <- data.frame(
    Pathway=names(scenarioresults), 
    meanAijNonComp=meanAijNonComp,
    meanAijComp=meanAijComp,
    Difference=meanAijNonComp-meanAijComp,
    pval=sapply(scenarioresults, function(results) results$pval),
    GenePegIDs=format(sapply(scenarios, function(scenario) paste(sapply(scenario, get_pegnum), collapse="\t")), justify="left")
  )
  rownames(all.results) <- NULL
  all.results <- all.results[order(all.results$Difference),]
  write.table(all.results, file = resultsFilename, sep = "\t", quote = FALSE, row.names = FALSE)
}
