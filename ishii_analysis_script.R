library(XML)
source("iMRM gene expression and kos.R")
source("wrAijContext.R")

model <- xmlInternalTreeParse("inputs/Carrera_iJO1366.sbml")
gene.associations <- get.gene.associations(model)

#################################################################
R_unsure_0.5 <- R  # defined in "iMRM gene expression and kos.R"; gives 0.5 when unsure of the right value
remove(R)  # can't accidentally refer to "R" now; must explicitly choose either R_unsure_0.5 or R_unsure_NA
R_unsure_NA <- function(snippet, aijs) {  # almost identical to R_unsure_0.5 except for handling of unknowns, which this returns as NAs
  words <- strsplit(snippet, "\\s+")[[1]]  # get a list of 'words' in snippet (separated by whitespace)
  words <- words[which(words!="")]   # remove any empty-string 'words'
  if(length(words)==0) return(NA)  # snippet was the empty string or just whitespace
  kind <- which.max(c("and" %in% words, "or" %in% words, TRUE))  
  # note that a snippet cannot have both 'and' and 'or'; this would be ambiguous
  # so kind==1 means this snippet contains 'and's; kind==2, 'or's; and kind==3, neither
  numbers <- unlist(sapply(words, function(word) {  # the unlist removes zero-length vectors
    if(word=="and" || word=="or") return(numeric(0))  # remove 'and' and 'or' words
    if(word=="Unknown") return(NA)  
    if(is.na(word) || word == "NA") return(NA)  # If NA was passed (either an actual NA or the string "NA"), return NA again
    numericword <- suppressWarnings(as.numeric(word))
    if(!is.na(numericword)) return(numericword)  # a number: R(s) for some other s
    return(aijs[word])  # by process of elimination, word must be a gene name
  }, USE.NAMES=FALSE))
  if(kind==1) {
    not.nas <- which(!is.na(numbers))
    if (length(not.nas) == 0)  return(NA)  # NA 'and' NA 'and' NA 'and'... = NA
    return(min(numbers[not.nas]))  # For any 'and' expression involving NA's, take the min of the non-NA values
  }
  else if(kind==2) {
    not.nas <- which(!is.na(numbers))
    if (length(not.nas) == 0)  return(NA)  # NA 'or' NA 'or' NA 'or'... = NA
    return(max(numbers[not.nas]))  # For any 'or' expression involving NA's, take the max of the non-NA values
  }
  else if(kind==3) {  # should be only one entry in numbers if kind==3
    return(numbers)
  }
}
#################################################################

# Get a yes/no response. Returns TRUE if the user entered Yes (or anything starting with upper or lowercase Y), or FALSE otherwise.
readYN <- function(prompt) {
  if(substr(prompt, nchar(prompt), nchar(prompt)) != " ") prompt <- paste0(prompt, " ")  # add trailing space if there isn't one
  string <- readline(paste0(prompt, "(y/n) "))
  toupper(substr(string, 1, 1)) == 'Y'
}

GPR.genebased.to.rxnbased.single.experiment <- function(values, reactions, evalFunc) {
  # values: vector of numeric values, with names being gene names
  # reactions: vector of names of reactions for which a rxn-based value is desired
  # evalFunc: passed to parseAssocStr; for now it's either R_unsure_NA or R_unsure_0.5
  # Returns a vector of rxn-based values
  sapply(gene.associations[reactions], parseAssocStr, evalFunc, aijs=values)
}

GPR.genebased.to.rxnbased.multiple.experiments <- function(values, reactions, evalFunc) {
  # values: matrix of numeric values, with rownames being gene names; every column will be converted by GPR
  # reactions: vector of names of reactions for which a rxn-based value is desired
  # evalFunc: passed to parseAssocStr; for now it's either R_unsure_NA or R_unsure_0.5
  # Returns a matrix of rxn-based values, with one row per reaction and columns corresponding to columns of values
  apply(values, 2, GPR.genebased.to.rxnbased.single.experiment, reactions, evalFunc)
}

plot.with.GPR.single.experiment <- function(values_by_gene, values_by_rxn, evalFunc, title, ...) {
  # values_by_gene: vector of numeric values, with names being gene names
  # values_by_rxn: vector of numeric values, with names being reaction names
  # evalFunc: passed to parseAssocStr; for now it's either R_unsure_NA or R_unsure_0.5
  # title: title to print on the plot
  # ... : more arguments to plot()
  # This function will first convert values_by_gene into reaction-based values
  #   using the rij GPR rules from the iMRM
  # then it will plot those vs values_by_rxn
  
  if(is.null(names(values_by_gene))) stop("values_by_gene must be named")
  if(is.null(names(values_by_rxn))) stop("values_by_rxn must be named")
  
  gene.values.by.rxn <- GPR.genebased.to.rxnbased.single.experiment(values_by_gene, names(values_by_rxn), evalFunc)
  correlation <- cor(gene.values.by.rxn, values_by_rxn, use="complete.obs")
  plot(gene.values.by.rxn, values_by_rxn, main=paste0(title,"\n(correlation coefficient ", round(correlation,3), ")"), ...)
  return(invisible(correlation))
}

plot.with.GPR.multiple.experiments <- function(values_by_gene, values_by_rxn, evalFunc, title, ..., makePlot=TRUE) {
  # values_by_gene: matrix of numeric values, with rownames being gene names and colnames being experiment names
  # values_by_rxn: matrix of numeric values, with rownames being reaction names and colnames being experiment names 
  # evalFunc: passed to parseAssocStr; for now it's either R_unsure_NA or R_unsure_0.5
  # title: title to print on the plot (ignored if makePlot=FALSE)
  # ... : more arguments to plot()
  # makePlot: whether to plot the results (or just return the correlation)
  
  if(is.null(rownames(values_by_gene))) stop("values_by_gene must have rownames")
  if(is.null(colnames(values_by_gene))) stop("values_by_gene must have colnames")
  if(is.null(rownames(values_by_rxn))) stop("values_by_rxn must have rownames")
  if(is.null(colnames(values_by_rxn))) stop("values_by_rxn must have colnames")
  
  values_by_rxn <- values_by_rxn[,colnames(values_by_gene)]  # make sure columns are in the same order
  gene.values.by.rxn <- GPR.genebased.to.rxnbased.multiple.experiments(values_by_gene, rownames(values_by_rxn), evalFunc)
  values_by_rxn <- as.vector(values_by_rxn)  # all values
  gene.values.by.rxn <- as.vector(gene.values.by.rxn)  # all values, in same order as values_by_rxn
  correlation <- cor(gene.values.by.rxn, values_by_rxn, use="complete.obs")
  if(makePlot) plot(gene.values.by.rxn, values_by_rxn, main=paste0(title,"\n(correlation coefficient ", round(correlation,3), ")"), ...)
  return(invisible(correlation))
}

test_current_aijs <- function() {
  fluxes <- as.matrix(read.csv("inputs/fluxomics_mmol_per_gDW_per_h.csv", header=T, row.names=1))
  ishiicontext <- get.wrapped.aij.context.from.ds("Ishii Ecoli")
  eijs <- get.expression.data.outliersNA(ishiicontext)
  aijs <- MultiMM(ishiicontext, 12)
  #for(experiment in colnames(fluxes)) {
  #  plot.with.GPR.single.experiment(eijs[,experiment], abs(fluxes[,experiment]), R_unsure_NA,
  #                title=paste("Abs(fluxes) vs eijs for experiment", experiment), xlab="Eij", ylab="Abs(Flux)")
  #  if(!readYN("Type 'y' to continue, 'n' to stop")) break
  #  plot.with.GPR.single.experiment(aijs[,experiment], abs(fluxes[,experiment]), R_unsure_0.5,
  #                title=paste("Abs(fluxes) vs aijs for experiment", experiment), xlab="Aij", ylab="Abs(Flux)")
  #  if(!readYN("Type 'y' to continue, 'n' to stop")) break
  #}
  plot.with.GPR.multiple.experiments(eijs, abs(fluxes), R_unsure_NA,
                                     title=paste("Abs(fluxes) vs eijs"), xlab="Eij", ylab="Abs(Flux)")
  readline()
  plot.with.GPR.multiple.experiments(aijs, abs(fluxes), R_unsure_NA,
                                     title=paste("Abs(fluxes) vs aijs"), xlab="Aij", ylab="Abs(Flux)")
  readline()
  plot.with.GPR.multiple.experiments(sqrt(eijs), sqrt(abs(fluxes)), R_unsure_NA,
                                     title=paste("sqrt(Abs(fluxes)) vs sqrt(eijs)"), xlab="sqrt(Eij)", ylab="sqrt(Abs(Flux))")
  readline()
  plot.with.GPR.multiple.experiments(aijs, sqrt(abs(fluxes)), R_unsure_NA,
                                     title=paste("sqrt(Abs(fluxes)) vs aijs"), xlab="Aij", ylab="sqrt(Abs(Flux))")
}

parameter.test <- function(numGibbsSamplerRuns=1000, numBurnedRuns=100, id=100) {
  # id: an integer, used for keeping apart save files. Multiple instances of this function can run concurrently AS LONG AS they have distinct id's. 
  fluxes <- as.matrix(read.csv("inputs/fluxomics_mmol_per_gDW_per_h.csv", header=T, row.names=1))
  removeishiifiles <- function() file.remove(list.files(path="Rsave", pattern=paste0(".*ishii",id,".*\\.Rdata$"), full.names=TRUE))
  mu1 <- 0.005
  sw <- 1e-3
  aij.correlations <- sapply(c(3e-6,1e-5,3e-5,1e-4), function(mu0.var) {
    sapply(seq(from=0.0005,to=0.004,by=0.0005), function(mu0) {
      removeishiifiles()
      cat("\n\nStarting run with mu0.var =", mu0.var, "and mu0 =", mu0, "\n\n\n")
      ishiicontext <- get.wrapped.aij.context.from.ds("Ishii Ecoli", id=id, mu0.mean=mu0, mu1.mean=mu1, mu0.var=mu0.var, wishart_scale=sw, numGibbsSamplerRuns=numGibbsSamplerRuns, numBurnedRuns=numBurnedRuns)
      aij.correlation <- plot.with.GPR.multiple.experiments(MultiMM(ishiicontext, 12), abs(fluxes), R_unsure_0.5, makePlot=FALSE)
      return(aij.correlation)
    })
  })
  removeishiifiles()
  save(aij.correlations, file=paste0("aij.correlations",id,".Rdata"))
  return(aij.correlations)
}

parameter.test.2 <- function(numGibbsSamplerRuns=1000, numBurnedRuns=100, id=1) {
  # id: an integer, used for keeping apart save files. Multiple instances of this function can run concurrently AS LONG AS they have distinct id's. 
  fluxes <- as.matrix(read.csv("inputs/fluxomics_mmol_per_gDW_per_h.csv", header=T, row.names=1))
  removeishiifiles <- function() file.remove(list.files(path="Rsave", pattern=paste0(".*ishii",id,".*\\.Rdata$"), full.names=TRUE))
  
  # Values to cycle through for each parameter
  mu0 <- c(0.002, 0.005)
  mu1.minus.mu0 <- c(0.001, 0.003, 0.006)
  mu0.var <- c(1e-7, 1e-5)
  sw <- c(3e-7, 3e-5, 1e-3)
  
  results <- do.call(rbind, lapply(mu0, function(mu0) {
    do.call(rbind, lapply(sw, function(sw) {
      do.call(rbind, lapply(mu0.var, function(mu0.var) {
        do.call(rbind, lapply(mu1.minus.mu0+mu0, function(mu1) {
          removeishiifiles()
          cat("\n\nStarting run with mu0 =", mu0, "; mu1 =", mu1, "; mu0.var =", mu0.var, "and sw =", sw, "\n\n\n")
          ishiicontext <- get.wrapped.aij.context.from.ds("Ishii Ecoli", id=id, mu0.mean=mu0, mu1.mean=mu1, mu0.var=mu0.var, wishart_scale=sw, numGibbsSamplerRuns=numGibbsSamplerRuns, numBurnedRuns=numBurnedRuns)
          aijs <- MultiMM(ishiicontext, 12)
          return(c(
            mu0 = mu0,
            mu1 = mu1,
            mu0.var = mu0.var,
            sw = sw,
            get.stats.vs.fluxes(aijs, fluxes)
          ))
        }))
      }))
    }))
  }))
  removeishiifiles()
  save(results, file=paste0("results",id,".Rdata"))
  return(results)
}

get.stats.vs.fluxes <- function(aijs, fluxes) {
  fluxes <- fluxes[,colnames(aijs)]  # make sure columns are in the same order
  aijs.by.rxn <- GPR.genebased.to.rxnbased.multiple.experiments(aijs, rownames(fluxes), R_unsure_NA)
  fluxes <- as.vector(abs(fluxes))  # all values; also take abs
  aijs.by.rxn <- as.vector(aijs.by.rxn)  # all values, in same order as fluxes
  cijs <- rep(NA, length(aijs.by.rxn))
  cijs[!is.na(aijs.by.rxn) & aijs.by.rxn>0.5 & fluxes>0] <- 1
  cijs[!is.na(aijs.by.rxn) & aijs.by.rxn>0.5 & fluxes==0] <- 0
  cijs[!is.na(aijs.by.rxn) & aijs.by.rxn<0.5 & fluxes>0] <- 0
  cijs[!is.na(aijs.by.rxn) & aijs.by.rxn<0.5 & fluxes==0] <- 1
  cijs[!is.na(aijs.by.rxn) & aijs.by.rxn==0.5] <- 0.5
  highconf <- !is.na(aijs.by.rxn) & (aijs.by.rxn >= 0.8 | aijs.by.rxn <= 0.2)
  medconf <- !is.na(aijs.by.rxn) & (aijs.by.rxn >= 0.6 | aijs.by.rxn <= 0.4) & !highconf
  lowconf <- !is.na(aijs.by.rxn) & !highconf & !medconf
  return(c(
    informativeness = mean(abs(aijs-0.5), na.rm=T),
    correlation.vs.fluxes = cor(aijs.by.rxn, fluxes, use="complete.obs"),
    correlation.vs.sqrt.fluxes = cor(aijs.by.rxn, sqrt(fluxes), use="complete.obs"),
    num.datapoints = sum(!is.na(aijs.by.rxn)),
    mean.cij = mean(cijs, na.rm=T),
    num.highconf = sum(highconf),
    mean.cij.highconf = mean(cijs[highconf]),
    num.medconf = sum(medconf),
    mean.cij.medconf = mean(cijs[medconf]),
    num.lowconf = sum(lowconf),
    mean.cij.lowconf = mean(cijs[lowconf])
  ))
}

get.lm.stats <- function(responsevec, explanatoryvec, explanatoryvec2=NULL) {
  if(is.null(explanatoryvec2)) fit <- lm(responsevec~explanatoryvec)
  else fit <- lm(responsevec~explanatoryvec+explanatoryvec2)
  coeffs <- summary(fit)$coefficients
  confint <- confint(fit)
  pval1 <- coeffs["explanatoryvec","Pr(>|t|)"]
  pval2 = if(is.null(explanatoryvec2)) NULL else coeffs["explanatoryvec2","Pr(>|t|)"]
  cat(sep="", coeffs["explanatoryvec","Estimate"], " (", confint["explanatoryvec",1], " - ", confint["explanatoryvec",2], ") ", if(pval1<0.001) "***" else if(pval1<0.01) "**" else if(pval1<0.05) "*" else "")
  if(!is.null(explanatoryvec2)) cat(sep="", "\n", coeffs["explanatoryvec2","Estimate"], " (", confint["explanatoryvec2",1], " - ", confint["explanatoryvec2",2], ") ", if(pval2<0.001) "***" else if(pval2<0.01) "**" else if(pval2<0.05) "*" else "")
}

standardize <- function(vec) {
  vec <- vec - mean(vec, na.rm=T)
  vec <- vec / sd(vec, na.rm=T)
  return(vec)
}
