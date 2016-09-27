source("Aij_validation.R")

run.standard.tests <- function(data.source, gold.source, outlierHandling="None") {
  # data.source: See notes at top of Aij_generation_input_data.R
  # gold.source: One of "matt", "iMRM", or "sim", indicating which calls should be used as the gold standard.
  #   "sim" should always be used if the expression data was simulated. 
  # outlierHandling: Passed directly to get.wrapped.aij.context.from.ds(); see notes there (in wrAijContext.R)
  
  cat("\n--- Beginning standard.tests for", data.source, if(outlierHandling=="None") "" else outlierHandling, "---\n")
  wrAijContext <- get.wrapped.aij.context.from.ds(data.source, outlierHandling=outlierHandling)
  cat("Got wrAijContext\n")

  aijList <- list(
    coinflip = coinflip(wrAijContext),
    MT = MT(wrAijContext),
    TT = TT(wrAijContext),
    RB = RB(wrAijContext),
    `UniMM(12)` = UniMM(wrAijContext, 12),
    `MultiMM(12)` = MultiMM(wrAijContext, 12),
    `UniMM(0)` = UniMM(wrAijContext, 0),
    `MultiMM(0)` = MultiMM(wrAijContext, 0),
    `UniMM(NBF)` = UniMM_NBF(wrAijContext),
    `MultiMM(NBF)` = MultiMM_NBF(wrAijContext)
  )
  cat("Got all aijs\n")
  
  valContext <- get.aij.validation.context(wrAijContext, aijList, gold.source)
  cat("Got valContext\n")

  # Note this function does not run ALL of the validation tests available in Aij_validation.R; only a collection of them
  
  gold_calls <- get.gold.calls(valContext)
  
  active.filter <- gold_calls & !is.na(gold_calls)
  inactive.filter <- !gold_calls & !is.na(gold_calls)
  always.on.filter <- matrix(apply(active.filter, 1, all), nr=nrow(active.filter), nc=ncol(active.filter))
  always.off.filter <- matrix(apply(inactive.filter, 1, all), nr=nrow(inactive.filter), nc=ncol(inactive.filter))
  nonconstant.filter <- !always.on.filter & !always.off.filter & !is.na(gold_calls)  # Either nonconstant, or an operon-mate has NA gold call
  
  onecomp.rownames <- intersect(get.valGenenames(valContext), get.1comp.genenames(valContext, 0))  # because get.1comp.genenames() includes genes from the original set, which are not necessarily in the validation set
  onecomp.filter <- matrix(FALSE, nr=nrow(gold_calls), nc=ncol(gold_calls), dimnames=dimnames(gold_calls))
  onecomp.filter[onecomp.rownames,] <- TRUE
  
  # operon_calls: Like gold_calls, but NA if the calls were inconsistent for the operon,experiment pair that gene belongs to
  # This code currently non-operational if gold_calls comes in already having NAs
  #if(any(is.na(gold_calls))) stop("NAs found in gold_calls")
  #operon_calls <- gold_calls
  #for(op in get.valOperons(valContext)) {
  #  for(colIndex in 1:ncol(gold_calls)) {
  #    op.ex.calls <- gold_calls[op, colIndex]
  #    if(any(op.ex.calls) && any(!op.ex.calls)) {
  #      operon_calls[op, colIndex] <- NA
  #    }
  #  }
  #}
  
  # These filters are dependent on operon_calls above
  #operon.on.filter <- operon_calls
  #operon.on.filter[is.na(operon.on.filter)] <- FALSE
  #operon.off.filter <- !operon_calls
  #operon.off.filter[is.na(operon.off.filter)] <- FALSE
  #operon.always.on.filter <- matrix(apply(operon.on.filter, 1, all), nr=nrow(operon.on.filter), nc=ncol(operon.on.filter))
  #operon.always.off.filter <- matrix(apply(operon.off.filter, 1, all), nr=nrow(operon.off.filter), nc=ncol(operon.off.filter))
  #operon.inconsistent.filter <- is.na(operon_calls)
  
  genes.from.onecomp.operons <- intersect(unlist(get.operons(valContext)[get.1comp.operons(valContext,0)]), get.valGenenames(valContext))  # NOT get.valOperons - because get.1comp.operons() gives indices in get.operons(valContext)
  operon.onecomp.filter <- matrix(FALSE, nr=nrow(gold_calls), nc=ncol(gold_calls), dimnames=dimnames(gold_calls))
  operon.onecomp.filter[genes.from.onecomp.operons,] <- TRUE
  
  all.boxplots <- function(type) {
    boxplot.deviations(valContext, byWhat="experiment", type=type)
    boxplot.deviations(valContext, byWhat="experiment", filter=active.filter, filterName="active", type=type)
    boxplot.deviations(valContext, byWhat="experiment", filter=inactive.filter, filterName="inactive", type=type)
    boxplot.deviations(valContext, byWhat="experiment", filter=always.on.filter, filterName="always_on", type=type)
    boxplot.deviations(valContext, byWhat="experiment", filter=always.off.filter, filterName="always_off", type=type)
    boxplot.deviations(valContext, byWhat="experiment", filter=nonconstant.filter, filterName="nonconstant", type=type)
    boxplot.deviations(valContext, byWhat="experiment", filter=onecomp.filter, filterName="onecomp", type=type)
    cat("Done generating", type, "boxplots\n")
  }
  
  #all.boxplots("abs")
  #all.boxplots("sq")
  #all.boxplots("log")
  
  myFilters <- list(
    overall = TRUE,
    #always_on = always.on.filter[,1],
    #always_off = always.off.filter[,1],
    nonconstant = nonconstant.filter[,1],
    constant = always.on.filter[,1] | always.off.filter[,1]
  )
  N <- length(get.expnames(valContext))
  sapply(1:length(myFilters), function(filterNum) hist.C.values.by.gene(filter=myFilters[[filterNum]], filterName=names(myFilters)[filterNum], valContext=valContext, A=-5/N, B=0, exp=FALSE))
  sapply(1:length(myFilters), function(filterNum) hist.C.values.by.operon(filter=myFilters[[filterNum]], filterName=names(myFilters)[filterNum], valContext=valContext, A=bquote(-exp(2*(3+2*p)/.(N))/p), B=quote(-1/p), exp=FALSE))
  sapply(1:length(myFilters), function(filterNum) hist.old.equiv.C.by.gene(filter=myFilters[[filterNum]], filterName=names(myFilters)[filterNum], valContext=valContext, BIC.bonus=0))
  sapply(1:length(myFilters), function(filterNum) hist.old.equiv.C.by.operon(filter=myFilters[[filterNum]], filterName=names(myFilters)[filterNum], valContext=valContext, BIC.bonus=0))
  sapply(1:length(myFilters), function(filterNum) hist.old.equiv.C.by.gene(filter=myFilters[[filterNum]], filterName=names(myFilters)[filterNum], valContext=valContext, BIC.bonus=12))
  sapply(1:length(myFilters), function(filterNum) hist.old.equiv.C.by.operon(filter=myFilters[[filterNum]], filterName=names(myFilters)[filterNum], valContext=valContext, BIC.bonus=12))
  
  sink(paste0("results/",get.valSourceStr(valContext),"/Cvalue_alignment.txt"))
  C.mean.sq.alignments <- sapply(names(myFilters), function(filterName) {
    filter <- myFilters[[filterName]]
    c(
      `UniCP(12)` = old.equiv.C.mean.square.alignment(valContext, 12, mult=F, filter=filter, filterName=filterName),
      `UniCP(0)` = old.equiv.C.mean.square.alignment(valContext, 0, mult=F, filter=filter, filterName=filterName),
      `UniCP(NBF)` = C.value.mean.square.alignment(valContext, A=-5/N, B=0, exp=FALSE, mult=F, filter=filter, filterName=filterName),
      `MultiCP(12)` = old.equiv.C.mean.square.alignment(valContext, 12, mult=T, filter=filter, filterName=filterName),
      `MultiCP(0)` = old.equiv.C.mean.square.alignment(valContext, 0, mult=T, filter=filter, filterName=filterName),
      `MultiCP(NBF)` = C.value.mean.square.alignment(valContext, A=bquote(-exp(2*(3+2*p)/.(N))/p), B=quote(-1/p), exp=FALSE, mult=T, filter=filter, filterName=filterName)
    )
  })
  print(C.mean.sq.alignments)
  for(filter in names(myFilters)) {
    png(paste0("results/", get.valSourceStr(valContext),"/C_mean_sq_alignments_",filter,".png"))
    par(mar=c(7,3,2,2)+0.1)
    barplot(C.mean.sq.alignments[,filter], ylab="Mean square alignment", ylim=c(0.4,1), xpd=FALSE, las=2, args.legend=list(x="topleft"))
    dev.off()
  }
  sink()
  
  sink(paste0("results/",get.valSourceStr(valContext),"/conftables.txt"))
  confidence.tables(valContext)
  confidence.tables(valContext, active.filter, "ACTIVE")
  confidence.tables(valContext, inactive.filter, "INACTIVE")
  confidence.tables(valContext, always.on.filter, "ALWAYS-ON")
  confidence.tables(valContext, always.off.filter, "ALWAYS-OFF")
  confidence.tables(valContext, nonconstant.filter, "NONCONSTANT")
  confidence.tables(valContext, onecomp.filter, "ONE-COMP(0)_GENES")
  confidence.tables(valContext, !onecomp.filter, "TWO-COMP(0)_GENES")
  confidence.tables(valContext, operon.onecomp.filter, "ONE-COMP(0)_OPERONS")
  confidence.tables(valContext, !operon.onecomp.filter, "TWO-COMP(0)_OPERONS")
  sink()
  cat("Done generating confidence tables\n")

  #sink(paste0("results/",get.valSourceStr(valContext),"/constables.txt"))
  #consistency.tables(valContext)
  #consistency.tables(valContext, operon.on.filter, "CONSISTENTLY ACTIVE")
  #consistency.tables(valContext, operon.off.filter, "CONSISTENTLY INACTIVE")
  #consistency.tables(valContext, operon.always.on.filter, "CONSISTENTLY ALWAYS-ON")
  #consistency.tables(valContext, operon.always.off.filter, "CONSISTENTLY ALWAYS-OFF")
  #consistency.tables(valContext, operon.inconsistent.filter, "INCONSISTENT")
  #consistency.tables(valContext, onecomp.filter, "ONE-COMP(0)_GENES")
  #consistency.tables(valContext, !onecomp.filter, "TWO-COMP(0)_GENES")
  #consistency.tables(valContext, operon.onecomp.filter, "ONE-COMP(0)_OPERONS")
  #consistency.tables(valContext, !operon.onecomp.filter, "TWO-COMP(0)_OPERONS")
  #sink()
  
  #totalValContext <- get.aij.validation.context(wrAijContext, aijList, "all-TRUE")
  #sink(paste0("results/",get.valSourceStr(valContext),"/total_constables.txt"))
  #consistency.tables(totalValContext, withCorrectness=FALSE)
  #consistency.tables(totalValContext, operon.on.filter, "CONSISTENTLY ACTIVE", withCorrectness=FALSE)
  #consistency.tables(totalValContext, operon.off.filter, "CONSISTENTLY INACTIVE", withCorrectness=FALSE)
  #consistency.tables(totalValContext, operon.always.on.filter, "CONSISTENTLY ALWAYS-ON", withCorrectness=FALSE)
  #consistency.tables(totalValContext, operon.always.off.filter, "CONSISTENTLY ALWAYS-OFF", withCorrectness=FALSE)
  #consistency.tables(totalValContext, operon.inconsistent.filter, "INCONSISTENT", withCorrectness=FALSE)
  #consistency.tables(totalValContext, onecomp.filter, "ONE-COMP(0)_GENES", withCorrectness=FALSE)
  #consistency.tables(totalValContext, !onecomp.filter, "TWO-COMP(0)_GENES", withCorrectness=FALSE)
  #consistency.tables(totalValContext, operon.onecomp.filter, "ONE-COMP(0)_OPERONS", withCorrectness=FALSE)
  #consistency.tables(totalValContext, !operon.onecomp.filter, "TWO-COMP(0)_OPERONS", withCorrectness=FALSE)
  #sink()
  
  cat("Done generating consistency tables\n")
  
  # Return the valContext in case the calling code has further use for it
  return(invisible(valContext))
  
}
