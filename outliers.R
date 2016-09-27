source("Aij_validation.R")
source("craig_utility_funcs.R")

wrContexts <- lapplyWithNames(c("None", "Pct 0.5", "Pct 1", "Pct 2", "SD 4", "SD 5", "Wins 4", "Wins 5"), function(outlierHandling) {
  get.wrapped.aij.context.from.ds("Hope Ecoli", outlierHandling=outlierHandling)
})
valContexts <- lapply(wrContexts, function(context) {
  get.aij.validation.context(context, gold="mattE", aijList=list(
    coinflip = coinflip(context),
    MT = MT(context),
    TT = TT(context),
    RB = RB(context),
    `UniMM(12)` = UniMM(context, 12),
    `MultiMM(12)` = MultiMM(context, 12),
    `UniMM(0)` = UniMM(context, 0),
    `MultiMM(0)` = MultiMM(context, 0),
    `UniMM(NBF)` = UniMM_NBF(context),
    `MultiMM(NBF)` = MultiMM_NBF(context)
  ))
})

results <- lapply(valContexts, function(context) {
  gold_calls <- get.gold.calls(context)
  active.filter <- gold_calls & !is.na(gold_calls)
  inactive.filter <- !gold_calls & !is.na(gold_calls)
  always.on.filter <- matrix(apply(active.filter, 1, all), nr=nrow(active.filter), nc=ncol(active.filter))
  always.off.filter <- matrix(apply(inactive.filter, 1, all), nr=nrow(inactive.filter), nc=ncol(inactive.filter))
  myFilters <- list(
    overall = TRUE,
    always_on = always.on.filter,
    always_off = always.off.filter,
    nonconstant = !always.on.filter & !always.off.filter & !is.na(gold_calls),  
    constant = always.on.filter | always.off.filter & !is.na(gold_calls)
  )
  lapply(myFilters, function(filter) all.mean.alignments(context, type="sq", filter=filter))
})

invisible(lapply(names(results[[1]]), function(filterName) {
  cat("\n", filterName, "\n\n")
  results.matrix <- t(sapply(results, `[[`, filterName))
  print(results.matrix)
  png(paste0("barplot_",filterName,".png"), width=800, height=500)
  par(mar=c(10,4,5,2)+0.1)
  barplot(results.matrix, legend.text=TRUE, beside=TRUE, col=rainbow(nrow(results.matrix)), 
          ylim=c(0.2, 0.8), xpd=FALSE, las=2, ylab="Mean square alignment", args.legend=list(x="topleft"),
          main=if(filterName=="overall") "Mean square alignment by outlier handling" else paste0("Mean square alignment by outlier handling\n",filterName," genes (as determined by FVA)"))
  dev.off()
  return(invisible(results.matrix))
}))
