source("craig_utility_funcs.R")

# This file contains a variety of useful functions for understanding what's going on in the aij algorithms, and delving into case studies
# Functions which take an argument "context" accept an aijContext or any of its subclasses

# Given pegnums (integers) get the Hope genename, and vice versa
hope.gene <- function(hope.pegnums) paste0("fig|83333.1.peg.", hope.pegnums)
hope.pegnum <- function(hope.genes) withoutNames(sapply(hope.genes, function(gene) if(grepl("peg", gene)) regmatches(gene, regexpr("[[:digit:]]+$", gene)) else NA))

# Convert fig id's to b-numbers
load("overallmap.Rdata")
fig_to_b <- function(genevec) overallmap$b[match(genevec, overallmap$fig)]
b_to_fig <- function(genevec) overallmap$fig[match(genevec, overallmap$b)]
com_to_fig <- function(genevec) overallmap$fig[match(genevec, overallmap$common.x)] 

# Get the index in get.operons(context) [NOT get.valOperons(context)] of the operon containing a given gene
get.operon.num <- function(context, gene) which(sapply(get.operons(context), function(op) gene %in% op))

# Basic histogram of expression data
hist.exp <- function(context, gene, title=paste("Expression data for",gene), ...) {
  # ... : Additional arguments to hist()
  hist(get.expression.data.outliersNA(context)[gene,],breaks=25,xlab="eij",main=title,...)
}

# Histogram of expression data, noting the Bayesian-1.0-style component determination
hist.exp.Bayesian1 <- function(context, gene, BIC.bonus, multivariate=TRUE) {
  # multivariate: If TRUE, MultiMM's determination, else UniMM's
  hist.exp(context, gene)
  one.comp <- 
    if(multivariate) get.operon.num(context, gene) %in% get.1comp.operons(context, BIC.bonus)
    else gene %in% get.1comp.genenames(context, BIC.bonus)
  one.comp.str <- 
    if(one.comp) paste0("1-component for Bayesian 1.0 (", if(multivariate) "MultiMM" else "UniMM", " bonus ", BIC.bonus, ")")
    else paste0("2-component for Bayesian 1.0 (", if(multivariate) "MultiMM" else "UniMM", " bonus ", BIC.bonus, ")")
  mtext(one.comp.str, 3, line=0, cex=1.25, font=1)
}

# Histogram of expression data, noting the Bayesian-2.0 "Cvalue"
hist.exp.Bayesian2 <- function(context, gene, multivariate=TRUE) {
  # multivariate: If TRUE, operon-level Cvalue, else the gene-level Cvalue
  hist.exp(context, gene)
  N <- length(get.expnames(context))
  Cvalue <- 
    if(multivariate) get.Cvalues.by.operon(context, A=bquote(-exp(2*(3+2*p)/.(N))/p), B=quote(-1/p), exponential=FALSE)[get.operon.num(context, gene)]
    else get.Cvalues.by.gene(context, A=-5/N, B=0, exponential=FALSE)[gene]
  mtext(paste(if(multivariate) "Multivariate" else "Univariate", "C =", Cvalue), 3, line=0, cex=1.25, font=1)
}

# Histogram of expression data, with the 2-component fitted distribution overlaid
hist.exp.with.overlay <- function(context, gene, multivariate=TRUE, title=paste("Expression data for",gene)) {
  # multivariate: If TRUE, the multivariate (MultiMM) 2-component fit is used, else the univariate (UniMM) one
  hist.exp(context, gene, title=title, prob=T)
  dnormOff <- 
    if(multivariate) function(x) (1-get.multi.pi(context)[gene])*dnorm(x, mean=get.multi.mu0(context)[gene], sd=sqrt(get.multi.var(context)[[get.operon.num(context, gene)]][gene,gene]))
    else function(x) (1-get.uni.pi(context)[gene])*dnorm(x, mean=get.uni.mu0(context)[gene], sd=sqrt(get.uni.var(context)[gene]))
  dnormOn <-
    if(multivariate) function(x) get.multi.pi(context)[gene]*dnorm(x, mean=get.multi.mu1(context)[gene], sd=sqrt(get.multi.var(context)[[get.operon.num(context, gene)]][gene,gene]))
    else function(x) get.uni.pi(context)[gene]*dnorm(x, mean=get.uni.mu1(context)[gene], sd=sqrt(get.uni.var(context)[gene]))
  curve(dnormOff, add=T, col="darkblue", lwd=3)
  curve(dnormOn, add=T, col="darkred", lwd=3)
  legend("topright", c("Inactive fit", "Active fit"), fill=c("darkblue", "darkred"))
}

# Basic gene vs. gene plot of expression data
# Returns the correlation between the two genes' expression data (invisibly)
exp.vs.exp <- function(context, gene1, gene2, ...) {
  # ... : Additional arguments to plot()
  expdata <- get.expression.data.outliersNA(context)
  plot(expdata[gene1,], expdata[gene2,], xlab=gene1, ylab=gene2, ...)
  return(invisible(cor(expdata[gene1,], expdata[gene2,], use="pairwise.complete.obs")))
}

# Plot of a gene's eijs vs its aijs
exp.vs.aij <- function(context, aijs, gene, ...) {
  # ... : Additional arguments to plot()
  plot(get.expression.data.outliersNA(context)[gene,], aijs[gene,], xlab="eij", ylab="aij", ...)
}

# Histogram of "C-values" (all either 0 or 1) produced by the Bayesian 1.0 Classification Procedures; returns said Cvalues (invisibly)
hist.Cvalues.Bayesian1 <- function(context, BIC.bonus, multivariate=TRUE, title=paste(if(multivariate) "Multivariate" else "Univariate", "Bayesian1 C-values with bonus", BIC.bonus)) {
  one.comp.genes <- if(multivariate) unlist(get.operons(context)[get.1comp.operons(context, BIC.bonus)])
                    else get.1comp.genenames(context, BIC.bonus)
  Cvalues <- as.integer(!(get.genenames(context) %in% one.comp.genes))
  names(Cvalues) <- get.genenames(context)
  hist(Cvalues, breaks=seq(from=0,to=1,by=0.1), ylim=c(0,length(get.genenames(context))), col="orange", main=title, xlab="C-value")
  return(invisible(Cvalues))
}

# Histogram of C-values produced by the Bayesian 2.0 Classification Procedures; returns said Cvalues (invisibly)
hist.Cvalues.Bayesian2 <- function(context, multivariate=TRUE, title=paste(if(multivariate) "Multivariate" else "Univariate", "Bayesian2 C-values")) {
  Cvalues <- if(multivariate) get.operon.Cvalues.by.gene(context, A=bquote(-exp(2*(3+2*p)/.(N))/p), B=quote(-1/p), exponential=FALSE) 
             else get.Cvalues.by.gene(context, A=-5/N, B=0, exponential=FALSE)
  hist(Cvalues, breaks=seq(from=0,to=1,by=0.1), ylim=c(0,length(get.genenames(context))), col="orange", main=title, xlab="C-value")
  return(invisible(Cvalues))
}
