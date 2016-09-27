source("wrAijContext.R")
source("simulation_from_fits.R")

# Accepts any matrix, e.g. aijs, aij calls, gold calls, whatever
remove_CEL_gz <- function(mat) {
  ncharColnames <- nchar(colnames(mat))  # just so we don't have to compute it over and over
  if(all(substr(colnames(mat), ncharColnames-6, ncharColnames) == ".CEL.gz")) {
    colnames(mat) <- substr(colnames(mat), 1, ncharColnames-7)  # Remove the .CEL.gz from each colname
  }
  return(mat)
}

get.aij.validation.context <- function(wrAijContext, aijList, gold.source) {
  # wrAijContext: A wrapped aij context (see wrAijContext.R)
  # aijList: Named list. Each component of the list is a matrix of aijs (rows=genes, columns=experiments). 
  #   All matrices in the list must have been generated from the wrAijContext. 
  #   Names of aijList will be used in various outputs to label the results.
  # gold.source: One of "mattE", "mattEA", "iMRM", or "sim", indicating which calls should be used as the gold standard.
  #   "sim" should always be used if the expression data was simulated. 
  #   Another option is "all-TRUE", which simply gives an "on" call for every value in the expression data
  #     (Intended for tests (e.g. consistency tables with correctness=FALSE) which ignore the calls, but want
  #     to guarantee a non-NA call for every value in the expression data to avoid dropping data)
  
  # Create gold_calls, a matrix of booleans
  if(gold.source == "iMRM") {
    load("results/AllCalls.Rdata")
    gold_calls <- allcalls==1
    remove(allcalls)
  } else if (gold.source == "mattE") {
    load("Rsave/mattcalls_4-14-16.Rdata")
    calls <- as.matrix(calls)
    gold_calls <- matrix(FALSE, nr=nrow(calls), nc=ncol(calls), dimnames=dimnames(calls))  # Start by making everything FALSE
    gold_calls[grepl("Essential",calls)] <- TRUE  # "Essential" or any conflict containing "Essential" is considered 'on'
    gold_calls[calls==""] <- NA  # no call made for these cells
    remove(calls)
  } else if (gold.source == "mattEA") {
    load("Rsave/mattcalls_4-14-16.Rdata")
    calls <- as.matrix(calls)
    gold_calls <- matrix(FALSE, nr=nrow(calls), nc=ncol(calls), dimnames=dimnames(calls))  # Start by making everything FALSE
    gold_calls[grepl("Essential",calls)|grepl("Active",calls)] <- TRUE  # "Essential", "Active", or any conflict containing either of these is considered 'on'
    gold_calls[calls==""] <- NA  # no call made for these cells
    remove(calls)
  } else if (gold.source == "sim") {
    first.space <- regexpr(" ", get.data.source(wrAijContext))  # character position of the first space in get.data.source(wrAijContext)
    type.of.sim <- substr(get.data.source(wrAijContext), 1, first.space-1)
    if(type.of.sim == "Sim-Uniform") type <- "uniform"
    else if(type.of.sim == "Sim-Fitted") type <- "fitted"
    else if(type.of.sim == "Sim-Unfitted") type <- "unfitted"
    else stop(paste("Unrecognized data.source:", get.data.source(wrAijContext)))
    if(type=="unfitted") sim.data.source <- get.unfitted.sim.data.source(get.data.source(wrAijContext))
    else sim.data.source <- substr(get.data.source(wrAijContext), first.space+1, nchar(get.data.source(wrAijContext)))  # everything after the first space
    gold_calls <- get.sim.data(sim.data.source, type)$gold_calls
  } else if (gold.source == "all-TRUE") {
    expdata <- get.expression.data.outliersNA(wrAijContext)
    gold_calls <- matrix(TRUE, nr=nrow(expdata), nc=ncol(expdata), dimnames=dimnames(expdata))
    remove(expdata)
  } else {
    stop(paste("Unrecognized gold.source:", gold.source))
  }
  
  aijList <- lapply(aijList, remove_CEL_gz)
  if(!all(sapply(aijList, function(aijs) dim(aijs) == dim(aijList[[1]])))) stop("Sets of aijs have different dimensions")
  if(!all(sapply(aijList, function(aijs) rownames(aijs) == rownames(aijList[[1]])))) stop("Rownames do not match up between the sets of aijs")
  if(!all(sapply(aijList, function(aijs) colnames(aijs) == colnames(aijList[[1]])))) {
    temp <- sapply(aijList, function(aijs) any(colnames(aijs)!=colnames(aijList[[1]])))
    cat("Offending aijs:", which(temp), "\n")
    stop("Colnames do not match up between the sets of aijs")
  }
  sapply(1:length(aijList), function(index) if(any(is.na(aijList[[index]]))) stop(paste("NA(s) found in", names(aijList)[index]))) 
  
  gold_calls <- remove_CEL_gz(gold_calls)
  genenames <- intersect(rownames(aijList[[1]]), rownames(gold_calls))
  expnames <- intersect(colnames(aijList[[1]]), colnames(gold_calls))
  
  # Sort all into the same order, remove any genes / experiments we don't have both aijs and gold calls for
  aijList <- lapply(aijList, function(aijs) aijs[genenames, expnames])
  gold_calls <- gold_calls[genenames, expnames] 
  
  operons <- lapply(get.operons(wrAijContext), function(operon) operon[operon %in% genenames])  # A copy of get.operons(wrAijContext), but without any operon entries that are not in genenames
  operons <- operons[sapply(operons, length)!=0]  # remove any resulting empty operons
  
  if(grepl("Hope Ecoli", get.data.source(wrAijContext))) {
    source("analyze_scenarios_complete.R")
    scenarios <- get.scenarios("inputs/83333.1.scenarios")
  } else {
    scenarios <- list()
  }
  scenarios <- lapply(scenarios, function(scenario) scenario[scenario %in% genenames])  # remove any entries not in genenames
  scenarios <- scenarios[sapply(scenarios, length)!=0]  # remove any resultin empty scenarios
  
  source.str <- paste0(get.source.str(wrAijContext), "_", gold.source)
  if(!file.exists(paste0("results/",source.str)))  dir.create(paste0("results/",source.str))
  
  cat("Commencing aij validation with", length(genenames), "genes and", length(expnames), "experiments\n")
  returnme <- list(wrAijContext=wrAijContext, aijList=aijList, gold_calls=gold_calls, valOperons=operons, valSourceStr=source.str, scenarios=scenarios)
  class(returnme) <- c("valContext", "list")  # inherit from list
  return(returnme)
}

# In addition to the standard 'context' methods (also implemented by aijContext and wrAijContext), valContext implements the following methods:
get.aijList <- function(context) UseMethod("get.aijList")
get.gold.calls <- function(context) UseMethod("get.gold.calls")
get.valGenenames <- function(context) UseMethod("get.valGenenames")  # gene names for only the genes being used for validation
get.valExpnames <- function(context) UseMethod("get.valExpnames")  # experiment names for only the experiments being used for validation
get.valOperons <- function(context) UseMethod("get.valOperons")  # operons for only the genes being used for validation
get.valSourceStr <- function(context) UseMethod("get.valSourceStr")  # source.str which also indicates the gold.source
get.scenarios <- function(context) UseMethod("get.scenarios")

### 'Public' S3 methods for class 'valContext' ###
# Methods described in aijContext.R
get.genenames.valContext <- function(context) get.genenames(context$wrAijContext)
get.expnames.valContext <- function(context) get.expnames(context$wrAijContext)
get.expression.data.outliersNA.valContext <- function(context) get.expression.data.outliersNA(context$wrAijContext)
get.expression.data.outliersInf.valContext <- function(context) get.expression.data.outliersInf(context$wrAijContext)
get.operons.valContext <- function(context) get.operons(context$wrAijContext)
get.sourcefile.valContext <- function(context) "valContext.R"
get.uni.raw.aijs.valContext <- function(context) get.uni.raw.aijs(context$wrAijContext)
get.multi.raw.aijs.valContext <- function(context) get.multi.raw.aijs(context$wrAijContext)
get.uni.mu0.valContext <- function(context) get.uni.mu0(context$wrAijContext)
get.multi.mu0.valContext <- function(context) get.multi.mu0(context$wrAijContext)
get.uni.mu1.valContext <- function(context) get.uni.mu1(context$wrAijContext)
get.multi.mu1.valContext <- function(context) get.multi.mu1(context$wrAijContext)
get.uni.var.valContext <- function(context) get.uni.var(context$wrAijContext)
get.multi.var.valContext <- function(context) get.multi.var(context$wrAijContext)
get.uni.pi.valContext <- function(context) get.uni.pi(context$wrAijContext)
get.multi.pi.valContext <- function(context) get.multi.pi(context$wrAijContext)
get.BICdiffs.by.gene.valContext <- get.BICdiffs.by.gene.wrAijContext
get.BICdiffs.by.operon.valContext <- get.BICdiffs.by.operon.wrAijContext
get.1comp.genenames.valContext <- get.1comp.genenames.wrAijContext
get.1comp.operons.valContext <- get.1comp.operons.wrAijContext
get.1comp.fits.valContext <- get.1comp.fits.wrAijContext
get.bayesfactors.by.gene.valContext <- get.bayesfactors.by.gene.wrAijContext
get.bayesfactors.by.operon.valContext <- get.bayesfactors.by.operon.wrAijContext
get.Cvalues.by.gene.valContext <- get.Cvalues.by.gene.wrAijContext
get.Cvalues.by.operon.valContext <- get.Cvalues.by.operon.wrAijContext
get.Cvalues.by.operon.oldmethod.valContext <- get.Cvalues.by.operon.oldmethod.wrAijContext
get.operon.Cvalues.by.gene.valContext <- get.operon.Cvalues.by.gene.wrAijContext
get.operon.Cvalues.by.gene.oldmethod.valContext <- get.operon.Cvalues.by.gene.oldmethod.wrAijContext

# Methods described in Aij_generation.R
coinflip.valContext <- coinflip.wrAijContext
MT.valContext <- MT.wrAijContext 
TT.valContext <- TT.wrAijContext
RB.valContext <- RB.wrAijContext
UniMM.valContext <- UniMM.wrAijContext
MultiMM.valContext <- MultiMM.wrAijContext
ConfMM.valContext <- ConfMM.wrAijContext
UniMM_NBF.valContext <- UniMM_NBF.wrAijContext
MultiMM_NBF.valContext <- MultiMM_NBF.wrAijContext
ConfMM_NBF.valContext <- ConfMM_NBF.wrAijContext

# Methods described in wrAijContext.R
get.data.source.valContext <- function(context) get.data.source(context$wrAijContext)
get.source.str.valContext <- function(context) get.source.str(context$wrAijContext)

# Methods introduced in this script
get.aijList.valContext <- function(context) context$aijList
get.gold.calls.valContext <- function(context) context$gold_calls
get.valGenenames.valContext <- function(context) rownames(context$gold_calls)
get.valExpnames.valContext <- function(context) colnames(context$gold_calls)
get.valOperons.valContext <- function(context) context$valOperons
get.valSourceStr.valContext <- function(context) context$valSourceStr
get.scenarios.valContext <- function(context) context$scenarios
