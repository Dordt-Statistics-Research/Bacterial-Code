source("Aij_generation.R")
source("Aij_generation_input_data.R")
library(mclust)

# This file defines "wrapped aij context" or "wrAijContext", a subclass of aijContext which has all the same methods (and a few more)
# but most importantly, has caching functionality enabled, for HUGE performance gains. In addition to internal usefulness (such as
# UniMM/MultiMM using caching for get.1comp.fits(), a huge performance gain), this also has great benefits for you, the end user.
# Whenever you ask for something computationally intensive done with a wrAijContext, the result is automatically saved, so that 
# whenever you make the same request in the future, it simply loads the saved result and gives that to you. No need to worry about 
# saving any results yourself, ever - it's automatically done. For example, once you calculate the MultiMM aijs for a given data set with given 
# parameters, outlier handling, BIC bonus, etc, you never have to calculate them again - whenever you make that specific MultiMM call,
# it won't do the long calculation but instead just load the result you calculated before and give that to you. 
# If you actually want things recalculated, you can try the fromSaved=FALSE arguments, or (easier in my opinion) you can just 
# delete the saved files in Rsave/ that you want recalculated. They're all named according to specific naming conventions so that 
# R knows which one to load in the future in response to a specific call with a specific wrAijContext. 

# To take advantage of caching functionality, make sure you are always using wrAijContext objects (which you get from one of the three
# 'constructors' below) and not normal aijContext objects (from get.aij.context).  
# Then just pass the wrAijContext objects to all of the same methods you've been using, which for convenience are all 
# listed towards the bottom of this file as well. 

get.wrapped.aij.context <- function(data.source.name, source.str.name, expression.data, operons, outlierHandling="None", useSavedAijContext=TRUE, ...) {
  # data.source.name and source.str.name are the names to use for the wrapped aij context
  # expression.data, operons, and outlierHandling are as for get.aij.context(), see notes there (in aijContext.R)
  # useSavedAijContext: If TRUE, and if there is a saved aij context with the same source.str.name and outlierHandling, 
  #   that aij context will be used and the 'expression.data' and 'operons' arguments to this function will be ignored.
  # ... : more arguments to get.aij.context(); ignored if the saved aij context is used, and furthermore these arguments
  #   are not (presently) taken into account when determining whether the saved aij context is valid
  
  if(!file.exists("Rsave")) dir.create("Rsave")
  
  source.str <- paste0(source.str.name,ifelse(outlierHandling=="None","",paste0("_",gsub(" ","_",outlierHandling))))
  savefile <- paste0("Rsave/aijContext_",source.str,".Rdata")
  if(file.exists(savefile) && useSavedAijContext) {
    load(savefile)
  } else {
    cat("Generating (or regenerating) aij context for", data.source.name, "with outlier handling", outlierHandling, "\n")
    aijContext <- get.aij.context(expression.data, operons, outlierHandling, ...)
    save(aijContext, file=savefile)
  }
  
  return(get.wrapped.aij.context.from.aijContext(aijContext, data.source.name, source.str))
}

get.wrapped.aij.context.from.ds <- function(data.source, operon.source="Microbes Online", outlierHandling="None", id="", ...) {
  # data.source: See notes at top of Aij_generation_input_data.R
  # operon.source: Passed directly to get.operons.from.ds() (in Aij_generation_input_data.R); see notes there
  #   Note that operon.source can simply be omitted for some data.sources
  # outlierHandling: Passed directly to get.wrapped.aij.context(), see notes there
  # id: an integer, used for keeping apart save files. Unless you specifically need this, keep it at the default "". 
  # ... : more arguments to get.wrapped.aij.context()
  get.wrapped.aij.context(data.source, get.source.str.from.ds(data.source, id), get.expression.data.from.ds(data.source), get.operons.from.ds(data.source, operon.source), outlierHandling=outlierHandling, ...)
}

# Shortcut to wrap an aijContext if you already have the aijContext
get.wrapped.aij.context.from.aijContext <- function(aijContext, data.source.name, source.str) {
  # USE WITH CAUTION!
  # In this case source.str should already be adjusted for outlier handling, and it goes without saying that this should match the setting originally used when the aijContext was generated
  if(class(aijContext)[1]!="aijContext") stop("The object you passed as the aijContext argument to this function, doesn't appear to be an aijContext")
  returnme <- list(aijContext=aijContext, data.source=data.source.name, source.str=source.str)
  class(returnme) <- c("wrAijContext", "list")
  return(returnme)
}

# In addition to the standard 'context' methods (also implemented by aijContext), wrAijContext implements the following methods:
get.data.source <- function(context) UseMethod("get.data.source")
get.source.str <- function(context) UseMethod("get.source.str")

# Returns a method for class 'wrAijContext' for f, which caches its results
withCaching <- function(f, wrapname) {
  # f must be an S3 generic function with first argument a context
  return(function(context, ..., fromSaved=TRUE) {
    savefile <- paste0("Rsave/",wrapname,"_",paste0(list(...),collapse="_"),"_",get.source.str(context),".Rdata")
    if(fromSaved && file.exists(savefile)) {
      load(savefile)
      return(saved)
    } else {
      cat(wrapname, "not found for", get.source.str(context), paste(list(...), collapse=" "), "; generating it\n")
      saved <- f(context, ...)
      save(saved, file=savefile)
      return(saved)
    }
  })
}

# Returns a method for class 'aijContext' for f, which temporarily wraps its context, calls the 'wrAijContext' method, and then deletes all signs of its presence. Much faster for some 'aijContext' methods. Watch for potential infinite recursion when using this (don't define the 'wrAijContext' method using this 'aijContext' method)
withTempWrapping <- function(f) {
  # f must be an S3 generic function with first argument a context
  return(function(context, ...) {
    # Wrap the aijContext so intermediate results such as get.1comp.fits() can be cached
    wrAijContext <- get.wrapped.aij.context.from.aijContext(context, "Temp-wrapped", "temp_wrapped")
    
    results <- f(wrAijContext, ...)
    
    # clear cache so the next run doesn't find any saved files
    for(cachefile in list.files(path="Rsave", pattern=".*temp_wrapped.*Rdata", full.names = TRUE)) file.remove(cachefile)
    
    return(results)
  })
}

### 'Public' S3 methods for class 'wrAijContext' ###
# Some methods aren't enough computation to be worth caching, especially if their savefiles 
#   would be large (and thus take a while to save/load); others benefit greatly from caching

# Methods described in aijContext.R
get.genenames.wrAijContext <- function(context) get.genenames(context$aijContext)
get.expnames.wrAijContext <- function(context) get.expnames(context$aijContext)
get.expression.data.outliersNA.wrAijContext <- function(context) get.expression.data.outliersNA(context$aijContext)
get.expression.data.outliersInf.wrAijContext <- function(context) get.expression.data.outliersInf(context$aijContext)
get.operons.wrAijContext <- function(context) get.operons(context$aijContext)
get.sourcefile.wrAijContext <- function(context) "wrAijContext.R"
get.uni.raw.aijs.wrAijContext <- function(context) get.uni.raw.aijs(context$aijContext)
get.multi.raw.aijs.wrAijContext <- function(context) get.multi.raw.aijs(context$aijContext)
get.uni.mu0.wrAijContext <- function(context) get.uni.mu0(context$aijContext)
get.multi.mu0.wrAijContext <- function(context) get.multi.mu0(context$aijContext)
get.uni.mu1.wrAijContext <- function(context) get.uni.mu1(context$aijContext)
get.multi.mu1.wrAijContext <- function(context) get.multi.mu1(context$aijContext)
get.uni.var.wrAijContext <- function(context) get.uni.var(context$aijContext)
get.multi.var.wrAijContext <- function(context) get.multi.var(context$aijContext)
get.uni.pi.wrAijContext <- function(context) get.uni.pi(context$aijContext)
get.multi.pi.wrAijContext <- function(context) get.multi.pi(context$aijContext)
get.BICdiffs.by.gene.wrAijContext <- withCaching(get.BICdiffs.by.gene.aijContext, "BICdiffs_by_gene")
get.BICdiffs.by.operon.wrAijContext <- withCaching(get.BICdiffs.by.operon.aijContext, "BICdiffs_by_operon")
get.1comp.genenames.wrAijContext <- get.1comp.genenames.aijContext
get.1comp.operons.wrAijContext <- get.1comp.operons.aijContext
get.1comp.fits.wrAijContext <- withCaching(get.1comp.fits.aijContext, "onecomp_fits")
get.bayesfactors.by.gene.wrAijContext <- withCaching(get.bayesfactors.by.gene.aijContext, "bayesfactors_gene")
get.bayesfactors.by.operon.wrAijContext <- withCaching(get.bayesfactors.by.operon.aijContext, "bayesfactors_operon")
get.Cvalues.by.gene.wrAijContext <- get.Cvalues.by.gene.aijContext
get.Cvalues.by.operon.wrAijContext <- get.Cvalues.by.operon.aijContext
get.Cvalues.by.operon.oldmethod.wrAijContext <- get.Cvalues.by.operon.oldmethod.aijContext
get.operon.Cvalues.by.gene.wrAijContext <- get.operon.Cvalues.by.gene.aijContext
get.operon.Cvalues.by.gene.oldmethod.wrAijContext <- get.operon.Cvalues.by.gene.oldmethod.aijContext

# Methods described in Aij_generation.R
coinflip.wrAijContext <- withCaching(coinflip.aijContext, "coinflip")
MT.wrAijContext <- MT.aijContext 
TT.wrAijContext <- TT.aijContext
RB.wrAijContext <- RB.aijContext
UniMM.wrAijContext <- withCaching(UniMM.aijContext, "UniMM")
MultiMM.wrAijContext <- withCaching(MultiMM.aijContext, "MultiMM")
ConfMM.wrAijContext <- ConfMM.aijContext
UniMM_NBF.wrAijContext <- withCaching(UniMM_NBF.aijContext, "UniMM_Bayesian2")
MultiMM_NBF.wrAijContext <- withCaching(MultiMM_NBF.aijContext, "MultiMM_Bayesian2")
ConfMM_NBF.wrAijContext <- ConfMM_NBF.aijContext

# Redefine ConfMM.aijContext and ConfMM_NBF.aijContext to be much more efficient, using the caching tools in this script
# ConfMM.wrAijContext and ConfMM_NBF.wrAijContext will remain the "former" ("real") versions of ConfMM, as defined above
ConfMM.aijContext <- withTempWrapping(ConfMM.wrAijContext)
ConfMM_NBF.aijContext <- withTempWrapping(ConfMM_NBF.wrAijContext)

# Methods introduced in this script
get.data.source.wrAijContext <- function(context) context$data.source
get.source.str.wrAijContext <- function(context) context$source.str
