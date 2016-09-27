source("valContext.R")

options(width=180)  # allow to print wide tables

dichotomous <- function(aijs) aijs>=0.5  # aijs of exactly 0.5 are given TRUE - i.e. considered as 'on' predictions
trichotomous <- function(aijs, bufSize) {
  # bufSize: size of the "unknown" buffer on either side of 0.5.  E.g. bufSize==0.1 means aijs between 0.4 and 0.6 (inclusive) are counted as unknown.
  #   Notably, bufSize==0 still considers aijs of exactly 0.5 to be 'unknown', and this situation does come up frequently in certain conservative aij-generation methods.
  returnme <- aijs
  returnme[aijs<=0.5+bufSize & aijs>=0.5-bufSize] <- NA
  return(returnme >= 0.5)
} 

only.1comp.operons <- function(context, aijs_or_calls, BIC.bonus) {
  # context: Any aij context (even an unwrapped aijContext)
  # BIC.bonus: See notes on get.1comp.operons (in aijContext.R)
  one.comp.operons <- get.1comp.operons(aijContext, BIC.bonus)
  genes.to.keep <- intersect(unlist(get.operons(context)[one.comp.operons]), rownames(aijs_or_calls))
  aijs_or_calls[genes.to.keep,]
}

only.1comp.genes <- function(context, aijs_or_calls, BIC.bonus) {
  # context: Any aij context (even an unwrapped aijContext)
  # BIC.bonus: See notes on get.1comp.genenames (in aijContext.R)
  genes.to.keep <- intersect(get.1comp.genenames(context, BIC.bonus), rownames(aijs_or_calls))
  aijs_or_calls[genes.to.keep,]
}

calibration.plots <- function(valContext, numBins) {
  # valContext: An aij validation context; see valContext.R
  
  single.calibration.plot <- function(aijs, methodName) {
    percents <- numeric(numBins)
    totals <- numeric(numBins)
    for(bin in 1:numBins) {
      bin.right.edge <- bin/numBins
      bin.left.edge <- bin.right.edge-1/numBins
      indices <- which(aijs >= bin.left.edge & aijs <= bin.right.edge)
      total.aijs.in.bin <- length(indices)
      on.aijs.in.bin <- sum(get.gold.calls(valContext)[indices]==1)
      percents[bin] <- as.numeric(on.aijs.in.bin)/as.numeric(total.aijs.in.bin)
      totals[bin] <- total.aijs.in.bin
    }
    
    bin.labels <- sapply(1:numBins, function(bin) paste((bin-1)/numBins, "to", bin/numBins))
    png(paste0("results/", get.valSourceStr(valContext), "/","calibration_plot_", methodName, ".png"))
    plot(1:numBins, percents, 
         xlab="Aij",ylab="% active",main=methodName,
         axes=F,
         ylim=c(0,1))
    axis(1, at=1:numBins, labels=bin.labels, las=3)
    axis(2, at=seq(from=0,t=1,by=0.1), las=1)
    abline(-1/numBins, 1/numBins)  # Lower extreme of each bin
    abline(0, 1/numBins)  # Upper extreme of each bin
    dev.off()
    return(list(percents=percents, totals=totals))
  }
  
  mapply(single.calibration.plot, get.aijList(valContext), names(get.aijList(valContext)))
  return(invisible(TRUE)) # success
  
}

get.deviations <- function(aijs, gold_calls, type="abs") {
  # type: Either "abs", "sq", or "log", giving the absolute value, the square, or the log base 10 of (1 - the absolute value) of the deviations, respectively
  # If "log" is selected, any aijs greater than 0.99 or less than 0.01 will also be replaced with 0.99 or 0.01 respectively
  if(type=="abs") return(abs(aijs - gold_calls))
  else if(type=="sq") return((aijs - gold_calls)^2)
  else if(type=="log") {
    aijs[aijs>0.9999] <- 0.9999
    aijs[aijs<0.0001] <- 0.0001
    return(log(1-abs(aijs - gold_calls), base=10))
  }
  else stop(paste("Invalid type:", type))
}

deviations.by.gene <- function(aijs, gold_calls, filter=TRUE, type="abs") {
  # filter: If not the default (TRUE), then a matrix the same size as aijs; any position in the matrix with "FALSE" will be filtered out of all calculations
  # type: Either "abs", "sq", or "log".  
  #   "abs" will give 1 minus the average (mean) absolute value of the deviations
  #   "sq" will give 1 minus the square root of the mean square deviation
  #   "log" will give e to the mean logloss, where the logloss is the log base 10 of (1 minus the absolute value of) the deviation
  if(any(is.na(filter))) stop("NAs present in filter")
  deviations <- get.deviations(aijs, gold_calls, type)
  deviations[!filter] <- NA
  meandevs <- apply(deviations, 1, mean, na.rm=T)
  if(type=="abs") {
    return(1-meandevs)
  } else if(type=="sq") {
    return(1-sqrt(meandevs))
  } else if(type=="log") {
    return(exp(meandevs))
  } else {
    stop(paste("Invalid type:", type))
  }
}

deviations.by.operon <- function(aijs, gold_calls, filter=TRUE, type="abs") {
  # filter: If not the default (TRUE), then a matrix the same size as aijs; any position in the matrix with "FALSE" will be filtered out of all calculations
  # type: Either "abs", "sq", or "log".  
  #   "abs" will give 1 minus the average (mean) absolute value of the deviations
  #   "sq" will give 1 minus the square root of the mean square deviation
  #   "log" will give e to the mean logloss, where the logloss is the log base 10 of (1 minus the absolute value of) the deviation
  if(any(is.na(filter))) stop("NAs present in filter")
  deviations <- get.deviations(aijs, gold_calls, type)
  deviations[!filter] <- NA
  meandevs <- sapply(get.valOperons(valContext), function(operon) mean(deviations[operon,], na.rm=T))
  if(type=="abs") {
    return(1-meandevs)
  } else if(type=="sq") {
    return(1-sqrt(meandevs))
  } else if(type=="log") {
    return(exp(meandevs))
  } else {
    stop(paste("Invalid type:", type))
  }
}

deviations.by.experiment <- function(aijs, gold_calls, filter=TRUE, type="abs") {
  # filter: If not the default (TRUE), then a matrix the same size as aijs; any position in the matrix with "FALSE" will be filtered out of all calculations
  # type: Either "abs", "sq", or "log".  
  #   "abs" will give 1 minus the average (mean) absolute value of the deviations
  #   "sq" will give 1 minus the square root of the mean square deviation
  #   "log" will give e to the mean logloss, where the logloss is the log base 10 of (1 minus the absolute value of) the deviation
  if(any(is.na(filter))) stop("NAs present in filter")
  deviations <- get.deviations(aijs, gold_calls, type)
  deviations[!filter] <- NA
  meandevs <- apply(deviations, 2, mean, na.rm=T)
  if(type=="abs") {
    return(1-meandevs)
  } else if(type=="sq") {
    return(1-sqrt(meandevs))
  } else if(type=="log") {
    return(exp(meandevs))
  } else {
    stop(paste("Invalid type:", type))
  }
}

hist.deviations <- function(valContext, byWhat="gene", filter=TRUE, filterName="overall", type="abs") {
  for(index in 1:length(get.aijList(valContext))) {
    aijs <- get.aijList(valContext)[[index]]
    deviations <- switch(byWhat, 
      gene=deviations.by.gene(aijs=aijs, gold_calls=get.gold.calls(valContext), filter=filter, type=type),
      operon=deviations.by.operon(aijs=aijs, gold_calls=get.gold.calls(valContext), filter=filter, type=type),
      experiment=deviations.by.experiment(aijs=aijs, gold_calls=get.gold.calls(valContext), filter=filter, type=type)
    )
    methodName <- names(get.aijList(valContext))[[index]]
    png(paste0("results/", get.valSourceStr(valContext),"/",methodName,"_",type,"_deviations_by_",byWhat,"_",filterName,".png"))
    hist(deviations, 
      xlab=paste(
        switch(type, abs="1 - Average absolute deviation", sq="1 - Square root of the mean square deviation", log="e to the mean logloss"),
        "for the", byWhat, "across all", ifelse(byWhat=="experiment", "genes", "experiments")), 
      ylab=paste0("Number of ",byWhat,"s"), 
      main=methodName)
    dev.off()
  }
}

boxplot.deviations <- function(valContext, byWhat="gene", filter=TRUE, filterName="overall", type="abs") {
  # valContext: An aij validation context; see valContext.R
  deviationsList <- do.call(c, lapply(get.aijList(valContext), function(aijs) {
    deviations <- switch(byWhat, 
      gene=deviations.by.gene(aijs=aijs, gold_calls=get.gold.calls(valContext), filter=filter, type=type),
      operon=deviations.by.operon(aijs=aijs, gold_calls=get.gold.calls(valContext), filter=filter, type=type),
      experiment=deviations.by.experiment(aijs=aijs, gold_calls=get.gold.calls(valContext), filter=filter, type=type)
    )
    deviations <- deviations[!is.nan(deviations)]  # NaNs result when an entire byWhat was filtered out (i.e. an entire gene, if byWhat="gene")
    if(length(deviations)==0) return(list())
    else return(list(deviations))
  }))
  if(length(deviationsList)==0) {
    cat("Skipping boxplot for", filterName, ": no data points pass filter\n")
  } else {
    png(paste0("results/", get.valSourceStr(valContext),"/","boxplot_",type,"_deviations_by_",byWhat,"_",filterName,".png"))
    par(mar=c(13,4,2,2)+0.1)
    boxplot(deviationsList, range=3, las=2)
    dev.off()
  }
}

mean.alignment <- function(aijs, gold_calls, filter=TRUE, type="abs") {
  # filter: If not the default (TRUE), then a matrix the same size as aijs; any position in the matrix with "FALSE" will be filtered out of all calculations
  # type: Either "abs", "sq", or "log".  
  #   "abs" will give 1 minus the average (mean) absolute value of the deviations
  #   "sq" will give 1 minus the square root of the mean square deviation
  #   "log" will give e to the mean logloss, where the logloss is the log base 10 of (1 minus the absolute value of) the deviation
  if(any(is.na(filter))) stop("NAs present in filter")
  deviations <- get.deviations(aijs, gold_calls, type)
  deviations[!filter] <- NA
  meandev <- mean(deviations, na.rm=T)
  if(type=="abs") {
    return(1-meandev)
  } else if(type=="sq") {
    return(1-sqrt(meandev))
  } else if(type=="log") {
    return(exp(meandev))
  } else {
    stop(paste("Invalid type:", type))
  }
}

# Returns a named vector with the result of mean.alignment for each set of aijs
all.mean.alignments <- function(valContext, filter=TRUE, type="abs") {
  sapply(get.aijList(valContext), mean.alignment, gold_calls=get.gold.calls(valContext), filter=filter, type=type)
}

confidence.tables <- function(valContext, filter=TRUE, filterName="OVERALL") {
  # valContext: An aij validation context; see valContext.R
  # filter: If not the default (TRUE), then a matrix the same size as a set of aijs; any position in the matrix with "FALSE" will be filtered out of all calculations
  # filterName: Only used to label the output
  
  if(any(is.na(filter))) stop("NAs present in filter")
  filter <- filter & !is.na(get.gold.calls(valContext))  # For this function, only consider gene-ex pairs with gold calls
  
  single.conftable <- function(aijs) {
    total <- length(aijs[filter])
    highconf <- (aijs >= 0.8 | aijs <= 0.2) & filter
    medconf <- (aijs >= 0.6 | aijs <= 0.4) & !highconf & filter
    lowconf <- !medconf & !highconf & filter
    results <- matrix(nr=4, nc=2, dimnames = list(c("Overall (all confidences)", "High confidence (0-0.2 or 0.8-1)", "Medium confidence (0.2-0.4 or 0.6-0.8)", "Low confidence (0.4-0.6)"), c("Count", "% Correct")))
    results[1,1] <- sprintf("%d (%.1f%%)", total, 100)
    results[2,1] <- sprintf("%d (%.1f%%)", sum(highconf), 100*sum(highconf)/total)
    results[3,1] <- sprintf("%d (%.1f%%)", sum(medconf), 100*sum(medconf)/total)
    results[4,1] <- sprintf("%d (%.1f%%)", sum(lowconf), 100*sum(lowconf)/total)
    correct <- (aijs>.5)==get.gold.calls(valContext) & aijs!=.5 & filter
    unk <- aijs==.5 & filter
    results[1,2] <- tryCatch(
      sprintf("%.1f%% (%d/%d)", 100*( sum(correct) + 0.5*sum(unk) )/total, floor(sum(correct) + 0.5*sum(unk)), total)
      , error=function(e) {
        cat("filterName =", filterName, "\n")
        cat("sum(correct) =", sum(correct), "\n")
        cat("sum(unk) =", sum(unk), "\n")
        cat("total =", total, "\n")
        "ERROR"
      }
    )
    results[2,2] <- sprintf("%.1f%% (%d/%d)", 100*sum(highconf & correct)/sum(highconf), sum(highconf & correct), sum(highconf))
    results[3,2] <- sprintf("%.1f%% (%d/%d)", 100*sum(medconf & correct)/sum(medconf), sum(medconf & correct), sum(medconf))
    results[4,2] <- sprintf("%.1f%% (%d/%d)", 100*( sum(lowconf & correct) + 0.5*sum(unk) )/sum(lowconf), floor(sum(lowconf & correct) + 0.5*sum(unk)), sum(lowconf))
    return(results)
  }
  
  cat(filterName, "\n")
  print(lapplyWithNames(get.aijList(valContext), single.conftable))
  cat("\n")
  
}

consistency.tables <- function(valContext, filter=TRUE, filterName="OVERALL", onWhat="operons", withCorrectness=TRUE) {
  # valContext: An aij validation context; see valContext.R
  # filter: If not the default (TRUE), then a matrix the same size as aijs. 
  #   Any position in the matrix with "FALSE" will be filtered out of all calculations,
  #   ALONG WITH the rest of its operon for that experiment
  # onWhat: Either "operons" or "scenarios"
  # withCorrectness: If TRUE, correctness statistics by consistency will also be computed, in addition to consistency statistics
  
  if(any(is.na(filter))) stop("NAs present in filter")
  if(withCorrectness) filter <- filter & !is.na(get.gold.calls(valContext))  # For this function, only consider gene-ex pairs with gold calls
  
  if(onWhat=="operons") {
    multi.gene.sets <- get.valOperons(valContext)[sapply(get.valOperons(valContext), length)>1]  # 'set' can be either an operon or a scenario
  } else if(onWhat=="scenarios") {
    multi.gene.sets <- get.scenarios(valContext)[sapply(get.scenarios(valContext), length)>1]
  } else {
    stop(paste("Unrecognized onWhat:", onWhat))
  }
  if(length(multi.gene.sets)==0) {
    cat("Omitting consistency tables for", filterName, "because there were no", onWhat, "with more than 1 gene\n\n")
    return()
  }
  
  single.constable <- function(aijs) {
    
    aijs[!filter] <- NA
    
    # Old rownames
    #result.rownames <- c(
    #  "Very consistent (all aijs >0.8 or all aijs <0.2)", 
    #  "Consistent (all aijs >0.6 or all aijs <0.4, but not VC)", 
    #  "Borderline (all aijs >0.4 or all aijs <0.6, but not VC or C)", 
    #  "Inconsistent (at least one aij >0.6 AND at least one aij <0.4, but not VI)", 
    #  "Very inconsistent (at least one aij >0.8 AND at least one aij <0.2)",
    #  "Leftover (fits into none of the above categories)"
    #)
    result.rownames <- c(
      "Very consistent (all aijs >0.8 or all aijs <0.2)",
      "Consistent (all aijs >0.4 or all aijs <0.6, but not VC)",
      "Ignore this line",
      "Inconsistent (at least one aij >0.6 AND at least one aij <0.4)",
      "Ignore this line",
      "Leftover (fits into none of the above categories)"
    )
    if(withCorrectness) result.colnames <- c("Total count", "% correct (when aijs and model cons)", "Count of cons (aijs, model)")
    else result.colnames <- c("Count in each category")
    
    cat("First set of Boolean matrices:\n"); print(system.time({
    
    # These functions take a vector of aijs for an operon,experiment (or scenario,experiment) pair.
    # Old definitions
    #is.vc <- function(vec) all(vec>0.8) | all(vec<0.2)
    #is.c <- function(vec) (all(vec>0.6) | all(vec<0.4)) & !is.vc(vec)
    #is.b <- function(vec) (all(vec>0.4) | all(vec<0.6)) & !is.vc(vec) & !is.c(vec)
    #is.i <- function(vec) (any(vec>0.6) & any(vec<0.4)) & !is.vi(vec)
    #is.vi <- function(vec) any(vec>0.8) & any(vec<0.2)
    
    # New definitions
    is.vc <- function(vec) all(vec>0.8) | all(vec<0.2)
    is.c <- function(vec) (all(vec>0.4) | all(vec<0.6)) & !is.vc(vec)
    is.b <- function(vec) FALSE
    is.i <- function(vec)  any(vec>0.6) & any(vec<0.4)
    is.vi <- function(vec) FALSE
      
    # Boolean matrices, each row is a multi-gene set and each experiment is a column
    # NA for entries that were filtered out
    vc <- t(sapply(multi.gene.sets, function(set) apply(aijs[set,], 2, is.vc)))
    c <- t(sapply(multi.gene.sets, function(set) apply(aijs[set,], 2, is.c)))
    b <- t(sapply(multi.gene.sets, function(set) apply(aijs[set,], 2, is.b)))
    i <- t(sapply(multi.gene.sets, function(set) apply(aijs[set,], 2, is.i)))
    vi <- t(sapply(multi.gene.sets, function(set) apply(aijs[set,], 2, is.vi)))
    
    l <- !vc & !c & !b & !i & !vi
  
    }, gcFirst=TRUE))
    
    
    cat("Second set of Boolean matrices:\n"); print(system.time({
    
    # Similar Boolean matrix, the overall prediction for the operon,experiment pair according to the aijs
    # NA for entries where the aijs were too inconsistent to determine a prediction, or filtered out
    aij.preds <- t(sapply(multi.gene.sets, function(set) apply(aijs[set,], 2, function(vec) {
      if(any(is.na(vec))) return(NA)
      if(all(vec>=0.5) && any(vec>0.5)) return(TRUE)
      if(all(vec<=0.5) && any(vec<0.5)) return(FALSE)
      return(NA)
    })))
    
    # Similar Boolean matrix, whether the aijs were consistent enough to make calls
    aijsconsistent <- !is.na(aij.preds)
      
    if(withCorrectness) {  
      
      # Similar Boolean matrix, what the model call was, or NA if the model was inconsistent
      model.calls <- t(sapply(multi.gene.sets, function(set) apply(get.gold.calls(valContext)[set,], 2, function(vec) {
        if(all(vec)) return(TRUE)
        if(all(!vec)) return(FALSE)
        return(NA)
      })))
  
      # Similar Boolean matrix, whether the model call matched the aij prediction
      # NA for entries where either were too inconsistent
      correct <- aij.preds==model.calls
      
      # Similar Boolean matrix, whether the model was consistent
      modelconsistent <- !is.na(model.calls)
      
      # Similar Boolean matrix, whether both were consistent enough to make calls
      bothconsistent <- modelconsistent & aijsconsistent
      
    }  
    
    }, gcFirst=TRUE))

    
    cat("Third set of Boolean matrices:\n"); print(system.time({
    
    # the number of operon,experiment pairs that were not filtered out
    total <- sum(!is.na(vc))  # all 6 Boolean matrices (vc-l) will have NAs in the same positions
    
    vc.total <- sum(vc, na.rm=T)
    c.total <- sum(c, na.rm=T)
    b.total <- sum(b, na.rm=T)
    i.total <- sum(i, na.rm=T)
    vi.total <- sum(vi, na.rm=T)
    l.total <- sum(l, na.rm=T)       # should be 0, if the bins above are meant to catch all possible situations
    
    vc.aijsconsistent <- sum(vc & aijsconsistent, na.rm=T)
    c.aijsconsistent <- sum(c & aijsconsistent, na.rm=T)
    b.aijsconsistent <- sum(b & aijsconsistent, na.rm=T)
    i.aijsconsistent <- sum(i & aijsconsistent, na.rm=T)
    vi.aijsconsistent <- sum(vi & aijsconsistent, na.rm=T)
    l.aijsconsistent <- sum(l & aijsconsistent, na.rm=T)
    
    if(withCorrectness) {
      vc.modelconsistent <- sum(vc & modelconsistent, na.rm=T)
      c.modelconsistent <- sum(c & modelconsistent, na.rm=T)
      b.modelconsistent <- sum(b & modelconsistent, na.rm=T)
      i.modelconsistent <- sum(i & modelconsistent, na.rm=T)
      vi.modelconsistent <- sum(vi & modelconsistent, na.rm=T)
      l.modelconsistent <- sum(l & modelconsistent, na.rm=T)
      
      vc.bothconsistent <- sum(vc & bothconsistent, na.rm=T)
      c.bothconsistent <- sum(c & bothconsistent, na.rm=T)
      b.bothconsistent <- sum(b & bothconsistent, na.rm=T)
      i.bothconsistent <- sum(i & bothconsistent, na.rm=T)
      vi.bothconsistent <- sum(vi & bothconsistent, na.rm=T)
      l.bothconsistent <- sum(l & bothconsistent, na.rm=T)
      
      vc.correct <- sum(vc & correct, na.rm=T)
      c.correct <- sum(c & correct, na.rm=T)
      b.correct <- sum(b & correct, na.rm=T)
      i.correct <- sum(i & correct, na.rm=T)    # should be 0
      vi.correct <- sum(vi & correct, na.rm=T)  # should be 0
      l.correct <- sum(l & correct, na.rm=T)    # should be 0, if the bins above are meant to catch all possible situations
    }
    
    }, gcFirst=TRUE))
    
    if(withCorrectness) results <- matrix(nr=6, nc=3, dimnames = list(result.rownames, result.colnames))
    else results <- matrix(nr=6, nc=1, dimnames = list(result.rownames, result.colnames))
    
    results[1,1] <- sprintf("%d/%d (%.1f%%)", vc.total, total, 100*vc.total/total)
    results[2,1] <- sprintf("%d/%d (%.1f%%)", c.total, total, 100*c.total/total)
    results[3,1] <- sprintf("%d/%d (%.1f%%)", b.total, total, 100*b.total/total)
    results[4,1] <- sprintf("%d/%d (%.1f%%)", i.total, total, 100*i.total/total)
    results[5,1] <- sprintf("%d/%d (%.1f%%)", vi.total, total, 100*vi.total/total)
    results[6,1] <- sprintf("%d/%d (%.1f%%)", l.total, total, 100*l.total/total)
    
    if(withCorrectness) {
      results[1,2] <- sprintf("%.1f%% (%d/%d)", 100*vc.correct/vc.bothconsistent, vc.correct, vc.bothconsistent)
      results[2,2] <- sprintf("%.1f%% (%d/%d)", 100*c.correct/c.bothconsistent, c.correct, c.bothconsistent)
      results[3,2] <- sprintf("%.1f%% (%d/%d)", 100*b.correct/b.bothconsistent, b.correct, b.bothconsistent)
      results[4,2] <- sprintf("%.1f%% (%d/%d)", 100*i.correct/i.bothconsistent, i.correct, i.bothconsistent)
      results[5,2] <- sprintf("%.1f%% (%d/%d)", 100*vi.correct/vi.bothconsistent, vi.correct, vi.bothconsistent)
      results[6,2] <- sprintf("%.1f%% (%d/%d)", 100*l.correct/l.bothconsistent, l.correct, l.bothconsistent)
      
      results[1,3] <- sprintf("%d, %d", vc.aijsconsistent, vc.modelconsistent)
      results[2,3] <- sprintf("%d, %d", c.aijsconsistent, c.modelconsistent)
      results[3,3] <- sprintf("%d, %d", b.aijsconsistent, b.modelconsistent)
      results[4,3] <- sprintf("%d, %d", i.aijsconsistent, i.modelconsistent)
      results[5,3] <- sprintf("%d, %d", vi.aijsconsistent, vi.modelconsistent)
      results[6,3] <- sprintf("%d, %d", l.aijsconsistent, l.modelconsistent)
    }
    
    return(results)
    
  }
  
  cat(filterName, "\n")
  print(lapplyWithNames(get.aijList(valContext), single.constable))
  cat("\n")
  
}

sensitivity <- function(aij_calls, gold_calls) {
  # aij_calls: a boolean matrix of what the aij thinks is active vs. inactive.  May contain NAs (don't knows).
  sum(aij_calls[gold_calls], na.rm = TRUE)/sum(!is.na(aij_calls[gold_calls]))
}

specificity <- function(aij_calls, gold_calls) {
  # aij_calls: a boolean matrix of what the aij thinks is active vs. inactive.  May contain NAs (don't knows).
  sum((!aij_calls)[!gold_calls], na.rm = TRUE)/sum(!is.na(!aij_calls[!gold_calls]))
}

sensitivity.plot <- function(valContext, pngFilename) {
  # valContext: An aij validation context; see valContext.R
  aijList <- get.aijList(valContext)
  bufSizes <- seq(from=0, to=0.5, by=0.01)
  png(pngFilename)
  plot(bufSizes, sapply(bufSizes, function(bufSize) sensitivity(trichotomous(aijList[[1]], bufSize), get.gold.calls(valContext))),
       xlim=c(0,0.5), ylim=c(0,1),
       xlab="Buffer size", ylab="Sensitivity", main="Sensitivity of various gene-calling methods", 
       pch=21, col=2, bg=2)
  if(length(aijList)>1) {
    for(index in 2:length(aijList)) {
      points(bufSizes, sapply(bufSizes, function(bufSize) sensitivity(trichotomous(aijList[[index]], bufSize), get.gold.calls(valContext))),
             pch=21, col=index+1, bg=index+1)
    }
  }
  legend("bottomright", names(aijList), pch=21, col=2:(length(aijList)+1), pt.bg=2:(length(aijList)+1))
  dev.off()
}

specificity.plot <- function(valContext, pngFilename) {
  # valContext: An aij validation context; see valContext.R
  aijList <- get.aijList(valContext)
  bufSizes <- seq(from=0, to=0.5, by=0.01)
  png(pngFilename)
  plot(bufSizes, sapply(bufSizes, function(bufSize) specificity(trichotomous(aijList[[1]], bufSize), get.gold.calls(valContext))),
       xlim=c(0,0.5), ylim=c(0,1),
       xlab="Buffer size", ylab="Specificity", main="Specificity of various gene-calling methods", 
         pch=21, col=2, bg=2)
  if(length(aijList)>1) {
    for(index in 2:length(aijList)) {
      points(bufSizes, sapply(bufSizes, function(bufSize) specificity(trichotomous(aijList[[index]], bufSize), get.gold.calls(valContext))),
             pch=21, col=index+1, bg=index+1)
    }
  }
  legend("bottomright", names(aijList), pch=21, col=2:(length(aijList)+1), pt.bg=2:(length(aijList)+1))
  dev.off()
}

roc.plot <- function(valContext, pngFilename) {
  # valContext: An aij validation context; see valContext.R
  aijThresholds <- seq(from=0, to=1, by=0.02)
  pointsOnCurve <- lapply(get.aijList(valContext), function(aijs) {
    sapply(aijThresholds, function(threshold) {
      calls <- aijs>threshold
      return(c(x=1-specificity(calls, get.gold.calls(valContext)), y=sensitivity(calls, get.gold.calls(valContext))))
    })
  })
  png(pngFilename)
  plot(x=c(0,0), type="n", xlim=c(0,1), ylim=c(0,1), xlab="False positive rate (1-specificity)", ylab="True positive rate (sensitivity)")
  abline(a=0, b=1)
  for(index in 1:length(pointsOnCurve)) {
    method <- pointsOnCurve[[index]]
    points(method["x",], method["y",], type="l", col=index+1)
  }
  legend("bottomright", names(get.aijList(valContext)), fill=2:(length(get.aijList(valContext))+1))
  dev.off()
}

# Returns a numeric vector (percent Consistent, percent Unknown, percent Inconsistent)
# percent Consistent: percentage of (operon, experiment) pairs (where the operon has more than one gene) where the calls are either all 'active' or all 'inactive'
# percent Unknown: percentage of (operon, experiment) pairs where the calls either include only 'active' and 'unknown', or only 'inactive' and 'unknown'
# percent Inconsistent: percentage of (operon, experiment) pairs where the calls include at least one 'active' AND at least one 'inactive'
operon.consistency <- function(valContext, aij_calls) {
  # valContext: An aij validation context; see valContext.R
  # aij_calls: a boolean matrix of what the aij thinks is active vs. inactive.  May contain NAs (don't knows).
  operonmatrices <- lapply(get.valOperons(valContext), function(operon) if(length(operon)<=1) return(numeric(0)) else return(aij_calls[operon,]))
  operonmatrices <- operonmatrices[sapply(operonmatrices, length)>0]
  if(length(operonmatrices)==0) return(c(NA, NA, NA))  # no multiple-gene operons in the data
  totalC <- 0
  totalU <- 0
  totalI <- 0
  for(operonmatrix in operonmatrices) {
    apply(operonmatrix, 2, function(operoncalls) {
      if(any(is.na(operoncalls))) {
        if(any(operoncalls, na.rm = T) && any(!operoncalls, na.rm = T)) {
          totalI <<- totalI + 1
        } else {
          totalU <<- totalU + 1
        }
      } else {
        if(any(operoncalls) && any(!operoncalls)) {
          totalI <<- totalI + 1
        } else {
          totalC <<- totalC + 1
        }
      }
    })
  }
  return(c(totalC, totalU, totalI)/sum(sapply(operonmatrices, function(mat) ncol(mat))))
}

operon.consistency.summary <- function(valContext, bufSizes=c(0, 0.1, 0.2)) {
  # valContext: An aij validation context; see valContext.R
  # bufSizes: Various parameters for trichotomous() to test
  
  printConsistency <- function(aijs, methodName) {
    cat("Operon consistency of", methodName, "\n")
    consistency <- operon.consistency(valContext, aijs)
    cat("Percent consistent:", consistency[1], "\n")
    cat("Percent unknown:", consistency[2], "\n")
    cat("Percent inconsistent:", consistency[3], "\n")
    cat("\n")
  }
  
  sapply(bufSizes, function(bufSize) mapply(printConsistency, lapplyWithNames(get.aijList(valContext), trichotomous, bufSize), paste(names(get.aijList(valContext)), "trichotomized by", bufSize)))
  
  return(invisible(TRUE)) # success
    
}

informativeness <- function(aijs) {
  return(mean(abs(aijs-0.5)))
}

hist.C.values.by.gene <- function(valContext, A, B, exponential, filter=TRUE, filterName="overall") {
  # filter, unlike in other functions, if not the default then should be just one boolean per gene in valContext
  png(paste0("results/", get.valSourceStr(valContext),"/hist_univariate_Cvalues_",filterName,".png"), width=436, height=294)
  hist(get.Cvalues.by.gene(valContext, A, B, exponential)[get.valGenenames(valContext)][filter], 
       breaks=seq(from=0,to=1,by=0.1), ylim=c(0,length(get.valGenenames(valContext)[filter])), col="orange",
       main=paste("Univariate C-values for",filterName,"genes"), xlab="C-value")
  dev.off()
}
hist.C.values.by.operon <- function(valContext, A, B, exponential, filter=TRUE, filterName="overall") {
  # filter, unlike in other functions, if not the default then should be just one boolean per gene in valContext
  png(paste0("results/", get.valSourceStr(valContext),"/hist_multivariate_Cvalues_",filterName,".png"), width=436, height=294)
  hist(get.operon.Cvalues.by.gene(valContext, A, B, exponential)[get.valGenenames(valContext)][filter], 
       breaks=seq(from=0,to=1,by=0.1), ylim=c(0,length(get.valGenenames(valContext)[filter])), col="orange",
       main=paste("Multivariate C-values for",filterName,"genes"), xlab="C-value")
  dev.off()
}
hist.old.equiv.C.by.gene <- function(valContext, BIC.bonus, filter=TRUE, filterName="overall") {
  # filter, unlike in other functions, if not the default then should be just one boolean per gene in valContext
  png(paste0("results/", get.valSourceStr(valContext),"/hist_univariate_Bayesian1_C_",BIC.bonus,"_",filterName,".png"), width=436, height=294)
  hist(as.integer(!(get.valGenenames(valContext)[filter] %in% get.1comp.genenames(valContext, BIC.bonus))), 
       breaks=seq(from=0,to=1,by=0.1), ylim=c(0,length(get.valGenenames(valContext)[filter])), col="orange",
       main=paste0("Univariate Bayesian1 C-values for ",filterName," genes, bonus=",BIC.bonus), xlab="C-value")
  dev.off()
}
hist.old.equiv.C.by.operon <- function(valContext, BIC.bonus, filter=TRUE, filterName="overall") {
  # filter, unlike in other functions, if not the default then should be just one boolean per gene in valContext
  genenames.in.1comp.operons <- unlist(get.operons(valContext)[get.1comp.operons(valContext, BIC.bonus)])
  png(paste0("results/", get.valSourceStr(valContext),"/hist_multivariate_Bayesian1_C_",BIC.bonus,"_",filterName,".png"), width=436, height=294)
  hist(as.integer(!(get.valGenenames(valContext)[filter] %in% genenames.in.1comp.operons)), 
       breaks=seq(from=0,to=1,by=0.1), ylim=c(0,length(get.valGenenames(valContext)[filter])), col="orange",
       main=paste0("Multivariate Bayesian1 C-values for ",filterName," genes, bonus=",BIC.bonus), xlab="C-value")
  dev.off()
}

C.value.mean.square.alignment <- function(valContext, A, B, exponential, multivariate=TRUE, filter=TRUE, filterName="overall") {
  # filter, unlike in other functions, if not the default then should be just one boolean per gene in valContext
  Cvalues <- if(multivariate) get.operon.Cvalues.by.gene(valContext, A, B, exponential)[get.valGenenames(valContext)][filter]
             else get.Cvalues.by.gene(valContext, A, B, exponential)[get.valGenenames(valContext)][filter]
  target.Cvalues <- get.target.Cvalues(valContext, filter)
  return(1-sqrt(mean((Cvalues-target.Cvalues)^2, na.rm=T)))
}
old.equiv.C.mean.square.alignment <- function(valContext, BIC.bonus, multivariate=TRUE, filter=TRUE, filterName="overall") {
  # filter, unlike in other functions, if not the default then should be just one boolean per gene in valContext
  genenames.in.1comp.operons <- unlist(get.operons(valContext)[get.1comp.operons(valContext, BIC.bonus)])
  Cvalues <- if(multivariate) as.integer(!(get.valGenenames(valContext)[filter] %in% genenames.in.1comp.operons))
             else as.integer(!(get.valGenenames(valContext)[filter] %in% get.1comp.genenames(valContext, BIC.bonus)))
  target.Cvalues <- get.target.Cvalues(valContext, filter)
  return(1-sqrt(mean((Cvalues-target.Cvalues)^2, na.rm=T)))
}

# Returns a vector of target Cvalues, named with genes. Some genes may have NA (unknown) target Cvalues. 
get.target.Cvalues <- function(valContext, filter=TRUE) {
  # filter, unlike in other functions, if not the default then should be just one boolean per gene in valContext
  gold.calls <- get.gold.calls(valContext)[filter,]
  constant <- apply(gold.calls, 1, all) | apply(!gold.calls, 1, all)  # either always-active or always-inactive
  return(as.integer(!constant))  # 1 for not-constant, 0 for constant
}
