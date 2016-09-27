source("iMRM.R")
source("expnameMaps.R")

lapplyWithNames <- function(...) sapply(..., simplify=FALSE, USE.NAMES=TRUE)

run <- function(modelName="ecoli", exp.to.run=character(0), GBP=TRUE) {
  # modelName: choose from "ecoli", "shewanella", "baby", "Carrera"
  # exp.to.run: an experiment name to run the iMRM for. character(0), the default,
  #   will be treated as "all experiments"
  # GBP: passed as the GBP argument to iMRM() or iMRM_batch() (see notes on those functions)
  
  if(modelName == "shewanella") {
    sbmlFile <- "inputs/iMR1_799_6-19-15.sbml"
    aijFile <- "aijs/multivar_aijs_shew_7-21-15_edited.txt"
    mediaFile <- "inputs/m3d_media_cpids.txt"
    celFile <- "inputs/S_oneidensis_cel_to_chip.txt"
    experimentsFile <- NA
    experimentsFileFormat <- NA
    commonGenenamesFile <- NA
  } else if (modelName == "ecoli") {
    sbmlFile <- "inputs/iJO1366_7-15-15.sbml"
    aijFile <- "aijs/ecoli_multimm_aijs_1-19-16.txt"
    mediaFile <- "inputs/media_ecoli_plus.txt"
    celFile <- "inputs/Ecoli_Cel_to_chip_fromClaire_edited.tab"
    experimentsFile <- "inputs/E_coli_v4_Build_6_6-23-15.experiment_feature_descriptions"
    experimentsFileFormat <- "Hope"
    commonGenenamesFile <- "inputs/E_coli_common_genenames_edited.txt"
  } else if (modelName == "baby") {
    sbmlFile <- "inputs/babyModel_edited.sbml"
    aijFile <- "aijs/baby_aijs_edited.txt"
    mediaFile <- "inputs/baby_media.txt"
    celFile <- "inputs/baby_cel_to_chip.txt"
    experimentsFile <- NA
    experimentsFileFormat <- NA
    commonGenenamesFile <- NA
  } else if (modelName == "Carrera") {
    sbmlFile <- "inputs/Carrera_iJO1366.sbml"
    aijFile <- "aijs/multivar_aijs_carrera_ecoli.txt"
    mediaFile <- NA
    celFile <- NA
    experimentsFile <- "inputs/Carrera_experiment_descriptions.txt"
    experimentsFileFormat <- "Carrera"
    commonGenenamesFile <- "inputs/Carrera_genenames.txt"
  } else {
    stop(paste("Unrecognized model name:", modelName))
  }
  
  aijdata <- loadaijdata(aijFile, exp.to.run)
  expnameMap <- get.expname.map(celFile)
  medias <- loadMediaData(mediaFile, aijdata, expnameMap)
  knockouts <- loadKnockoutData(experimentsFile, experimentsFileFormat, commonGenenamesFile, expnameMap)
  
  # Indicate "no knockouts" for all experiments with no knockouts
  noKnockExps <- colnames(aijdata)[!(colnames(aijdata) %in% names(knockouts))]
  knockouts[noKnockExps] <- rep(list(character(0)), length(noKnockExps))
  
  # For these two models, the compound names in the media files are missing "_e0" tags
  if(modelName=="ecoli" || modelName=="shewanella") {
    medias <- lapplyWithNames(medias, function(media) {
      names(media) <- paste0(names(media), "_e0")
      return(media)
    })
  }
  
  # This model needs standardization of gene names
  if(modelName=="shewanella") {
    rownames(aijdata) <- ifFigConvertToKBase(rownames(aijdata))
    knockouts <- ifFigConvertToKBase(knockouts)
  }
  
  results <- iMRM_batch(sbmlFile, aijdata, medias, knockouts, 0.1, 0.5, GBP=GBP)
  
  cat("results:"); str(results)
  
  sapply(1:ncol(aijdata), function(expnum) {
    filename <- paste0("results/",colnames(aijdata)[expnum],"_",format(Sys.time(), format="%Y-%m-%d_%H-%M-%S"),".txt")
    write.table(results[[expnum]], filename, sep="\t\t", quote=F, col.names=T)
    cat("Results written to", filename, "\n")
  })
  
  return(invisible(results))
}

loadaijdata <- function(aijFile, experiments=character(0)) {
  # Returns aij data for the specified experiment(s), in the form of a matrix
  # experiments: a string vector of experiment names (stripped CEL format); character(0) will be treated as "all experiments"
  
  aijdata <- as.matrix(read.table(aijFile))
  expnames <- colnames(aijdata)
  if(all(substr(expnames, nchar(expnames)-6, nchar(expnames)) == ".CEL.gz")) {
    colnames(aijdata) <- substr(expnames, 1, nchar(expnames)-7)  # Remove the .CEL.gz from each colname
  }
  if(length(experiments)==0) return(aijdata)
  return(aijdata[,experiments, drop=FALSE])  # keep it as a matrix, even if there is only one experiment
}

loadMediaData <- function(mediaFile, aijdata, expnameMap) {
  
  mediadata <- read.table(mediaFile, quote="", header=F, fill=T, stringsAsFactors=F, row.names=NULL, flush=T)
  mediadata <- mediadata[mediadata[[2]]!="is",]  # Drop lines of the form "expname is Complete", if any
  breaks <- which(mediadata[[1]]=="Adding")  # these line numbers start new experiments
  breaks <- c(breaks, nrow(mediadata)+1)  # add a "break" after the last line
  medias <- lapply(1:(length(breaks)-1), function(index) {
    mediasub <- mediadata[(breaks[index]+1):(breaks[index+1]-1),c(1,3)]  # the relevant snippet of mediadata (a data frame)
    mediavector <- as.vector(mediasub[[2]])  # the concentration values
    names(mediavector) <- mediasub[[1]]  # the actual compounds
    return(mediavector)  # return a named vector of concentrations, where the compounds are the names
  })
  breaks <- head(breaks,-1)  # remove the phony "break" again
  names(medias) <- mediadata[breaks,2]  # medias is now a named list of mediavectors; names(medias) are the experiment names (chip format)

  CELnames.in.medias <- unlist(get_all_CELnames(names(medias), expnameMap=expnameMap))  # all CELnames which have data in medias
  presentFilter <- colnames(aijdata) %in% CELnames.in.medias  # whether each experiment name from aijdata is present in the media file
  experimentsPresent <- colnames(aijdata)[presentFilter]
  experimentsNotPresent <- colnames(aijdata)[!presentFilter]
  if(all(!presentFilter)) warning("None of the experiments have media data, so all will be run as Complete") 

  medias <- medias[get_chipname(expnameMap, experimentsPresent)]  # Select the appropriate columns to match the columns in aijdata (with duplication) (and sort into the order they appear in aijdata)
  names(medias) <- experimentsPresent

  # All experiments not present in the media file will be treated as "complete" experiments
  medias[experimentsNotPresent] <- NA
  
  return(medias)
  
}

loadKnockoutData <- function(experimentsFile, experimentsFileFormat, commonGenenamesFile, expnameMap) {
  # Returns a named list, where the names of the list are experiment names (CEL format) 
  # and the elements of the list are character vectors of gene names that were knocked out in the experiment
  
  if(is.na(experimentsFile)) return(list())
  
  source("genenameMaps.R")
  genenameMap <- getGenenameMap(commonGenenamesFile)
  
  temp <- read.delim(experimentsFile, fill=T, header=T)
  
  if(experimentsFileFormat=="Hope") {
    
    knockoutrownums <- which(temp$value=="knockout")+1  # rows with "value"=="knockout" are immediately followed by the rows where "value"== the knocked out gene
    knockoutsChip <- strsplit(as.character(temp$value[knockoutrownums]), "\\s+and\\s+")
    knockoutsChip <- lapply(knockoutsChip, function(commonNames) {
      pegNames <- as.vector(genenameMap$fig[match(commonNames, genenameMap$common)])
      if(any(is.na(pegNames))) {
        stop(paste("Failed to find pegNames for some or all of these genes:", paste(commonNames, collapse=", ")))
      }
      return(pegNames)
    })
    names(knockoutsChip) <- temp$experiment_name[knockoutrownums]
    
    # At this point knockoutsChip contains exactly the data we want to return, in the right format, except that it uses chip names for experiments
    
    knockouts <- do.call(c, lapply(1:length(knockoutsChip), function(index) {
      CELnames <- unlist(get_all_CELnames(expnameMap, names(knockoutsChip)[index]))
      knockoutsThisChip <- lapply(CELnames, function(name) knockoutsChip[[index]])  # give each CELname the same knockouts: the ones for the current chip name
      names(knockoutsThisChip) <- CELnames
      return(knockoutsThisChip)
    }))

    return(knockouts)
    
  } else if (experimentsFileFormat=="Carrera") {
    
    perturbs <- strsplit(as.character(temp$Type.of.perturbation), split="\\s+")
    genes <- strsplit(as.character(temp$Gene.perturbated), split="\\s+")
    knockouts <- lapply(1:nrow(temp), function(rownum) {
      return(genes[[rownum]][perturbs[[rownum]]=="D"])  # Return the gene if the corresponding position in perturbs is "D"
    })
    names(knockouts) <- make.names(temp$CEL.file.name)
    return(knockouts)
    
  } else {
    stop(paste("Unrecognized experimentsFileFormat:",experimentsFileFormat))
  }
  
}

loadSupplementData <- function(experimentsFile, experimentsFileFormat) {
  if(experimentsFileFormat=="Hope") {
    stop("Hope-format experiment file not supported for loadSupplementData")
  } else if(experimentsFileFormat=="Carrera") {
    
    
    
  } else {
    stop(paste("Unrecognized experimentsFileFormat:",experimentsFileFormat))
  }
}

ifKBaseConvertToFig <- function(genenames) {
  # Takes a vector of (possibly KBase-format) genenames.
  # Returns a vector of genenames in the same order, with the KBase ones
  #   converted to fig format.
  
  f2p <- read.delim("inputs/f2p.txt", header = F, colClasses = "character")
  figs <- f2p[,3] 
  KBase <- f2p[,2]
  
  return(sapply(genenames, function(genename) {
    index <- match(genename, KBase)
    if(is.na(index)) {
      return(genename)
    } else {
      return(figs[index])
    }
  }))

}

ifFigConvertToKBase <- function(genenames) {
  # Takes a vector of (possibly fig-format) genenames.
  # Returns a vector of genenames in the same order, with the fig-format ones
  #   converted to KBase format.
  
  f2p <- read.delim("inputs/f2p.txt", header = F, colClasses = "character")
  figs <- f2p[,3] 
  KBase <- f2p[,2]
  
  return(sapply(genenames, function(genename) {
    index <- match(genename, figs)
    if(is.na(index)) {
      return(genename)
    } else {
      return(KBase[index])
    }
  }))
  
}
