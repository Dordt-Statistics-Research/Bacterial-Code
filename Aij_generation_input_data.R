source("genenameMaps.R")

# Many of the functions in this document take an argument 'data.source'.
# This argument can be one of "Hope Ecoli", "Hope Shewanella", "Carrera Ecoli", "Ishii Ecoli", "Ishii Non-std Ecoli", or "Hope Ecoli Compendium [X]" where [X] identifies a compendium
#   or one of several options for simulated data based on each of those real data sources:
#     "Sim-Uniform " prepended to a data.source, which uses uniform random mixing parameters 
#     "Sim-Fitted " prepended to a data.source, which uses the fitted mixing parameters from the Gibbs sampler
#     "Sim-Unfitted x y" where x and y are numbers bounding the mean gaps in the simulated data
#       unlike the above simulated options, this simulation option isn't based on real data fits at all
#   Examples: "Sim-Uniform Hope Ecoli" or "Sim-Fitted Hope Shewanella" or "Sim-Unfitted 2 3" etc
#   More info about simulated data in the file simulation_from_fits.R

get.expression.data.from.ds <- function(data.source) {
  # data.source: See notes at top of this document
  
  if(data.source == "Hope Ecoli") {
    
    load("inputs/83333.1.Rdata")
    return(as.matrix(temp))
    
  } else if(data.source == "Hope Shewanella") {
    
    load("inputs/211586.9.Rdata")
    return(as.matrix(temp))
    
  } else if(data.source == "Carrera Ecoli") {
    
    data <- read.delim(file = "inputs/Carrera_EcoMAC.txt", header = F, row.names = NULL, colClasses="numeric")
    exptnames <- read.table("inputs/Carrera_exptnames.txt", header=F, row.names=1)[[1]]
    genenames <- read.delim("inputs/Carrera_genenames.txt", header=F, row.names=1, col.names=c("rowname","b","common"))
    colnames(data) <- exptnames
    rownames(data) <- genenames$b
    return(as.matrix(data))
    
  } else if(data.source == "Ishii Ecoli") {
    
    return(as.matrix(read.csv("inputs/transcriptomics_std_nmol_per_gDW.csv", row.names=1)))
    
  } else if(data.source == "Ishii Non-std Ecoli") {
    
    return(as.matrix(read.csv("inputs/transcriptomics_nmol_per_gDW.csv", row.names=1)))
    
  } else if(substr(data.source, 1, 22)=="Hope Ecoli Compendium ") {
    
    compendiumName <- substr(data.source, 23, nchar(data.source))
    data <- read.delim("inputs/83333.1.compendia.txt", header = FALSE, sep = "\t", row.names = NULL)
    return(get.expression.data.from.ds("Hope Ecoli")[,data[,1][data[,2]==compendiumName]])  # only experiments in that compendium
    
  } else if(substr(data.source, 1, 12)=="Sim-Uniform ") {
    
    source("simulation_from_fits.R")
    sim.data.source <- substr(data.source, 13, nchar(data.source))  # data.source for the desired simulated data
    return(as.matrix(get.sim.data(sim.data.source, type="uniform")$expression.data))
    
  } else if(substr(data.source, 1, 11)=="Sim-Fitted ") {
    
    source("simulation_from_fits.R")
    sim.data.source <- substr(data.source, 12, nchar(data.source))  # data.source for the desired simulated data
    return(as.matrix(get.sim.data(sim.data.source, type="fitted")$expression.data))
    
  } else if(substr(data.source, 1, 13)=="Sim-Unfitted ") {
    
    source("simulation_from_fits.R")
    sim.data.source <- get.unfitted.sim.data.source(data.source)
    return(as.matrix(get.sim.data(sim.data.source, type="unfitted")$expression.data))
    
  } else {
    
    stop(paste("Unrecognized data.source:", data.source))
    
  }
  
}

get.operons.from.ds <- function(data.source, operon.source) {
  # data.source: See notes at top of this document
  # operon.source: Only relevant for data.sources "Hope Ecoli" or "Hope Ecoli Compendium [X]"
  #   in which case it should be either "Microbes Online", "RegulonDB All", or "RegulonDB Strong/Confirmed"
  #   For other data.sources, operon.source is ignored
  
  if(data.source == "Hope Ecoli") {
    
    if(operon.source == "Microbes Online") {
      operondata <- read.csv(file = "inputs/83333.1.edited.operons", header = F, col.names = 1:28, fill=T, stringsAsFactors=F)
    } else if(substr(operon.source, 1, 10) == "RegulonDB ") {
      regulondb <- read.delim("inputs/RegulonDB_operons_raw.txt", header = FALSE, row.names = 1, col.names = c("name","startpos","endpos","strand","size","genes","evidence","confidence"), fill=T, comment.char = "#")
      regulondb <- regulondb[regulondb$size>1,]  # per email from AAB on 5/19/16, use only operons of size > 1 in RegulonDB
      if(operon.source == "RegulonDB Strong/Confirmed") regulondb <- regulondb[regulondb$confidence %in% c("Strong", "Confirmed"),]  # Use only "Strong" or "Confirmed" operons; note this excludes e.g. araBAD (at least at the time of this writing)
      else if(operon.source != "RegulonDB All") stop(paste("Unrecognized operon.source:", operon.source))
      operons <- strsplit(as.character(regulondb$genes), ",")
      load("overallmap.Rdata")
      missingcount <- 0
      operons <- lapply(operons, function(operon) {
        operon <- sapply(operon, function(gene) {
          rownum <- match(gene, overallmap$common.x)
          if(is.na(rownum)) rownum <- match(gene, overallmap$common.y)
          if(is.na(rownum)) {missingcount <<- missingcount+1; return(NA)}
          return(overallmap$fig[rownum])
        }, USE.NAMES=FALSE)
        return(operon[!is.na(operon)])
      })
      if(missingcount>1) warning(paste(missingcount, "genes in RegulonDB were not found in genename mappings"))
      return(operons[sapply(operons,length)>0])
    } else {
      stop(paste("Unrecognized operon.source:", operon.source))
    }
    
  } else if(data.source == "Hope Shewanella") {
    
    operondata <- read.csv(file = "inputs/211586.9.7-17-15.edited.operons", header = F, col.names = 1:28, fill=T, stringsAsFactors=F)
  
  } else if(data.source == "Carrera Ecoli") {
    
    operondata <- as.matrix(read.csv(file = "inputs/83333.1.edited.operons", header = F, col.names = 1:28, fill=T, stringsAsFactors = F))
    genenameMap <- getGenenameMap("inputs/E_coli_common_genenames_edited.txt")
    operondata <- as.matrix(apply(operondata, 2, function(column) {
      sapply(column, function(figname) {
        bnumber <- genenameMap[which(genenameMap$fig==figname),"b"]
        if(length(bnumber)==0) return(NA)
        return(bnumber)
      }, USE.NAMES = FALSE)
    }))
    
  } else if(substr(data.source, 1, 6) == "Ishii ") {
    
    return(get.operons.from.ds("Hope Ecoli", "Microbes Online"))
    
  } else if(substr(data.source, 1, 22)=="Hope Ecoli Compendium ") {
    
    return(get.operons.from.ds("Hope Ecoli"))
    
  } else if(substr(data.source, 1, 12)=="Sim-Uniform ") {
    
    source("simulation_from_fits.R")
    sim.data.source <- substr(data.source, 13, nchar(data.source))  # data.source for the desired simulated data
    return(get.sim.data(sim.data.source, type="uniform")$operons)  
    
  } else if(substr(data.source, 1, 11)=="Sim-Fitted ") {
    
    source("simulation_from_fits.R")
    sim.data.source <- substr(data.source, 12, nchar(data.source))  # data.source for the desired simulated data
    return(get.sim.data(sim.data.source, type="fitted")$operons)  
    
  } else if(substr(data.source, 1, 13)=="Sim-Unfitted ") {
    
    source("simulation_from_fits.R")
    sim.data.source <- get.unfitted.sim.data.source(data.source)
    return(get.sim.data(sim.data.source, type="unfitted")$operons)
    
  } else {
    
    stop(paste("Unrecognized data.source:", data.source))
    
  }
  
  return(lapply(1:nrow(operondata), function(i) {
    genes <- unlist(operondata[i, !is.na(operondata[i,]) & operondata[i,]!=""])
    names(genes) <- NULL
    return(genes)
  }))
  
}

get.source.str.from.ds <- function(data.source, id="") {
  # data.source: See notes at top of this document
  # id: an integer, used for keeping apart save files. Unless you specifically need this, keep it at the default "". 
  
  if(data.source == "Hope Ecoli") {
    return(paste0("hope",id,"_ecoli"))
  } else if(data.source == "Hope Shewanella") {
    return(paste0("hope",id,"_shew"))
  } else if(data.source == "Carrera Ecoli") {
    return(paste0("carrera",id,"_ecoli"))
  } else if(data.source == "Ishii Ecoli") {
    return(paste0("ishii",id,"_ecoli"))
  } else if(data.source == "Ishii Non-std Ecoli") {
    return(paste0("ishii",id,"_nonstd_ecoli"))
  } else if(substr(data.source, 1, 22)=="Hope Ecoli Compendium ") {
    compendiumName <- substr(data.source, 23, nchar(data.source))
    return(paste0("hope",id,"_ecoli_comp", compendiumName))
  } else if(substr(data.source, 1, 12)=="Sim-Uniform ") {
    return(paste0("sim_uniform_", get.source.str.from.ds(substr(data.source, 13, nchar(data.source)), id)))
  } else if(substr(data.source, 1, 11)=="Sim-Fitted ") {
    return(paste0("sim_fitted_", get.source.str.from.ds(substr(data.source, 12, nchar(data.source)), id)))
  } else if(substr(data.source, 1, 13)=="Sim-Unfitted ") {
    sim.data.source <- get.unfitted.sim.data.source(data.source)
    return(paste0("sim_unfitted",id,"_",do.call(paste,c(sim.data.source,list(sep="_",collapse="_")))))
  } else {
    stop(paste("Unrecognized data.source:", data.source))
  }
}

get.unfitted.sim.data.source <- function(data.source) {
  second.space <- gregexpr(" ", data.source)[[1]][2]  # index of the second space in data.source
  return(list(numGenes=5000, numExps=1000, gene.centers.mean=10, gene.centers.sd=2, gene.sd.min=0.5, gene.sd.max=2, 
                          meangap.min=as.numeric(substr(data.source, 14, second.space-1)),  # 14 because "Sim-Unfitted " takes 13 characters
                          meangap.max=as.numeric(substr(data.source, second.space+1, nchar(data.source))), 
                          mixing.min=0.2, mixing.max=0.8, onecompPct=0.25))
}
