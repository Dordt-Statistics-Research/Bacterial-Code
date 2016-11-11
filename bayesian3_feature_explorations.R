setwd("/home/jason/Documents/School/College/Fall Junior/STAT RESEARCH/MyWD/Bacterial-Code")
source("expnameMaps.R")
source("Aij_utility_funcs.R")
source("E_coli_v4_Build_6_6-23-15_parser.R")
source("/home/jason/Documents/School/College/Fall Junior/STAT RESEARCH/MyWD/Bacterial-Code/expnameMaps.R")
cel_chip_map = get.expname.map("/home/jason/Documents/School/College/Fall Junior/STAT RESEARCH/MyWD/Bacterial-Code/inputs/Ecoli_Cel_to_chip_fromClaire_edited.tab")
load("/home/jason/Documents/School/College/Fall Junior/STAT RESEARCH/MyWD/Bacterial-Code/bayesian3_exploratory_inputs/aijResults.Rdata")
multiaijs <- aijResults$aijs

#Returns a character vector of chip names for experiments where all functions in FUNs return true given the properties of the experiment
#FUNs should be a function or list of functions each with five feature parameters (name, value, units, type, url) that returns a boolean
get.by.features <- function(FUNs) {
  if (typeof(FUNs)=="closure") {
    #gets the subset of the data based FUN then takes the union with an empty vector to generate a set (i.e. to remove duplicates)
    #returns the set of names of experiments as a character vector
    return(union(subset(feature_data, FUNs(feature_name, value, feature_units, feature_type, feature_url))[,1],c()))
  } else { #(typeof(FUNs)=="list")
    set = c()
    for (i in 1:length(FUNs)) {
      if (i==1) set <- get.by.features(FUNs[[i]])
      else set <- intersect(set, get.by.features(FUNs[[i]]))
    }
    return(set)
  }
}

#Given a matrix of lists, a matrix of the same dimensions will be returned with each entry 
#representing the length of the list in the same index of the parameter matrix
list_to_length <- function(list_matrix) {
  nrows <- dim(list_matrix)[1]
  ncols <- dim(list_matrix)[2]
  matrix <- matrix(nrow=nrows, ncol=ncols)
  for (r in 1:nrows){
    for (c in 1:(ncols)){
      matrix[r,c] <- length(list_matrix[[r,c]])
    }
  }
  return(matrix)
}

#Prints a table to the console with a fisher p-value
print.table.and.fisher <- function(aijs,rowFUNs,colFUNs,rowfun_names,colfun_names,default_name="Other") {
  list_matrix <- get.table(aijs,rowFUNs,colFUNs)
  matrix <- list_to_length(list_matrix)
  
  nrows <- length(rowFUNs)+1
  ncols <- length(colFUNs)+1

  p <- fisher.test(matrix,hybrid=TRUE,simulate.p.value=TRUE)$p
  
  frame <- as.data.frame(matrix)
  frame[,ncols+1] <- 0 #add a column of 0's
  frame[nrows+1,] <- 0 #add a row of 0's

  for (r in 1:nrows){
    frame[r,ncols+1] <- sum(frame[r,])
  }
  for (c in 1:(ncols+1)){
    frame[nrows+1,c] <- sum(frame[,c])
  }
  rownames(frame) <- c(rowfun_names, default_name, "Totals")
  colnames(frame) <- c(colfun_names, default_name, "Totals")
  print(frame)
  writeLines(paste0("Fisher p-value: ",p,"\n"))
  return(table)
}

#Given a data.frame (I think a matrix should also work) of aijs where the column names are CEL names and row names are ped IDs 
#and two lists of boolean functions, a matrix of size (length(rowFUNs)+1, length(colFUNs)+1) will be returned.
#Index [r,c] is a list of CEL names of experiments where both rowFUNs[[r]] and colFUNs[[c]] return true.
#where the 'extra functions' for rows and columns is true for all experiments for which no other function was true and false otherwise.
#To get desired behavior, it is the caller's responsibility to choose mutually exclusive functions within rowFUNs and within colFUNs
#rowFUNs and colFUNs should be lists of functions of two parameters (aijs, cel) that returns a boolean
#aijs will be a named numeric vector and cel will be the CEL name for the experiment (the CEL name is unique to each experiment)
get.table <- function(aijs,rowFUNs,colFUNs) {
  colnames <- colnames(aijs)
  rownames <- rownames(aijs)
  m <- matrix(list(), nrow=length(rowFUNs)+1, ncol=length(colFUNs)+1)
  row_solutions <- rep(list(c()),length(rowFUNs)+1)
  col_solutions <- rep(list(c()),length(colFUNs)+1)
  all_row_solutions <- c()
  all_col_solutions <- c()
  for (col in 1:length(aijs[1,])) {
    experiment_name <- colnames[col]
    experiment_column <- aijs[,col]
    names(experiment_column) <- rownames
    for (rowfun_index in 1:length(rowFUNs)) {
      if (rowFUNs[[rowfun_index]](experiment_column, experiment_name)) {
        row_solutions[[rowfun_index]] <- c(row_solutions[[rowfun_index]], experiment_name)
        all_row_solutions <- c(all_row_solutions, experiment_name)
      }
    }
    for (colfun_index in 1:length(colFUNs)) {
      if (colFUNs[[colfun_index]](experiment_column, experiment_name)) {
        col_solutions[[colfun_index]] <- c(col_solutions[[colfun_index]], experiment_name)
        all_col_solutions <- c(all_col_solutions, experiment_name)
      }
    }
  }
  row_solutions[[length(rowFUNs)+1]] <- setdiff(colnames, all_row_solutions)
  col_solutions[[length(colFUNs)+1]] <- setdiff(colnames, all_col_solutions)
  for (r in 1:(length(rowFUNs)+1)) {
    for (c in 1:(length(colFUNs)+1)) {
      m[r,c] <- list(intersect(row_solutions[[r]], col_solutions[[c]]))
    }
  }
  return(m)
}

# Some general forms for classification functions
in.chips <- function(chip_set, aijs, cel) {
  chip <- get_chipname(cel_chip_map, gsub("\\.CEL\\.gz","",cel))
  return(chip %in% chip_set);
}

active.gene <- function(gene, aijs, cel, threshold=0.8) {
  return(aijs[com_to_fig(gene)] >= threshold)
}

inactive.gene <- function(gene, aijs, cel, threshold=0.2) {
  return(aijs[com_to_fig(gene)] <= threshold)
}

activity.gene <- function(gene, aijs, cel, thresholds=c(0.2,0.8)) {
  return((thresholds[1] < aijs[com_to_fig(gene)]) & (aijs[com_to_fig(gene)] <= thresholds[2]))
}

conflict.tf.tg <- function(tf, tg, aijs, cel, threshold_lo=0.2, threshold_hi=0.8) {
  return(aijs[com_to_fig(tf)] >= threshold_hi & aijs[com_to_fig(tg)] <= threshold_lo)
}

# active functions for specific genes
active.araC <- function(aijs, cel) {return(active.gene("araC", aijs, cel))}
active.araB <- function(aijs, cel) {return(active.gene("araB", aijs, cel))}
active.araA <- function(aijs, cel) {return(active.gene("araA", aijs, cel))}
active.araD <- function(aijs, cel) {return(active.gene("araD", aijs, cel))}

# inactive functions for specific genes
inactive.araC <- function(aijs, cel) {return(inactive.gene("araC", aijs, cel))}
inactive.araB <- function(aijs, cel) {return(inactive.gene("araB", aijs, cel))}
inactive.araA <- function(aijs, cel) {return(inactive.gene("araA", aijs, cel))}
inactive.araD <- function(aijs, cel) {return(inactive.gene("araD", aijs, cel))}

# mid activity functions for specific genes
activity.4.6.araC <- function(aijs, cel) {return(activity.gene("araC", aijs, cel, c(0.4,0.6)))}

# functions based on feature constraints
chips.strain.bw25113 <- get.by.features(function(nam,val,uts,typ,url) {nam=="strain" & val=="BW25113"})
chips.perturbation_gene.ccdB <- get.by.features(function(nam,val,uts,typ,url) {nam=="perturbation_gene" & val=="ccdB"})
chips.arabinose <- get.by.features(function(nam,val,uts,typ,url) {nam=="arabinose"})
chips.arabinose.high <- get.by.features(function(nam,val,uts,typ,url) {nam=="arabinose" & (val=="6.66" | val=="6.6609")})
chips.arabinose.mid <- get.by.features(function(nam,val,uts,typ,url) {nam=="arabinose" & val=="8.3261"})
chips.arabinose.low <- get.by.features(function(nam,val,uts,typ,url) {nam=="arabinose" & (val=="13.32" | val=="16.6522")})

chips.strain.bw25113.perturbation_gene.ccdB <- intersect(chips.strain.bw25113, chips.perturbation_gene.ccdB)
chips.arabinose.perturbation_gene.ccdB <- intersect(chips.arabinose, chips.perturbation_gene.ccdB)
    
feature.strain.bw25113 <- function(aijs, cel) {return(in.chips(chips.strain.bw25113, aijs, cel))}
feature.strain.bw25113.perturbation_gene.ccdB <- function(aijs, cel) {return(in.chips(chips.strain.bw25113.perturbation_gene.ccdB, aijs, cel))}
feature.arabinose <- function(aijs,cel) {return(in.chips(chips.arabinose, aijs, cel))}
feature.arabinose.perturbation_gene.ccdB <- function(aijs,cel) {return(in.chips(chips.arabinose.perturbation_gene.ccdB, aijs, cel))}
feature.arabinose.high <- function(aijs,cel) {return(in.chips(chips.arabinose.high, aijs, cel))}
feature.arabinose.mid <- function(aijs,cel) {return(in.chips(chips.arabinose.mid, aijs, cel))}
feature.arabinose.low <- function(aijs,cel) {return(in.chips(chips.arabinose.low, aijs, cel))}

# conflict functions for tf and tg
conflict.araC.araA <- function(aijs, cel) {return(conflict.tf.tg("araC","araA",aijs,cel))}

# test functions for constant boolean value
fun.true <- function(aijs, cel){return(TRUE)}
fun.false <- function(aijs, cel){return(FALSE)}


test.feature.strain.bw25113.activty.araB.araA.araD.rhaD.rhaA.rhaB.lacZ.hsdR.rph <- function() {
  genes <- c("araB", "araA", "araD", "rhaD", "rhaA", "rhaB", "lacZ", "hsdR", "rph")
  for (gene in genes) {
    thr = 0.2
    print.table.and.fisher(multiaijs, 
                           c(
                             function(aijs, cel) {return(inactive.gene(gene, aijs, cel, threshold=thr))}
                           ), 
                           c(
                             feature.strain.bw25113
                           ),
                           c(paste0(gene, " is inactive (aij<=", thr, ")")), 
                           c("In BW25113 Strain")
    )
  }
  for (gene in genes) {
    thr = 0.8
    print.table.and.fisher(multiaijs, 
                           c(
                             function(aijs, cel) {return(active.gene(gene, aijs, cel, threshold=thr))}
                           ), 
                           c(
                             feature.strain.bw25113
                           ),
                           c(paste0(gene, " is active (aij>=", thr, ")")), 
                           c("In BW25113 Strain")
    )
  }
  for (gene in genes) {
    thr = 0.003
    print.table.and.fisher(multiaijs, 
                           c(
                             function(aijs, cel) {return(inactive.gene(gene, aijs, cel, threshold=thr))}
                           ), 
                           c(
                             feature.strain.bw25113
                           ),
                           c(paste0(gene, " is inactive (aij<=", thr, ")")), 
                           c("In BW25113 Strain")
    )
  }
  for (gene in genes) {
    thr = 0.997
    print.table.and.fisher(multiaijs, 
                           c(
                             function(aijs, cel) {return(active.gene(gene, aijs, cel, threshold=thr))}
                           ), 
                           c(
                             feature.strain.bw25113
                           ),
                           c(paste0(gene, " is active (aij>=", thr, ")")), 
                           c("In BW25113 Strain")
    )
  }
}

test.feature.strain.bw25113.activty.2way.araB.araA.araD.rhaD.rhaA.rhaB.lacZ.hsdR.rph <- function() {
  genes <- c("araB", "araA", "araD", "rhaD", "rhaA", "rhaB", "lacZ", "hsdR", "rph")
  for (gene in genes) {
    thr_act = 0.997
    thr_inact = 0.003
    print.table.and.fisher(multiaijs, 
                           c(
                             function(aijs, cel) {return(active.gene(gene, aijs, cel, threshold=thr_act))},
                             function(aijs, cel) {return(inactive.gene(gene, aijs, cel, threshold=thr_inact))}
                           ), 
                           c(
                             feature.strain.bw25113
                           ),
                           c(
                             paste0(gene, " is active (aij>=", thr_act, ")"),
                             paste0(gene, " is inactive (aij<=", thr_inact, ")")
                             ), 
                           c("In BW25113 Strain")
    )
  }
}

test.feature.arabinose.perturbation_gene.ccdB.activity.araC <- function() {
  print.table.and.fisher(multiaijs, 
                         c(
                           function(aijs, cel) {return(active.gene("araC", aijs, cel, threshold=0.8))}
                           ), 
                         c(
                           feature.arabinose.perturbation_gene.ccdB
                           ),
                         c("araC is active (aij>=0.8)"), 
                         c("Arabinose Present and ccdB Perturbation Gene")
                         )
  print.table.and.fisher(multiaijs, 
                         c(
                           function(aijs, cel) {return(inactive.gene("araC", aijs, cel, threshold=0.2))}
                         ), 
                         c(
                           feature.arabinose.perturbation_gene.ccdB
                         ),
                         c("araC is active (aij<=0.2)"), 
                         c("Arabinose Present and ccdB Perturbation Gene")
  )
  print.table.and.fisher(multiaijs, 
                         c(
                           function(aijs, cel) {return(active.gene("araC", aijs, cel, threshold=0.9))},
                           function(aijs, cel) {return(inactive.gene("araC", aijs, cel, threshold=0.15))}
                         ), 
                         c(
                           feature.arabinose.perturbation_gene.ccdB
                         ),
                         c("araC is active (aij>=0.9)", "araC is inactive (aij<=0.15)"), 
                         c("Arabinose Present and ccdB Perturbation Gene")
  )
}

plot.distributions <- function(aijs, gene, FUN, condition) {
  cel_names <- colnames(aijs)
  experiments <- aijs[com_to_fig(gene),]
  bools <- vector(length=length(cel_names))
  for (i in 1:length(cel_names)) {
    if (FUN(aijs[,i], cel_names[i])) {
      bools[i] <- TRUE
    }
  }
  true <- experiments[bools]
  false <- experiments[!bools]
  minimum <- min(experiments)
  maximum <- max(experiments)
  
  layout(matrix(c(1,2,3,3),ncol=2,byrow=TRUE),heights=c(0.8,0.2))
  plot(density(true), main=paste("Activity Density for",gene), xlim=c(0,1), col="darkred")
  lines(density(false), col="darkblue")
  boxplot(true, false, main=paste("Activity Levels for",gene), col=c("darkred","darkblue"))
  oldmar<-par(mar=c(1,1,1,1))
  plot.new() 
  legend(0.45,1,c(condition, "Other"),fill=c("darkred","darkblue")) 
  par(oldmar)
  
}

plot.feature.strain.bw25113.inactivty.araB.araA.araD.rhaD.rhaA.rhaB.lacZ.hsdR.rph <- function() {
  for (gene in c("araB", "araA", "araD", "rhaD", "rhaA", "rhaB", "lacZ", "hsdR", "rph")) {
    plot.distributions(multiaijs, gene, feature.strain.bw25113, "BW25113 strain")
  }
}

plot.feature.arabinose.3levels <- function() {
  print.table.and.fisher(multiaijs,
                         c(
                           function(aijs, cel) {return(active.gene("araC", aijs, cel, threshold=0.9))},
                           function(aijs, cel) {return(inactive.gene("araC", aijs, cel, threshold=0.15))}
                         ),
                         c(
                           feature.arabinose.high,
                           feature.arabinose.mid,
                           feature.arabinose.low
                         ),
                         c("araC is active (aij>=0.9)", "araC is inactive (aij<=0.15)"),
                         c("High Arabinose Content","Medium Arabinose Content","Low Arabinose Content")
  )
}




# MAIN CODE TO RUN
    # test.feature.arabinose.perturbation_gene.ccdB.activity.araC()
    # test.feature.strain.bw25113.activty.2way.araB.araA.araD.rhaD.rhaA.rhaB.lacZ.hsdR.rph()
    # test.feature.strain.bw25113.activty.araB.araA.araD.rhaD.rhaA.rhaB.lacZ.hsdR.rph()
    # plot.distributions(multiaijs, "araB", feature.strain.bw25113)
    # plot.feature.strain.bw25113.inactivty.araB.araA.araD.rhaD.rhaA.rhaB.lacZ.hsdR.rph()
    plot.feature.arabinose.3levels()


