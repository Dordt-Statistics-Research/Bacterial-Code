source("Aij_utility_funcs.R")
source("valContext.R")
context <- get.wrapped.aij.context.from.ds("Hope Ecoli")
eijs <- get.expression.data.outliersNA(context)
aijs <- MultiMM_NBF(context)
load("nonCompleteExpts_hope_ecoli.Rdata")
source("expnameMaps.R")
expnameMap <- get.expname.map("inputs/Ecoli_Cel_to_chip_fromClaire_edited.tab")

is.complete <- function(expname) !(get_chipname(expnameMap, expname) %in% nonCompleteExpts)

myf <- function(numIn, numAc) fisher.test(matrix(c(numIn, numAc, length(offexps)-numIn, length(onexps)-numAc), nrow=2, ncol=2))$p

efd <- read.delim("inputs/E_coli_v4_Build_6_6-23-15.experiment_feature_descriptions", fill=T, header=T)

get.counts <- function(feature) {
  inactive <- ifelse(length(offexps)==0, 0, sum(sapply(offexps, function(offexp) {
    feature %in% efd[efd$experiment_name==get_chipname(expnameMap, offexp),]$feature_name
  })))
  active <- ifelse(length(onexps)==0, 0, sum(sapply(onexps, function(onexp) {
    feature %in% efd[efd$experiment_name==get_chipname(expnameMap, onexp),]$feature_name
  })))
  unknown <- ifelse(length(unkexps)==0, 0, sum(sapply(unkexps, function(unkexp) {
    feature %in% efd[efd$experiment_name==get_chipname(expnameMap, unkexp),]$feature_name
  })))
  return(c(inactive, active, unknown))
}

feature.p <- function(feature) {
  counts <- get.counts(feature)
  return(myf(counts[1],counts[2]))
}

contingency.table <- function(feature) {
  print(feature.p(feature))
  counts <- get.counts(feature)
  numOn <- length(onexps); numOff <- length(offexps); numUnk <- length(unkexps)
  numTotal <- numOn + numOff + numUnk
  matrix(c(counts[1], counts[2], counts[3], sum(counts), 
         numOff-counts[1], numOn-counts[2], numUnk-counts[3], numTotal-sum(counts),
         numOff, numOn, numUnk, numTotal),
         nrow=4, ncol=3, dimnames = list(c("Inactive", "Active", "Unknown", "Total"), c("Y", "N", "Total")))
}

all.tests <- function(pegnum, aijs) {
  offexps <<- colnames(remove_CEL_gz(eijs))[aijs[hope.gene(pegnum),]<0.2]
  onexps <<- colnames(remove_CEL_gz(eijs))[aijs[hope.gene(pegnum),]>0.8]
  unkexps <<- colnames(remove_CEL_gz(eijs))[aijs[hope.gene(pegnum),]>=0.2 & aijs[hope.gene(pegnum),]<=0.8]
  
  if(length(offexps)==0) {
    cat("No inactive experiments\n")
  } else {
    numComplete_In <- sum(sapply(offexps, is.complete))
    cat("Inactive experiments which are complete:", numComplete_In, "of", length(offexps), "\n")
  }
  if(length(onexps)==0) {
    cat("No active experiments\n")
  } else {
    numComplete_Ac <- sum(sapply(onexps, is.complete))
    cat("Active experiments which are complete:", numComplete_Ac, "of", length(onexps), "\n")
  }
  if(length(unkexps)==0) {
    cat("No unknown experiments\n")
  } else {
    numComplete_Unk <- sum(sapply(unkexps, is.complete))
    cat("Unknown experiments which are complete:", numComplete_Unk, "of", length(unkexps), "\n")
  }
  
  all.features <- unique(efd$feature_name)
  all.feature.ps <<- sapply(all.features, feature.p)
  names(all.feature.ps) <<- all.features
  complete.p <- myf(numComplete_In, numComplete_Ac)
  names(complete.p) <- "Complete"
  all.feature.ps <<- c(all.feature.ps, complete.p)
  all.feature.ps <<- sort(all.feature.ps)
}
