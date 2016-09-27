verbose <- TRUE

# Gives the hypothetical reaction call for a set of gene calls.
# Gene calls of 1 or 0 indicate gene is 'on' or 'off' respectively.
# Gene calls of any number other than 1 or 0 should only be passed to this function if the reaction does not depend on those genes.
# Gene calls of 'NA' are treated as 'no information known about gene.'
# The return value is 1 or 0 if the given information is enough to call the reaction 'on' or 'off' respectively.
# The return value is NA if the given information is not enough to call the reaction (i.e. the reaction could still go either way based on the eventual value of the NAs.)
# The return value is 0.5 if the reaction could go either way regardless of the amount of information present. (This could happen, for instance, if the reaction's gene association is 'Unknown', or 'gene1 and (Unknown or gene2)' where gene1 is called 1 and gene2 is called 0.)
get.rxn.call <- function(snippet, gene.calls) {
  words <- strsplit(snippet, "\\s+")[[1]]  # get a list of 'words' in snippet (separated by whitespace)
  words <- words[which(words!="")]   # remove any empty-string 'words'
  if(length(words)==0) return(0.5)  # use 0.5 as the rij for an empty string or just whitespace
  kind <- which.max(c("and" %in% words, "or" %in% words, TRUE))  
  # note that a snippet cannot have both 'and' and 'or'; this would be ambiguous
  # so kind==1 means this snippet contains 'and's; kind==2, 'or's; and kind==3, neither
  numbers <- unlist(sapply(words, function(word) {  # the unlist removes zero-length vectors
    if(word=="and" || word=="or") return(numeric(0))  # remove 'and' and 'or' words
    if(word=="Unknown") return(0.5)  # use 0.5 as the call for Unknown gene
    if(word=="NA") return(NA)  
    numericword <- suppressWarnings(as.numeric(word))
    if(!is.na(numericword)) return(numericword)  # a previously calculated expression score
    return(gene.calls[word])  # word must be a gene name
  }, USE.NAMES=FALSE))
  if(kind==1) {
    if(0 %in% numbers)  return(0)  # If there's any 0, it will force the works to 0 regardless of what else is present (i.e. 0 'and' 1 'and' NA 'and' 0.5 = 0).
    if(any(is.na(numbers)))  return(NA)  # No 0 but any NAs gives an NA - we don't know if the NA could be a 0, and if it was, that would change our answer.
    return(min(numbers))  # No 0 or NA, we just take the minimum. 
  }
  else if(kind==2) {
    if(1 %in% numbers)  return(1)  # If there's any 1, it will force the works to 1 regardless of what else is present (i.e. 0 'or' 1 'or' NA 'or' 0.5 = 1). 
    if(any(is.na(numbers)))  return(NA)  # No 1 but any NAs gives an NA - we don't know if the NA could be a 1, and if it was, that would change our answer.
    return(max(numbers))  # No 1 or NA, we just take the maximum. 
  }
  else if(kind==3) {  # should be only one entry in numbers if kind==3
    return(numbers)
  }
}

# Returns a named vector (names are gene names) where the value for each gene is
#   1 if the gene must be on to achieve this set of reaction calls
#   0 if the gene must be off to achieve this set of reaction calls
#   0.5 if the set of reaction calls could be achieved with the gene in either state
# Function errors out if the set of reaction calls cannot be achieved with any set of gene calls
# NOTE: This function makes the assumption that a gene involved in different reactions is in the same state for all of those reactions. This assumption was determined to be problematic at a meeting on 8-12-15.
calculate.gene.calls <- function(model, rxncalls) {
  # rxncalls: A named vector of the reaction calls (0s, 1s, and possibly 0.5s) (names = reaction ids)
  
  gene.associations <- get.gene.associations(model)
  
  # rxnGeneList: A named list (names are reaction ids); each element is a vector of the genenames on which that reaction depends.
  # The list only contains reactions which still can still contribute useful information to determine gene calls. Other reactions are removed as we go. 
  rxnGeneList <- lapply(gene.associations, getAllGenes) 
  names(rxnGeneList) <- names(gene.associations)
  
  # rxnsByGene: A named list (names are genenames); each element is a vector of the reaction ids which depend on that gene.
  # CURRENTLY rxnsByGene is UNUSED
  genenames <- unique(unlist(rxnGeneList))
  rxnsByGene <- sapply(genenames, function(genename) {
    names(rxnGeneList)[which(sapply(rxnGeneList, function(geneList) genename %in% geneList))]
  }, simplify=FALSE, USE.NAMES=TRUE)  # leave it a named list
  
  # Initialize gene.calls
  gene.calls <- rep(NA,length(genenames))
  names(gene.calls) <- genenames
  
  somethingChanged <- TRUE
  while(somethingChanged) {
    somethingChanged <- FALSE
    rxnsToRemove <- c()  # list of reactions that no longer can contribute any useful information regarding gene calls
    for(index in 1:length(rxnGeneList)) {
      geneList <- rxnGeneList[[index]]
      rxnid <- names(rxnGeneList)[index]
      if(verbose) cat("Processing reaction",rxnid,"with association string",gene.associations[rxnid], "\n")
      if(verbose) cat("Gene calls",paste(gene.calls[geneList]), "\n")
      if(verbose) cat("Reaction call",rxncalls[rxnid], "\n")
      uncalledGenes <- geneList[which(is.na(gene.calls[geneList]))]
      numUncalledGenes <- length(uncalledGenes)
      if(numUncalledGenes==0) {
        
        if(verbose) cat("No uncalled genes.\n")
        rij <- parseAssocStr(gene.associations[rxnid], get.rxn.call, gene.calls)
        if(rij != rxncalls[rxnid] && rij != 0.5) stop(paste("Contradictory reaction call: reaction", rxnid))
        else rxnsToRemove <- c(rxnsToRemove, index)  # Reaction has no more useful information for determining gene calls
        
      } else {
        
        for(gene in uncalledGenes) {
          
          if(verbose) cat("Processing uncalled gene",gene, "\n")
          
          # Calculate the rij with the gene called 0 and with the gene called 1.
          # If the result is a number, then that gene call forces the reaction to have that rij, regardless of which way other uncalled genes are called.
          # If the result is NA, then that reaction could still take on multiple rijs even with that gene call. 
          gene.calls[gene] <- 0
          rijWith0 <- parseAssocStr(gene.associations[rxnid], get.rxn.call, gene.calls)
          gene.calls[gene] <- 1
          rijWith1 <- parseAssocStr(gene.associations[rxnid], get.rxn.call, gene.calls)
          gene.calls[gene] <- NA
          
          if(is.na(rijWith0) && is.na(rijWith1)) {
            if(verbose) cat("Both rijs NA\n")  # do nothing - gene does not (yet) fully determine reaction call either way
          } else if(is.na(rijWith0) && rijWith1==rxncalls[rxnid]) {
            if(verbose) cat("Gene=1 would be correct, but Gene=0 gives NA\n")  # still do nothing. Although calling the gene 1 would cause the reaction to line up with its call, it's not necessarily necessary to call this gene 1 to get that result, and perhaps a different reaction needs this gene 0
          } else if(is.na(rijWith1) && rijWith0==rxncalls[rxnid]) {
            if(verbose) cat("Gene=0 would be correct, but Gene=1 gives NA\n")  # dual of the above case
          } else if(is.na(rijWith0) && rijWith1==0.5) {
            if(verbose) cat("Gene=1 gives 0.5, but we need to see how Gene=0 resolves\n")  # how the NA resolves still makes a difference in this case, so can't do anything
          } else if(is.na(rijWith1) && rijWith0==0.5) {
            if(verbose) cat("Gene=0 gives 0.5, but we need to see how Gene=1 resolves\n")  # dual of the above case
          } else if(is.na(rijWith0)) { # then, by process of elimination, rijWith1 is opposite of the reaction call
            if(verbose) cat("Gene=1 gives wrong, so must set Gene=0\n")
            gene.calls[gene] <- 0  # Set for real (the change sticks) - we have determined that this gene cannot be called 1, so we must call it 0.
            somethingChanged <- TRUE
          } else if(is.na(rijWith1)) { # then, by process of elimination, rijWith0 is opposite of the reaction call
            if(verbose) cat("Gene=0 gives wrong, so must set Gene=1\n")
            gene.calls[gene] <- 1  # dual of the above case
            somethingChanged <- TRUE
          } else {  # neither rij is NA
            if(rijWith0==rijWith1) {
              if(rijWith0==rxncalls[rxnid] || rijWith0==0.5) {
                # With the actual calls made so far, reaction no longer depends on this gene.  (E.g., maybe the reaction is associated with "gene1 and (gene2 or gene3)" and we already have gene3=1, then the reaction no longer depends on gene2.)
                if(verbose) cat("Both rijs are either correct or 0.5. Reaction does not depend on this gene.\n")
                rxnGeneList[index] <- geneList[which(geneList!=gene)]  # Remove this gene from the geneList for this reaction
                rxnsByGene[gene] <- rxnsByGene[gene][which(rxnsByGene[gene]!=rxnid)]  # Remove this reaction from the rxnsByGene for this gene
                somethingChanged <- TRUE
              } else {
                stop(paste("Impossible to get reaction to satisfy its call: reaction", rxnid))
              }
            } else if (rijWith0==0.5) {
              if(rijWith1==rxncalls[rxnid]) {
                # This reaction cannot be used to call the gene. If the gene were called 1, this reaction would be satisfied; but if it were called 0, it could also be satisfied.
                # So, the reaction no longer depends on this gene.
                if(verbose) cat("Gene=1 gives correct, Gene=0 gives 0.5. Reaction can't decide this gene.\n")
                rxnGeneList[index] <- geneList[which(geneList!=gene)]  # Remove this gene from the geneList for this reaction
                rxnsByGene[gene] <- rxnsByGene[gene][which(rxnsByGene[gene]!=rxnid)]  # Remove this reaction from the rxnsByGene for this gene
                somethingChanged <- TRUE
              } else {
                if(verbose) cat("Gene=1 gives wrong, but Gene=0 gives 0.5; must set Gene=0\n")
                gene.calls[gene] <- 0  # Set for real (the change sticks) - we have determined that this gene cannot be called 1, so we must call it 0.
                somethingChanged <- TRUE
              }
            } else if (rijWith1==0.5) {
              # dual of the rijWith0==0.5 case
              if(rijWith0==rxncalls[rxnid]) {
                if(verbose) cat("Gene=0 gives correct, Gene=1 gives 0.5. Reaction can't decide this gene.\n")
                rxnGeneList[index] <- geneList[which(geneList!=gene)]
                rxnsByGene[gene] <- rxnsByGene[gene][which(rxnsByGene[gene]!=rxnid)]
                somethingChanged <- TRUE
              } else {
                if(verbose) cat("Gene=0 gives wrong, but Gene=1 gives 0.5; must set Gene=1\n")
                gene.calls[gene] <- 1
                somethingChanged <- TRUE
              }
            } else {
              # Both rijs are either 0 or 1, and do not match
              if (rijWith0==rxncalls[rxnid]) {
                if(verbose) cat("Gene=0 is required\n")
                gene.calls[gene] <- 0  # Set for real (the change sticks) - this reaction requires this gene to be called 0 in order to match its reaction call
                somethingChanged <- TRUE
              } else if (rijWith1==rxncalls[rxnid]) {
                if(verbose) cat("Gene=1 is required\n")
                gene.calls[gene] <- 1  # Likewise with 1
                somethingChanged <- TRUE
              }
            }
          }
        }
        
      }
    }
    
    if(length(rxnsToRemove) > 0) {
      if(verbose) {
        cat("Removing the following reactions which have no information to contribute:\n")
        for(rxn in rxnsToRemove) {
          cat("Reaction",rxn,"with association",gene.associations[rxn],"and gene calls",gene.calls[[rxnGeneList[rxn]]], "\n")
        }
      }
      rxnGeneList <- rxnGeneList[-rxnsToRemove]
      somethingChanged <- TRUE
    }
    
  }
  
  if(length(rxnGeneList)==0) {
    # We finished; all reactions have contributed all the information they can
    gene.calls[which(is.na(gene.calls))] <- 0.5  # Any remaining uncalled genes could go either way without affecting the outcome, so assign them all 0.5
    return(gene.calls)
  } else {
    stop(paste("Got stuck; need better coding. Contact Craig Disselkoen for help."))
  }
  
}

# Returns a named vector (names are gene names) where the value for each gene is
#   1 if the gene must be on to achieve this set of reaction calls
#   0 if the gene state cannot be determined from the reaction calls
calculate.simple.gene.calls <- function(model, rxncalls) {
  # rxncalls: A named vector of the reaction calls (0s, 1s, and possibly 0.5s) (names = reaction ids)
  
  gene.associations <- get.gene.associations(model)
  
  # rxnGeneList: A named list (names are reaction ids); each element is a vector of the genenames on which that reaction depends.
  rxnGeneList <- lapply(gene.associations, getAllGenes) 
  names(rxnGeneList) <- names(gene.associations)
  genenames <- unique(unlist(rxnGeneList))
  # We're only going to work with reactions that were called 'on'
  rxnGeneList <- rxnGeneList[rxncalls[names(rxnGeneList)]==1]
  gene.associations <- gene.associations[names(rxnGeneList)]
  
  # Initialize gene.calls
  gene.calls <- rep(0,length(genenames))
  names(gene.calls) <- genenames
  
  # dep.genes: For each reaction, a character string with a space-separated list of all genes that have to for sure be on, given that reaction is on 
  # For example, the association string "gene1 and gene2" would give "gene1 gene2", as both genes have to be on since the reaction is on
  # But the association string "gene1 or gene2" would give "", as neither gene necessarily has to be on
  # The association string "gene1 and (gene2 or (gene3 and gene4 and gene5))" would give "gene1", as gene1 has to be on, but none of the others necessarily have to be.
  dep.genes <- sapply(gene.associations, parseAssocStr, function(snippet, moreargs=NA) {
    words <- strsplit(snippet, "\\s+")[[1]]  # get a list of 'words' in snippet (separated by whitespace)
    words <- words[which(words!="")]   # remove any empty-string 'words'
    if(length(words)==0) return("")  # No genes present, just the empty string or whitespace
    kind <- which.max(c("and" %in% words, "or" %in% words, TRUE))  
    # note that a snippet cannot have both 'and' and 'or'; this would be ambiguous
    # so kind==1 means this snippet contains 'and's; kind==2, 'or's; and kind==3, neither
    genes <- words[words!="and" & words!="or" & words!="Unknown"]   # remove 'and's, 'or's, and 'Unknown's
    if(kind==1) return(paste(genes, collapse = " "))  # an 'and' requires all genes to be on
    else if(kind==2) return("")  # an 'or' does not require any of the involved genes to be on
    else if(kind==3) return(paste(genes, collapse = " "))  # No 'and' and no 'or', just pass the previous results on.  (We could possibly be looking at several genes from a previous recursion, or just one gene on which the reaction depends.)
  })
  
  # convert the space-separated lists into actual lists
  dep.genes <- sapply(dep.genes, function(sp.sep.list) strsplit(sp.sep.list, "\\s+")[[1]] )
  
  gene.calls[unique(unlist(dep.genes))] <- 1
  return(gene.calls)
  
}

outputGeneCalls <- function(cplex, experiment.name) {
  fluxes <- getVarValues(cplex)
  non.rev.zeroes <- rep(0, length(non.rev.rxns))
  names(non.rev.zeroes) <- paste("RF_",non.rev.rxnIds,sep="")
  fluxes <- c(fluxes, non.rev.zeroes)  # Add a 0 reverse flux for non-reversible reactions
  rxncalls <- 1*( fluxes[paste("FF_",rxnIds,sep="")]>0 | fluxes[paste("RF_",rxnIds,sep="")]>0 )
  names(rxncalls) <- rxnIds
  gene.calls <- calculate.simple.gene.calls(model, rxncalls)
  filename <- paste("results/GeneCalls/",experiment.name,".txt",sep="")
  write.table(format(gene.calls), filename, sep="\t\t", quote=F, col.names=F)
  cat("Gene calls written to", filename, "\n")
}
