source("iMRM gene expression and kos.R")

buildCPLEXprob <- function(model, GBP) {
  
  # Open a blank CPLEX problem
  cplex <- openProbCPLEX(pname = "iMRM", ptrtypeENV = "cplex_env", ptrtypePROB = "cplex_prob")
  
  # Add empty rows for Constraint 1, one for each compound - we'll fill in the data in a moment
  newRowsCPLEX(cplex$env, cplex$lp, length(compounds), rhs=NULL, sense=NULL, rnames=paste0("Constraint1_", compoundIds))
  
  # Add FF and RF variables, and data for Constraint 1 rows
  mats <- compute.stoichiometric.mats(reactions)
  addColsCPLEX(cplex$env, cplex$lp, length(reactions), length(mats$matval), NULL, mats$matbeg, mats$matind, mats$matval, lb=rep(0, length(reactions)), ub=rep(1000, length(reactions)), cnames=paste0("FF_", rxnIds))
  rev.mats <- compute.stoichiometric.mats(rev.rxns)
  addColsCPLEX(cplex$env, cplex$lp, length(rev.rxns), length(rev.mats$matval), NULL, rev.mats$matbeg, rev.mats$matind, -rev.mats$matval, lb=rep(0, length(rev.rxns)), ub=rep(1000, length(rev.rxns)), cnames=paste0("RF_", rev.rxnIds))
  
  # Add a beta, delta, and z variable for each reaction
  newColsCPLEX(cplex$env, cplex$lp, 3*length(reactions), lb=NULL, ub=rep(1, 3*length(reactions)), xctype=rep(CPX_BINARY, 3*length(reactions)), cnames=c(paste0("beta_", rxnIds), paste0("delta_", rxnIds), paste0("z_", rxnIds)))
  
  if(GBP) {
    # Add a b and d variable for each metabolic gene
    newColsCPLEX(cplex$env, cplex$lp, 2*length(metabolic.genes), lb=NULL, ub=rep(1, 2*length(metabolic.genes)), xctype=rep(CPX_BINARY, 2*length(metabolic.genes)), cnames=c(paste0("b_", metabolic.genes), paste0("d_", metabolic.genes)))
  }
  
  # Some lists of column indices which we'll need
  rev.ff.col.indices        <<- sapply(rev.rxnIds,     function(name) getColIndexCPLEX(cplex$env, cplex$lp, paste0("FF_",   name)))
  non.rev.ff.col.indices    <<- sapply(non.rev.rxnIds, function(name) getColIndexCPLEX(cplex$env, cplex$lp, paste0("FF_",   name)))
  rev.rf.col.indices        <<- sapply(rev.rxnIds,     function(name) getColIndexCPLEX(cplex$env, cplex$lp, paste0("RF_",   name)))
  rev.beta.col.indices      <<- sapply(rev.rxnIds,     function(name) getColIndexCPLEX(cplex$env, cplex$lp, paste0("beta_", name)))
  non.rev.beta.col.indices  <<- sapply(non.rev.rxnIds, function(name) getColIndexCPLEX(cplex$env, cplex$lp, paste0("beta_", name)))
  rev.delta.col.indices     <<- sapply(rev.rxnIds,     function(name) getColIndexCPLEX(cplex$env, cplex$lp, paste0("delta_",name)))
  non.rev.delta.col.indices <<- sapply(non.rev.rxnIds, function(name) getColIndexCPLEX(cplex$env, cplex$lp, paste0("delta_",name)))
  rev.z.col.indices         <<- sapply(rev.rxnIds,     function(name) getColIndexCPLEX(cplex$env, cplex$lp, paste0("z_",    name)))
  non.rev.z.col.indices     <<- sapply(non.rev.rxnIds, function(name) getColIndexCPLEX(cplex$env, cplex$lp, paste0("z_",    name)))
  ff.col.indices    <<- c(rev.ff.col.indices,    non.rev.ff.col.indices)
  beta.col.indices  <<- c(rev.beta.col.indices,  non.rev.beta.col.indices)
  delta.col.indices <<- c(rev.delta.col.indices, non.rev.delta.col.indices)
  z.col.indices     <<- c(rev.z.col.indices,     non.rev.z.col.indices)
  ex.ff.col.indices <<- sapply(ex.rxnIds, function(name) getColIndexCPLEX(cplex$env, cplex$lp, paste0("FF_", name)))
  ex.rf.col.indices <<- sapply(ex.rxnIds, function(name) getColIndexCPLEX(cplex$env, cplex$lp, paste0("RF_", name)))
  biomass.col.index <<- getColIndexCPLEX(cplex$env, cplex$lp, paste0("FF_", biomass.rxn))
  if(GBP) {
    b.col.indices <<- sapply(metabolic.genes, function(name) getColIndexCPLEX(cplex$env, cplex$lp, paste0("b_", name)))
    d.col.indices <<- sapply(metabolic.genes, function(name) getColIndexCPLEX(cplex$env, cplex$lp, paste0("d_", name)))
  }

  # Add rows for Constraint 2
  # first for non-reversible reactions
  addCPLEXconstraint(cplex, list(non.rev.ff.col.indices, non.rev.beta.col.indices), c(-1,1000), 'G', 0, paste0("Constraint2_", non.rev.rxnIds))
  # now for reversible reactions
  addCPLEXconstraint(cplex, list(rev.ff.col.indices, rev.rf.col.indices, rev.beta.col.indices), c(-1,-1,1000), 'G', 0, paste0("Constraint2_", rev.rxnIds))

  # Add rows for Constraint 3
  # first for non-reversible reactions
  addCPLEXconstraint(cplex, list(non.rev.ff.col.indices, non.rev.delta.col.indices), c(1000,1), 'G', 1, paste0("Constraint3_", non.rev.rxnIds))
  # now for reversible reactions
  addCPLEXconstraint(cplex, list(rev.ff.col.indices, rev.rf.col.indices, rev.delta.col.indices), c(1000,1000,1), 'G', 1, paste0("Constraint3_", rev.rxnIds))

  # Add rows for Constraint 4a
  addCPLEXconstraint(cplex, list(rev.ff.col.indices, rev.z.col.indices), c(1,-1000), 'L', 0, paste0("Constraint4a_", rev.rxnIds))

  # Add rows for Constraint 4b
  addCPLEXconstraint(cplex, list(rev.rf.col.indices, rev.z.col.indices), c(1,1000), 'L', 1000, paste0("Constraint4b_", rev.rxnIds))

  if(GBP) {
    # Add rows for Constraint 5a
    addCPLEXconstraint(cplex, list(beta.col.indices[rxnIds], delta.col.indices[rxnIds]), c(1,1), 'E', 1, paste0("Constraint5a_", rxnIds))

    # Add rows for Constraint 5b
    addCPLEXconstraint(cplex, list(b.col.indices, d.col.indices), c(1,1), 'E', 1, paste0("Constraint5b_", metabolic.genes))

    gene.associations <- get.gene.associations(model)
    sigma.s.col.indices <- mapply(parseAssocStr,  # make several calls to parseAssocStr()
      assocStr = gene.associations,  # each parseAssocStr() call takes an element of gene.associations for its assocStr argument
      MoreArgs=list(evalFunc=psi, cplex=cplex))  # the MoreArgs argument to mapply(), which in this case specifies two arguments to be included in
        # every invocation of parseAssocStr: an evalFunc argument and a cplex argument (both of these are part of the ... argument to parseAssocStr() as well)
    names(sigma.s.col.indices) <- names(gene.associations)  # names are now reaction ids
    
    # Add rows for R
    addCPLEXconstraint(cplex, list(beta.col.indices[rxnIds], sigma.s.col.indices[rxnIds]), c(1,-1), 'L', 0, paste0("ConstraintR_", rxnIds))
    
    # Add rows for Q
    for(gene in metabolic.genes) {
      S <- rxnIds[sapply(gene.associations, function(assoc) gene %in% get.genenames(assoc))]
      addCPLEXconstraint(cplex, as.list(c(b.col.indices[gene], beta.col.indices[S])), c(1, rep(-1, length(S))), 'L', 0, paste0("ConstraintQ_", gene))
    }
    
  }
  
  # We've left the objective coefficients all 0 for now, but before returning we will set the objective sense to "minimize" because that is what it will be throughout
  setObjDirCPLEX(cplex$env, cplex$lp, CPX_MIN)
  
  # and set CPLEX parameters which we'll use throughout
  setIntParmCPLEX(cplex$env, CPXPARAM_MIP_Pool_Capacity, 0)  # turn off all features of the solution pool.  We only care about one (the optimal) solution.
  setIntParmCPLEX(cplex$env, CPXPARAM_Parallel, -1)  # Non-deterministic parallel algorithms allowed.  Gains efficiency and speed at the cost of not guaranteeing to get exactly the same solution every run for the same inputs.
  setDblParmCPLEX(cplex$env, CPXPARAM_TimeLimit, 300)  # Quit after 5 minutes (per call to optimizer) and return best solution found so far.
  
  # use an extra field in 'cplex' to track which type of problem was initialized
  cplex$GBP <- GBP
  
  return(cplex)
  
}

addCPLEXconstraint <- function(cplex, indices, coeffs, sense, rhs, rnames) {
  # Suppose we wish to add the constraint 2a + 3b <= 6
  # indices would be a list, where the first component is the indices of all the 'a' variables,
  #   and likewise the second component is the list of all the 'b' variables
  # coeffs would be c(2,3)
  # sense would be 'L' for <=  (the other options are 'G' for >= and 'E' for ==)
  # rhs would be 6
  # rnames specifies the name of each row, should be a vector the same length as either component of indices
  matind <- as.vector(do.call(rbind, indices))
  numRowsToAdd <- length(indices[[1]])
  addRowsCPLEX(cplex$env, cplex$lp, 0, numRowsToAdd, length(indices)*numRowsToAdd, seq(from=0,by=length(indices),length=numRowsToAdd), matind, rep(coeffs,numRowsToAdd), sense=rep(sense,numRowsToAdd), rhs=rep(rhs,numRowsToAdd), rnames=rnames)
}

compute.stoichiometric.mats <- function(rxns.to.add) {
  # List of lists of reaction coefficients
  coefflists <- sapply(rxns.to.add, function(rxn.node) {
    species <- getRxnSpecies(rxn.node)
    reactantCoeffs <- as.numeric(sapply(species$reactants, xmlGetAttr, "stoichiometry", default=1))  # some SBML models (e.g. Carrera) have compounds with no "stoichiometry" attribute; in these models, this means an implied 1
    productCoeffs <- as.numeric(sapply(species$products, xmlGetAttr, "stoichiometry", default=1))
    return(c( (-1)*reactantCoeffs, productCoeffs))
  })
  matval <- unlist(coefflists)
  lengths <- sapply(coefflists,length)
  matbeg <- c(0, cumsum(lengths)[-length(lengths)])
  
  # List of lists. "rxnspecies" has one entry per reaction, that is a list of species IDs (strings)
  rxnspecies <- sapply(rxns.to.add, getAllRxnSpeciesIds)
  matind <- match(unlist(rxnspecies), compoundIds)-1
  
  return(list(matbeg=matbeg,matind=matind,matval=matval))
}

# returns two lists, one of XML speciesReference nodes for reactants, the other the same for products
getRxnSpecies <- function(rxn.node) {
  subDoc <- xmlDoc(rxn.node)
  reactants <- getNodeSet(subDoc, "//x:listOfReactants/x:speciesReference", c(x="http://www.sbml.org/sbml/level2"))
  products <- suppressWarnings(getNodeSet(subDoc, "//x:listOfProducts/x:speciesReference[substring(@species,string-length(@species)-1)!='_b']", c(x="http://www.sbml.org/sbml/level2")))
  free(subDoc)  # avoid memory leak
  return(list(reactants=reactants, products=products))
}

getAllRxnSpeciesIds <- function(rxn.node) {
  species <- unlist(getRxnSpecies(rxn.node))  # we don't care whether they're reactants or products, just that they're in the right order
  return(sapply(species, xmlGetAttr, "species"))
}

numSigmasCreated <- 0  # global, used by the below function so it knows what sigma to create next
# Creates a new sigma variable in the MILP and returns its column index
create.new.sigma <- function(cplex) {
  numSigmasCreated <<- numSigmasCreated + 1
  sigma.name <- paste0("sigma_", numSigmasCreated)
  newColsCPLEX(cplex$env, cplex$lp, 1, lb=NULL, ub=1, xctype=CPX_BINARY, cnames=sigma.name)
  return(getColIndexCPLEX(cplex$env, cplex$lp, sigma.name))
}

# psi, unlike how it's strictly defined in "iMRM method description.docx", has side-effects.
# Rather than returning an ordered triple (sigma, C, V), it simply returns sigma (namely, 
#   the CPLEX column index of sigma) and adds constraints C and variables V to the MILP.
# Also, here psi works only on "snippets" containing no parentheses.  parseAssocStr() handles
#   the parentheses and the recursion. 
psi <- function(snippet, cplex) {
  
  sigma.s <- create.new.sigma(cplex)
  
  words <- strsplit(snippet, "\\s+")[[1]]  # get a list of 'words' in snippet (separated by whitespace)
  words <- words[which(words!="")]   # remove any empty-string 'words'
  if(length(words)!=0) {  # if it was 0, snippet was the empty string or only whitespace: don't add anything more to the MILP
    
    kind <- which.max(c("and" %in% words, "or" %in% words, TRUE))  
    # note that a snippet cannot have both 'and' and 'or'; this would be ambiguous
    # so kind==1 means this snippet contains 'and's; kind==2, 'or's; and kind==3, neither
    
    colIndices <- unlist(sapply(words, function(word) {  # the unlist removes zero-length vectors
      if(word=="and" || word=="or") return(numeric(0))  # remove 'and' and 'or' words
      numericword <- suppressWarnings(as.numeric(word))
      if(!is.na(numericword)) return(numericword)  # a number: the sigma column index for some other sigma (from a previous evaluation of psi)
      if(word=="Unknown") {
        return(create.new.sigma(cplex))  # sigma column index for a new sigma, which will be unconstrained (this abbreviates an entire recursive call to psi with the argument "Unknown")
      }
      return(b.col.indices[match(word, metabolic.genes)])  # by process of elimination, word must be a gene name. Return its 'b' column index. This also abbreviates an entire recursive call to psi, which would create a new sigma, constrain it equal to this gene's 'b', and return the new sigma. We simply use the gene's 'b' instead. 
    }, USE.NAMES=FALSE))

    if(kind==1) {
      # In iMRM method description.docx, we treated only the case of "a and b" and assumed if there were
      #   more gene-strings, we could recurse down to pairs.  Here we generalize the constraints to treat "a and b and c and ..." directly.
      indices <- as.list(c(sigma.s, colIndices))
      coeffs1 <- c(length(colIndices), rep(-1, length(colIndices)))
      coeffs2 <- c(-1, rep(1, length(colIndices)))
      addCPLEXconstraint(cplex, indices, coeffs1, 'L', 0, rnames=NULL)  # no names for these constraints
      addCPLEXconstraint(cplex, indices, coeffs2, 'L', length(colIndices)-1, rnames=NULL)
    }
    else if(kind==2) {
      indices <- as.list(c(sigma.s, colIndices))
      coeffs1 <- c(1, rep(-1, length(colIndices)))
      coeffs2 <- c(-length(colIndices), rep(1, length(colIndices)))
      addCPLEXconstraint(cplex, indices, coeffs1, 'L', 0, rnames=NULL)
      addCPLEXconstraint(cplex, indices, coeffs2, 'L', 0, rnames=NULL)
    }
    else if(kind==3) {  # should be only one entry in numbers if kind==3
      addCPLEXconstraint(cplex, list(sigma.s, colIndices), c(1, -1), 'E', 0, rnames=NULL)
    }
  
  }
  
  return(sigma.s)
}
