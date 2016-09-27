# This code implements the method described in "iMRM method description.docx" version 1-15-16.
# See that document for much more explanation, details, definitions, notation, etc.

library(XML)
library(cplexAPI)

source("iMRM CPLEX.R")
source("iMRM media.R")
source("iMRM gene expression and kos.R")

iMRM <- function(SBML, aijs, media, knockouts, w, a, GBP=TRUE, useCarreraMediaMethod=FALSE) {
  # SBML: filename of the SBML model
  # aijs: Named vector of aijs for the given experimental condition, with names 
  #   being gene names.  At minimum, every metabolic gene must be included,
  #   but 'aijs' may contain aijs for more genes, which will be ignored. 
  # media: Named vector of concentrations, with names being compound names, 
  #   indicating the media components present in the given experimental condition.
  #   Any compounds not included in 'media' will be assumed to be not available, 
  #   that is, have concentration 0.  Alternately, media==NA will be taken to 
  #   indicate that the media in the given experimental condition is "complete". 
  # knockouts: vector of gene names which are knocked-out in the given experimental
  #   condition.  May be a vector of length 0, to indicate no knockouts. 
  # w: the "omega" parameter to the MILP (see "iMRM method description.docx")
  # a: the "alpha" parameter to the MILP (see "iMRM method description.docx")
  # GBP: whether to use the GBP iMRM (TRUE) or the RBP iMRM (FALSE)
  # useCarreraMediaMethod: setting this to TRUE forces an alternate implementation
  #   of media conditions that attempts to mimic the handling of media conditions in
  #   the Carrera et al paper.  In this case, 'media' is interpreted completely 
  #   differently: instead of a vector of concentrations, 'media' should be a vector
  #   of one or more compound names to be considered as "supplemental nutrients" in
  #   the given experimental condition.
  #   useCarreraMediaMethod=TRUE is not recommended for general use. 
  # Gene names in 'aijs' and 'knockouts' must exactly match the gene names present in
  #   the SBML model.  Likewise, compound names in 'media' must exactly match the 
  #   compound names present in the SBML model.
  # This function returns the values of all MILP variables as a named vector. 
  
  ### STEP 1 OF PROCEDURE ###
  model <- xmlInternalTreeParse(SBML)
  defineGlobals(model)  # does some preprocessing, defining some lists we'll need later
  cplex <- buildCPLEXprob(model, GBP)  # sets up the MILP

  ### STEPS 2 AND 3 OF PROCEDURE ###
  solveMILP(cplex, model, aijs, 1, 0.5, 1, media, knockouts, useCarreraMediaMethod)
  ComputedMaxBiomass = -getObjValCPLEX(cplex$env, cplex$lp)  # the objective coefficient was -1, so the negative of the objective value will be the biomass flux
  cat("The maximum biomass ignoring expression data is", ComputedMaxBiomass, "\n")
  
  ### STEP 4 OF PROCEDURE ###
  if(ComputedMaxBiomass==0) {
    cat("Since this is 0, we cannot continue, as the objective function is undefined.\n")
    return(NA)
  }
  solveMILP(cplex, model, aijs, w, a, ComputedMaxBiomass, media, knockouts, useCarreraMediaMethod)

  varValues <- getVarValues(cplex)
  closeProbCPLEX(cplex)
  return(varValues)
  
}

iMRM_batch <- function(SBML, aijs, medias, knockouts, w, a, GBP=TRUE, useCarreraMediaMethod=FALSE) {
  # Arguments are the same as the arguments to iMRM() with the following exceptions:
  # aijs should be a matrix where each row is a gene and each column is an experimental
  #   condition, named appropriately. 
  # medias should be a named list, with names being experimental conditions, and each 
  #   element of the list being the 'media' argument to iMRM() for that experiment.
  # knockouts should likewise be a list, with names being experimental conditions, and
  #   each element of the list being the 'knockouts' argument to iMRM() for that experiment.
  # Names of experimental conditions in 'aijs', 'medias', and 'knockouts' should all be
  #   identical.
  
  # Naive implementation - this is conceptually what's going on
  #return(lapply(colnames(aijs), function(experiment.name) {
  #  iMRM(SBML, aijs[,experiment.name], medias[[experiment.name]], knockouts[[experiment.name]], w, a, GBP, useCarreraMediaMethod)
  #}))
  
  # Real implementation - simply a computational optimization of the naive implementation
  
  # Do Step 1 just once, for all experiments
  model <- xmlInternalTreeParse(SBML)
  defineGlobals(model)  # does some preprocessing, defining some lists we'll need later
  cplex <- buildCPLEXprob(model, GBP)  # sets up the MILP
  
  # Then do Steps 2-4 once for each experiment and collect the results
  lapply(colnames(aijs), function(experiment.name) {
    
    # Normal arguments to iMRM()
    aijs <- aijs[,experiment.name]
    media <- medias[[experiment.name]]
    knockouts <- knockouts[[experiment.name]]
    
    ### STEPS 2 AND 3 OF PROCEDURE ###
    solveMILP(cplex, model, aijs, 1, 0.5, 1, media, knockouts, useCarreraMediaMethod)
    ComputedMaxBiomass = -getObjValCPLEX(cplex$env, cplex$lp)  # the objective coefficient was -1, so the negative of the objective value will be the biomass flux
    cat("The maximum biomass ignoring expression data is", ComputedMaxBiomass, "\n")
    
    ### STEP 4 OF PROCEDURE ###
    if(ComputedMaxBiomass==0) {
      cat("Since this is 0, we cannot continue, as the objective function is undefined.\n")
      return(NA)
    }
    solveMILP(cplex, model, aijs, w, a, ComputedMaxBiomass, media, knockouts, useCarreraMediaMethod)
    
    varValues <- getVarValues(cplex)
    closeProbCPLEX(cplex)
    return(varValues)
    
  })
}

# Create data objects, lists etc needed for processing
defineGlobals <- function(model) {
  # List of the XML reaction nodes for non-reversible reactions
  non.rev.rxns <<- getNodeSet(model, "//x:reaction[@reversible='false']", c(x="http://www.sbml.org/sbml/level2"))
  # List of all XML reaction nodes
  reactions <<- getNodeSet(model, "//x:reaction", c(x="http://www.sbml.org/sbml/level2"))
  # List of the XML reaction nodes for reversible reactions
  rev.rxns <<- setdiff(reactions, non.rev.rxns)  # We do it this way because in some SBML documents (e.g. Carrera), reversible reactions are not marked reversible=true but rather have no 'reversible' attribute
  # List of reaction names for each of those three
  rev.rxnIds <<- sapply(rev.rxns, xmlGetAttr, "id")
  non.rev.rxnIds <<- sapply(non.rev.rxns, xmlGetAttr, "id")
  rxnIds <<- sapply(reactions, xmlGetAttr, "id")
  # List of all exchange reactions, and list of their ids
  ex.rxn.indices <- which(grepl("EX_", rxnIds))
  ex.rxns <<- reactions[ex.rxn.indices]
  ex.rxnIds <<- rxnIds[ex.rxn.indices]
  # The biomass reaction (id)
  biomass.rxn <<- getNodeSet(model, "//x:parameter[@id='OBJECTIVE_COEFFICIENT' and @value!='0']", c(x="http://www.sbml.org/sbml/level2"), fun=function(parameternode) {
    return(xmlGetAttr(xmlParent(xmlParent(xmlParent(parameternode))),"id"))  # return the 'id' attribute of the great-grandparent of the parameter node (which is the reaction node)
  })[[1]]
  
  # List of all XML species nodes, and list of their ids
  compounds <<- getNodeSet(model, "//x:species[substring(@id,string-length(@id)-1)!='_b']", c(x="http://www.sbml.org/sbml/level2"))
  compoundIds <<- sapply(compounds, xmlGetAttr, "id")
  
  # List of all gene names present in the SBML
  metabolic.genes <<- unique(unlist(lapply(get.gene.associations(model), get.genenames)))
  
  return(invisible(TRUE))  # success
}

solveMILP <- function(cplex, model, aijs, w, a, maxBiomass, media, knockouts, useCarreraMediaMethod) {

  # Set media constraints
  if(useCarreraMediaMethod) {
    setMedia_CarreraMethod(cplex, media)
  } else {
    setMedia(cplex, media)
  }
  
  # Set objective function and knockout constraints
  if(any(is.na(knockouts))) stop("NA found in knockouts")
  badknockouts <- knockouts[!(knockouts %in% metabolic.genes)]
  if(length(badknockouts)>0) stop(paste("Invalid knockout(s):", paste(badknockouts, collapse=", ")))
  aijs[knockouts] <- NA
  if(cplex$GBP) {
    set.obj.GBP(cplex, model, aijs, w, a, maxBiomass)
    set.knockouts.GBP(cplex, knockouts)
  } else {
    set.obj.and.knockouts.RBP(cplex, model, aijs, knockouts, w, a, maxBiomass)
  }

  # Write the LP (for debugging)
  lpFilename <- paste0("lps/",format(Sys.time(), format="%Y-%m-%d_%H-%M-%S"),".lp")
  writeProbCPLEX(cplex$env, cplex$lp, lpFilename)
  cat("Solving MILP; LP written to", lpFilename, "\n")
  
  # Tell CPLEX to solve the MILP
  mipoptCPLEX(cplex$env, cplex$lp)
  
}

getVarValues <- function(cplex) {
  ncols <- getNumColsCPLEX(cplex$env, cplex$lp)
  values <- getProbVarCPLEX(cplex$env, cplex$lp, 0, ncols-1)
  names(values) <- getColNameCPLEX(cplex$env, cplex$lp, 0, ncols-1)
  return(values)
}
