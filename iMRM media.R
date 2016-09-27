setMedia <- function(cplex, media) {
  # As always, see comments in "iMRM method description.docx" for explanation and details
  
  if(is.na(media)) {
    
    # treat as complete: all compounds with exchange reactions are available
    chgBndsCPLEX(cplex$env, cplex$lp, length(ex.rxnIds), ex.rf.col.indices, rep("U",length(ex.rxnIds)), rep(1000,length(ex.rxnIds)))
    
  } else {
    
    # Start by turning off all compounds
    chgBndsCPLEX(cplex$env, cplex$lp, length(ex.rxnIds), ex.rf.col.indices, rep("U",length(ex.rxnIds)), rep(0,length(ex.rxnIds)))
    
    # Now turn on the compounds that are present.
    col.indices <- sapply(names(media), function(name) {
      format1name <- paste0("EX_",name)  # naming convention in SBML documents we've received from Matt or exported from KBase
      if(format1name %in% ex.rxnIds) return(getColIndexCPLEX(cplex$env, cplex$lp, paste("RF_",format1name,sep="")))
      format2name <- paste0("R_EX_",substr(name,3,nchar(name)-2),"_LPAREN_e_RPAREN_")  # naming convention in Carrera SBML document
      if(format2name %in% ex.rxnIds) return(getColIndexCPLEX(cplex$env, cplex$lp, paste("RF_",format2name,sep="")))
      warning(paste("Couldn't find an exchange reaction for this compound:",name))
    })
    chgBndsCPLEX(cplex$env, cplex$lp, length(media), col.indices, rep("U",length(media)), media)
    
  }
  
}

setMedia_CarreraMethod <- function(cplex, supplNutrient) {
  # Unlike the setMedia function, this uses the data in the SBML itself to determine the media
  # supplNutrient: a single compound ID (or vector length 0) which is to be allowed up to 100 reverse flux. 
  #   This exactly mimics the behavior of the Carrera code; but it would be trivial to generalize to multiple supplemental nutrients
  
  if(length(supplNutrient)>1) stop("Only one supplNutrient is allowed (per experiment); this is to mimic the behavior of the Carrera code")
  else if(length(supplNutrient)==1) {
    format1name <- paste0("EX_",supplNutrient)  # naming convention in SBML documents we've received from Matt or exported from KBase
    format2name <- paste0("R_EX_",substr(name,3,nchar(name)-2),"_LPAREN_e_RPAREN_")  # naming convention in Carrera SBML document
    if(format1name %in% ex.rxnIds) { 
      supplRxn <- ex.rxns[which(sapply(ex.rxns, xmlGetAttr, "id")==format1name)]
    } else if(format2name %in% ex.rxnIds) {
      supplRxn <- ex.rxns[which(sapply(ex.rxns, xmlGetAttr, "id")==format2name)]
    } else {
      stop(paste("Couldn't find an exchange reaction for supplNutrient",supplNutrient))
    }
  } else {
    supplRxn <- NA
  }
  
  lower.bounds <- as.numeric(sapply(ex.rxns, function(rxn) {
    if(!is.na(supplRxn) && rxn==supplRxn) return(-100)
    subDoc <- xmlDoc(rxn)
    lower.bound <- xmlGetAttr(getNodeSet(subDoc, "//parameter[@id='LOWER_BOUND']")[[1]], "value")
    free(subDoc)  # avoid memory leak
    return(lower.bound)
  }))
  
  chgBndsCPLEX(cplex$env, cplex$lp, length(ex.rxns), ex.rf.col.indices, rep("U",length(ex.rxns)), -lower.bounds)  # set the upper bound of the reverse flux to the negative of the lower bound of the reversible reaction
}
