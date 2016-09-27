set.obj.and.knockouts.RBP <- function(cplex, model, aijs, knockouts, w, a, maxBiomass) {
  
  r.values.by.rxn <- sapply(get.gene.associations(model), parseAssocStr, R, aijs=aijs)
  gammas <- get.gammas(r.values.by.rxn)
  
  # sort the gamma vectors to be consistent with the ordering of reactions that we've been using
  gammas$off <- gammas$off[rxnIds]
  gammas$on <- gammas$on[rxnIds]
  
  # Any reaction with NA for its R value has been effectively knocked out due to gene knockouts.
  # Knocked-out reactions will have beta bounded to 0; this is very different from the "soft constraint" represented by an R value of 0.
  # First we remove any of these hard constraints which were added previously
  chgBndsCPLEX(cplex$env, cplex$lp, length(beta.col.indices), beta.col.indices, rep("U",length(beta.col.indices)), rep(1,length(beta.col.indices)))
  # Now we add the hard constraints for knocked-out reactions
  whichnas <- which(is.na(r.values.by.rxn)) 
  chgBndsCPLEX(cplex$env, cplex$lp, length(whichnas), beta.col.indices[whichnas], rep("U",length(whichnas)), rep(0,length(whichnas)))
  
  maxPenalty <- 0.5*length(reactions)
  
  ncols <- getNumColsCPLEX(cplex$env, cplex$lp)
  objCoeffs <- rep(0, ncols)  # start with a vector of all zeroes
  beta_objCoeffs <- (1-w)/maxPenalty*(1-a)*gammas$on
  delta_objCoeffs <- (1-w)/maxPenalty*a*gammas$off
  objCoeffs[beta.col.indices+1] <- beta_objCoeffs   # adding 1 to the beta.col.indices vector compensates for CPLEX returning zero-based indices
  objCoeffs[delta.col.indices+1] <- delta_objCoeffs
  objCoeffs[biomass.col.index+1] <- -w/maxBiomass
  
  chgObjCPLEX(cplex$env, cplex$lp, ncols, 0:(ncols-1), objCoeffs)
}

# Returns a named vector of strings; names are reaction ids, strings are the gene association strings
get.gene.associations <- function(model) {
  gene.association.nodes <- getNodeSet(model, "//x:reaction/x:notes/html:p[1]", c(x="http://www.sbml.org/sbml/level2", html="http://www.w3.org/1999/xhtml"))
  if(length(gene.association.nodes)==0) {
    gene.association.nodes <- getNodeSet(model, "//x:reaction/x:notes/html:body/html:p[1]", c(x="http://www.sbml.org/sbml/level2", html="http://www.w3.org/1999/xhtml"))  # an alternate format in some SBML files
    ids <- sapply(gene.association.nodes, function(node) xmlGetAttr(xmlParent(xmlParent(xmlParent(node))), "id"))  # id-getting in the alternate format
  } else {
    ids <- sapply(gene.association.nodes, function(node) xmlGetAttr(xmlParent(xmlParent(node)), "id"))  # id-getting in the original format
  }
  texts <- sapply(gene.association.nodes, xmlValue)
  gene.associations <- substr(texts, 18, nchar(texts))  # strip the 'GENE_ASSOCIATION:'
  names(gene.associations) <- ids
  return(gene.associations)
}

get.genenames <- function(gene.association) {
  assocStrNoParens <- parseAssocStr(gene.association, function(snippet, moreargs) snippet)
  words <- strsplit(assocStrNoParens, "\\s+")[[1]]  # get a list of 'words' (separated by whitespace)
  return(words[which(words!="" & words!="and" & words!="or" & words!="Unknown")])  # Remove things which are not genenames
}

withoutNames <- function(anything) {
  names(anything) <- NULL
  return(anything)
}

parseAssocStr <- function(assocStr, evalFunc, ...) {
  # Parses and evaluates assocStr, using evalFunc on every 'snippet' (i.e. string fragment between parentheses). 
  # Every snippet is replaced with the result of evalFunc(snippet), and this process continues recursively. 
  # ... contains additional arguments to be passed to evalFunc at every invocation
  
  innermost <- regmatches(assocStr,regexpr("\\([^()]*\\)",assocStr))  # find the first expression surrounded by a left and right paren but containing no parens.
  if (length(innermost)==0) return(withoutNames(evalFunc(assocStr, ...)))  # if this does not exist, then evaluate the expression plain and return.
  innermost <- substr(innermost,2,nchar(innermost)-1)  # strip the parens off each end
  if (length(innermost)==nchar(assocStr)-2) return(withoutNames(evalFunc(innermost, ...)))  # if the entire assocStr was surrounded by parentheses, just evaluate it without the parentheses and return that
  value <- as.character(evalFunc(innermost, ...))  # evaluate the sub-expression, put the result back in a string
  remainder <- regmatches(assocStr,regexpr("\\([^()]*\\)",assocStr), invert=T)[[1]]  # get everything else left from the string
  return(parseAssocStr(paste(remainder[[1]], value, remainder[[2]], sep=" "), evalFunc, ...))  # paste it all together and repeat
}

R <- function(snippet, aijs) {
  words <- strsplit(snippet, "\\s+")[[1]]  # get a list of 'words' in snippet (separated by whitespace)
  words <- words[which(words!="")]   # remove any empty-string 'words'
  if(length(words)==0) return(0.5)  # snippet was the empty string or just whitespace
  kind <- which.max(c("and" %in% words, "or" %in% words, TRUE))  
  # note that a snippet cannot have both 'and' and 'or'; this would be ambiguous
  # so kind==1 means this snippet contains 'and's; kind==2, 'or's; and kind==3, neither
  numbers <- unlist(sapply(words, function(word) {  # the unlist removes zero-length vectors
    if(word=="and" || word=="or") return(numeric(0))  # remove 'and' and 'or' words
    if(word=="Unknown") return(0.5)  
    if(is.na(word) || word == "NA") return(NA)  # If NA was passed (either an actual NA or the string "NA"), return NA again
    numericword <- suppressWarnings(as.numeric(word))
    if(!is.na(numericword)) return(numericword)  # a number: R(s) for some other s
    return(aijs[word])  # by process of elimination, word must be a gene name
  }, USE.NAMES=FALSE))
  if(kind==1) {
    return(min(numbers))  # If one of 'numbers' was NA, then one of the genes was knocked out; this correctly returns NA, as a knocked-out gene 'and' any other gene(s) is a knocked-out expression
  }
  else if(kind==2) {
    not.nas <- which(!is.na(numbers))
    if (length(not.nas) == 0)  return(NA)  # NA 'or' NA 'or' NA 'or'... = NA
    return(max(numbers[not.nas]))  # For any 'or' expression involving NA's (knocked-out genes), take the max of the non-NA values
  }
  else if(kind==3) {  # should be only one entry in numbers if kind==3
    return(numbers)
  }
}

get.gammas <- function(values) {
  list(
    off = sapply(values, function(val) {
      if(is.na(val)) return(0)
      if(val<0.5) return(0)
      return(val - 0.5)
    }),
    on = sapply(values, function(val) {
      if(is.na(val)) return(0)
      if(val>0.5) return(0)
      return(0.5 - val)
    })
  )
}

set.obj.GBP <- function(cplex, model, aijs, w, a, maxBiomass) {
  gammas <- get.gammas(aijs)
  
  # sort the gamma vectors to be consistent with the ordering of genes that we've been using, and drop any non-metabolic genes
  gammas$off <- gammas$off[metabolic.genes]
  gammas$on <- gammas$on[metabolic.genes]
  
  maxPenalty <- 0.5*length(metabolic.genes)
  
  ncols <- getNumColsCPLEX(cplex$env, cplex$lp)
  objCoeffs <- rep(0, ncols)  # start with a vector of all zeroes
  b_objCoeffs <- (1-w)/maxPenalty*(1-a)*gammas$on
  d_objCoeffs <- (1-w)/maxPenalty*a*gammas$off
  objCoeffs[b.col.indices+1] <- b_objCoeffs   # adding 1 to the b.col.indices vector compensates for CPLEX returning zero-based indices
  objCoeffs[d.col.indices+1] <- d_objCoeffs
  objCoeffs[biomass.col.index+1] <- -w/maxBiomass
  
  chgObjCPLEX(cplex$env, cplex$lp, ncols, 0:(ncols-1), objCoeffs)
}

set.knockouts.GBP <- function(cplex, knockouts) {
  # First remove any knockouts set previously
  chgBndsCPLEX(cplex$env, cplex$lp, length(b.col.indices), b.col.indices, rep("U",length(b.col.indices)), rep(1,length(b.col.indices)))
  # Now, set the new knockouts
  knockout.b.indices <- b.col.indices[match(knockouts, metabolic.genes)]
  chgBndsCPLEX(cplex$env, cplex$lp, length(knockouts), knockout.b.indices, rep("U",length(knockouts)), rep(0,length(knockouts)))
}