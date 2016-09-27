# Some generic utility functions I like. 


# Returns a version of df, with the selected colnames ('cols') at the front
moveColsToFront <- function(df, cols) df[,c(cols, setdiff(colnames(df),cols))]


# Returns a version of whatever it is passed, with names removed
withoutNames <- function(anything) {
  names(anything) <- NULL
  return(anything)
}


# From any named object x, this removes the elements with the given names
removeByName <- function(x, namesToRemove) x[setdiff(names(x), namesToRemove)]
# If you have indices (numbers), say numsToRemove, you can just do x[-numsToRemove]
# but you can't do x[-namesToRemove] because you can't negate a character vector.
# Therefore this function, which performs the intent of x[-namesToRemove].
# Note that if you try to remove names which aren't part of names(x), the function
#   will simply return x - no warning - which is both a feature and a pitfall


# For any named object x, return a version with some of 
#   its names (specified in character vector oldNames) replaced
#   with the names in the character vector newNames
renameByName <- function(x, oldNames, newNames) {
  if(length(oldNames)!=length(newNames)) stop("oldNames and newNames must be the same length")
  colNums <- match(oldNames, names(x))
  if(any(is.na(colNums))) stop("Some or all of oldNames were not names of the object")
  names(x)[colNums] <- newNames
  return(x)
}
# As with removeByName, if you had indices (numbers) of the columns to rename, this
# would be easy; but you can't do something like names(x)[oldNames] <- newNames, because
# names(x) isn't itself named.  Therefore, this function.  Unlike removeByName, this 
# function does throw an error if oldNames is not a subset of names(x). 


lapplyWithNames <- function(...) sapply(..., simplify=FALSE, USE.NAMES=TRUE)  
# Normal lapply() gives a result with the same names as whatever it was applied to.
# I.e. if you lapply() over a named vector, the result will have the same names as the
# vector you started with.  This is good and intuitive.
# What I also want, is that if I lapply() over a (un-named) character vector, the
# character vector itself is used as names for the result.  E.g. if I lapply() over
# 'phases', then I want the names of the result to be the phases. 
# sapply() has this behavior by default, but sometimes I don't want the simplifying
# behavior of sapply().  So, lapplyWithNames(): the naming of sapply() with the 
# non-simplifying of lapply().  Also note that when used on a named vector of any kind,
# lapplyWithNames() behaves identically to lapply(). 
