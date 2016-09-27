# Returns a "expnameMap" object suitable to pass to the other functions in this script
get.expname.map <- function(mapFile) {
  # tested with mapFile=="inputs/Ecoli_Cel_to_chip_fromClaire_edited.tab"
  
  if(is.na(mapFile)) return(NA)
  
  expnameMap <- read.table(mapFile, header=T, quote="", colClasses="character")  # First column, CEL name; second column, chip name
  expnameMap[[1]] <- substr(expnameMap[[1]], 1, nchar(expnameMap[[1]])-4)  # Remove the .CEL
  expnameMap[[2]] <- substr(expnameMap[[2]], 1, nchar(expnameMap[[2]])-3)  # Remove the replicate number
  
  return(expnameMap)
}

# Returns the *first* CELname matching each chipname
get_first_CELname <- function(expnameMap, chipnames) {
  # chipnames: A character vector of chipnames (or a single chipname)
  expnameMap[match(chipnames,expnameMap[[2]]), 1]
}

# Returns a list with one component for each chipname
# Each component is a vector of all CELnames which the chipname maps to
get_all_CELnames <- function(expnameMap, chipnames) {
  lapply(chipnames, function(chipname) {
    matchrows <- which(expnameMap[[2]]==chipname)
    expnameMap[matchrows, 1]
  })
}

get_chipname <- function(expnameMap, CELnames) {
  # CELnames: A character vector of CELnames (or a single CELname)
  expnameMap[match(CELnames,expnameMap[[1]]), 2]
}