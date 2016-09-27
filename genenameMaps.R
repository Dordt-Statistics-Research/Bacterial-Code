# Returns a dataframe with three columns "fig", "b", and "common"
# Names in each row correspond to each other, i.e. refer to the same gene.
# Any names not found are NA entries in the dataframe.
# Notably, this function works with both "Carrera_genenames.txt" and "E_coli_common_genenames_edited.txt", despite the differences between those two files.
getGenenameMap <- function(genenameFile) {
  map <- read.delim(genenameFile, fill=T, header=F, row.names = NULL, stringsAsFactors = F)
  nrows <- nrow(map)
  newmap <- data.frame(fig=rep(NA,nrows), b=rep(NA,nrows), common=rep(NA,nrows), row.names=1:nrows)
  newmap$fig <- sapply(map[[1]], function(name) if(substr(name,1,3)=="fig") name else NA)
  newmap$b <- apply(map, 1, function(row) {
    result <- regmatches(row, regexpr("^\\s*b\\S+", row))
    if(length(result)==0) return(NA)
    return(regmatches(result[[1]], regexpr("\\S+", result[[1]])))  # strip whitespace from result[[1]]
  })
  newmap$common <- apply(map, 1, function(row) {
    result <- regmatches(row, regexpr("(?<=:)\\S+", row, perl=TRUE))
    if(length(result)==0) result <- regmatches(row, regexpr("^[^[:digit:][:space:]]\\S*", row))  # grabs all fields that don't start with a digit or whitespace
    result <- result[!(result %in% newmap$fig) & !(result %in% newmap$b)]  # remove any fields previously identified as a fig or b name
    if(length(result)==0) return(NA)
    return(result[[1]])
  })
  
  reportDuplicates <- function(vec) {
    vec <- vec[which(!is.na(vec))]
    if(anyDuplicated(vec)) {
      duplicates <- vec[which(duplicated(vec))]
      stop(paste("Duplicates of the following found in common genes file:", paste(duplicates, collapse=" ")))
    }
  }
  
  reportDuplicates(newmap$fig)
  reportDuplicates(newmap$b)
  reportDuplicates(newmap$common)
  
  return(newmap)
}
