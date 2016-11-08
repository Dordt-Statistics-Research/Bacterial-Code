
processFileTRN = function(filepath, map) {
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0) {
      break
    }
    if (nchar(line)==0 || grepl("# RNA", line)) {
      #SKIP
    }else if (grepl("# TF", line)) {
      tail = substring(line,8,nchar(line))
      tail = paste0(tolower(substr(tail,1,1)), substr(tail,2,nchar(tail)))
      s = strsplit(tail,": ")
      bid = s[[1]][2]
      name = s[[1]][1]
      stored = list("bid"=bid, "name"=name)
      map[[bid]] <- stored
      map[[name]] <- map[[bid]]
    } else {
      s = strsplit(line,"\t")
      bid = s[[1]][2]
      name = s[[1]][3]
      stored = list("bid"=bid, "name"=name)
      map[[bid]] <- stored
      map[[name]] <- map[[bid]]
    }
  }
  close(con)
}


processFileAnno = function(filepath, map) {
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    s = strsplit(line,"\t")
    bid = s[[1]][2]
    peg = s[[1]][1]
    if (exists(bid, envir=map, inherits=FALSE)) {
      map[[bid]][["peg"]] <- peg
      map[[peg]] <- map[[bid]]
      map[[map[[bid]][["name"]]]] <- map[[bid]]
    }
  }
  close(con)
}

generate_conversion_map = function(trn_file, peg_file, map) {
  processFileTRN(trn_file, map)
  processFileAnno(peg_file, map)
}

generate_better_trn = function(trn_file, map) {
  con = file(trn_file, "r")
  current_tf = ""
  tfs_bid = character(0)
  tgs_bid = character(0)
  tfs_peg = character(0)
  tgs_peg = character(0)
  tfs_name = character(0)
  tgs_name = character(0)
  skip_rna = FALSE
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0) {
      break
    }
    if (nchar(line)==0) {
      #SKIP
    } else if (grepl("# RNA", line)){
      skip_rna = TRUE
    } else if (grepl("# TF", line)) {
      skip_rna = FALSE
      tail = substring(line,8,nchar(line))
      tail = paste0(tolower(substr(tail,1,1)), substr(tail,2,nchar(tail)))
      s = strsplit(tail,": ")
      bid = s[[1]][2]
      current_tf = bid
    } else if (!skip_rna){
      s = strsplit(line,"\t")
      bid = s[[1]][2]
      tfs_bid <- c(tfs_bid, current_tf)
      tgs_bid <- c(tgs_bid, bid)
      temp = map[[current_tf]][["peg"]]
      tfs_peg <- c(tfs_peg, ifelse(is.null(temp), "", temp))
      temp = map[[bid]][["peg"]]
      tgs_peg <- c(tgs_peg, ifelse(is.null(temp), "", temp))
      temp = map[[current_tf]][["name"]]
      tfs_name <- c(tfs_name, ifelse(is.null(temp), "", temp))
      temp = map[[bid]][["name"]]
      tgs_name <- c(tgs_name, ifelse(is.null(temp), "", temp))
    }
  }
  close(con)
  data = data.frame(tgs_bid, tfs_bid, tgs_peg, tfs_peg, tgs_name, tfs_name)
  return(data)
}

setwd("/home/jason/Documents/School/College/Fall Junior/STAT RESEARCH/MyWD/")
name_map <- new.env(hash=T, parent=emptyenv())
generate_conversion_map("TRN.txt", "83333.1.anno", name_map)
save(name_map, file = "NamingConversionMap.Rda")
trn_struct = generate_better_trn("TRN.txt", name_map)
write.csv(trn_struct, file="TRN.csv",row.names=FALSE, na="")