# 
# generate_feature_map = function(feature_file, map) {
#   con = file(feature_file, "r")
#   exp_name = ""
#   while ( TRUE ) {
#     line = readLines(con, n = 1)
#     if ( length(line) == 0) {
#       break
#     }
#     if (nchar(line)==0) {
#       #SKIP
#     } else if (grepl("# TF", line)) {
#       skip_rna = FALSE
#       tail = substring(line,8,nchar(line))
#       tail = paste0(tolower(substr(tail,1,1)), substr(tail,2,nchar(tail)))
#       s = strsplit(tail,": ")
#       bid = s[[1]][2]
#       current_tf = bid
#     } else {
#       s = strsplit(line,"\t")
#       bid = s[[1]][2]
#       tfs_bid <- c(tfs_bid, current_tf)
#       tgs_bid <- c(tgs_bid, bid)
#       temp = map[[current_tf]][["peg"]]
#       tfs_peg <- c(tfs_peg, ifelse(is.null(temp), "", temp))
#       temp = map[[bid]][["peg"]]
#       tgs_peg <- c(tgs_peg, ifelse(is.null(temp), "", temp))
#       temp = map[[current_tf]][["name"]]
#       tfs_name <- c(tfs_name, ifelse(is.null(temp), "", temp))
#       temp = map[[bid]][["name"]]
#       tgs_name <- c(tgs_name, ifelse(is.null(temp), "", temp))
#     }
#   }
#   close(con)
#   data = data.frame(tgs_bid, tfs_bid, tgs_peg, tfs_peg, tgs_name, tfs_name)
#   return(data)
# }
# 
# 
# setwd("/home/jason/Documents/School/College/Fall Junior/STAT RESEARCH/MyWD/")
# exp_feature_map <- new.env(hash=T, parent=emptyenv())
# generate_feature_map("E_coli_v4_Build_6_6-23-15.experiment_feature_descriptions", exp_feature_map)
# save(exp_feature_map, file = "ExperimentFeatureMap.Rda")

setwd("/home/jason/Documents/School/College/Fall Junior/STAT RESEARCH/MyWD/Bacterial-Code/inputs")
feature_data <- read.table("E_coli_v4_Build_6_6-23-15.experiment_feature_descriptions", sep=",",header=TRUE)

