source("expnameMaps.R")
expnameMap <- get.expname.map("inputs/Ecoli_Cel_to_chip_fromClaire_edited.tab")

rawcalls <- read.delim("inputs/ecoli_revised_activity_4-14-16.txt")
rownames(rawcalls) <- rawcalls$gene_id
rawcalls$gene_id <- NULL

celnames <- unlist(get_all_CELnames(expnameMap, colnames(rawcalls)))
calls <- rawcalls[,get_chipname(expnameMap, celnames)]  # duplicate experiments as necessary, and sort into the order of celnames
colnames(calls) <- celnames  # convert to CEL names (they were already in exactly the right order)

save(calls, file="Rsave/mattcalls_4-14-16.Rdata")