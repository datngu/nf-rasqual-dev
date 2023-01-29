#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./slipt_bed.R attac_filtered.txt fold_size out_prefix\n'



args = commandArgs(trailingOnly = TRUE)

if(length(args) < 3 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax)
  quit()
}

in_bed = args[1]
fold_size = as.integer(args[2])
out_prefix = args[3]

# in_bed = "/Users/datn/github/nf-rasqual-dev/data/atlantic_salmon_v3.1_trf.bed"
# fold_size = 10
# out_prefix = "atac_bed"

data = read.table(in_bed)
data = as.data.frame(data)
data = data[,c(2,3,4)]
pick = data[,1] %in% as.character(1:29)
data = data[pick,]


data$fold_id = c(1:nrow(data))
set.seed(2021)

id  = sample(1:nrow(data))
data = data[id,]
step = nrow(data) %/% fold_size + 1 
id = rep(1:fold_size, step)
data$fold_id = id[1:nrow(data)]
for(i in 1:fold_size){
  out = paste0(out_prefix, i, ".txt")
  tem = data[data$fold_id == i,]
  tem = tem[,-4]
  write.table(tem, file = out, sep = "\t", quote = F, row.names = F, col.names = F)
}

