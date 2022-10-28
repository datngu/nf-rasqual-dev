#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./this_script.R meta_csv'


args = commandArgs(trailingOnly = TRUE)

if(length(args) < 1 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax)
  quit()
}

meta_fn = args[1]

require(data.table)
# meta_fn = "data/meta/Brain.csv"

# meta
meta = fread(meta_fn, header = T)
meta = as.data.frame(meta)

dir.create("tem")

for(i in 1:nrow(meta)){
    genotype_id = meta$genotype_id[i]
    meta_tem = meta[-i,]
    out_meta = paste0("tem/", genotype_id, ".csv", sep = "")
    fwrite(meta_tem, out_meta, sep = "\t")
}






