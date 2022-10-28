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


for(i in 1:nrow(meta)){
    genotype_id = meta$genotype_id[i]
    meta_tem = meta[-i,]
    out_meta = paste0(genotype_id, "_meta.csv", sep = "")
    fwrite(meta_tem, out_meta, sep = "\t")
}





