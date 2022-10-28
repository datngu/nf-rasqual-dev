#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./this_script.R meta_csv atac_count_filtered_txt'


require(data.table)

args = commandArgs(trailingOnly = TRUE)

if(length(args) < 1 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax)
  quit()
}

meta_fn = args[1]
atac_fn = args[2]



require(data.table)
# meta_fn = "data/meta/Brain.csv"
# atac_fn = "data/atac_consensus_peak_featureCounts_filtered.txt"



# meta
meta = fread(meta_fn, header = T)
meta = as.data.frame(meta)

# atac seq
atac = fread(atac_fn, header = T)
atac = as.data.frame(atac)


for(i in 1:nrow(meta)){
    genotype_id = meta$genotype_id[i]
    meta_tem = meta[-i,]
    atac_tem = atac[, -(which(names(atac) == genotype_id))]

    out_atac = paste0(genotype_id, "_atac_count.txt", sep = "")
    fwrite(atac_tem, out_atac, sep = "\t")
}






