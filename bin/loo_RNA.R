#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./this_script.R meta_csv rna_count_filtered_txt'


require(data.table)

args = commandArgs(trailingOnly = TRUE)

if(length(args) < 1 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax)
  quit()
}

meta_fn = args[1]
rna_fn = args[2]



require(data.table)
# meta_fn = "data/meta/Brain.csv"
# rna_fn = "data/rna_gene_level_count_salmon_filtered.txt"

# meta
meta = fread(meta_fn, header = T)
meta = as.data.frame(meta)

# rna seq
rna = fread(rna_fn, header = T)
rna = as.data.frame(rna)


for(i in 1:nrow(meta)){
    genotype_id = meta$genotype_id[i]
    meta_tem = meta[-i,]
    rna_tem = rna[, -(which(names(rna) == genotype_id))]

    out_rna = paste0(genotype_id, "_rna_count.txt", sep = "")
    fwrite(rna_tem, out_rna, sep = "\t")
}






