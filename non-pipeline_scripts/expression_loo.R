#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./this_script.R meta_csv feature_count_atac_txt salmon_merged_rna_txt'


require(data.table)

args = commandArgs(trailingOnly = TRUE)

if(length(args) < 3 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax)
  quit()
}

meta_fn = args[1]
atac_fn = args[2]
rna_fn = args[3]


require(data.table)
# meta_fn = "data/meta/brain.csv"
# atac_fn = "data/atac_consensus_peak_featureCounts.txt"
# rna_fn = "data/rna_gene_level_count_salmon.txt"


# meta
meta_out_dir = dirname(meta_fn)
meta_base_name = basename(meta_fn)

meta = fread(meta_fn)
meta = as.data.frame(meta)

# atac seq
atac_out_dir = dirname(atac_fn)
atac_base_name = basename(atac_fn)

atac = fread(atac_fn, skip = "Geneid")
atac = as.data.frame(atac)

# rna seq
rna_out_dir = dirname(rna_fn)
rna_base_name = basename(rna_fn)

rna = fread(rna_fn, skip = "gene_id")
rna = as.data.frame(rna)

for(i in 1:nrow(meta)){
    genotype_id = meta$genotype_id[i]
    atac_id = meta$atac_count_id[i]
    rna_id = meta$rna_count_id[i]


    meta_tem = meta[-i,]
    atac_tem = atac[, -(which(names(atac) == atac_id))]
    rna_tem = rna[, -(which(names(rna) == rna_id))]

    out_meta = paste0(meta_out_dir, "/", genotype_id, "_", meta_base_name, sep = "")
    fwrite(meta_tem, out_meta, sep = "\t")
    
    out_atac = paste0(atac_out_dir, "/", genotype_id, "_", atac_base_name, sep = "")
    fwrite(atac_tem, out_atac, sep = "\t")
    
    out_rna = paste0(rna_out_dir, "/", genotype_id, "_", rna_base_name, sep = "")
    fwrite(rna_tem, out_rna, sep = "\t")
}






