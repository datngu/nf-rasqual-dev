#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./ATAC_filtering.R in_count out_count exp_prop fpkm_cutoff meta_csv'


args = commandArgs(trailingOnly = TRUE)

if(length(args) < 4 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax)
  quit()
}

in_count = args[1]
out_fn = args[2]
exp_prop = as.numeric(args[3])
fpkm_cutoff = as.numeric(args[4])
meta_csv = args[5]


compute_size_factor <- function(counts){

  #Calculate library sizes
  library_size = colSums(counts)
  size_factors = library_size/mean(library_size) #Standardise
  size_matrix = matrix(rep(size_factors, nrow(counts)), nrow = nrow(counts), byrow = TRUE)
  rownames(size_matrix) = rownames(counts)
  colnames(size_matrix) = colnames(counts)
  return(size_matrix)
}

compute_fpkm <- function(counts){
  Y = counts
  K = compute_size_factor(Y)
  n=ncol(Y)
  fpkm=t(t(Y/K+1)/apply(Y/K,2,sum))*1e6
  return(fpkm)
}

require(data.table)

# singularity run /mnt/SCRATCH/ngda/nf-rasqual/shared_dir/singularity/ndatth-rasqual-v0.0.0.img
# in_count = "data/atac_consensus_peak_featureCounts.txt"
# out_fn = 
# exp_prop = 0.5
# fpkm_cutoff = 0.5
# meta_csv = "data/meta/brain.csv"


meta = fread(meta_csv)
meta = as.data.frame(meta)
count = fread(in_count, skip = "Geneid", header = T)
count = as.data.frame(count)

count1 = count[,c(1:6)]
count2 = count[,-c(1:6)]
count2 = count2[, meta$atac_count_id]
names(count2) = meta$genotype_id

fpkm = compute_fpkm(count2)
fpkm2 = fpkm > fpkm_cutoff
pick = rowSums(fpkm2)/ncol(fpkm2) >= exp_prop

count3 = cbind(count1, count2)
count3 = count3[pick,]

fwrite(count3, file = out_fn, sep = "\t")