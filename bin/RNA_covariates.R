#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./RNA_covariates.R meta_csv salmon_gene_level_count_filtered_txt phenotype_PCs'



args = commandArgs(trailingOnly = TRUE)

if(length(args) < 3 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax)
  quit()
}

meta_fn = args[1]
count_fn = args[2]
phenotype_PCs = as.integer(args[3])



# ############
# setwd("/Users/datn/github/nf-rasqual-dev/data")
# # input
#  meta_fn = "meta/Brain.csv"
#  count_fn = "rna_gene_level_count_salmon_filtered.txt"
#  phenotype_PCs = 2

# ############
# # output:
#  out_fn = rna.covs_all_chrom.bin - rna.covs_all_chrom.txt

require(rasqualTools)
require(data.table)



set.seed(2022)


PCA_Covariates <- function(counts, size_factors, n_PCs = 2) {
  # author Natsuhiko Kumasaka
  #Map parameters to Natsuhiko's variables
  Y = counts
  K = size_factors
  n=ncol(Y)
  
  # fpm calculation
  fpkm=t(t(Y/K+1)/apply(Y/K,2,sum))*1e6 #  /len*1e9
  
  ## author Dat T Nguyen
  p = prcomp(fpkm)
  pca = p$rotation[,1:n_PCs]
  # Write covariates
  return(pca)
}




meta = fread(meta_fn, sep = ",")
meta = as.data.frame(meta)
count = fread(count_fn, skip = "gene_id", header = T, sep = "\t")
count = as.data.frame(count)


gene_id = count$gene_id

## count maxtrix processing
count2 = count[,-c(1:5)]
row.names(count2) = gene_id

## size factor with no GC correction
# comput size factors
size_factors = rasqualCalculateSampleOffsets(count2, gc_correct = FALSE)

## covariates
covs = PCA_Covariates(count2, size_factors, phenotype_PCs)
covs = cbind(meta[,-c(1:6)], covs)

#fwrite(covs, file = out_fn, sep = "\t")
saveRasqualMatrices(list(rna = covs), ".", file_suffix = "covs_all_chrom")
