#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./RNA_covariates.R meta_csv salmon_gene_level_count_filtered_txt genotype_vcf phenotype_PCs'



args = commandArgs(trailingOnly = TRUE)

if(length(args) < 4 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax)
  quit()
}

meta_fn = args[1]
count_fn = args[2]
geno_fn = args[3]
phenotype_PCs = as.integer(args[4])



# # ############
# setwd("/mnt/SCRATCH/ngda/nf-rasqual")
# # input
#  meta_fn = "data/meta/brain.csv"
#  count_fn = "/mnt/SCRATCH/ngda/nf-rasqual/results/RNA_filtering_expression/rna_gene_level_count_salmon_filtered.txt"
#  geno_fn = "data/genotype.vcf.gz"
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
genotype = fread(geno_fn, skip = "CHROM", sep = "\t")
genotype = as.data.frame(genotype)

# ordering count data by genotype data
geno_id = colnames(genotype)[-c(1:9)]
meta = meta[meta$genotype_id %in% geno_id,]
od = match(geno_id, meta$genotype_id)
meta = meta[od,]
rownames(meta) = meta$genotype_id

#atac_peaks = paste(count$Geneid, count$Chr, count$Start, count$End, count$Strand, count$Length, sep = ":")
gene_id = count$gene_id

## count maxtrix processing
count2 = count[,-c(1:5)]
row.names(count2) = gene_id
count2 = count2[ , meta$rna_count_id]
colnames(count2) = meta$genotype_id

## size factor with no GC correction
# comput size factors
size_factors = rasqualCalculateSampleOffsets(count2, gc_correct = FALSE)

## covariates
covs = PCA_Covariates(count2, size_factors, phenotype_PCs)
covs = cbind(meta[,-c(1:6)], covs)

#fwrite(covs, file = out_fn, sep = "\t")
saveRasqualMatrices(list(rna = covs), ".", file_suffix = "covs_all_chrom")
