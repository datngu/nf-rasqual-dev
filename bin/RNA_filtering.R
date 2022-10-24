#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./RNA_filtering.R in_count out_count gene_info exp_prop fpkm_cutoff meta_csv'


args = commandArgs(trailingOnly = TRUE)

if(length(args) < 5 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax)
  quit()
}

in_count = args[1]
out_fn = args[2]
gene_info = args[3]
exp_prop = as.numeric(args[4])
fpkm_cutoff = as.numeric(args[5])
meta_csv = args[6]

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


# in_count = "/mnt/users/ngda/ngs_data/atlantic_salmon/brain/rna_gene_level_count_salmon.txt"
# out_fn = "rna_gene_level_count_salmon_filtered.txt"
# gene_info = "/mnt/SCRATCH/ngda/nf-rasqual/gene_info/gene_info.txt"
# exp_prop = 0.5
# fpkm_cutoff = 0.5



count = fread(in_count, header = T)
count = as.data.frame(count)
meta = fread(meta_csv)
meta = as.data.frame(meta)

count1 = count[,c(1,2)]
count2 = count[,-c(1,2)]
count2 = count2[, meta$rna_count_id]
names(count2) = meta$genotype_id

fpkm = compute_fpkm(count2)
fpkm2 = fpkm > fpkm_cutoff
pick = rowSums(fpkm2)/ncol(fpkm2) >= exp_prop

count3 = cbind(count1, count2)
count3 = count3[pick,]
count3 = count3[,-2]

info = fread(gene_info, header = F)
names(info) = c("gene_id", "chr", "strand", "exon_starts", "exon_ends")
info = info[info$gene_id %in% count3$gene_id,]

res = merge(info, count3, by.x ="gene_id", by.y ="gene_id")


fwrite(res, file = out_fn, sep = "\t")