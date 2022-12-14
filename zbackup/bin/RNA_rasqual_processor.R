#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./RNA_rasqual_processor.R meta_csv feature_count_txt genotype_vcf genome.fa cis_window n_core'

# defaut outputs:
# rna.covs.bin - rna.covs.txt
# rna.exp.bin - rna.exp.txt
# rna.size_factors.bin - rna.size_factors.txt
# snp_counts.tsv



args = commandArgs(trailingOnly = TRUE)

if(length(args) < 5 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax)
  quit()
}

cpu = NA
meta_fn = args[1]
count_fn = args[2]
geno_fn = args[3]
genome_fn = args[4]
cis_window = as.integer(args[5])
if(length(args) >= 6) cpu = as.integer(args[6])
if (is.na(cpu)) cpu=1
cpu = as.integer(cpu)



#setwd("/Users/datn/github/nf-rasqual/data")

#############

# setwd("/mnt/SCRATCH/ngda/nf-rasqual")
# # input
#  cis_window = 10000
#  meta_fn = "data/meta/brain.csv"
#  count_fn = "/mnt/SCRATCH/ngda/nf-rasqual/results/RNA_split_chrom/8_count.txt"
#  geno_fn = "data/genotype.vcf.gz"
#  genome_fn = "/mnt/users/ngda/genomes/atlantic_salmon/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa"

#############
# output:
# rna.exp.bin - rna.exp.txt
# rna.size_factors.bin - rna.size_factors.txt
#############


require(rasqualTools)
require(GenomicFeatures)
require(Biostrings)
require(data.table)
require(foreach)
require(doParallel)
registerDoParallel(cores=cpu)


set.seed(2022)


get_GC <- function(genome, feature_info){
  #for(i in 1:nrow(feature_info)){
  seqs = foreach(i = 1:nrow(feature_info),.combine = c) %dopar% {
    x = ""
    row = feature_info[i,]
    chr = row$chr
    s = unlist(strsplit(row$exon_starts, ","))
    s = as.integer(s)
    e = unlist(strsplit(row$exon_ends, ","))
    e = as.integer(e)
    x = substring(genome[chr], s, e)
    x = paste(x, collapse = "")
  }
  seqs = DNAStringSet(seqs)
  alf <- Biostrings::alphabetFrequency(seqs, as.prob=TRUE)
  gc = rowSums(alf[,c("G", "C"), drop=FALSE])
  return(gc)
}

# HANDLING UNKNOWN ERROR WHEN COMPUTE OFFSETS
get_offset <- function(count2, GC) {
  x = try(rasqualCalculateSampleOffsets(count2, GC))
  if(class(x) == "matrix"){
    res = x
  }else{
    res = rasqualCalculateSampleOffsets(count2, gc_correct = FALSE)
  }
  return(res)
}


meta = fread(meta_fn, sep = ",")
meta = as.data.frame(meta)
count = fread(count_fn, skip = "gene_id", header = T, sep = "\t")
count = as.data.frame(count)
genotype = fread(geno_fn, skip = "CHROM", sep = "\t")
genotype = as.data.frame(genotype)

# process genome for GC counting
genome = readDNAStringSet(genome_fn)
tem = strsplit(names(genome), " ")
names(genome) = do.call("rbind", tem)[,1]


# ordering count data by genotype data
geno_id = colnames(genotype)[-c(1:9)]
meta = meta[meta$genotype_id %in% geno_id,]
od = match(geno_id, meta$genotype_id)
meta = meta[od,]
rownames(meta) = meta$genotype_id

gene_id = count$gene_id


## counting snps overlapping atac peaks
gene_info = data.frame(gene_id = gene_id)
gene_info$chr = as.character(count$chr)
gene_info$strand = 1
gene_info$strand[which(count$strand == "-")] = -1
gene_info$strand = as.integer(gene_info$strand)
gene_info$exon_starts = count$exon_starts
gene_info$exon_ends = count$exon_ends

snp_info = genotype[,c(1:3)]
colnames(snp_info) = c("chr", "pos", "snp_id")
#snp_info$chr = gsub("ssa0", "", snp_info$chr, fixed = TRUE)
#snp_info$chr = gsub("ssa", "", snp_info$chr, fixed = TRUE)

snp_counts = countSnpsOverlapingExons(gene_info, snp_info, cis_window = cis_window)
fwrite(snp_counts, file = "snp_counts.tsv", sep = "\t")



## save count maxtrix
count2 = count[,-c(1:5)]
row.names(count2) = gene_id
count2 = count2[ ,meta$rna_count_id]
colnames(count2) = meta$genotype_id
saveRasqualMatrices(list( rna = count2), ".", file_suffix = "exp")

## size factor
# GC counting
gc = get_GC(genome, gene_info)

# # fix inf values
# pick = gc > 0 & gc <1
# gc[!pick] = 0.5

gc_percentage = gc*100
GC = data.frame(gene_id = row.names(count2), percentage_gc_content = gc_percentage)
# comput size factors
size_factors = get_offset(count2, GC)

# size_factors = rasqualCalculateSampleOffsets(count2, gc_correct = FALSE)
# size_factors = rasqualCalculateSampleOffsets(count2, GC)

saveRasqualMatrices(list(rna = size_factors), ".", file_suffix = "size_factors")




