#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./ATAC_rasqual_processor.R feature_count_txt genotype_vcf genome.fa cis_window n_core'

# defaut outputs:
# atac.covs.bin - atac.covs.txt
# atac.exp.bin - atac.exp.txt
# atac.size_factors.bin - atac.size_factors.txt
# snp_counts.tsv



args = commandArgs(trailingOnly = TRUE)

if(length(args) < 4 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax)
  quit()
}

cpu = NA
count_fn = args[1]
geno_fn = args[2]
genome_fn = args[3]
cis_window = as.integer(args[4])
if(length(args) >= 5) cpu = as.integer(args[5])
if (is.na(cpu)) cpu=1
cpu = as.integer(cpu)



#setwd("/Users/datn/github/nf-rasqual-dev/data")

#############

# input
#  cis_window = 10000
#  count_fn = "atac_consensus_peak_featureCounts_filtered.txt"
#  geno_fn = "genotype.vcf.gz"
#  genome_fn = "/Users/datn/GENOMES/atlantic_salmon/Salmo_salar.Ssal_v3.1.dna.toplevel.fa"

#############
# output:
# atac.exp.bin - atac.exp.txt
# atac.size_factors.bin - atac.size_factors.txt
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


count = fread(count_fn, skip = "Geneid", header = T, sep = "\t")
count = as.data.frame(count)
# count = count[1:100,]

genotype = fread(geno_fn, skip = "CHROM", sep = "\t")
genotype = as.data.frame(genotype)

# process genome for GC counting
genome = readDNAStringSet(genome_fn)
tem = strsplit(names(genome), " ")
names(genome) = do.call("rbind", tem)[,1]


atac_peaks = paste(count$Geneid, count$Chr, count$Start, count$End, count$Strand, count$Length, sep = ":")

## counting snps overlapping atac peaks
peak_info = data.frame(gene_id = atac_peaks)
peak_info$chr = as.character(count$Chr)
peak_info$strand = 1
peak_info$strand[which(count$Strand == "-")] = -1
peak_info$strand = as.integer(peak_info$strand)
peak_info$exon_starts = as.character(count$Start)
peak_info$exon_ends = as.character(count$End)


snp_info = genotype[,c(1:3)]
colnames(snp_info) = c("chr", "pos", "snp_id")

# counting SNPs
snp_counts = countSnpsOverlapingExons(peak_info, snp_info, cis_window = cis_window)
fwrite(snp_counts, file = "snp_counts.tsv", sep = "\t")



## save count maxtrix
count2 = count[,-c(1:6)]
row.names(count2) = atac_peaks
saveRasqualMatrices(list( atac = count2), ".", file_suffix = "exp")

## size factor
# GC counting
gc = get_GC(genome, peak_info)
gc_percentage = gc*100
GC = data.frame(gene_id = row.names(count2), percentage_gc_content = gc_percentage)

# comput size factors
size_factors = get_offset(count2, GC)
saveRasqualMatrices(list(atac = size_factors), ".", file_suffix = "size_factors")




