#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./this_script.R meta_csv vcf_file'


require(data.table)

args = commandArgs(trailingOnly = TRUE)

if(length(args) < 1 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax)
  quit()
}

meta_fn = args[1]
vcf_fn = args[2]


require(data.table)
# meta_fn = "data/meta/Brain.csv"
# vcf_fn = "data/genotype.vcf.gz"

# meta
meta = fread(meta_fn, header = T)
meta = as.data.frame(meta)

for(i in 1:nrow(meta)){
    genotype_id = meta$genotype_id[i]
    # bcftools view  genotype.vcf.gz -s ^${sample} -o ${sample}_genotype.vcf.gz -Oz

    out_vcf = paste0(genotype_id, "_loo.vcf.gz", sep = "")
    cmd = paste0("bcftools view  ", vcf_fn, " -s ^", genotype_id, " -Oz -o ", out_vcf)
    try(system(cmd))
    cmd2 = paste0("bcftools index -t ", out_vcf)
    try(system(cmd2))
}






