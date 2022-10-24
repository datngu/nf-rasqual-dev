#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
args = commandArgs(trailingOnly = TRUE)

meta_path = out_file = args[1]
# meta_path = "data/meta/brain.csv"
meta = read.csv(meta_path, header = TRUE)
dir.create("copied_files")

for( i in 1: nrow(meta)){
    # bam
    new_bam = paste0("copied_files/", meta$genotype_id[i], ".bam")
    old_bam = meta$atac_bam_id[i]
    cmd = paste0("cp -L ", old_bam, " ", new_bam)
    print(cmd)
    system(cmd)
    # bai
    new_bai = paste0("copied_files/", meta$genotype_id[i], ".bam.bai")
    old_bai = paste0(meta$atac_bam_id[i], ".bai")
    cmd2 = paste0("cp -L ", old_bai, " ", new_bai)
    print(cmd2)
    system(cmd2)
}


