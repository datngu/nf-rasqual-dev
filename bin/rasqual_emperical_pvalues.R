#!/usr/bin/env Rscript

# Author: Dat T Nguyen <ndat@utexas.edu>
# Date: 30 Sep 2022

require(data.table)
require(qvalue)

options(stringsAsFactors=FALSE)
syntax='Usage:
              ./this_script.R out_file rasqual_result rasqual_permute_1 rasqual_permute_2 rasqual_permute_3 ...'

args = commandArgs(trailingOnly = TRUE)

if(length(args) < 4 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax, useSource = TRUE)
  quit()
}

out_file = args[1]
rasqual_result = args[2]
in_files = args[-c(1,2)]


# out_file="qvalue_rasqual.txt"
# rasqual_result = "/mnt/SCRATCH/ngda/nf-rasqual/ATAC_results_rasqual/all_chromosome_rasqual_lead_snp.txt"
# fdr = 0.05
# in_files <- list.files(pattern = "^permute_")


cat("\nout file: ", out_file, "\n")
cat("\nrasqual_result: ", rasqual_result, "\n")
x = paste(in_files, collapse = " ,")
cat("\nrasqual permutation list:\n", x, "\n")


get_chi_square <- function(rasqual_res_path){
  d = fread(rasqual_res_path, header = F, fill = T)
  c = d$V11
  return(c)
}


df = fread(rasqual_result, header = F, fill = T)
#df = as.data.frame(df)
res_cq = df[,c(1,2)]
colnames(res_cq) = c("feature_id", "snp_id")
res_cq$chi_square = df$V11

for(file in in_files){
  col = gsub("_all_chromosome_rasqual_lead_snp.txt", "", file)
  res_cq[,col] = get_chi_square(file)
}
#res2 = res
res_cq = res_cq[!res_cq$snp_id == "SKIPPED", ]


x1 = as.numeric(res_cq[[3]])
x2 = as.matrix(res_cq[, -c(1:3)])
#emp = empPvals(x1, x2, pool = FALSE)
emp = empPvals(x1, x2, pool = T)
emp_fdr = qvalue(emp)
q = emp_fdr$qvalues
# sum(q<0.05)
out_data = data.frame(empPvals = emp, qvalues = q)
out_data = cbind(out_data, res_cq)

fwrite(out_data, file = out_file, sep = "\t")
