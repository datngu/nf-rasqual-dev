#!/usr/bin/env python

import sys
#import gzip

in_fn = sys.argv[1]
out_fn = sys.argv[2]

print(in_fn)
print(out_fn)

#in_fn = "/Users/datn/GENOMES/atlantic_salmon/Salmo_salar.Ssal_v3.1.106.gtf"
#out_fn = "gene.txt"



chr_dict = dict()
strand_dict = dict()
exon_stat_dict = dict()
exon_end_dict = dict()

f = open(in_fn,'rt')
for l in f:
    if l[0] == "#":
        continue
    tem = l.split("\t")
    if tem[2] == "gene":
        tem2 = tem[8].split(";")[0]
        gene_ID = tem2.split(" ")[1]
        chr_dict[gene_ID] = tem[0]
        strand_dict[gene_ID] = tem[6]
        exon_stat_dict[gene_ID] = ""
        exon_end_dict[gene_ID] = ""
    if tem[2] == "exon":
        tem3 = tem[8].split(";")[0]
        exon_gene_ID = tem3.split(" ")[1]
        exon_stat_dict[exon_gene_ID] = exon_stat_dict[exon_gene_ID] + "," + tem[3]
        exon_end_dict[exon_gene_ID] = exon_end_dict[exon_gene_ID] + "," + tem[4]

f.close()

out_file = open(out_fn, 'w')
for gene in chr_dict.keys():
    res = gene.strip('"') + "\t" + chr_dict[gene] + "\t" + strand_dict[gene] + "\t" + exon_stat_dict[gene][1:] + "\t" + exon_end_dict[gene][1:] + "\n"
    out_file.writelines(res)

out_file.close()