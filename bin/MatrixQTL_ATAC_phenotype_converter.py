#!/usr/bin/env python

import sys
import argparse

parser = argparse.ArgumentParser(description = "Convert count file to phenotype data input of MatrixQTL", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--count", help = "ATAC count file input")
parser.add_argument("--out_phenotype", help = "phenotype expression output")
parser.add_argument("--out_phenotype_position", help = "phenotype position output")

args = parser.parse_args()

fi = open(args.count, 'rt')
fout_phen = open(args.out_phenotype, 'w')
fout_pos = open(args.out_phenotype_position, 'w')

for l in fi:
    tokens = l.split("\t")
    if tokens[0] == "Geneid":
        tokens[0] = "ID"
        l1 = "\t".join(tokens[0:1] + tokens[6:])
        #l1 = l1 + "\n"
        fout_phen.writelines(l1)
        l2 =  "\t".join(["gene_id", "chrom_probe", "s1", "s2"])
        l2 = l2 + "\n"
        fout_pos.writelines(l2)
    else:
        gene_id = ":".join(tokens[0:6])
        l1 = "\t".join([gene_id] + tokens[0:1] + tokens[5:])
        #l1 = l1 + "\n"
        fout_phen.writelines(l1)
        
        chr = tokens[1]
        stat = tokens[2]
        end = tokens[3]
        l2 = "\t".join([ gene_id, chr, stat, end])
        l2 = l2 + "\n"
        fout_pos.writelines(l2)

fi.close()
fout_phen.close()
fout_pos.close()

# python3 ./bin/MatrixQTL_ATAC_phenotype_converter.py --count data/atac_consensus_peak_featureCounts_filtered.txt --out_phenotype tes1 --out_phenotype_position tes2