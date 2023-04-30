#!/usr/bin/env python

import sys
import gzip
import argparse

parser = argparse.ArgumentParser(description = "Extract rasqual lead SNP - matrixQTL formated", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--normial", help = "normial pval input - matrixQTL formated")
parser.add_argument("--eigenMT", help = "out put of eigenMT")



args = parser.parse_args()

fi_BF = open(args.eigenMT, 'rt')
fi_norm = open(args.normial, 'rt')

lead_snps = dict()

for l in fi_BF:
    if l[0:3] == "snp": continue

    tokens = l.split("\t")
    key = tokens[1]
    val = tokens[0] + "__" + tokens[1]
    lead_snps[key] = val


for l in fi_norm:
    if l[0:3] == "snp":
        print(l[:-1])
    else:
        tokens = l.split("\t")
        key = tokens[1]
        val = tokens[0] + "__" + tokens[1]
        if lead_snps.get(key, 'NA') == val:
            print(l[:-1])





# python3 extract_lead_SNPs_nomial_pval.py --normial res_collected/$tis/obs_atac_rasqual_normial_pval.txt --eigenMT res_collected/$tis/obs_atac_rasqual_eigenMT_pval.txt > test.txt

# python3 extract.py --normial res_collected/$tis/obs_atac_rasqual_normial_pval.txt --eigenMT res_collected/$tis/obs_atac_rasqual_eigenMT_pval.txt > test.txt