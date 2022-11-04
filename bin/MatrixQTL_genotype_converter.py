#!/usr/bin/env python

import sys
import gzip
import argparse

parser = argparse.ArgumentParser(description = "Convert count file to phenotype data input of MatrixQTL", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--vcf", help = "vcf.gz file input")
parser.add_argument("--out_genotype", help = "genotype output")
parser.add_argument("--out_genotype_position", help = "genotype position output")



args = parser.parse_args()

fi = gzip.open(args.vcf, 'rt')
fout_gen = open(args.out_genotype, 'w')
fout_pos = open(args.out_genotype_position, 'w')

for l in fi:
    if l[0:2] == "##": continue

    tokens = l.split("\t")
    if tokens[0] == "#CHROM":
        l1 = "\t".join(tokens[2:3] + tokens[9:])
        #l1 = l1 + "\n"
        fout_gen.writelines(l1)
        l2 =  "\t".join(["snp", "chr_snp", "pos"])
        l2 = l2 + "\n"
        fout_pos.writelines(l2)
    else:
        genotype_count = []
        for t in tokens[9:]:
            #print(t[0:3], "\n")
            if(t[0:3] == "1|0"): tem = "1"
            if(t[0:3] == "0|1"): tem = "1"
            if(t[0:3] == "1|1"): tem = "2"
            if(t[0:3] == "0|0"): tem = "0"
            genotype_count.append(tem)

        l1 = "\t".join(tokens[2:3] + genotype_count)
        l1 = l1 + "\n"
        fout_gen.writelines(l1)
        snp = tokens[2]
        chr = tokens[0]
        pos = tokens[1]
        l2 = "\t".join([ snp, chr, pos])
        l2 = l2 + "\n"
        fout_pos.writelines(l2)

fi.close()
fout_gen.close()
fout_pos.close()

# python3 ./bin/MatrixQTL_genotype_converter.py --vcf ./data/genotype.vcf.gz --out_genotype tes1 --out_genotype_position tes2
