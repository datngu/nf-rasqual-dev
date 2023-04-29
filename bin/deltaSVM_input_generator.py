#!/usr/bin/env python

import sys
import gzip
import argparse
import pysam
import pandas as pd
import numpy as np


parser = argparse.ArgumentParser(description = "Generating input fasta for deltaSVM from genome sequence and vcf.gz variant file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--vcf", help = "vcf.gz file input (SNP only)", required = True)
parser.add_argument("--genome", help = "genome file input", required = True)
#parser.add_argument("--bed", help = "bed file input")
parser.add_argument("--out", default = "out", help = "output prefix, outputs are: <out>_ref.fa and <out>_alt.fa")
parser.add_argument('--kmer', type=int, default = 10, help = 'kmer length of the gkmSVM weight')

args = parser.parse_args()

vcf_fn = args.vcf
genome_fn = args.genome
#bed_fn = args.bed
out = args.out
kmer = args.kmer

# genome_fn = "/Users/datn/GENOMES/atlantic_salmon/Salmo_salar.Ssal_v3.1.dna.toplevel.fa"
# vcf_fn = "data/genotype.vcf.gz"
# bed_fn = "atac1.txt"
# out = "out"
# kmer = 10

# bin/deltaSVM_input_generator.py --genome /Users/datn/GENOMES/atlantic_salmon/Salmo_salar.Ssal_v3.1.dna.toplevel.fa --vcf data/genotype.vcf.gz --kmer 10 --out test


# reading genome fasta
ref = pysam.FastaFile(genome_fn)
ref_fa = out + "_ref.fa"
alt_fa = out + "_alt.fa"


# reading vcf files
vcf = {"chr":[], "pos": [], "id" : [], "ref": [], "alt": []}
v = gzip.open(vcf_fn, "rt")
for l in v:
    if l[0] != "#":
        tokens = l.split()
        # consider SNP only
        if( len(tokens[3]) == 1 and len(tokens[4]) == 1 ):
            vcf["chr"].append(str(tokens[0]))
            vcf["pos"].append(int(tokens[1]))
            vcf["id"].append(str(tokens[2]))
            vcf["ref"].append(str(tokens[3]))
            vcf["alt"].append(str(tokens[4]))

v.close()



f_ref = open(ref_fa, 'w')
f_alt = open(alt_fa, 'w')

for i, chr in enumerate(vcf["chr"]):
    pos = vcf["pos"][i] -1
    stat = pos - (kmer - 1)
    end = pos + kmer
    seq_ref = ref.fetch(chr, stat, end)
    # handle extreme cases
    if len(seq_ref) == kmer*2-1:
        seq_alt = list(seq_ref)
        seq_alt[kmer-1] = vcf["alt"][i]
        seq_alt = ''.join(seq_alt)
        id = vcf["id"][i]
        f_ref.writelines(">" + id + "\n")
        f_ref.writelines(seq_ref + "\n")
        f_alt.writelines(">" + id + "\n")
        f_alt.writelines(seq_alt + "\n")


f_ref.close()    
f_alt.close()




