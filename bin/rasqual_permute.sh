#!/bin/bash
vcf=$1
y=$2
k=$3
x=$4
x_txt=$5
meta=$6
cpu=$7
permute=$8
chr=$9


for i in $(seq 1 $permute)
do
   rasqual_permute.R vcf=$vcf y=$y k=$k x=$x x_txt=$x_txt meta=$meta out=${chr}_permute_${i}_rasqual_lead_snp.txt cpu=$cpu
done
