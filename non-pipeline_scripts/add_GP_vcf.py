#!/usr/bin/env python

import os
import sys
import gzip


in_vcf_GP = sys.argv[1]
in_vcf_phased = sys.argv[2]
out_vcf = sys.argv[3]

print(in_vcf_GP)
print(in_vcf_phased)
print(out_vcf)

#########################
# in_vcf_phased = 'ssa29.vcf.gz'
# in_vcf_GP = 'ssa29_tem.vcf.gz'
# out_vcf = 'out.vcf'


fi_phased = gzip.open(in_vcf_phased,'rt')
genotype_dic = {}

for l in fi_phased:
    if l[0:2] == "##":
        continue
    elif l[0:2] == "#C":
        l = l[0:-1]
        header = l.split("\t")
    else:
        l = l[0:-1]
        row = l.split("\t")
        for x in range(9, len(row)):
            key = row[2] + ":" + header[x]
            genotype_dic[key] = row[x]

fi_phased.close()       



fi_gp = gzip.open(in_vcf_GP,'rt')
fi_out = open(out_vcf, 'w')

for l in fi_gp:
    if l[0:2] == "##":
        fi_out.writelines(l)
    elif l[0:2] == "#C":
        fi_out.writelines(l)
        l = l[0:-1]
        header = l.split("\t")
    else:
        l = l[0:-1]
        row = l.split("\t")
        for x in range(9, len(row)):
            key = row[2] + ":" + header[x]
            tem = row[x][0:3]
            row[x] = row[x].replace(tem, genotype_dic.get(key, "0|0"))
        res = "\t".join(row)
        res = res + "\n"
        fi_out.writelines(res)

fi_gp.close()  
fi_out.close()

