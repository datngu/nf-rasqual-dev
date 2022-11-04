#!/bin/bash
#SBATCH --ntasks=16
#SBATCH --nodes=1                
#SBATCH --job-name=vcf   
#SBATCH --mem=32G                
#SBATCH --partition=smallmem     
#SBATCH --mail-user=nguyen.thanh.dat@nmbu.no
#SBATCH --mail-type=ALL


module load BCFtools/1.10.2-GCC-8.3.0
module load SAMtools/1.11-GCC-10.2.0

beagle4=/mnt/users/ngda/sofware/beagle.27Jan18.7e1.jar
beagle5=/mnt/users/ngda/sofware/beagle.22Jul22.46e.jar

# singularity run /mnt/users/ngda/proj/paper1/nf-rasqual/work/singularity/ndatth-rasqual-v0.0.0.img

cd /mnt/users/ngda/ngs_data/atlantic_salmon/wgs

for i in {1..29}
do
    bcftools index -t -f vcf_cp_raw/ssa${i}.DP10.GQ10.MS0.7.recode.vcf.gz
done


## phasing
mkdir vcf_phased
for i in {1..29}
do
    # estimating genotypes
    java -jar $beagle4 gtgl=vcf_cp_raw/ssa${i}.DP10.GQ10.MS0.7.recode.vcf.gz nthreads=16 niterations=20 gprobs=true out=vcf_phased/ssa${i}_tem
    # phasing
    java -jar $beagle4 gt=vcf_phased/ssa${i}_tem.vcf.gz nthreads=16 niterations=20 gprobs=true out=vcf_phased/ssa${i}
done


## add GP
for i in {1..29}
do
    /mnt/users/ngda/sofware/add_GP_vcf.py vcf_phased/ssa${i}_tem.vcf.gz vcf_phased/ssa${i}.vcf.gz vcf_phased/ssa${i}_added_GP.vcf
    bgzip -f vcf_phased/ssa${i}_added_GP.vcf
done

for i in {1..29}
do
    #echo vcf_phased/ssa${i}_added_GP.vcf.gz >> file_list.txt
done

bcftools concat -n -f file_list.txt -Oz -o all_chr_added_GP.vcf.gz

# fai=/mnt/users/ngda/genomes/atlantic_salmon/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa.fai

# bcftools reheader --fai $fai all_chr_added_GP.vcf.gz | bcftools sort -Oz > sorted_all_chr_added_GP.vcf.gz
# bftools index -t sorted_all_chr_added_GP.vcf.gz

# bftools index -t sorted_all_chr_added_GP.vcf.gz

zcat all_chr_added_GP.vcf.gz | sed 's/ssa0//g' | sed 's/ssa//g' | bgzip > processed_all_chrom.vcf.gz

bftools index -t processed_all_chrom.vcf.gz