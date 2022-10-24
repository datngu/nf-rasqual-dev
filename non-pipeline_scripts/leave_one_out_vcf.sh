#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1                
#SBATCH --job-name=vcf_loo   
#SBATCH --mem=4G                
#SBATCH --partition=gpu
#SBATCH --mail-user=nguyen.thanh.dat@nmbu.no
#SBATCH --mail-type=ALL



module load BCFtools/1.10.2-GCC-8.3.0

cd /mnt/ScratchProjects/Aqua-Faang/dat_projects/paper1/data/salmon/WGS

cp /mnt/users/ngda/ngs_data/atlantic_salmon/wgs/processed_all_chrom.vcf.gz genotype.vcf.gz

sample_list="A14_AF A15_AF A24_AF A19_AF A20_AF A22_AF J15_AF J16_AF J18_A J8_AF J9_AF J10_AF"

for sample in $sample_list
do
    bcftools view  genotype.vcf.gz -s ^${sample} -o ${sample}_genotype.vcf.gz -Oz
done



