#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1                
#SBATCH --job-name=test   
#SBATCH --mem=4G                
#SBATCH --partition=gpu
#SBATCH --mail-user=nguyen.thanh.dat@nmbu.no
#SBATCH --mail-type=ALL


module load BCFtools/1.10.2-GCC-8.3.0
module load git/2.23.0-GCCcore-9.3.0-nodocs
module load Nextflow/21.03
module load singularity/rpm

# export PATH=$PATH:/Users/datn/Downloads/plink_mac_20220402
# git clone https://github.com/datngu/nf-rasqual-dev.git
#cd /mnt/SCRATCH/ngda/nf-rasqual-dev
git pull

#
genome=/mnt/users/ngda/genomes/atlantic_salmon/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa
annotation=/mnt/users/ngda/genomes/atlantic_salmon/Salmo_salar.Ssal_v3.1.106.gtf


# ATAC seq input
atac_bam=/mnt/users/ngda/ngs_data/atlantic_salmon/brain/atac_bam/*{.bam,.bai}
atac_count=/mnt/users/ngda/ngs_data/atlantic_salmon/brain/atac_consensus_peak_featureCounts.txt

# RNA input
rna_bam=/mnt/users/ngda/ngs_data/atlantic_salmon/brain/rna_bam/*{.bam,.bai}
rna_count=/mnt/users/ngda/ngs_data/atlantic_salmon/brain/rna_gene_level_count_salmon.txt

# SNP genotype input - phased - added GP
genotype=/mnt/users/ngda/ngs_data/atlantic_salmon/wgs/processed_all_chrom.vcf.gz

##
outdir=results

export NXF_SINGULARITY_CACHEDIR=/mnt/users/ngda/sofware/singularity
##
nextflow run main.nf -resume --genome $genome --annotation $annotation --atac_bam $atac_bam --atac_count $atac_count --rna_bam $rna_bam --rna_count $rna_count --genotype $genotype --outdir $outdir