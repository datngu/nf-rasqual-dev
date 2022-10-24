#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1                
#SBATCH --job-name=nf-rasqual   
#SBATCH --mem=4G                
#SBATCH --partition=smallmem     
#SBATCH --mail-user=nguyen.thanh.dat@nmbu.no
#SBATCH --mail-type=ALL


module load BCFtools/1.10.2-GCC-8.3.0
module load git/2.23.0-GCCcore-9.3.0-nodocs
module load Nextflow/21.03
module load singularity/rpm


cd /mnt/SCRATCH/ngda/paper1/nf-rasqual

git pull

atac_bam=/mnt/users/ngda/ngs_data/atlantic_salmon/atac/brain/bam/*{.bam,.bai}
atac_count=/mnt/users/ngda/ngs_data/atlantic_salmon/atac/brain/consensus_peaks.mLb.clN.featureCounts.txt
genotype=/mnt/users/ngda/ngs_data/atlantic_salmon/wgs/processed_all_chrom.vcf.gz

outdir=results

nextflow run main.nf -resume --atac_bam $atac_bam --atac_count $atac_count --genotype $genotype --outdir $outdir