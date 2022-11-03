# aqua faang dir
cd /mnt/project/Aqua-Faang/seq_data/AtlanticSalmon/ATAC




###############################
# working dir
#module load BCFtools/1.10.2-GCC-8.3.0
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


########### debugs
# singularity build ./ndatth-rasqual-v0.0.0.img docker://ndatth/rasqual:v0.0.0

# singularity run /mnt/SCRATCH/ngda/nf-rasqual/shared_dir/singularity/ndatth-rasqual-v0.0.0.img


RASQUALDIR=/rasqual
/mnt/users/ngda/proj/paper1/nf-rasqual/bin/createASVCF_fixed_path.sh paired_end bam_list.txt processed_all_chrom.vcf.gz out.vcf.gz atac
