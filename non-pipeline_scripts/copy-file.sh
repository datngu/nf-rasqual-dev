#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1                
#SBATCH --job-name=cpfiles   
#SBATCH --mem=4G                
#SBATCH --partition=smallmem     
#SBATCH --mail-user=nguyen.thanh.dat@nmbu.no
#SBATCH --mail-type=ALL




cd /mnt/ScratchProjects/Aqua-Faang/dat_projects/paper1/data/salmon


# # 1. Brain
# tissue=Brain
# cd /mnt/ScratchProjects/Aqua-Faang/dat_projects/paper1/data/salmon

# home_dir=$PWD/${tissue}
# mkdir $home_dir

# cp /mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/ATAC/${tissue}/old_results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.featureCounts.txt ${home_dir}/atac_consensus_peak_featureCounts.txt

# cp /mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/RNA/${tissue}/results/star_salmon/salmon.merged.gene_counts.tsv ${home_dir}/rna_gene_level_count_salmon.txt


# ## ATAC_seq
# mkdir ${home_dir}/atac_bam

# source_data=/mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/ATAC/${tissue}/old_results/bwa/mergedLibrary

# cd $source_data
# for fn in *bam*
# do
# 	cp $fn ${home_dir}/atac_bam/$fn
# done

# ## RNA seq
# mkdir ${home_dir}/rna_bam

# source_data=/mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/RNA/${tissue}/results/star_salmon

# cd $source_data
# for fn in *bam*
# do
# 	cp $fn ${home_dir}/rna_bam/$fn
# done



# 2. Gill

tissue=Gill
cd /mnt/ScratchProjects/Aqua-Faang/dat_projects/paper1/data/salmon

home_dir=$PWD/${tissue}
mkdir $home_dir

cp /mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/ATAC/${tissue}/old_results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.featureCounts.txt ${home_dir}/atac_consensus_peak_featureCounts.txt

cp /mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/RNA/${tissue}/results/star_salmon/salmon.merged.gene_counts.tsv ${home_dir}/rna_gene_level_count_salmon.txt


## ATAC_seq
mkdir ${home_dir}/atac_bam

source_data=/mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/ATAC/${tissue}/old_results/bwa/mergedLibrary

cd $source_data
for fn in *bam*
do
	cp $fn ${home_dir}/atac_bam/$fn
done

## RNA seq
mkdir ${home_dir}/rna_bam

source_data=/mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/RNA/${tissue}/results/star_salmon

cd $source_data
for fn in *bam*
do
	cp $fn ${home_dir}/rna_bam/$fn
done



# 3. Gonad

tissue=Gonad
cd /mnt/ScratchProjects/Aqua-Faang/dat_projects/paper1/data/salmon

home_dir=$PWD/${tissue}
mkdir $home_dir

cp /mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/ATAC/${tissue}/old_results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.featureCounts.txt ${home_dir}/atac_consensus_peak_featureCounts.txt

cp /mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/RNA/${tissue}/results/star_salmon/salmon.merged.gene_counts.tsv ${home_dir}/rna_gene_level_count_salmon.txt


## ATAC_seq
mkdir ${home_dir}/atac_bam

source_data=/mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/ATAC/${tissue}/old_results/bwa/mergedLibrary

cd $source_data
for fn in *bam*
do
	cp $fn ${home_dir}/atac_bam/$fn
done

## RNA seq
mkdir ${home_dir}/rna_bam

source_data=/mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/RNA/${tissue}/results/star_salmon

cd $source_data
for fn in *bam*
do
	cp $fn ${home_dir}/rna_bam/$fn
done




# 4. Liver

tissue=Liver
cd /mnt/ScratchProjects/Aqua-Faang/dat_projects/paper1/data/salmon

home_dir=$PWD/${tissue}
mkdir $home_dir

cp /mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/ATAC/${tissue}/old_results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.featureCounts.txt ${home_dir}/atac_consensus_peak_featureCounts.txt

cp /mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/RNA/${tissue}/results/star_salmon/salmon.merged.gene_counts.tsv ${home_dir}/rna_gene_level_count_salmon.txt


## ATAC_seq
mkdir ${home_dir}/atac_bam

source_data=/mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/ATAC/${tissue}/old_results/bwa/mergedLibrary

cd $source_data
for fn in *bam*
do
	cp $fn ${home_dir}/atac_bam/$fn
done

## RNA seq
mkdir ${home_dir}/rna_bam

source_data=/mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/RNA/${tissue}/results/star_salmon

cd $source_data
for fn in *bam*
do
	cp $fn ${home_dir}/rna_bam/$fn
done


# 5. Muscle

tissue=Muscle
cd /mnt/ScratchProjects/Aqua-Faang/dat_projects/paper1/data/salmon


home_dir=$PWD/${tissue}
mkdir $home_dir

cp /mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/ATAC/${tissue}/old_results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.featureCounts.txt ${home_dir}/atac_consensus_peak_featureCounts.txt

cp /mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/RNA/${tissue}/results/star_salmon/salmon.merged.gene_counts.tsv ${home_dir}/rna_gene_level_count_salmon.txt


## ATAC_seq
mkdir ${home_dir}/atac_bam

source_data=/mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/ATAC/${tissue}/old_results/bwa/mergedLibrary

cd $source_data
for fn in *bam*
do
	cp $fn ${home_dir}/atac_bam/$fn
done

## RNA seq
mkdir ${home_dir}/rna_bam

source_data=/mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/RNA/${tissue}/results/star_salmon

cd $source_data
for fn in *bam*
do
	cp $fn ${home_dir}/rna_bam/$fn
done

