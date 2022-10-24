# working dir
/mnt/users/ngda/proj/paper1
# attac dir
#cd /mnt/project/Aqua-Faang/seq_data/AtlanticSalmon/ATAC

############ WGS
cd /mnt/users/ngda/ngs_data/atlantic_salmon/wgs



########### ATAC BAM
## tissues code:
	# tissueCode <- c(
	# 	"1" = "Brain",
	# 	"2" = "Gill",
	# 	"3" = "Gonad",
	# 	"4" = "HeadKidney",
	# 	"5" = "Liver",
	# 	"6" = "Muscle",
	# 	"7" = "SuppDistalIntestine",
	# 	"8" = "Sperm"
	# )
#


cd /mnt/ScratchProjects/Aqua-Faang/dat_projects/paper1/data/salmon


# 1. brain
home_dir=$PWD/brain
mkdir $home_dir

cp /mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/ATAC/Brain/old_results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.featureCounts.txt ${home_dir}/atac_consensus_peak_featureCounts.txt

cp /mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/RNA/Brain/results/star_salmon/salmon.merged.gene_counts.tsv ${home_dir}/rna_gene_level_count_salmon.txt


#cp /mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/RNA/Brain/results/star_salmon/salmon.merged.gene_tpm.tsv ${home_dir}/rna_gene_level_tpm_salmon.txt


## ATAC_seq
mkdir ${home_dir}/atac_bam

source_data=/mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/ATAC/Brain/old_results/bwa/mergedLibrary

cd $source_data
for fn in *bam*
do
	cp $fn ${home_dir}/atac_bam/$fn
done

## RNA seq
mkdir ${home_dir}/rna_bam

source_data=/mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/RNA/Brain/results/star_salmon

cd $source_data
for fn in *bam*
do
	cp $fn ${home_dir}/rna_bam/$fn
done



# 2. gill


home_dir=$PWD/gill
mkdir $home_dir

cp /mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/ATAC/Gill/old_results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.featureCounts.txt ${home_dir}/atac_consensus_peak_featureCounts.txt

cp /mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/RNA/Gill/results/star_salmon/salmon.merged.gene_counts.tsv ${home_dir}/rna_gene_level_count_salmon.txt



## ATAC_seq
mkdir ${home_dir}/atac_bam

source_data=/mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/ATAC/Gill/old_results/bwa/mergedLibrary

cd $source_data
for fn in *bam*
do
	cp $fn ${home_dir}/atac_bam/$fn
done

## RNA seq
mkdir ${home_dir}/rna_bam

source_data=/mnt/project/Aqua-Faang/seq_results/AtlanticSalmon/BodyMap/RNA/Gill/results/star_salmon

cd $source_data
for fn in *bam*
do
	cp $fn ${home_dir}/rna_bam/$fn
done


