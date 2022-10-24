#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1                
#SBATCH --job-name=exp_loo   
#SBATCH --mem=4G                
#SBATCH --partition=gpu
#SBATCH --mail-user=nguyen.thanh.dat@nmbu.no
#SBATCH --mail-type=ALL




expression_loo=/mnt/users/ngda/sofware/nf-rasqual/non-pipeline_scripts/expression_loo.R
data=/mnt/ScratchProjects/Aqua-Faang/dat_projects/paper1/data/salmon
meta_dir=/net/fs-2/scale/OrionStore/ScratchProjects/Aqua-Faang/dat_projects/paper1/data/salmon/meta

# 1. Brain

tissue=Brain

meta=${meta_dir}/${tissue}.csv
home_dir=${data}/${tissue}

atac=${home_dir}/atac_consensus_peak_featureCounts.txt
rna=${home_dir}/rna_gene_level_count_salmon.txt

Rscript $expression_loo $meta $atac $rna


# 2. Gill

tissue=Gill

meta=${meta_dir}/${tissue}.csv
home_dir=${data}/${tissue}

atac=${home_dir}/atac_consensus_peak_featureCounts.txt
rna=${home_dir}/rna_gene_level_count_salmon.txt

Rscript $expression_loo $meta $atac $rna



# 3. Gonad

tissue=Gonad

meta=${meta_dir}/${tissue}.csv
home_dir=${data}/${tissue}

atac=${home_dir}/atac_consensus_peak_featureCounts.txt
rna=${home_dir}/rna_gene_level_count_salmon.txt

Rscript $expression_loo $meta $atac $rna


# 4. Liver

tissue=Liver

meta=${meta_dir}/${tissue}.csv
home_dir=${data}/${tissue}

atac=${home_dir}/atac_consensus_peak_featureCounts.txt
rna=${home_dir}/rna_gene_level_count_salmon.txt

Rscript $expression_loo $meta $atac $rna


# 5. Muscle

tissue=Muscle

meta=${meta_dir}/${tissue}.csv
home_dir=${data}/${tissue}

atac=${home_dir}/atac_consensus_peak_featureCounts.txt
rna=${home_dir}/rna_gene_level_count_salmon.txt

Rscript $expression_loo $meta $atac $rna
