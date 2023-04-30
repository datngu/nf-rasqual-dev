cd /mnt/ScratchProjects/Aqua-Faang/dat_projects/paper1/results
rm -rf res_collected
mkdir res_collected

tissue_list="Brain Gill Gonad Liver Muscle"

for tis in $tissue_list
do 
    mkdir res_collected/$tis
    
    cp $tis/ATAC_eigenMT_results_merged/ALL_eigenMT_results.txt res_collected/$tis/obs_atac_rasqual_eigenMT_pval.txt
    cp $tis/ATAC_eigenMT_results_merged_permute/ALL_eigenMT_results.txt res_collected/$tis/nul_atac_rasqual_eigenMT_pval.txt
    
    cp $tis/RNA_eigenMT_results_merged/ALL_eigenMT_results.txt res_collected/$tis/obs_rna_rasqual_eigenMT_pval.txt
    cp $tis/RNA_eigenMT_results_merged_permute/ALL_eigenMT_results.txt res_collected/$tis/nul_rna_rasqual_eigenMT_pval.txt

    cat $tis/ATAC_results_rasqual_eigenMT_processed/*txt >> res_collected/$tis/obs_atac_rasqual_normial_pval.txt
    
    python3 extract.py --normial res_collected/$tis/obs_atac_rasqual_normial_pval.txt --eigenMT res_collected/$tis/obs_atac_rasqual_eigenMT_pval.txt > res_collected/$tis/obs_atac_rasqual_normial_pval_lead_snp.txt

    cat $tis/ATAC_results_rasqual_eigenMT_processed_permute/*txt >> res_collected/$tis/nul_atac_rasqual_normial_pval.txt
    
    python3 extract.py --normial res_collected/$tis/nul_atac_rasqual_normial_pval.txt --eigenMT res_collected/$tis/nul_atac_rasqual_eigenMT_pval.txt > res_collected/$tis/nul_atac_rasqual_normial_pval_lead_snp.txt
    
    cat $tis/RNA_results_rasqual_eigenMT_processed/*txt >> res_collected/$tis/obs_rna_rasqual_normial_pval.txt
    python3 extract.py --normial res_collected/$tis/obs_rna_rasqual_normial_pval.txt --eigenMT res_collected/$tis/obs_rna_rasqual_eigenMT_pval.txt > res_collected/$tis/obs_rna_rasqual_normial_pval_lead_snp.txt

    cat $tis/RNA_results_rasqual_eigenMT_processed_permute/*txt >> res_collected/$tis/nul_rna_rasqual_normial_pval.txt

    python3 extract.py --normial res_collected/$tis/nul_rna_rasqual_normial_pval.txt --eigenMT res_collected/$tis/nul_rna_rasqual_eigenMT_pval.txt > res_collected/$tis/nul_rna_rasqual_normial_pval_lead_snp.txt

done

rm res_collected/*/*rasqual_normial_pval.txt

# get -r /mnt/ScratchProjects/Aqua-Faang/dat_projects/paper1/results/res_collected /Users/datn/DATA_ANALYSES/aqua_fang_QTL