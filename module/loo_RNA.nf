#!/usr/bin/env nextflow

// include { foo } from './module/loo_RNA'


nextflow.enable.dsl=2

process LOO_rna {

    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/loo_rna", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    path meta
    path rna_count_filtered

    output:
    path "*_rna_count.txt"

    script:
    """
    loo_RNA.R $meta $rna_count_filtered
    """
}


process LOO_rna_vcf {

    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/loo_rna_vcf", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    path meta
    path vcf

    output:
    path("*_loo.vcf.*")


    script:
    """
    loo_VCF.R $meta "processed.vcf.gz"
    """
}

process LOO_RNA_PROCESS_covariates {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/loo_RNA_covariates", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    val ID
    path meta
    path rna_count

    output:
    tuple path("${ID}_rna.covs_all_chrom.bin"), path("${ID}_rna.covs_all_chrom.txt")

    script:
    """
    RNA_covariates.R ${ID}.csv ${ID}_rna_count.txt $params.phenotype_PCs

    mv rna.covs_all_chrom.bin ${ID}_rna.covs_all_chrom.bin
    mv rna.covs_all_chrom.txt ${ID}_rna.covs_all_chrom.txt
    """
}



process LOO_RNA_SPLIT_chromosome {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/loo_RNA_split_chrom", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    tuple val(ID), val (chr)
    path in_vcf
    path in_exp

    output:
    tuple path("${ID}_${chr}.vcf.gz"), path("${ID}_${chr}.vcf.gz.tbi"), path("${ID}_${chr}_count.txt")

    script:
    """
    awk 'NR==1{print }' ${ID}_rna_count.txt > ${ID}_${chr}_count.txt
    awk -v chr=$chr '{ if (\$2 == $chr) { print } }' ${ID}_rna_count.txt >> ${ID}_${chr}_count.txt

    bcftools view ${ID}_loo.vcf.gz --regions ${chr} -Oz -o ${ID}_${chr}.vcf.gz
    bcftools index -t ${ID}_${chr}.vcf.gz
    """
}




process LOO_RNA_PREPROCESS_rasqual {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/loo_RNA_qtl_input", mode: 'symlink', overwrite: true
    memory '64 GB'
    cpus 8

    input:
    tuple val(ID), val (chr)
    path split_chrom
    path genome

    output:
    tuple path("${ID}_${chr}_rna.exp.bin"), path("${ID}_${chr}_rna.exp.txt"), path("${ID}_${chr}_rna.size_factors.bin"), path("${ID}_${chr}_rna.size_factors.txt"), path("${ID}_${chr}_snp_counts.tsv")


    script:
    """

    RNA_rasqual_processor.R ${ID}_${chr}_count.txt ${ID}_${chr}.vcf.gz $genome $params.eqtl_window ${task.cpus}
    ## rename files
    mv rna.exp.bin ${ID}_${chr}_rna.exp.bin
    mv rna.exp.txt ${ID}_${chr}_rna.exp.txt
    mv rna.size_factors.bin ${ID}_${chr}_rna.size_factors.bin
    mv rna.size_factors.txt ${ID}_${chr}_rna.size_factors.txt
    mv snp_counts.tsv ${ID}_${chr}_snp_counts.tsv

    """
}


process LOO_RNA_RUN_rasqual_eigenMT {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/loo_RNA_results_rasqual_eigenMT", mode: 'symlink', overwrite: true
    memory '64 GB'
    cpus 16

    input:
    tuple val(ID), val (chr)
    path preproces_data
    path split_chrom
    path covariates

    output:
    tuple val("${ID}"), path("${ID}_${chr}_rasqual_all_snp.txt")


    script:
    """
    echo \$HOSTNAME
    rasqual_eigenMT.R vcf=${ID}_${chr}.vcf.gz y=${ID}_${chr}_rna.exp.bin k=${ID}_${chr}_rna.size_factors.bin x=${ID}_rna.covs_all_chrom.bin x_txt=${ID}_rna.covs_all_chrom.txt meta=${ID}_${chr}_snp_counts.tsv out=${ID}_${chr}_rasqual_all_snp.txt cpu=${task.cpus}
    """
}




process LOO_RNA_rasqual_TO_eigenMT {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/loo_RNA_results_rasqual_eigenMT_processed", mode: 'symlink', overwrite: true
    memory '64 GB'
    cpus 16

    input:
    tuple val(ID), val (chr)
    path rasqual_res

    output:
    path("${ID}_${chr}_formated_EigenMT.txt")


    script:
    """
    echo \$HOSTNAME
    rasqualToEigenMT.py --rasqualOut ${ID}_${chr}_rasqual_all_snp.txt > ${ID}_${chr}_formated_EigenMT.txt
    """
}



process LOO_RNA_eigenMT_process_input {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/loo_RNA_eigenMT_process_input", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    tuple val(ID), val (chr)
    path split_chrom

    output:
    tuple path("${ID}_${chr}_genotype.txt"), path("${ID}_${chr}_genotype_position.txt"), path("${ID}_${chr}_phenotype_position.txt")

    script:
    """
    MatrixQTL_genotype_converter.py --vcf ${ID}_${chr}.vcf.gz --out_genotype ${ID}_${chr}_genotype.txt --out_genotype_position ${ID}_${chr}_genotype_position.txt

    MatrixQTL_RNA_phenotype_converter.py --count ${ID}_${chr}_count.txt --out_phenotype ${ID}_${chr}_phenotype.txt --out_phenotype_position ${ID}_${chr}_phenotype_position.txt

    """
}



process LOO_RNA_eigenMT {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/loo_RNA_eigenMT_results", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    tuple val(ID), val (chr)
    path rasqual_eigenMT
    path all_input

    output:
    tuple val("${ID}"), path("${ID}_${chr}_eigenMT_results.txt")

    script:
    """
    
    eigenMT.py --CHROM ${chr} \
	    --QTL ${ID}_${chr}_formated_EigenMT.txt \
	    --GEN ${ID}_${chr}_genotype.txt \
	    --GENPOS ${ID}_${chr}_genotype_position.txt \
	    --PHEPOS ${ID}_${chr}_phenotype_position.txt \
        --cis_dist ${params.eqtl_window} \
	    --OUT ${ID}_${chr}_eigenMT_results.txt

    """
}



process LOO_RNA_MERGE_eigenMT {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/loo_RNA_eigenMT_results_merged", mode: 'copy', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    val max_chr
    tuple val(ID), path(rasqual_results)

    output:
    path("${ID}_ALL_eigenMT_results.txt")


    script:
    """
    
    cat ${ID}_1_eigenMT_results.txt > ${ID}_ALL_eigenMT_results.txt
    ## exclude header
    for chr in \$(seq 2 $max_chr)
    do
        awk 'NR!=1' ${ID}_\${chr}_eigenMT_results.txt >> ${ID}_ALL_eigenMT_results.txt
    done
    
    """
}
