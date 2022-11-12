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
