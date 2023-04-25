#!/usr/bin/env nextflow

// include { foo } from './module/deltaSVM'


nextflow.enable.dsl=2


process EXTERNAL_LD_SPLIT_chromosome {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/EXTERNAL_LD_SPLIT_chromosome", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    val chr
    path in_vcf

    output:
    tuple path("${chr}.vcf.gz"), path("${chr}.vcf.gz.tbi")

    script:
    """
    bcftools index -t ${in_vcf}

    bcftools view ${in_vcf} --regions $chr -Oz -o ${chr}.vcf.gz
    bcftools index -t ${chr}.vcf.gz
    """
}


process EXTERNAL_LD_eigenMT_process_input {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/EXTERNAL_LD_eigenMT_process_input", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    val chr
    path split_chrom

    output:
    tuple path("${chr}_genotype.txt"), path("${chr}_genotype_position.txt")

    script:
    """
    MatrixQTL_genotype_converter.py --vcf ${chr}.vcf.gz --out_genotype ${chr}_genotype.txt --out_genotype_position ${chr}_genotype_position.txt
    
    """
}


process EXTERNAL_LD_ATAC_eigenMT_process_input {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/EXTERNAL_LD_ATAC_eigenMT_process_input", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    val chr
    path split_chrom

    output:
    path ("${chr}_phenotype_position.txt")

    script:
    """
    MatrixQTL_ATAC_phenotype_converter.py --count ${chr}_count.txt --out_phenotype ${chr}_phenotype.txt --out_phenotype_position ${chr}_phenotype_position.txt

    """
}



process EXTERNAL_LD_ATAC_RUN_rasqual_eigenMT {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_results_rasqual_eigenMT", mode: 'symlink', overwrite: true
    memory '64 GB'
    cpus 16

    input:
    val chr
    path preproces_data
    path split_chrom
    path covariates

    output:
    path("${chr}_rasqual_all_snp.txt")


    script:
    """
    echo \$HOSTNAME
    rasqual_eigenMT.R vcf=${chr}.vcf.gz y=${chr}_atac.exp.bin k=${chr}_atac.size_factors.bin x=atac.covs_all_chrom.bin x_txt=atac.covs_all_chrom.txt meta=${chr}_snp_counts.tsv out=${chr}_rasqual_all_snp.txt cpu=${task.cpus}
    """
}



process EXTERNAL_LD_ATAC_eigenMT {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/EXTERNAL_LD_ATAC_eigenMT_results", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    val chr
    path rasqual_eigenMT
    path genotype_input
    path phenotype_input

    output:
    path("${chr}_eigenMT_results.txt")

    script:
    """

    eigenMT.py --CHROM ${chr} \
	    --QTL ${chr}_formated_EigenMT.txt \
	    --GEN ${chr}_genotype.txt \
	    --GENPOS ${chr}_genotype_position.txt \
	    --PHEPOS ${chr}_phenotype_position.txt \
        --cis_dist ${params.atac_window} \
	    --OUT ${chr}_eigenMT_results.txt \
        --external

    """
}



process EXTERNAL_LD_ATAC_MERGE_eigenMT {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/EXTERNAL_LD_ATAC_eigenMT_results_merged", mode: 'copy', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    val max_chr
    path results

    output:
    path("ALL_eigenMT_results.txt")


    script:
    """
    
    cat 1_eigenMT_results.txt > ALL_eigenMT_results.txt
    ## exclude header
    for chr in \$(seq 2 $max_chr)
    do
        awk 'NR!=1' \${chr}_eigenMT_results.txt >> ALL_eigenMT_results.txt
    done
    """
}

