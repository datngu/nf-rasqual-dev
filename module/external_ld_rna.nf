#!/usr/bin/env nextflow

// include { foo } from './module/external_ld_rna'


nextflow.enable.dsl=2




process EXTERNAL_LD_RNA_eigenMT_process_input {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/EXTERNAL_LD_RNA_eigenMT_process_input", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    val chr
    path split_chrom

    output:
    path ("${chr}_phenotype_position.txt")

    script:
    """
    MatrixQTL_RNA_phenotype_converter.py --count ${chr}_count.txt --out_phenotype ${chr}_phenotype.txt --out_phenotype_position ${chr}_phenotype_position.txt

    """
}





process EXTERNAL_LD_RNA_eigenMT{
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/EXTERNAL_LD_RNA_eigenMT_results", mode: 'symlink', overwrite: true
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
        --cis_dist ${params.eqtl_window} \
	    --OUT ${chr}_eigenMT_results.txt \
        --external

    """
}



process EXTERNAL_LD_RNA_MERGE_eigenMT {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/EXTERNAL_LD_RNA_eigenMT_results_merged", mode: 'copy', overwrite: true
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



process EXTERNAL_LD_RNA_eigenMT_permute {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/EXTERNAL_LD_RNA_eigenMT_permute", mode: 'symlink', overwrite: true
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
        --cis_dist ${params.eqtl_window} \
	    --OUT ${chr}_eigenMT_results.txt \
        --external

    """
}





process EXTERNAL_LD_RNA_MERGE_eigenMT_permute {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/EXTERNAL_LD_RNA_eigenMT_merged_permute", mode: 'copy', overwrite: true
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