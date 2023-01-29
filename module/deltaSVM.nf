#!/usr/bin/env nextflow

// include { foo } from './module/deltaSVM'


nextflow.enable.dsl=2


process ATAC_deltaSVM_slipt_bed {

    container 'ndatth/delta-svm:v0.0.0'
    publishDir "${params.outdir}/deltaSVM", mode: 'symlink', overwrite: true
    memory '16 GB'

    input:
    path atac_filtered
    
    output:
    path "atac_*.txt"

    script:
    """
    slipt_bed.R atac_filtered $params.deltaSVM_folds atac_
    """
}


process ATAC_deltaSVM_gen_null_seqs {

    container 'ndatth/delta-svm:v0.0.0'
    publishDir "${params.outdir}/deltaSVM", mode: 'symlink', overwrite: true
    memory '16 GB'

    input:
    path atac_filtered
    
    output:
    path "atac_*.txt"

    script:
    """
    slipt_bed.R $atac_filtered $params.deltaSVM_folds atac_
    """
}