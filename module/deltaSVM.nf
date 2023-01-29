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
    slipt_bed.R $atac_filtered $params.deltaSVM_folds atac_
    """
}


process ATAC_deltaSVM_gen_null_seqs {

    container 'ndatth/delta-svm:v0.0.0'
    publishDir "${params.outdir}/deltaSVM", mode: 'symlink', overwrite: true
    memory '64 GB'
    cpus 16

    input:
    path trf_bed
    path list_atac_bed
    
    output:
    path "*.fa"

    script:
    """
    for i in \$(seq 1 $params.deltaSVM_folds)
    do
        genNullSeqs_tfr.R bed=atac_\$i.txt trf=$trf_bed bsgenome="BSgenome.Salmo.Salar.Ensembl.106" xfold=1 out_prefix=\$i batchsize=20000 &
    done
    wait
    """
}