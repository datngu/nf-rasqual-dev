#!/usr/bin/env nextflow
/*
========================================================================================
                          nf-rasqual
========================================================================================
                RASQUAL Analysis Pipeline with nextflow.
                https://github.com/datngu/nf-rasqual
                Author: Dat T Nguyen
                Contact: ndat<at>utexas.edu
----------------------------------------------------------------------------------------
*/





/*
 Define the default parameters
*/ 
params.genome          = "$baseDir/data/ref/genome.fa"
params.annotation      = "$baseDir/data/ref/annotation.gtf"
params.atac_bam        = "$baseDir/data/atac_bam/*.bam"
params.atac_count      = "$baseDir/data/atac_consensus_peak_featureCounts.txt"
params.rna_bam         = "$baseDir/data/rna_bam/*.bam"
params.rna_count       = "$baseDir/data/rna_gene_level_count_salmon.txt"
params.genotype        = "$baseDir/data/genotype.vcf.gz"
params.ld_genotype     = "$baseDir/data/ld_genotype.vcf.gz"
params.meta            = "$baseDir/data/meta/Brain.csv"
params.trf_bed         = "$baseDir/data/atlantic_salmon_v3.1_trf.bed"
params.outdir          = "results"
params.trace_dir       = "trace_dir"

// running options
params.chrom           = 1..29 
params.permute         = 1..12
params.phenotype_PCs   = 2 
params.exp_prop        = 0.5
params.fpkm_cutoff     = 0.5
params.maf             = 0.05
params.fdr             = 0.1
params.atac_window     = 10000
params.eqtl_window     = 500000
params.deltaSVM_folds  = 10

// pipeline options

params.atac_qtl          = true
params.eqtl_qtl          = true
external_ld              = false
params.deltaSVM          = true

// control for FDR
params.loo               = false
params.eigenMT_fdr       = true



log.info """\
================================================================
                        nf-rasqual
================================================================
    genome              : $params.genome
    annotation          : $params.annotation
    atac_bam            : $params.atac_bam
    atac_count          : $params.atac_count
    rna_bam             : $params.rna_bam
    rna_count           : $params.rna_count
    genotype            : $params.genotype
    ld_genotype         : $params.ld_genotype 
    meta                : $params.meta
    outdir              : $params.outdir
    trace_dir           : $params.trace_dir

    
    fpkm_cutoff         : $params.fpkm_cutoff
    chrom               : $params.chrom
    permute             : $params.permute
    maf                 : $params.maf
    fdr                 : $params.fdr
    eqtl_window         : $params.eqtl_window
    atac_window         : $params.atac_window
    phenotype_PCs       : $params.phenotype_PCs

    atac_qtl            : $params.atac_qtl
    eqtl_qtl            : $params.eqtl_qtl
    external_ld         : $params.external_ld
    deltaSVM            : $params.deltaSVM
================================================================
"""

nextflow.enable.dsl=2

// include module loo_RNA

include { LOO_rna_vcf;  LOO_rna; LOO_RNA_PROCESS_covariates; LOO_RNA_SPLIT_chromosome; LOO_RNA_PREPROCESS_rasqual; LOO_RNA_RUN_rasqual_eigenMT; LOO_RNA_rasqual_TO_eigenMT; LOO_RNA_eigenMT_process_input; LOO_RNA_eigenMT; LOO_RNA_MERGE_eigenMT} from './module/loo_RNA'

include { LOO_atac_vcf;  LOO_atac; LOO_ATAC_PROCESS_covariates; LOO_ATAC_SPLIT_chromosome; LOO_ATAC_PREPROCESS_rasqual; LOO_ATAC_RUN_rasqual_eigenMT; LOO_ATAC_rasqual_TO_eigenMT; LOO_ATAC_eigenMT_process_input; LOO_ATAC_eigenMT; LOO_ATAC_MERGE_eigenMT} from './module/loo_ATAC'


include { ATAC_deltaSVM_slipt_bed; ATAC_deltaSVM_gen_null_seqs; ATAC_deltaSVM_train } from './module/deltaSVM'

include { EXTERNAL_LD_SPLIT_chromosome; EXTERNAL_LD_eigenMT_process_input; EXTERNAL_LD_ATAC_eigenMT_process_input; EXTERNAL_LD_ATAC_eigenMT; EXTERNAL_LD_ATAC_MERGE_eigenMT; EXTERNAL_LD_ATAC_eigenMT_permute; EXTERNAL_LD_ATAC_MERGE_eigenMT_permute } from './module/external_ld_atac'


//include {  } from './module/external_ld'


workflow {

    // channel general processing
    chrom_list_ch = channel.from(params.chrom)
    permute_ch = channel.from(params.permute)
    VCF_filtering(params.genotype, params.meta)


    // Leave one out implementation
    if(params.loo){
        LOO_meta_csv(params.meta)
        ID_ch = LOO_meta_csv.out.map{ infile -> tuple( infile.baseName) }.flatten()
    }

    // external LD
    if(params.external_ld){
        EXTERNAL_LD_SPLIT_chromosome(chrom_list_ch, params.ld_genotype)
    }


    // ATAC QTL
    if( params.atac_qtl ){
        atac_bam_ch = channel.fromPath( params.atac_bam, checkIfExists: true )
        ATAC_BAM_rename(params.meta, atac_bam_ch.collect())
        ATAC_ADD_AS_vcf(VCF_filtering.out, ATAC_BAM_rename.out)

        ATAC_FILTERING_expression(params.atac_count, params.meta)
        ATAC_PROCESS_covariates(params.meta, ATAC_FILTERING_expression.out)
        ATAC_SPLIT_chromosome(chrom_list_ch, ATAC_ADD_AS_vcf.out, ATAC_FILTERING_expression.out )
        ATAC_PREPROCESS_rasqual(chrom_list_ch, ATAC_SPLIT_chromosome.out.collect(), params.genome)
       
        // standard run
        ATAC_RUN_rasqual_eigenMT(chrom_list_ch, ATAC_PREPROCESS_rasqual.out.collect(), ATAC_SPLIT_chromosome.out.collect(), ATAC_PROCESS_covariates.out)
        ATAC_rasqual_TO_eigenMT(chrom_list_ch, ATAC_RUN_rasqual_eigenMT.out.collect())
        // permutation run
        ATAC_RUN_rasqual_eigenMT_permute(chrom_list_ch, ATAC_PREPROCESS_rasqual.out.collect(), ATAC_SPLIT_chromosome.out.collect(), ATAC_PROCESS_covariates.out)
        ATAC_rasqual_TO_eigenMT_permute(chrom_list_ch, ATAC_RUN_rasqual_eigenMT_permute.out.collect())        
        
        // eigenMT
        if(params.external_ld){
            EXTERNAL_LD_eigenMT_process_input(chrom_list_ch, EXTERNAL_LD_SPLIT_chromosome.out.collect())
            EXTERNAL_LD_ATAC_eigenMT_process_input(chrom_list_ch, ATAC_SPLIT_chromosome.out.collect())

            EXTERNAL_LD_ATAC_eigenMT(chrom_list_ch, ATAC_rasqual_TO_eigenMT.out.collect(), EXTERNAL_LD_eigenMT_process_input.out.collect(), EXTERNAL_LD_ATAC_eigenMT_process_input.out.collect())

            EXTERNAL_LD_ATAC_MERGE_eigenMT(chrom_list_ch.max(), EXTERNAL_LD_ATAC_eigenMT.out.collect())

            EXTERNAL_LD_ATAC_eigenMT_permute(chrom_list_ch, ATAC_rasqual_TO_eigenMT_permute.out.collect(), EXTERNAL_LD_eigenMT_process_input.out.collect(), EXTERNAL_LD_ATAC_eigenMT_process_input.out.collect())

            EXTERNAL_LD_ATAC_MERGE_eigenMT_permute(chrom_list_ch.max(), EXTERNAL_LD_ATAC_eigenMT_permute.out.collect())

        }else{

            ATAC_eigenMT_process_input(chrom_list_ch, ATAC_SPLIT_chromosome.out.collect())
            ATAC_eigenMT(chrom_list_ch, ATAC_rasqual_TO_eigenMT.out.collect(), ATAC_eigenMT_process_input.out.collect())
            ATAC_MERGE_eigenMT(chrom_list_ch.max(), ATAC_eigenMT.out.collect())
        
            ATAC_eigenMT_permute(chrom_list_ch, ATAC_rasqual_TO_eigenMT_permute.out.collect(), ATAC_eigenMT_process_input.out.collect())
            ATAC_MERGE_eigenMT_permute(chrom_list_ch.max(), ATAC_eigenMT_permute.out.collect())
        }



        // deltaSVM
        if(params.deltaSVM){

            ATAC_deltaSVM_slipt_bed(ATAC_FILTERING_expression.out)
            ATAC_deltaSVM_gen_null_seqs(params.trf_bed, ATAC_deltaSVM_slipt_bed.out)
            ATAC_deltaSVM_train(ATAC_deltaSVM_gen_null_seqs.out)

        }    
    }

    if( params.eqtl_qtl ){
        rna_bam_ch = channel.fromPath( params.rna_bam, checkIfExists: true )
        RNA_BAM_rename(params.meta, rna_bam_ch.collect())
        RNA_ADD_AS_vcf(VCF_filtering.out, RNA_BAM_rename.out)
        GTF_GENE_INFO_parser(params.annotation)
    
        RNA_FILTERING_expression(params.rna_count, GTF_GENE_INFO_parser.out, params.meta)
        RNA_PROCESS_covariates(params.meta, RNA_FILTERING_expression.out)
        RNA_SPLIT_chromosome(chrom_list_ch, RNA_ADD_AS_vcf.out, RNA_FILTERING_expression.out )
        RNA_PREPROCESS_rasqual(chrom_list_ch, RNA_SPLIT_chromosome.out.collect(), params.genome)
        
        // eigenMT part
        RNA_RUN_rasqual_eigenMT(chrom_list_ch, RNA_PREPROCESS_rasqual.out.collect(), RNA_SPLIT_chromosome.out.collect(), RNA_PROCESS_covariates.out)

        RNA_rasqual_TO_eigenMT(chrom_list_ch, RNA_RUN_rasqual_eigenMT.out.collect())
            
        RNA_eigenMT_process_input(chrom_list_ch, RNA_SPLIT_chromosome.out.collect())

        RNA_eigenMT(chrom_list_ch, RNA_rasqual_TO_eigenMT.out.collect(), RNA_eigenMT_process_input.out.collect())

        RNA_MERGE_eigenMT(chrom_list_ch.max(), RNA_eigenMT.out.collect())

        // eigenMT permutation - for FDR checking

        RNA_RUN_rasqual_eigenMT_permute(chrom_list_ch, RNA_PREPROCESS_rasqual.out.collect(), RNA_SPLIT_chromosome.out.collect(), RNA_PROCESS_covariates.out)

        RNA_rasqual_TO_eigenMT_permute(chrom_list_ch, RNA_RUN_rasqual_eigenMT_permute.out.collect())

        RNA_eigenMT_permute(chrom_list_ch, RNA_rasqual_TO_eigenMT_permute.out.collect(), RNA_eigenMT_process_input.out.collect())

        RNA_MERGE_eigenMT_permute(chrom_list_ch.max(), RNA_eigenMT_permute.out.collect())

    }
}


// vcf filtering

// GENOTYPE PROCESSING
process VCF_filtering { 
    
    publishDir "${params.trace_dir}/vcf_filtering", mode: 'symlink', overwrite: true
    container 'ndatth/rasqual:v0.0.0'
    memory '8 GB'
    
    input:
    path in_vcf
    path meta
 
    output:
    path("genotype_filtered.vcf.gz")
 
    script:
    """
    # in_vcf=genotype.vcf.gz
    # meta=meta/brain.csv
    bcftools index -t $in_vcf
    grep -v ^genotype_id ${meta} | cut -d , -f 1 > genotype_sample.txt
    bcftools view $in_vcf -S genotype_sample.txt | sed 's/chr//g' | bgzip > genotype_filtered.vcf.gz
    """
}



// rename BAM

process ATAC_BAM_rename {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_bam_dir", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    path meta
    path bamfiles

    output:
    path "copied_files/*{.bam,.bai}"

    script:
    """
    ATAC_rename_bam.R ${meta}
    """
}

process RNA_BAM_rename {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_bam_dir", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    path meta
    path bamfiles

    output:
    path "copied_files/*{.bam,.bai}"

    script:
    """
    RNA_rename_bam.R ${meta}
    """
}




// add allel specific inforation

process ATAC_ADD_AS_vcf {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_AS_vcf", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    path in_vcf
    path bamfiles

    output:
    tuple path("processed.vcf.gz"), path("processed.vcf.gz.tbi")

    script:
    """
    ls \$PWD/*bam > bam_list.txt
    zcat $in_vcf | sed 's/ssa0//g' | sed 's/ssa//g' | bgzip > tem.vcf.gz
    bcftools index -t tem.vcf.gz
    createASVCF_fixed_path.sh paired_end bam_list.txt tem.vcf.gz processed.vcf.gz rna
    bcftools index -t processed.vcf.gz
    rm tem.vcf.gz tem.vcf.gz.tbi
    """
}


process RNA_ADD_AS_vcf {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_AS_vcf", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    path in_vcf
    path bamfiles

    output:
    tuple path("processed.vcf.gz"), path("processed.vcf.gz.tbi")

    script:
    """
    ls \$PWD/*bam > bam_list.txt
    zcat $in_vcf | sed 's/ssa0//g' | sed 's/ssa//g' | bgzip > tem.vcf.gz
    bcftools index -t tem.vcf.gz
    createASVCF_fixed_path.sh paired_end bam_list.txt tem.vcf.gz processed.vcf.gz atac
    bcftools index -t processed.vcf.gz
    rm tem.vcf.gz tem.vcf.gz.tbi
    """
}

// expression filtering


process ATAC_FILTERING_expression {

    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_filtering_expression", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    path atac_count
    path meta

    output:
    path "atac_consensus_peak_featureCounts_filtered.txt"

    script:
    """
    ATAC_filtering.R $atac_count atac_consensus_peak_featureCounts_filtered.txt $params.exp_prop $params.fpkm_cutoff $meta
    """
}



process GTF_GENE_INFO_parser {

    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/gene_info", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    path gtf

    output:
    path "gene_info.txt"

    script:
    """
    gtf_gene_parser.py $gtf gene_info.txt
    """
}


process RNA_FILTERING_expression {

    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_filtering_expression", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    path rna_count
    path gene_info
    path meta

    output:
    path "rna_gene_level_count_salmon_filtered.txt"

    script:
    """
    RNA_filtering.R $rna_count rna_gene_level_count_salmon_filtered.txt $gene_info $params.exp_prop $params.fpkm_cutoff $meta
    """
}




// PCA

process ATAC_PROCESS_covariates {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_covariates", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    path meta
    path atac_count

    output:
    tuple path("atac.covs_all_chrom.bin"), path("atac.covs_all_chrom.txt")

    script:
    """
    ATAC_covariates.R $meta $atac_count $params.phenotype_PCs
    """
}


process RNA_PROCESS_covariates {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_covariates", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    path meta
    path rna_count_filtered

    output:
    tuple path("rna.covs_all_chrom.bin"), path("rna.covs_all_chrom.txt")

    script:
    """
    RNA_covariates.R $meta $rna_count_filtered $params.phenotype_PCs
    """
}



// slipt chomosome

process ATAC_SPLIT_chromosome {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_split_chrom", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    val chr
    path in_vcf
    path in_exp

    output:
    tuple path("${chr}.vcf.gz"), path("${chr}.vcf.gz.tbi"), path("${chr}_count.txt")

    script:
    """
    awk 'NR==1{print }' $in_exp > ${chr}_count.txt
    awk -v chr=$chr '{ if (\$2 == $chr) { print } }' $in_exp >> ${chr}_count.txt

    bcftools view processed.vcf.gz --regions $chr -Oz -o ${chr}.vcf.gz
    bcftools index -t ${chr}.vcf.gz
    """
}


process RNA_SPLIT_chromosome {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_split_chrom", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    val chr
    path in_vcf
    path in_exp

    output:
    tuple path("${chr}.vcf.gz"), path("${chr}.vcf.gz.tbi"), path("${chr}_count.txt")

    script:
    """
    awk 'NR==1{print }' $in_exp > ${chr}_count.txt
    awk -v chr=$chr '{ if (\$2 == $chr) { print } }' $in_exp >> ${chr}_count.txt
    
    bcftools view processed.vcf.gz --regions $chr -Oz -o ${chr}.vcf.gz
    bcftools index -t ${chr}.vcf.gz
    """
}










// proprocessing


process ATAC_PREPROCESS_rasqual {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_qtl_input", mode: 'symlink', overwrite: true
    memory '64 GB'
    cpus 8

    input:
    val chr
    path split_chrom
    path genome

    output:
    tuple path("${chr}_atac.exp.bin"), path("${chr}_atac.exp.txt"), path("${chr}_atac.size_factors.bin"), path("${chr}_atac.size_factors.txt"), path("${chr}_snp_counts.tsv")


    script:
    """
    ATAC_rasqual_processor.R ${chr}_count.txt ${chr}.vcf.gz $genome $params.atac_window ${task.cpus}
    ## rename files
    mv atac.exp.bin ${chr}_atac.exp.bin
    mv atac.exp.txt ${chr}_atac.exp.txt
    mv atac.size_factors.bin ${chr}_atac.size_factors.bin
    mv atac.size_factors.txt ${chr}_atac.size_factors.txt
    mv snp_counts.tsv ${chr}_snp_counts.tsv
    """
}



process RNA_PREPROCESS_rasqual {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_qtl_input", mode: 'symlink', overwrite: true
    memory '64 GB'
    cpus 8

    input:
    val chr
    path split_chrom
    path genome

    output:
    tuple path("${chr}_rna.exp.bin"), path("${chr}_rna.exp.txt"), path("${chr}_rna.size_factors.bin"), path("${chr}_rna.size_factors.txt"), path("${chr}_snp_counts.tsv")


    script:
    """
    RNA_rasqual_processor.R ${chr}_count.txt ${chr}.vcf.gz $genome $params.eqtl_window ${task.cpus}
    ## rename files
    mv rna.exp.bin ${chr}_rna.exp.bin
    mv rna.exp.txt ${chr}_rna.exp.txt
    mv rna.size_factors.bin ${chr}_rna.size_factors.bin
    mv rna.size_factors.txt ${chr}_rna.size_factors.txt
    mv snp_counts.tsv ${chr}_snp_counts.tsv
    """
}



// QTL mapping with rasqual

process ATAC_RUN_rasqual {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_results_rasqual", mode: 'symlink', overwrite: true
    memory '64 GB'
    cpus 16

    input:
    val chr
    path preproces_data
    path split_chrom
    path covariates

    output:
    path("${chr}_rasqual_lead_snp.txt")


    script:
    """
    echo \$HOSTNAME
    rasqual.R vcf=${chr}.vcf.gz y=${chr}_atac.exp.bin k=${chr}_atac.size_factors.bin x=atac.covs_all_chrom.bin x_txt=atac.covs_all_chrom.txt meta=${chr}_snp_counts.tsv out=${chr}_rasqual_lead_snp.txt cpu=${task.cpus}
    """
}


process RNA_RUN_rasqual {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_results_rasqual", mode: 'symlink', overwrite: true
    memory '64 GB'
    cpus 16

    input:
    val chr
    path preproces_data
    path split_chrom
    path covariates

    output:
    path("${chr}_rasqual_lead_snp.txt")


    script:
    """
    echo \$HOSTNAME
    rasqual.R vcf=${chr}.vcf.gz y=${chr}_rna.exp.bin k=${chr}_rna.size_factors.bin x=rna.covs_all_chrom.bin x_txt=rna.covs_all_chrom.txt meta=${chr}_snp_counts.tsv out=${chr}_rasqual_lead_snp.txt cpu=${task.cpus}
    """
}





// merge rasqual results


process ATAC_MERGE_rasqual {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_results_rasqual_merged", mode: 'symlink', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    val max_chr
    path rasqual_results

    output:
    path("all_chromosome_rasqual_lead_snp.txt")


    script:
    """
    for chr in \$(seq 1 $max_chr)
    do
        cat \${chr}_rasqual_lead_snp.txt >> all_chromosome_rasqual_lead_snp.txt
    done
    """
}


process RNA_MERGE_rasqual {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_results_rasqual_merged", mode: 'symlink', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    val max_chr
    path rasqual_results

    output:
    path("all_chromosome_rasqual_lead_snp.txt")


    script:
    """
    for chr in \$(seq 1 $max_chr)
    do
        cat \${chr}_rasqual_lead_snp.txt >> all_chromosome_rasqual_lead_snp.txt
    done
    """
}






// run rasqual permulation

process ATAC_RUN_rasqual_permutation {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_results_rasqual_permutaion", mode: 'symlink', overwrite: true
    memory '64 GB'
    cpus 16

    input:
    tuple val(chr), val(permute_flag)
    path preproces_data
    path split_chrom
    path covariates


    output:
    tuple val("${permute_flag}"), path("${chr}_permute_${permute_flag}_rasqual_lead_snp.txt")


    script:
    """

    echo \$HOSTNAME
    rasqual_permute.R vcf=${chr}.vcf.gz y=${chr}_atac.exp.bin k=${chr}_atac.size_factors.bin x=atac.covs_all_chrom.bin x_txt=atac.covs_all_chrom.txt meta=${chr}_snp_counts.tsv out=${chr}_permute_${permute_flag}_rasqual_lead_snp.txt cpu=${task.cpus}
    
    """
}



process RNA_RUN_rasqual_permutation {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_results_rasqual_permutaion", mode: 'symlink', overwrite: true
    memory '64 GB'
    cpus 16

    input:
    tuple val(chr), val(permute_flag)
    path preproces_data
    path split_chrom
    path covariates


    output:
    tuple val("${permute_flag}"), path("${chr}_permute_${permute_flag}_rasqual_lead_snp.txt")

    script:
    """

    echo \$HOSTNAME
    rasqual_permute.R vcf=${chr}.vcf.gz y=${chr}_rna.exp.bin k=${chr}_rna.size_factors.bin x=rna.covs_all_chrom.bin x_txt=rna.covs_all_chrom.txt meta=${chr}_snp_counts.tsv out=${chr}_permute_${permute_flag}_rasqual_lead_snp.txt cpu=${task.cpus}


    """
}


// merge rasqual permutation results

process ATAC_MERGE_rasqual_permutation {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_results_rasqual_permutaion_merged", mode: 'symlink', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    val max_chr
    tuple val("permute_flag"), path(rasqual_results)
    

    output:
    path("permute_${permute_flag}_all_chromosome_rasqual_lead_snp.txt")


    script:
    """
    for chr in \$(seq 1 $max_chr)
    do
        cat \${chr}_permute_${permute_flag}_rasqual_lead_snp.txt >> permute_${permute_flag}_all_chromosome_rasqual_lead_snp.txt
    done
    """
}



process RNA_MERGE_rasqual_permutation {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_results_rasqual_permutaion_merged", mode: 'symlink', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    val max_chr
    tuple val("permute_flag"), path(rasqual_results)

    output:
    path("permute_${permute_flag}_all_chromosome_rasqual_lead_snp.txt")


    script:
    """
    for chr in \$(seq 1 $max_chr)
    do
        cat \${chr}_permute_${permute_flag}_rasqual_lead_snp.txt >> permute_${permute_flag}_all_chromosome_rasqual_lead_snp.txt
    done
    """
}


// compute rasqual emperical pvalues


process ATAC_COMPUTE_rasqual_emperical_pvalues {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_results_emperical_pvalues", mode: 'symlink', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    path merged_results
    path permuation_merged_results

    output:
    path("rasqual_emperical_pvalues.txt")


    script:
    """
    rasqual_emperical_pvalues.R rasqual_emperical_pvalues.txt $merged_results $permuation_merged_results
    """
}


process RNA_COMPUTE_rasqual_emperical_pvalues {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_results_emperical_pvalues", mode: 'symlink', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    path merged_results
    path permuation_merged_results

    output:
    path("rasqual_emperical_pvalues.txt")


    script:
    """
    rasqual_emperical_pvalues.R rasqual_emperical_pvalues.txt $merged_results $permuation_merged_results
    """
}





// eigenMT rasqual runnning


process ATAC_eigenMT_process_input {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_eigenMT_process_input", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    val chr
    path split_chrom

    output:
    tuple path("${chr}_genotype.txt"), path("${chr}_genotype_position.txt"), path("${chr}_phenotype_position.txt")

    script:
    """
    MatrixQTL_genotype_converter.py --vcf ${chr}.vcf.gz --out_genotype ${chr}_genotype.txt --out_genotype_position ${chr}_genotype_position.txt

    MatrixQTL_ATAC_phenotype_converter.py --count ${chr}_count.txt --out_phenotype ${chr}_phenotype.txt --out_phenotype_position ${chr}_phenotype_position.txt

    """
}



process RNA_eigenMT_process_input {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_eigenMT_process_input", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    val chr
    path split_chrom

    output:
    tuple path("${chr}_genotype.txt"), path("${chr}_genotype_position.txt"), path("${chr}_phenotype_position.txt")

    script:
    """
    MatrixQTL_genotype_converter.py --vcf ${chr}.vcf.gz --out_genotype ${chr}_genotype.txt --out_genotype_position ${chr}_genotype_position.txt

    MatrixQTL_RNA_phenotype_converter.py --count ${chr}_count.txt --out_phenotype ${chr}_phenotype.txt --out_phenotype_position ${chr}_phenotype_position.txt

    """
}



// run rasqual eigenMT

process ATAC_RUN_rasqual_eigenMT {
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


process RNA_RUN_rasqual_eigenMT {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_results_rasqual_eigenMT", mode: 'symlink', overwrite: true
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
    rasqual_eigenMT.R vcf=${chr}.vcf.gz y=${chr}_rna.exp.bin k=${chr}_rna.size_factors.bin x=rna.covs_all_chrom.bin x_txt=rna.covs_all_chrom.txt meta=${chr}_snp_counts.tsv out=${chr}_rasqual_all_snp.txt cpu=${task.cpus}
    """
}



// permuatation sample


process ATAC_RUN_rasqual_eigenMT_permute {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_results_rasqual_eigenMT_permute", mode: 'symlink', overwrite: true
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
    rasqual_permute_eigenMT.R vcf=${chr}.vcf.gz y=${chr}_atac.exp.bin k=${chr}_atac.size_factors.bin x=atac.covs_all_chrom.bin x_txt=atac.covs_all_chrom.txt meta=${chr}_snp_counts.tsv out=${chr}_rasqual_all_snp.txt cpu=${task.cpus}
    """
}


process RNA_RUN_rasqual_eigenMT_permute {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_results_rasqual_eigenMT_permute", mode: 'symlink', overwrite: true
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
    rasqual_permute_eigenMT.R vcf=${chr}.vcf.gz y=${chr}_rna.exp.bin k=${chr}_rna.size_factors.bin x=rna.covs_all_chrom.bin x_txt=rna.covs_all_chrom.txt meta=${chr}_snp_counts.tsv out=${chr}_rasqual_all_snp.txt cpu=${task.cpus}
    """
}


// convert to eigenMT input


process ATAC_rasqual_TO_eigenMT {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_results_rasqual_eigenMT_processed", mode: 'symlink', overwrite: true
    memory '64 GB'
    cpus 16

    input:
    val chr
    path rasqual_res

    output:
    path("${chr}_formated_EigenMT.txt")


    script:
    """
    echo \$HOSTNAME
    rasqualToEigenMT.py --rasqualOut ${chr}_rasqual_all_snp.txt > ${chr}_formated_EigenMT.txt
    """
}


process RNA_rasqual_TO_eigenMT {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_results_rasqual_eigenMT_processed", mode: 'symlink', overwrite: true
    memory '64 GB'
    cpus 16

    input:
    val chr
    path rasqual_res

    output:
    path("${chr}_formated_EigenMT.txt")


    script:
    """
    echo \$HOSTNAME
    rasqualToEigenMT.py --rasqualOut ${chr}_rasqual_all_snp.txt > ${chr}_formated_EigenMT.txt
    """
}


// permuate


process ATAC_rasqual_TO_eigenMT_permute {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_results_rasqual_eigenMT_processed_permute", mode: 'symlink', overwrite: true
    memory '64 GB'
    cpus 16

    input:
    val chr
    path rasqual_res

    output:
    path("${chr}_formated_EigenMT.txt")


    script:
    """
    echo \$HOSTNAME
    rasqualToEigenMT.py --rasqualOut ${chr}_rasqual_all_snp.txt > ${chr}_formated_EigenMT.txt
    """
}


process RNA_rasqual_TO_eigenMT_permute {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_results_rasqual_eigenMT_processed_permute", mode: 'symlink', overwrite: true
    memory '64 GB'
    cpus 16

    input:
    val chr
    path rasqual_res

    output:
    path("${chr}_formated_EigenMT.txt")


    script:
    """
    echo \$HOSTNAME
    rasqualToEigenMT.py --rasqualOut ${chr}_rasqual_all_snp.txt > ${chr}_formated_EigenMT.txt
    """
}



// FDR correction with eigenMT

process ATAC_eigenMT {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_eigenMT_results", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    val chr
    path rasqual_eigenMT
    path all_input

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
	    --OUT ${chr}_eigenMT_results.txt

    """
}



process RNA_eigenMT {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_eigenMT_results", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    val chr
    path rasqual_eigenMT
    path all_input

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
	    --OUT ${chr}_eigenMT_results.txt

    """
}


// permuate


process ATAC_eigenMT_permute {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_eigenMT_results_permute", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    val chr
    path rasqual_eigenMT
    path all_input

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
	    --OUT ${chr}_eigenMT_results.txt

    """
}



process RNA_eigenMT_permute {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_eigenMT_results_permute", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    val chr
    path rasqual_eigenMT
    path all_input

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
	    --OUT ${chr}_eigenMT_results.txt

    """
}


// merge rasqual results


process ATAC_MERGE_eigenMT {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_eigenMT_results_merged", mode: 'copy', overwrite: true
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


process RNA_MERGE_eigenMT {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_eigenMT_results_merged", mode: 'copy', overwrite: true
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


// permuate


process ATAC_MERGE_eigenMT_permute {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/ATAC_eigenMT_results_merged_permute", mode: 'copy', overwrite: true
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


process RNA_MERGE_eigenMT_permute {
    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/RNA_eigenMT_results_merged_permute", mode: 'copy', overwrite: true
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

// LEAVE ONE OUT ANALYSES

process LOO_meta_csv {

    container 'ndatth/rasqual:v0.0.0'
    publishDir "${params.outdir}/loo_meta", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    path meta

    output:
    path "tem/*"

    script:
    """
    loo_meta.R $meta
    """
}


// external LD reference for eigenMT

