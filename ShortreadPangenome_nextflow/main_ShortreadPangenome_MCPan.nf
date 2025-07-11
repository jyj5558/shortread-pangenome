#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//where are the inputs located
//path to CSV file
params.csv = "" // e.g., "/scratch/negishi/jeon96/shortread_pan/params.csv"

//option to change from default directory to path of choice
params.savepath = "" // e.g., "/scratch/negishi/jeon96/shortread_pan/"

//user email information
params.email = "" //your email address; e.g., "jeon96@purdue.edu"

/*
========================================================================================
========================================================================================
//A TEMPLATE FOR A NEXTFLOW WORKFLOW FOR USE ON PURDUE RCAC CLUSTER
//DESCRIPTION: WHAT DOES THIS WORKFLOW DO?
//currently it reads an input CSV, processes genomic data, maps sequences, distinguishes mappable and unmappable reads,
//then finally output short-read-based linear pangenome and finally short-read-based Minigraph-Cactus graph pangenome
// DO NOT forget to store MASURCA config example file in your home directory in advance
// If any step fails, make sure whether the last result files are corrupted or not
// !! PLEASE LOOK AT EVERY SCRIPT THOROUGHLY BEFORE RUNNING THE NEXTFLOW PIPELINE AND MODIFY THEM ACCORDINGLY !!"
========================================================================================
========================================================================================
*/

/*
========================================================================================
========================================================================================
//PROCESSES FOR THE PIPELINE
//PIPELINE EXECUTION CONTROL IN WORKFLOW AT THE END
========================================================================================
    SETUP A NEEDED DIRECTORY (e.g., /scratch/negishi/jeon96/shortread_pangenome/) 
    AND STORE NEEDED SCRIPTS IN IT
========================================================================================
*/

process step1{
    tag "$step1"
    clusterOptions '--job-name=Step1 -n 32 -N 1 -t 5-00:00:00 -A fnrdewoody --mail-user $params.email --mail-type END,FAIL' 
    //errorStrategy 'ignore' #don't want to ignore errors for this process

    input:
    tuple val(sra), val(db), val(ref), val(linpan), val(contig), val(mcpan)

    output:
    tuple val(sra), val(db), val(ref), val(linpan), val(contig), val(mcpan)

    script:
    """
    cd /scratch/negishi/allen715/shortread_pan/
    bash Step1_SequencePreprocessing_nf.sh ${sra} ${db} ${params.N}
    """
}

process step2_1{
    tag "$step2_1"
    clusterOptions '--job-name=Step2.1 -n 64 -N 1 -t 1-00:00:00 -A highmem --mail-user $params.email --mail-type END,FAIL'
    //errorStrategy 'ignore' #don't want to ignore errors for this process

    input:
    tuple val(sra), val(db), val(ref), val(linpan), val(contig), val(mcpan)

    output:
    tuple val(sra), val(db), val(ref), val(linpan), val(contig)

    script:
    """
    bash Step2.1_FromAssemblingToLinPan_nf.sh ${sra} ${ref} ${linpan} ${contig} ${params.GENOME} ${params.N} ${params.APP} ${params.MASURCA} ${params.GENOME_SIZE} ${params.PREFIX}
    """
}

process step2_2{
    tag "$step2_2"
    clusterOptions '--job-name=Step2.2 -n 32 -N 1 -t 1-00:00:00 -A highmem --mail-user $params.email --mail-type END,FAIL'
    errorStrategy 'retry'
    maxRetries 5

    input:
    tuple val(sra), val(db), val(ref), val(linpan), val(contig)

    output:
    tuple val(sra), val(db), val(ref), val(linpan), val(contig)

    script:
    """
    bash Step2.2_FromAssemblingToLinPan_nf.sh ${sra} ${ref} ${linpan} ${contig} ${params.GENOME} ${params.N} ${params.APP} ${params.MASURCA} ${params.GENOME_SIZE} ${params.PREFIX} ${params.MASURCA_CONFIG_MAPPED} ${params.MASURCA_CONFIG_UNMAPPED}
    """
}

process step2_3{
    tag "$step2_3"
    clusterOptions '--job-name=Step2.3 -n 64 -N 1 -t 1-00:00:00 -A highmem --mail-user $params.email --mail-type END,FAIL'
    //errorStrategy 'ignore' #don't want to ignore errors for this process

    input:
    tuple val(sra), val(db), val(ref), val(linpan), val(contig)

    output:
    tuple val(sra), val(db), val(ref), val(linpan), val(contig)

    script:
    """
    bash Step2.3_FromAssemblingToLinPan_nf.sh ${sra} ${ref} ${linpan} ${contig} ${params.GENOME} ${params.N} ${params.APP} ${params.MASURCA} ${params.GENOME_SIZE} ${params.PREFIX}
    """
}

process step2_4{
    tag "$step2_4"
    clusterOptions '--job-name=Step2.4 -n 64 -N 1 -t 5-00:00:00 -A fnrdewoody --mail-user $params.email --mail-type END,FAIL'
    errorStrategy 'retry'
    maxRetries 5

    input:
    tuple val(sra), val(db), val(ref), val(linpan), val(contig)

    output:
    tuple val(sra), val(db), val(ref), val(linpan), val(contig)

    script:
    """
    bash Step2.4_FromAssemblingToLinPan_nf.sh ${sra} ${ref} ${linpan} ${contig} ${params.GENOME} ${params.N} ${params.APP} ${params.MASURCA} ${params.GENOME_SIZE} ${params.PREFIX} ${params.MASURCA_CONFIG_UNMAPPED2}
    """
}

process step2_5{
    tag "$step2_5"
    clusterOptions '--job-name=Step2.5 -n 64 -N 1 -t 5-00:00:00 -A fnrdewoody --mail-user $params.email --mail-type END,FAIL'
    //errorStrategy 'ignore' #don't want to ignore errors for this process

    input:
    tuple val(sra), val(db), val(ref), val(linpan), val(contig)

    output:
    tuple val(sra), val(db), val(ref), val(linpan), val(contig)

    script:
    """
    bash Step2.5_FromAssemblingToLinPan_nf.sh ${sra} ${ref} ${linpan} ${contig} ${params.GENOME} ${params.N} ${params.APP} ${params.MASURCA} ${params.GENOME_SIZE} ${params.PREFIX} ${params.MASURCA_CONFIG_UNMAPPED2}
    """
}

process step2_6{
    tag "$step2_6"
    clusterOptions '--job-name=Step2.6 -n 64 -N 1 -t 1-00:00:00 -A highmem --mail-user $params.email --mail-type END,FAIL'
    //errorStrategy 'ignore' #don't want to ignore errors for this process

    input:
    tuple val(sra), val(db), val(ref), val(linpan), val(contig)

    output:
    tuple val(sra), val(db), val(ref), val(linpan), val(contig)

    script:
    """
    bash Step2.6_FromAssemblingToLinPan_nf.sh ${sra} ${ref} ${linpan} ${contig} ${params.GENOME} ${params.N} ${params.APP} ${params.MASURCA} ${params.GENOME_SIZE} ${params.PREFIX} ${params.MASURCA_CONFIG_UNMAPPED2}
    """
}

process step3{
    tag "$step3"
    clusterOptions '--job-name=Step3 -n 32 -N 1 -t 5-00:00:00 -A fnrdewoody --mail-user $params.email --mail-type END,FAIL'
    //errorStrategy 'ignore' #don't want to ignore errors for this process

    input:
    tuple val(sra), val(db), val(ref), val(linpan), val(contig), val(mcpan)

    output:
    tuple val(sra), val(db), val(ref), val(linpan), val(contig), val(mcpan)

    script:
    """
    bash Step3.A_FromRetrievingToMCPan_nf.sh ${sra} ${ref} ${linpan} ${contig} ${mcpan} ${params.N} ${params.PREFIX} ${params.GENOME} ${params.RefInd} ${params.AltInd} ${params.AltHap} ${params.APP} ${params.CACTUS}
    """
}


/*
========================================================================================
========================================================================================
PIPELINE EXECUTION CONTROL
*/
workflow{
    //input channel for a CSV that could contain parameters for Theta workflow
    print "Welcome to this CSV workflow to construct a short-read-based pangenome"
    print 'You will be working on the following files: '+params.csv
    Channel.fromPath(params.csv) \
        | splitCsv(header:true) \
        | map { row-> [row.sra, row.db, row.ref, row.linpan, row.contig, row.mcpan] } \
        | view() \
        | step1 \
        | step2_1 \
        | step2_2 \
        | step2_3 \
        | step2_4 \
        | step2_5 \
        | step2_6 \
        | step3 
}
