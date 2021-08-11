#!/usr/bin/env nextflow
// Structure based on https://github.com/epi2me-labs/wf-template/blob/master/main.nf

nextflow.enable.dsl = 2

process merge_reads {
    input: 
        file "run_directory"
    output:
        file "merged_reads/${params.output}.fastq"
    shell:
    """
    mkdir -p merged_reads
    fastcat -x run_directory > merged_reads/${params.output}.fastq
    """
}

process basecall_reads {
    input:
        file "run_directory"
    output:
        file "basecalled"
        file "basecalled/sequencing_summary.txt"
    shell:
    """
    $params.guppy -r --num_callers 8 --gpu_runners_per_device 4 --chunks_per_runner 512 -c dna_r9.4.1_450bps_hac.cfg -i run_directory -x "cuda:0 cuda:1" -s basecalled
    """
}

process align_reads {
    cpus params.threads
    publishDir "${params.output}_results", mode: 'copy'
    input:
        file "${params.output}.fastq"
    output:
        file "${params.output}.sorted.bam"
    shell:
    """
    minimap2 -ax map-ont -t $task.cpus $params.reference ${params.output}.fastq | samtools sort -T tmp > ${params.output}.sorted.bam
    """
}

process index_bam {
    input:
        file "${params.output}.sorted.bam"
    output:
        file "${params.output}.sorted.bam.bai"
    shell:
    """
    samtools index ${params.output}.sorted.bam
    """
}

process nanopolish_index {
    cpus params.threads
    input:
        file "${params.output}.fastq"
        file "sequencing_summary.txt"
        file "run_directory"
    output:
        file "${params.output}.fastq.index"
        file "${params.output}.fastq.index.gzi"
        file "${params.output}.fastq.index.fai"
        file "${params.output}.fastq.index.readdb"
    shell:
    """
    $params.nanopolish index -s sequencing_summary.txt -d run_directory ${params.output}.fastq
    """

}

process nanopolish_call_methylation {
    cpus params.threads
    publishDir "${params.output}_results", mode: 'copy'

    input:
        file "${params.output}.fastq"
        file "${params.output}.sorted.bam"
        file "${params.output}.sorted.bam.bai"
        file "run_directory"
        file "${params.output}.fastq.index"
        file "${params.output}.fastq.index.gzi"
        file "${params.output}.fastq.index.fai"
        file "${params.output}.fastq.index.readdb"
    output:
        file "${params.output}.modifications.sorted.bam"
    shell:
    """
    $params.nanopolish call-methylation -b ${params.output}.sorted.bam -r ${params.output}.fastq -g $params.reference --modbam-output ${params.output}.modifications.sorted.bam -t $task.cpus
    """
}

process run_nanoplot_bam {
    publishDir "${params.output}_results", mode: 'copy'
    input:
        file "${params.output}.sorted.bam"
    output:
        file "nanoplot_bam"
    shell:
    """
    NanoPlot --bam ${params.output}.sorted.bam -o nanoplot_bam
    """
}

workflow pipeline {
    take:
        run_directory
    main:

        if(params.rebasecall) {
            basecall_output = basecall_reads(run_directory)
            basecalled_directory = basecall_output[0]
            sequencing_summary = basecall_output[1]
        } else {
            basecalled_directory = run_directory
            sequencing_summary = file("${params.run_directory}/**_sequencing_summary.txt", type: 'file')
        }

        merged_fastq = merge_reads(basecalled_directory)
        np_index = nanopolish_index(merged_fastq, sequencing_summary, run_directory)

        bam = align_reads(merged_fastq)
        bai = index_bam(bam)
        modbam = nanopolish_call_methylation(merged_fastq, bam, bai, run_directory, np_index)
        nanoplot_bam = run_nanoplot_bam(bam)
}

workflow {

    if(!params.run_directory) {
        log.info """
Error, no run directory provided
        """
        exit 1
    }

    // if the user does not provide a sample name, create one based on the input directory
    if(!params.sample_name) {
        rd = file(params.run_directory)
        params.output = rd.simpleName
    } else {
        params.output = params.sample_name
    }

    reads = Channel.fromPath(params.run_directory, type: 'dir', checkIfExists: true)
    results = pipeline(reads)
    // output(results)
}
