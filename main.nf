#!/usr/bin/env nextflow
// Structure based on https://github.com/epi2me-labs/wf-template/blob/master/main.nf

nextflow.enable.dsl = 2

process merge_reads {
    cpus 1
    memory '1 G'

    input: 
        val sample_name
        file basecall_directory
    output:
        file "merged_reads/${sample_name}.fastq"
    shell:
    """
    mkdir -p merged_reads
    fastcat -x ${basecall_directory} > merged_reads/${sample_name}.fastq
    """
}

process basecall_reads {
    label 'GPU'

    memory '16 G'
    cpus params.threads
    time '3d'
    queue 'gpu.q'

    input:
        tuple val(sample_name), file(run_directory)
    output:
        val sample_name, emit: sample_name
        path run_directory, emit: run_directory
        path "basecalled", emit: basecalled_directory
        path "basecalled/sequencing_summary.txt", emit: sequencing_summary
    shell:
    """
    $params.guppy -r --num_callers $task.cpus --gpu_runners_per_device 4 --chunks_per_runner 512 -c dna_r9.4.1_450bps_hac.cfg -i $run_directory -x "cuda:0 cuda:1" -s basecalled
    """
}

process align_reads {
    cpus params.threads
    memory '32 G'

    publishDir "${sample_name}_results", mode: 'copy'
    input:
        val sample_name
        file "${sample_name}.fastq"
    output:
        path "${sample_name}.sorted.bam", emit: bam
        path "${sample_name}.sorted.bam.bai", emit: bai
    shell:
    """
    minimap2 -ax map-ont -t $task.cpus $params.reference ${sample_name}.fastq | samtools sort -T tmp > ${sample_name}.sorted.bam
    samtools index ${sample_name}.sorted.bam
    """
}

process nanopolish_index {
    time '1d'
    memory '24 G'

    input:
        val sample_name
        file "${sample_name}.fastq"
        file "sequencing_summary.txt"
        file "run_directory"
    output:
        tuple file("${sample_name}.fastq.index"), file("${sample_name}.fastq.index.gzi"), file("${sample_name}.fastq.index.fai"), file("${sample_name}.fastq.index.readdb")
    shell:
    """
    $params.nanopolish index -s sequencing_summary.txt -d run_directory ${sample_name}.fastq
    """

}

process nanopolish_call_methylation {
    cpus params.threads
    memory '32 G'
    time '3d'
    
    publishDir "${sample_name}_results", mode: 'copy'

    input:
        each chr
        val sample_name
        file "${sample_name}.fastq"
        file "${sample_name}.sorted.bam"
        file "${sample_name}.sorted.bam.bai"
        file "run_directory"
        tuple file("${sample_name}.fastq.index"), file("${sample_name}.fastq.index.gzi"), file("${sample_name}.fastq.index.fai"), file("${sample_name}.fastq.index.readdb")
    output:
        tuple( val(sample_name), path("${sample_name}.modifications.${chr}.sorted.bam"), emit: modbam)
    shell:
    """
    $params.nanopolish call-methylation -b ${sample_name}.sorted.bam -r ${sample_name}.fastq -g $params.reference --modbam-output ${sample_name}.modifications.${chr}.sorted.bam -t $task.cpus -w $chr
    """
}

process merge_bam {
    cpus params.threads
    memory '32 G'
    time '1d'

    publishDir "${sample_name}_results", mode: 'copy'

    input:
        tuple( val(sample_name), val(bams) )
    output:
        val sample_name, emit: sample_name
        path "${sample_name}.modifications.sorted.bam", emit: bam
    shell:
    """
    samtools merge -o ${sample_name}.modifications.sorted.bam ${bams.join(' ')}
    """
}

process calculate_read_modification_frequency {
    cpus params.threads
    memory '32 G'
    time '1d'
    
    publishDir "${sample_name}_results", mode: 'copy'

    input:
        val sample_name
        file "${sample_name}.modifications.sorted.bam"
    output:
        val sample_name, emit: sample_name
        path "${sample_name}.read_modifications.tsv", emit: readmod
    shell:
    """
    $params.mbtools read-frequency ${sample_name}.modifications.sorted.bam > ${sample_name}.read_modifications.tsv
    """
}

process calculate_fragmentation {
    cpus params.threads
    memory '32 G'
    time '1d'
    
    publishDir "${sample_name}_results", mode: 'copy'

    input:
        val sample_name
        file "${sample_name}.read_modifications.tsv"
    output:
        file "${sample_name}.fragmentation.ratios.bed"
    shell:
    """
    python $params.fragmentation -o ${sample_name}.fragmentation.ratios.bed -s 100 151 -l 151 221 -b 500000 ${sample_name}.read_modifications.tsv
    """
}

workflow pipeline {
    take:
        input
    main:
        
        if(params.rebasecall) {
            basecall_output = basecall_reads(input)
            //basecalled_directory = basecall_output[0]
            //sequencing_summary = basecall_output[1]
        } else {
            basecalled_directory = run_directory
            sequencing_summary = file("${params.run_directory}/**_sequencing_summary.txt", type: 'file')
        }
        
        merge_fastq_output = merge_reads(basecall_output.sample_name, basecall_output.basecalled_directory)
        index_output = nanopolish_index(basecall_output.sample_name, merge_fastq_output, basecall_output.sequencing_summary, basecall_output.run_directory)
        align_output = align_reads(basecall_output.sample_name, merge_fastq_output)

        chrs = Channel.from(1 .. 22).map { "chr" + it }
        chrs_sex = Channel.of("chrX", "chrY")
        chrs = chrs.concat(chrs_sex)

        modbam_output = nanopolish_call_methylation(chrs, basecall_output.sample_name, merge_fastq_output, align_output.bam, align_output.bai, basecall_output.run_directory, index_output)
        merge_bam_output = merge_bam(modbam_output.modbam.groupTuple())
        read_frequency_output = calculate_read_modification_frequency(merge_bam_output.sample_name, merge_bam_output.bam)
        fragmentation_output = calculate_fragmentation(read_frequency_output.sample_name, read_frequency_output.readmod)
}

workflow {

    runs = Channel.fromPath( 'data/*', type: 'dir')

    // determine sample name for each input
    input = runs.map { [it.simpleName, it] }
    input.view()

    pipeline(input)
}
