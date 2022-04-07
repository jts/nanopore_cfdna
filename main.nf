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
    $params.guppy -r --num_callers $task.cpus --gpu_runners_per_device 4 --chunks_per_runner 512 -c res_dna_r9.4.1_e8.1_hac_v033.cfg -d $params.reriomodels -i $run_directory -x "cuda:0 cuda:1" -s basecalled
    """
}
process align_reads {
    cpus params.threads
    memory '32 G'

    publishDir "${sample_name}_results", mode: 'symlink'
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
process bam_stats {
    memory '32 G'
    time '1d'

    publishDir "${sample_name}_results", mode: 'symlink'

    input:
        val sample_name
        file "${sample_name}.sorted.bam"
        file "${sample_name}.sorted.bam.bai"
    output:
        path "${sample_name}.bamstats.tsv", emit: stats
    shell:
    """
    stats_from_bam ${sample_name}.sorted.bam -o ${sample_name}.bamstats.tsv
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
    time '7d'
    
    input:
        each chr
        val sample_name
        file "${sample_name}.fastq"
        file "${sample_name}.sorted.bam"
        file "${sample_name}.sorted.bam.bai"
        file "run_directory"
        tuple file("${sample_name}.fastq.index"),
            file("${sample_name}.fastq.index.gzi"),
            file("${sample_name}.fastq.index.fai"),
            file("${sample_name}.fastq.index.readdb")
    output:
        tuple( val(sample_name), path("${sample_name}.${chr}.modifications.sorted.bam"), emit: modbam)
    shell:
    """
    $params.nanopolish call-methylation -b ${sample_name}*.bam -r ${sample_name}.fastq -g $params.reference -w ${chr} --modbam-output ${sample_name}.${chr}.modifications.sorted.bam -t $task.cpus
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
        path "${sample_name}.merged.modifications.sorted.bam", emit: modbam
    shell:
    """
    samtools merge -o ${sample_name}.merged.modifications.sorted.bam ${bams.join(' ')} --threads $params.threads
    """
}
process read_modification_frequency {
    cpus params.threads
    memory '32 G'
    time '1d'
    
    publishDir "${sample_name}_results", mode: 'symlink'

    input:
        val sample_name
        file "${sample_name}.merged.modifications.sorted.bam"
    output:
        val sample_name , emit: sample_name
        path "${sample_name}.read_modifications.tsv", emit: modtsv
    shell:
    """
        $params.mbtools read-frequency  ${sample_name}.merged.modifications.sorted.bam  > ${sample_name}.read_modifications.tsv 
    """
}
process clean_bams {
    input:
        val sample_name
    shell:
    """
    rm \$(find -type f | grep modifications.sorted.bam | grep chr | grep ${sample_name})
    """
}
process region_modification_frequency {
    memory '32 G'
    time '7d'
    publishDir "${sample_name}_results", mode: 'symlink'
    input:
        val sample_name
        file "${sample_name}.merged.modifications.sorted.bam"
    output:
        val sample_name, emit: sample_name
        path "${sample_name}.region_modifications.tsv", emit: modtsv
    shell:
        """
        $params.mbtools region-frequency  ${sample_name}.merged.modifications.sorted.bam -r ${projectDir}/atlases/cfDNAme.bed  > ${sample_name}.region_modifications.tsv 
        """
}
process fragmentation {
    cpus params.threads
    memory '32 G'
    time '1d'

    publishDir "${sample_name}_results", mode: 'symlink'
    
    input:
        val sample_name
        file "${sample_name}.read_modifications.tsv"
    output:
        file "${sample_name}.fragmentation.ratios.tsv"
    shell:
    """
    ${projectDir}/scripts/fragmentation.py -o ${sample_name}.fragmentation.ratios.tsv -s 100 151 -l 151 221 -b 5000000 ${sample_name}.read_modifications.tsv
    """
}
process plot_fragmentation {
    publishDir "plots", mode: 'copy'

    input:
        val fragmentomes
    output:
        path "fragmentome.pdf"
    shell:
    """
    ${projectDir}/scripts/plot_fragmentome.r ${fragmentomes.join(" ")}
    """
}
process deconvolve {
    cpus params.threads
    memory '32 G'
    time '1d'

    publishDir "plots", mode: 'copy'
    
    input:
        val modbams
    output:
        path "deconv_output.tsv"
    shell:
    """
    ${projectDir}/scripts/nanomix.py --atlas ${projectDir}/atlases/cfDNAme.agg.tsv ${modbams.join(" ")} > "deconv_output.tsv"
    """
}

process plot_deconvolution {

    publishDir "plots", mode: 'copy'

    input:
        val deconv_output
    output:
        file "deconv_output.png"
    shell:
    """
    ${projectDir}/scripts/plot_deconv.py -name deconv_output.png $deconv_output
    """
}
workflow pipeline {
    take:
        input
    main:
        basecall_output = basecall_reads(input)
        merge_fastq_output = merge_reads(
            basecall_output.sample_name,
            basecall_output.basecalled_directory
        )
        index_output = nanopolish_index(
            basecall_output.sample_name,
            merge_fastq_output,
            basecall_output.sequencing_summary,
            basecall_output.run_directory
        )
        align_output = align_reads(
            basecall_output.sample_name,
            merge_fastq_output
        )
        bam_stats(
            basecall_output.sample_name,
            align_output.bam,
            align_output.bai
        )
        chrs = Channel.from(1 .. 22).map { "chr" + it }
        chrs_sex = Channel.of("chrX", "chrY")
        chrs = chrs.concat(chrs_sex)
        modbams = nanopolish_call_methylation(chrs,
                        basecall_output.sample_name,
                        merge_fastq_output,
                        align_output.bam,
                        align_output.bai,
                        basecall_output.run_directory,
                        index_output
                        )
        modbam_output = merge_bam(modbams.modbam.groupTuple())
        read_frequency_output = read_modification_frequency(
            modbam_output.sample_name,
            modbam_output.modbam
        )
        fragmentation_output = fragmentation(
            read_frequency_output.sample_name,
            read_frequency_output.modtsv
        )
        plot_fragmentation(fragmentation_output.collect())

        region_frequency_output = region_modification_frequency(
            modbam_output.sample_name,
            modbam_output.modbam,
            )
        deconv_output = deconvolve(region_frequency_output.modtsv.collect())
        plot_deconvolution(deconv_output)
}

workflow {
    runs = Channel.fromPath( 'data/*', type: 'dir')

    // determine sample name for each input
    input = runs.map { [it.simpleName, it] }
    pipeline(input)
}
