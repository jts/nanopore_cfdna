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

process get_nanopolish_input {
    input:
        tuple val(sample_name), file("run_directory")
    output:
        tuple val(sample_name),
            path("**/fastq/${sample_name}.fastq.gz"),
            path("**/mapped/${sample_name}*.bam"),
            path("**/mapped/${sample_name}*.bam.bai"),
            path("**/fastq/${sample_name}.fastq.gz.index"),
            path("**/fastq/${sample_name}.fastq.gz.index.gzi"),
            path("**/fastq/${sample_name}.fastq.gz.index.fai"),
            path("**/fastq/${sample_name}.fastq.gz.index.readdb")
    shell:
    """
    ls $run_directory/fastq/${sample_name}*
    ls $run_directory/mapped/${sample_name}*.bam* 
    """
}
process get_mergebam_output {
    publishDir "${sample_name}_results", mode: 'symlink'

    input:
        tuple val(sample_name), file("run_directory")
    output:
        val sample_name, emit: sample_name
        path "${sample_name}.merged.modifications.sorted.bam", emit: modbam
    shell:
    """
    ln -s \$(find $workDir -type f | grep ${sample_name}.merged.modifications.sorted.bam | tail -n1) ${sample_name}.merged.modifications.sorted.bam
    """
}
process get_nanopolish_output {

    input:
        each chr
        tuple val(sample_name), file("run_directory")
    output:
        tuple( val(sample_name), path("${sample_name}.${chr}.modifications.sorted.bam"), emit: modbam)
    shell:
    """
    ln -s \$(find $launchDir -type f | grep ${sample_name}.${chr}.modifications.sorted.bam) ${sample_name}.${chr}.modifications.sorted.bam
    """
}
process nanopolish_call_methylation {
    cpus params.threads
    memory '32 G'
    time '7d'
    
    input:
        each chr
        tuple val(sample_name),
            path("${sample_name}.fastq.gz"),
            path("${sample_name}*.bam"),
            path("${sample_name}*.bam.bai"),
            path("${sample_name}.fastq.gz.index"),
            path("${sample_name}.fastq.gz.index.gzi"),
            path("${sample_name}.fastq.gz.index.fai"),
            path("${sample_name}.fastq.gz.index.readdb")
    output:
        tuple( val(sample_name), path("${sample_name}.${chr}.modifications.sorted.bam"), emit: modbam)
    shell:
    """
    $params.nanopolish call-methylation -b ${sample_name}*.bam -r ${sample_name}.fastq.gz -g $params.reference -w ${chr} --modbam-output ${sample_name}.${chr}.modifications.sorted.bam -t $task.cpus
    """
}
process merge_bam {
    cpus params.threads
    memory '32 G'
    time '1d'
    
    publishDir "${sample_name}_results", mode: 'symlink'

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
process clean_bams {
    input:
        val sample_name
    shell:
    """
    rm \$(find -type f | grep modifications.sorted.bam | grep chr | grep ${sample_name})
    """
}
process downsample_reference_frequency {
    memory '32 G'
    time '2d'

    input:
        val sample_name
        file "${sample_name}.modifications.sorted.bam"
        val coverage
        each target_coverage
    output:
        tuple (val(sample_name), path("${sample_name}.C*C.reference_modifications.tsv"), emit: refmod)
    shell:
    if (target_coverage < coverage as float)
        """
        $params.mbtools reference-frequency ${sample_name}.modifications.sorted.bam > ${sample_name}.C${target_coverage}C.reference_modifications.tsv -D ${target_coverage} -C ${coverage} 
        """
    else
        """
        f=\$(find $workDir -type f -name "${sample_name}.C${coverage as float}C.reference_modifications.tsv")
        if [ -z \$f ]
        then
            $params.mbtools reference-frequency ${sample_name}.modifications.sorted.bam > ${sample_name}.C${coverage as float}C.reference_modifications.tsv
        else
            touch ${sample_name}.C${target_coverage}C.reference_modifications.tsv
        fi
        """
}
process get_ref_freq_output {

    input:
        val sample_name
        val coverage
        each target_coverage
    output:
        tuple (val(sample_name), path("${sample_name}.C*C.reference_modifications.tsv"), emit: refmod)
    shell:
    if (target_coverage < coverage as float)
        """
        ln -s \$(ls \$(find $workDir -type f -name "${sample_name}.C${target_coverage}C.reference_modifications.tsv") -t | head -n1) ${sample_name}.C${target_coverage}C.reference_modifications.tsv
        """
    else
        """
        ln -s \$(ls \$(find $workDir -type f -name "${sample_name}.C${coverage as float}C.reference_modifications.tsv") -t | head -n1) ${sample_name}.C${coverage as float}C.reference_modifications.tsv
        """
}
process calculate_coverage {
    memory '8 G'
    time '1d'

    input:
        val sample_name
        file "${sample_name}.sorted.bam"
    output:
        val sample_name, emit: sample_name
        path "${sample_name}.sorted.bam", emit: modbam
        stdout emit: sample_coverage
    shell:
    """
    samtools coverage ${sample_name}.sorted.bam | head -n25 | tail -n24 | awk -F'\\t' '{sum+=\$7} END {print(sum/NR) }'
    """
}
process get_cpgs {
    memory '32 G'
    time '1d'
    publishDir "${sample_name}_results", mode: 'symlink'

    input:
        tuple(val(sample_name), val(refmods))

    output:
        val sample_name, emit: sample_name
        path "${sample_name}.cpgs.csv", emit: csv

    shell:
    """
    ${projectDir}/scripts/get_cpgs.py ${refmods.join(' ')} -o ${sample_name}.cpgs.csv
    """
}
process deconvolve {
    cpus params.threads
    memory '32 G'
    time '1d'

    publishDir "${sample_name}_results", mode: 'symlink'
    
    input:
        val sample_name 
        file "${sample_name}.cpgs.csv"
    output:
        tuple( path ("${sample_name}_llse.cpgs_deconv_output.csv"), 
               path("${sample_name}_nnls.cpgs_deconv_output.csv")) 
        /*file "${sample_name}.cpgs_deconv_plot.png*/
    shell:
    """
    ${params.deconvolve} --input ${sample_name}.cpgs.csv > "${sample_name}_llse.cpgs_deconv_output.csv"
    ${params.deconvolve} --input ${sample_name}.cpgs.csv --model nnls > "${sample_name}_nnls.cpgs_deconv_output.csv"
    """
}
process plot_accuracy {

    publishDir "plots", mode: 'symlink'

    input:
        val deconv_output
    output:
        file "*coverageVaccuracy.png"
        file "*coverageVaccuracy.tsv"
    shell:
    """
    ${projectDir}/scripts/plot_accuracy.r ${deconv_output.join(' ')}
    """
}
workflow pipeline {
    take:
        input
    main:
        chrs = Channel.from(1 .. 22).map { "chr" + it }
        chrs_sex = Channel.of("chrX", "chrY")
        chrs = chrs.concat(chrs_sex)
        cvrg = Channel.from([0.1, 0.5, 1, 2, 5, 10, 50])
        if (params.call_nanopolish) {
            nanopolish_input = get_nanopolish_input(input)
            modbams = nanopolish_call_methylation(chrs,
                            nanopolish_input
                            )
            modbam_output = merge_bam(modbams.modbam.groupTuple())
	    }
        else if (params.merge_bam) {
            modbams = get_nanopolish_output(chrs, input) 
            modbam_output = merge_bam(modbams.modbam.groupTuple())
        }
        else {
            modbam_output = get_mergebam_output(input)
            coverage = calculate_coverage(
			    modbam_output.sample_name,
			    modbam_output.modbam
            )
            if (params.ref_freq){
            reference_frequency_output = downsample_reference_frequency(
                coverage.sample_name,
                coverage.modbam,
                coverage.sample_coverage,
                cvrg
                )
            }
            else {
                reference_frequency_output = get_ref_freq_output(coverage.sample_name,
                                                                 coverage.sample_coverage,
                                                                 cvrg)
            }
        }
        cpgs = get_cpgs(reference_frequency_output.groupTuple())
        deconv_output = deconvolve(cpgs.sample_name, cpgs.csv)
        plot_accuracy(deconv_output)
        if (params.clean_bams) {
            clean_bams(modbam_output.sample_name)
        }
    emit:
        deconv_output
}

workflow {
    runs = Channel.fromPath( 'data/MMinden_*', type: 'dir')

    // determine sample name for each input
    input = runs.map { [it.simpleName, it] }
    deconv_output = pipeline(input)
}
