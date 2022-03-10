#!/usr/bin/env nextflow
// Structure based on https://github.com/epi2me-labs/wf-template/blob/master/main.nf

nextflow.enable.dsl = 2

process get_mergebam_output {
    publishDir "${sample_name}_results", mode: 'symlink'

    input:
        tuple val(sample_name), file("run_directory")
    output:
        val sample_name, emit: sample_name
        path "${sample_name}.merged.modifications.sorted.bam", emit: modbam
    shell:
    """
    ln -s \$(find ${params.merge_bam_dir} -name ${sample_name}.merged.modifications.sorted.bam) ${sample_name}.merged.modifications.sorted.bam
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

    publishDir "refmods", mode: "copy"

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
        val modbam
    output:
        val sample_name, emit: sample_name
        val modbam, emit: modbam
        stdout emit: sample_coverage
    shell:
    """
    samtools coverage ${modbam} | head -n25 | tail -n24 | awk -F'\\t' '{sum+=\$7} END {print(sum/NR) }'
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
    $projectDir/scripts/nanomix.py --model llse ${sample_name}.cpgs.csv > "${sample_name}_llse.cpgs_deconv_output.csv"
    $projectDir/scripts/nanomix.py --model nnls ${sample_name}.cpgs.csv > "${sample_name}_nnls.cpgs_deconv_output.csv"
    """
}
process plot_accuracy {

    publishDir "plots", mode: 'symlink'

    input:
        val deconv_output
    output:
        file "*coverageV*.png"
        file "*coverageV*.tsv"
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
        cvrg = Channel.from([20])
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
        cpgs = get_cpgs(reference_frequency_output.groupTuple())
        deconv_output = deconvolve(cpgs.sample_name, cpgs.csv)
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

    plot_accuracy(deconv_output.collect())

}
