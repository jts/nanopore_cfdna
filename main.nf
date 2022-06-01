#!/usr/bin/env nextflow
// Structure based on https://github.com/epi2me-labs/wf-template/blob/master/main.nf

nextflow.enable.dsl = 2

process downsample_modification_freq {
    memory '32 G'
    time '7d'

    publishDir "mbtools", mode: "symlink"

    input:
        val sample_name
        file "${sample_name}.modifications.sorted.bam"
        val coverage
        each target_coverage
    output:
        tuple( val(sample_name), path("${sample_name}_${target_coverage}.modifications.tsv"), emit: modbam )
    shell:
    if (target_coverage < coverage as float)
        """
        $params.mbtools region-frequency  ${sample_name}.modifications.sorted.bam -r ${projectDir}/atlases/meth_atlas.bed -D ${target_coverage} -C ${coverage as float} > ${sample_name}_${target_coverage}.modifications.tsv 
        """
    else
        """
        touch ${sample_name}_${target_coverage}.modifications.tsv
        """
}
process modification_freq {
    memory '32 G'
    time '7d'

    publishDir "mbtools", mode: "symlink"

    input:
        val sample_name
        file "${sample_name}.modifications.sorted.bam"
    output:
        tuple( val(sample_name), path("${sample_name}_fullCoverage.modifications.tsv"), emit: modbam )
    shell:
        """
        $params.mbtools region-frequency  ${sample_name}.modifications.sorted.bam -r ${projectDir}/atlases/meth_atlas.bed > ${sample_name}_fullCoverage.modifications.tsv 
        """
}
process get_modification_freq_output {

    input:
        val sample_name
        val coverage
        each target_coverage
    output:
        tuple (val(sample_name), path("${sample_name}.C*C.modifications.tsv"), emit: refmod)
    shell:
    if (target_coverage < coverage as float)
        """
        ln -s \$(ls \$(find $workDir -type f -name "${sample_name}.C${target_coverage}C.modifications.tsv") -t | head -n1) ${sample_name}.C${target_coverage}C.modifications.tsv
        """
    else
        """
        ln -s \$(ls \$(find $workDir -type f -name "${sample_name}.C${coverage as float}C.modifications.tsv") -t | head -n1) ${sample_name}.C${coverage as float}C.modifications.tsv
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
process deconvolve {
    cpus params.threads
    memory '32 G'
    time '1d'

    publishDir "${sample_name}_results", mode: 'symlink'
    
    input:
        tuple(val(sample_name), val(modbams))
    output:
        tuple(path("${sample_name}-llse.deconv_output.csv"), path("${sample_name}-nnls.deconv_output.csv"))
    shell:
    """
    ${projectDir}/scripts/nanomix.py --atlas ${projectDir}/atlases/meth_atlas.tsv ${modbams.join(" ")} --fill > "${sample_name}-llse.deconv_output.csv"
    ${projectDir}/scripts/nanomix.py --atlas ${projectDir}/atlases/meth_atlas.tsv --model nnls ${modbams.join(" ")} --fill > "${sample_name}-nnls.deconv_output.csv"
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
        cvrg = Channel.from([0.1,0.25,0.5,0.75,1,2,3,4,5,6,7,8,9,10,20])
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
            /*coverage = calculate_coverage(*/
				/*modbam_output.sample_name,*/
				/*modbam_output.modbam*/
            /*)*/
            if (params.run_mbtools){
            modification_freq_output = modification_freq(
                modbam_output.sample_name,
                modbam_output.modbam,
                )
            }
            else {
                modification_freq_output = get_modification_freq_output(coverage.sample_name,
                                                                 coverage.sample_coverage,
                                                                 cvrg)
            }
        }
        deconv_output = deconvolve(modification_freq_output.modbam.groupTuple())
        if (params.clean_bams) {
            clean_bams(modbam_output.sample_name)
        }
    emit:
        deconv_output
}

workflow {
    runs = Channel.fromPath( 'data/*')

    // determine sample name for each input
    input = runs.map { [it.simpleName, it] }
    deconv_output = pipeline(input)

    plot_accuracy(deconv_output.collect())

}
