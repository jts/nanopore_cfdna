#!/usr/bin/env nextflow
// Structure based on https://github.com/epi2me-labs/wf-template/blob/master/main.nf

nextflow.enable.dsl = 2
process fragmentation {
    cpus params.threads
    memory '32 G'
    time '1d'

    publishDir "fragmentation", mode: 'copy'
    
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
    publishDir "fragmentation", mode: 'copy'

    input:
        val fragmentomes
    output:
        file "*fragmentome.pdf"
    shell:

    """
    ${projectDir}/scripts/plot_fragmentome.r ${fragmentomes.join(" ")}
    """
}
process plot_deconvolution {

    publishDir "deconvolution", mode: 'copy'

    input:
        val deconv_outputs
    output:
        file "*.png"
    shell:
    """
    for f in ${deconv_outputs.join(" ")}
        do name=\$(echo \$f | sed s/.tsv/.png/g) 
        ${projectDir}/scripts/plot_deconv.py \$f 
    done
    """
}
process region_modification_frequency {
    memory '32 G'
    time '7d'
    publishDir "modifications/read_frequency", mode: 'copy'
    cache 'lenient'
    input:
        tuple val(sample_name), path("${sample_name}.bam")
        each atlas
    output:
        tuple(val(sample_name), val(atlas), path("${sample_name}.${atlas}Atlas.read_modifications.tsv"), emit: tsv)
    shell:
        """
        bedtools intersect -abam ${sample_name}.bam -b ${projectDir}/atlases/${atlas}Atlas.bed > ${sample_name}.${atlas}.intersect.bam
        $params.mbtools read-frequency  ${sample_name}.${atlas}.intersect.bam --cpg  -g ~/simpsonlab/data/references/GRCh38_no_alt_analysis_set.GCA_000001405.15.fna \
                | cut -f2,3,4,8,9 > ${sample_name}.${atlas}Atlas.read_modifications.tsv 
        """
}
process deconvolute {
    cpus params.threads
    memory '32 G'
    time '1d'

    label 'parallel'

    publishDir "deconvolution", mode: 'copy'
    
    input:
        tuple(val(sample_name), val(atlas), val(tsv))
    output:
        tuple val(atlas),
        path("${atlas}Atlas.llse.deconv_output.tsv"),
        path("${atlas}Atlas.llse_em.deconv_output.tsv")
    shell:
    """
    ~/jbroadbent/code/cfdna/nanomix_em/target/release/nanomix_em deconvolute \
            ${tsv.join(" ")} \
            --atlas ${projectDir}/atlases/${atlas}Atlas.tsv \
            --p01 0.0909 --p11 0.949\
            --max_iter 300\
            --min_proportion 0.015\
            --stop_thresh 0.001 > "${atlas}Atlas.llse_em.deconv_output.tsv"
    """
}

process plot_deconvolution {

    publishDir "deconvolution", mode: 'copy'

    input:
        tuple val(atlas),
        path("${atlas}Atlas.llse.deconv_output.tsv"),
        path("${atlas}Atlas.llse_em.deconv_output.tsv")
    output:
        file "*.png"
    shell:
    """
    ${projectDir}/scripts/plot_deconv.py "${atlas}Atlas.llse.deconv_output.tsv" -name ${atlas}Atlas.llse.deconv_output.png
    ${projectDir}/scripts/plot_deconv.py "${atlas}Atlas.llse_em.deconv_output.tsv" -name ${atlas}Atlas.llse_em.deconv_output.png
    """
}
workflow pipeline {
    take:
        input
    main:
        basecall_output = basecall_reads(input)
        modbam_output = merge_bam(basecall_output)

        atlas = Channel.from("berman", "loyfer250")
        region_frequency_output = region_modification_frequency(modbam_output, atlas)
        deconv_output = deconvolute(region_frequency_output.groupTuple(by: 1))
        plot_deconvolution(deconv_output.groupTuple(by: [1, 2] ))
}

workflow {
    input = Channel.fromPath('data/*', type: 'dir')
        .map {[it.simpleName, it, ]}
    pipeline(input)
}
