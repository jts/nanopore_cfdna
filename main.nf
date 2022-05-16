#!/usr/bin/env nextflow
// Structure based on https://github.com/epi2me-labs/wf-template/blob/master/main.nf

nextflow.enable.dsl = 2

process basecall_reads {
    label 'GPU'

    memory '16 G'
    cpus params.threads
    time '3d'
    queue 'gpu.q'
    cache 'lenient'

    publishDir "basecalls/${sample_name}", mode: "symlink", overwrite: false

    input:
        tuple val(sample_name), path(run_directory)
    output:
        val sample_name, emit: sample_name
        path run_directory, emit: run_directory
        path "basecalled", emit: basecalled_directory
    shell:
    """
    $params.guppy -r --num_callers $task.cpus --gpu_runners_per_device 4 --chunks_per_runner 512 -c $params.guppy_model -i $run_directory -x "cuda:0 cuda:1" -s basecalled
    """
}
process get_basecall_reads {
    cache 'lenient'
    input:
        tuple val(sample_name), path(run_directory)
    output:
        val sample_name
        path run_directory
        path "basecalled"
    shell:
    """
    ln -s ${launchDir}/basecalls/${sample_name}/basecalled basecalled
    """
}
process merge_reads {
    cpus 1
    memory '1 G'
    
    publishDir "basecalls/${sample_name}", mode: "symlink"
    cache 'lenient'

    input: 
        val sample_name
        path run_directory
        path basecall_directory
    output:
        val sample_name
        path run_directory
        path "basecalled"
        path "merged_reads/${sample_name}.fastq"
    shell:
    """
    mkdir -p merged_reads
    fastcat -x ${basecall_directory} > merged_reads/${sample_name}.fastq
    """
}
process align_reads {
    cpus params.threads
    memory '32 G'
    time '2d'
    cache 'lenient'

    publishDir "alignments", mode: 'symlink'
    input:
        val sample_name
        path run_directory
        path "basecalled"
        path "${sample_name}.fastq"
    output:
        tuple val(sample_name),
            path("${sample_name}.sorted.bam"),
            path("${sample_name}.sorted.bam.bai")
    shell:
    """
    minimap2 -ax map-ont -t $task.cpus $params.reference ${sample_name}.fastq | samtools sort -T tmp > ${sample_name}.sorted.bam
    samtools index ${sample_name}.sorted.bam
    """
}
process bam_stats {
    memory '32 G'
    time '1d'
    cache 'lenient'

    publishDir "alignments", mode: 'symlink'

    input:
        tuple val(sample_name),
            path("${sample_name}.sorted.bam"),
            path("${sample_name}.sorted.bam.bai")
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
    cache 'lenient'

    input:
        val sample_name
        path run_directory
        path basecalled
        path "${sample_name}.fastq"
    output:
        tuple val(sample_name),
            path(run_directory),
            path("${sample_name}.fastq"),
            path("${sample_name}.fastq.index"),
            path("${sample_name}.fastq.index.gzi"),
            path("${sample_name}.fastq.index.fai"),
            path("${sample_name}.fastq.index.readdb")
    shell:
    """
    $params.nanopolish index -s ${basecalled}/sequencing_summary.txt -d ${run_directory}/ ${sample_name}.fastq
    """
}

process nanopolish_call_methylation {
    cpus params.threads
    memory '32 G'
    time '7d'
    cache 'lenient'
    
    input:
        each chr
        tuple val(sample_name),
            path("${sample_name}.sorted.bam"),
            path("${sample_name}.sorted.bam.bai"),
            path(run_directory),
            path("${sample_name}.fastq"),
            path("${sample_name}.fastq.index"),
            path("${sample_name}.fastq.index.gzi"),
            path("${sample_name}.fastq.index.fai"),
            path("${sample_name}.fastq.index.readdb")
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
    
    publishDir "modifications/bams", mode: 'symlink', overwrite: true

    input:
        tuple( val(sample_name), val(bams) )
    output:
        val sample_name
        path "${sample_name}.merged.modifications.sorted.bam"
    shell:
    """
    samtools merge -o ${sample_name}.merged.modifications.sorted.bam ${bams.join(' ')} --threads $params.threads
    """
}
process read_modification_frequency {
    cpus params.threads
    memory '32 G'
    time '4d'
    
    publishDir "modifications/read_frequency", mode: 'symlink', overwrite: true

    input:
        val sample_name
        file "${sample_name}.merged.modifications.sorted.bam"
    output:
        val sample_name
        path "${sample_name}.read_modifications.tsv"
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
process region_modification_frequency {
    memory '32 G'
    time '7d'
    publishDir "modifications/region_frequency", mode: 'copy'
    input:
        val sample_name
        file "${sample_name}.bam"
        each atlas
    output:
        tuple(val(atlas), path("${sample_name}.${atlas}Atlas.region_modifications.tsv"), emit: tsv)
    shell:
        """
        $params.mbtools region-frequency  ${sample_name}.bam -r ${projectDir}/atlases/${atlas}Atlas.bed --cpg  --reference-genome ~/simpsonlab/data/references/GRCh38_no_alt_analysis_set.GCA_000001405.15.fna > ${sample_name}.${atlas}Atlas.region_modifications.tsv 
        """
}
process deconvolve {
    cpus params.threads
    memory '32 G'
    time '1d'

    publishDir "deconvolution", mode: 'copy'
    
    input:
        tuple(val(atlas), val(tsv))
    output:
        tuple(path("${atlas}Atlas.llse.deconv_output.tsv"),
        path("${atlas}Atlas.nnls.deconv_output.tsv"))
    shell:
    """
    ${projectDir}/scripts/nanomix.py --model llse --atlas ${projectDir}/atlases/${atlas}Atlas.tsv ${tsv.join(" ")}  > "${atlas}Atlas.llse.deconv_output.tsv"
    ${projectDir}/scripts/nanomix.py --model nnls --atlas ${projectDir}/atlases/${atlas}Atlas.tsv ${tsv.join(" ")} > "${atlas}Atlas.nnls.deconv_output.tsv"
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
workflow pipeline {
    take:
        input
    main:
        if (params.basecall){
            basecall_output = basecall_reads(input)
        }else{
            basecall_output = get_basecall_reads(input)
        }
        merge_fastq_output = merge_reads(basecall_output)   
        index_output = nanopolish_index(merge_fastq_output)
        align_output = align_reads(merge_fastq_output)
        bam_stats(align_output)

        chrs = Channel.from(1 .. 22).map { "chr" + it }
        chrs_sex = Channel.of("chrX", "chrY")
        chrs = chrs.concat(chrs_sex)
        modbams = nanopolish_call_methylation(chrs, align_output.join(index_output))
        modbam_output = merge_bam(modbams.modbam.groupTuple())
        read_frequency_output = read_modification_frequency(modbam_output)
        fragmentation_output = fragmentation(read_frequency_output)
        plot_fragmentation(fragmentation_output.collect())

        atlas = Channel.from("Berman", "moss", "cheng", "chengOrig")
        region_frequency_output = region_modification_frequency(modbam_output, atlas)
        deconv_output = deconvolve(region_frequency_output.tsv.groupTuple())
        plot_deconvolution(deconv_output)
}

workflow {
    input = Channel.fromPath( 'data/*', type: 'dir')
        .map { [it.simpleName, it] }
    pipeline(input)
}
