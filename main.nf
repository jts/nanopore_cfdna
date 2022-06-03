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
    $params.guppy -r --num_callers $task.cpus --gpu_runners_per_device 4 --chunks_per_runner 512 -c $params.guppy_model -i $run_directory -x "cuda:0 cuda:1" -s basecalled
    """
}
process align_reads {
    cpus params.threads
    memory '32 G'

    publishDir "alignments", mode: 'symlink'
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

    publishDir "alignments", mode: 'symlink'

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
    
    publishDir "modifications/bams", mode: 'copy'

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
    
    publishDir "modifications/read_frequency", mode: 'symlink'

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
process fragmentation {
    cpus params.threads
    memory '32 G'
    time '1d'

    publishDir "fragmentation", mode: 'symlink'
    
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
process compute_coverage {
    input:
        tuple val(sample_name), val(modbam)
    output:
        tuple val(sample_name),
            val(modbam),
            stdout
    shell:
    """
    summary=~/jbroadbent/projects/methylation_deconvolution/pipeline_outputs/MLS/mls_dataset_summary.tsv
    grep -w $sample_name \$summary | cut -f4 
    """
    /*"""*/
    /*samtools coverage ${modbam} | head -n25 | tail -n24 | awk -F'\\t' '{sum+=\$7} END {print(sum/NR) }'*/
    /*"""*/
}
process region_modification_frequency {
    memory '32 G'
    time '7d'
    publishDir "modifications/region_frequency", mode: 'copy'
    input:
        tuple val(sample_name), file("${sample_name}.bam"), val(orig_cvrg)
        each atlas
        each cvrg
    output:
        tuple(val(sample_name), val(atlas), path("${sample_name}.${atlas}Atlas.${cvrg}.region_modifications.tsv"), emit: tsv)
    shell:
    """
    $params.mbtools region-frequency  ${sample_name}.bam -r ${projectDir}/atlases/${atlas}Atlas.bed --cpg  --reference-genome ~/simpsonlab/data/references/GRCh38_no_alt_analysis_set.GCA_000001405.15.fna -s ${(cvrg as float)/(orig_cvrg as float)} > ${sample_name}.${atlas}Atlas.${cvrg}.region_modifications.tsv 
    """
}
process deconvolve {
    cpus params.threads
    memory '32 G'
    time '1d'

    publishDir "deconvolution", mode: 'copy'
    
    input:
        tuple(val(sample_name), val(atlas), val(tsv))
        each model
    output:
        tuple val(sample_name),
            val(atlas),
            val(model),
            path("${sample_name}.${atlas}Atlas.${model}.deconv_output.tsv")
    shell:
    """
    $params.nanomix --model ${model} --atlas ${projectDir}/atlases/${atlas}Atlas.tsv ${tsv.join(" ")}  > "${sample_name}.${atlas}Atlas.${model}.deconv_output.tsv"
    """
}

process plot_deconvolution {

    publishDir "deconvolution", mode: 'copy'

    input:
        tuple val(sample_name),
            val(atlas),
            val(model),
            path("${sample_name}.${atlas}Atlas.${model}.deconv_output.tsv")
    output:
        file "*.png"
    shell:
    """
    ${projectDir}/scripts/plot_deconv.py "${sample_name}.${atlas}Atlas.${model}.deconv_output.tsv" -name ${sample_name}.${atlas}Atlas.${model}.deconv_output.png
    """
}
process plot_accuracy_by_atlas {
    publishDir "deconvolution_loss/${atlas}", mode: 'copy'

    input:
        tuple val(sample_name),
            val(atlas),
            val(model),
            val(tsv)
    output:
        path "deconvolution_loss.*${atlas}.png"
    shell:
    """
    ${projectDir}/scripts/plot_accuracy.r ${tsv.join(" ")} ${atlas}
    """
}
process plot_accuracy_by_model {
    publishDir "deconvolution_loss/${model}", mode: 'copy'

    input:
        tuple val(sample_name),
            val(atlas),
            val(model),
            val(tsv)
    output:
        path "deconvolution_loss.*${model}.png"
    shell:
    """
    ${projectDir}/scripts/plot_accuracy.r ${tsv.join(" ")} ${model}
    """
}
workflow pipeline {
    take:
        input
    main:
        cvrg = Channel.from([0.1,0.25,0.5,0.75,1,2,3,4,5,6,7,8,9,10,20])
        atlas = Channel.from("loyfer25","loyfer50","loyfer100","loyfer150","loyfer200","loyfer250")
        model = Channel.from("llse", "nnls")
        orig_cvrg = compute_coverage(input)
        region_frequency_output = region_modification_frequency(orig_cvrg, atlas, cvrg)
        deconv_output = deconvolve(region_frequency_output.tsv.groupTuple(by: [0,1]), model)
        plot_deconvolution(deconv_output)
    emit:
        deconv_output

}

workflow {
    input = Channel.fromPath('data/*.bam', type: 'file')
        .map {[it.simpleName, it]}
    deconv_output = pipeline(input)
    plot_accuracy_by_model(deconv_output.groupTuple(by: 2))
    plot_accuracy_by_atlas(deconv_output.groupTuple(by: 1))
}
