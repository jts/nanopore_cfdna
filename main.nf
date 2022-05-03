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
process region_modification_frequency {
    memory '32 G'
    time '7d'
    publishDir "modifications/region_frequency", mode: 'copy'
    input:
        tuple val(sample_name), file("${sample_name}.bam")
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
        ${projectDir}/scripts/plot_deconv_berman.py \$f -s HU10,HU12,HU11,bc05,bc02,bc04,bc03,bc09,bc08,bc01,S1,bc11,bc10
    done
    """
}
workflow pipeline {
    take:
        input
    main:
        atlas = Channel.from("Berman", "moss", "cheng", "chengOrig")
        region_frequency_output = region_modification_frequency(input, atlas)
        deconv_output = deconvolve(region_frequency_output.tsv.groupTuple())
        plot_deconvolution(deconv_output)
}

workflow {
    input = Channel.fromPath('data/*.bam', type: 'file')
        .map {[it.simpleName, it, ]}
    pipeline(input)
}
