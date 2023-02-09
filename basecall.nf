
nextflow.enable.dsl = 2
process BASECALL{
    label 'GPU'

    memory '16 G'
    cpus params.threads
    time '12d'
    queue 'gpu.q'
    cache 'lenient'

    publishDir "basecalls/${sample_name}", mode: "symlink", overwrite: false

    input:
        tuple val(sample_name), path(run_directory)
    output:
        tuple val(sample_name), path("basecalled")
    shell:
    """
    $params.guppy -r --num_callers $task.cpus --gpu_runners_per_device 4 --chunks_per_runner 512 -c $params.guppy_model -i $run_directory -x "cuda:0 cuda:1" -s basecalled --bam_out
    """
}
process BAM_STATS {
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
process MERGE_BAM{
    cpus 8
    memory '32 G'
    time '1d'
    
    publishDir "modifications/bams", mode: 'symlink', overwrite: true

    input:
        tuple val(sample_name), path(basecalled)
    output:
        tuple val(sample_name), path("${sample_name}.bam")
    shell:
    """
    samtools merge -o ${sample_name}.bam ${basecalled}/pass/*.bam --threads $task.cpus
    """
}
