#!/usr/bin/env nextflow
// Structure based on https://github.com/epi2me-labs/wf-template/blob/master/main.nf
import groovy.json.JsonBuilder
nextflow.enable.dsl = 2
process BASECALL {
    debug true
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
    $params.guppy --recursive --compress_fastq --align_ref $params.reference  --num_callers $task.cpus --gpu_runners_per_device 5 --chunks_per_runner 512 -c $params.guppy_nanomix_model -i $run_directory -x "cuda:0 cuda:1" -s basecalled --bam_out
    """
}
process MERGE_BAM {
    cpus params.threads
    time '3d'
    memory { 8.GB * task.attempt }
    cache 'lenient'
    maxRetries 3
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    
    publishDir "${launchDir}/bams", mode: 'symlink', overwrite: true

    input:
        tuple val(sample_name), path(basecalled)
    output:
        tuple val(sample_name), path("${sample_name}.bam")
    shell:
    """
    samtools merge -o ${sample_name}.bam ${basecalled}/pass/*.bam --threads $task.cpus
    """
}
process MBTOOLS {
    publishDir "${launchDir}/methylomes", mode: 'symlink'
    memory { 32.GB * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    tag "$name"

    input:
        tuple val(name), val(bam)
    output:
        tuple path("*.methylome.tsv"), val(name)
    shell:
    """
    mbtools read-region-frequency \
        -r $params.atlas \
        --cpg \
        -g $params.reference \
        -f $params.read_end_filter \
        $bam > ${name}.methylome.tsv
    """
}
process DECONVOLUTE {
    publishDir "${launchDir}/mixture_proportions", mode: 'symlink'
    cpus params.threads

    tag "$name"

    input:
        tuple val(methylome), val(name)
    output:
        path("*sigma.tsv")
    shell:
    """
    nanomix deconvolute \
        -a $params.atlas \
        -@ $task.cpus \
        -nt $params.num_trials \
        --concentration $params.concentration \
        --model $params.nanomix_model \
        --nnls_init \
        -p11 $params.p11 \
        -p01 $params.p01 \
        -n 298 \
        -t 0.00001 \
        -p $params.min_proportion \
        $methylome > ${name}.sigma.tsv
    """
}
process PLOT {
    publishDir "${launchDir}/plots", mode: 'symlink'

    input:
        val(sigmas)
    output:
        path "*.png"
    shell:
    """
    nanomix plot -o deconvolution.png ${sigmas.join(' ')} 
    """
}
def check_for_bam = { sample_name, run_directory ->
    def bam = file("${launchDir}/bams/${sample_name}.bam")
    if (bam.exists()) {
        println "Found BAM file for ${sample_name} at ${bam}"
        return [sample_name, bam, true]
    } else {
        println "No BAM file found for ${sample_name} at ${bam}"
        return [sample_name, run_directory, false]
    }
}

workflow basecall {
    take:
        input
    main:
        BASECALL(input)
        MERGE_BAM(BASECALL.out)
    emit:
        MERGE_BAM.out
}
workflow deconvolute {
    take:
        modbam
    main:
        MBTOOLS(modbam)
        DECONVOLUTE(MBTOOLS.out)
        PLOT(DECONVOLUTE.out.collect())
}

// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
    }
    // Check if BAM files have been generated for input samples
    // If not, run basecalling
    Channel.fromPath('data/*', type: 'dir')
        .map {[it.simpleName, it, ]}
        .map(check_for_bam)
        .branch {
            modbam: it[2]
            basecall: !it[2]
        }
        .set { input }
    input.basecall
        .map { [it[0], it[1]] }
        .set { basecall_input }
    input.modbam
        .map { [it[0], it[1]] }
        .set { modbam_input }

    basecall(basecall_input)
        .set { basecalled }
    deconvolute( modbam_input.concat(basecalled) )

    /*software_versions = getVersions()*/
    /*workflow_params = getParams()*/
    /*report = makeReport(*/
        /*metadata, software_versions.collect(), workflow_params*/
        /*)*/
}
if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
    }

    workflow.onError {
        Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
    }
}

