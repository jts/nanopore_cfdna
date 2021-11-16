#!/usr/bin/env nextflow



process nanopolish_shard {
    cpus params.threads
    memory '32 G'
    time '3d'
    
    publishDir "${sample_name}_results/shards", mode: 'copy'

    input:
        val sample_name
        val shard_name
        file fastq
        file bam
        file bai
        tuple file(nano_index, nano_index.gz, nano_index.fai, nano_index.readdb) 
    output:
        val sample_name, emit: sample_name
        path "${shard_name}.modifications.sorted.bam", emit: modbam
    shell:
    """
    export HDF5_PLUGIN_PATH=/.mounts/labs/simpsonlab/users/jsimpson/code/cfdna_pipeline/etc/
    $params.nanopolish call-methylation -b ${sample_name}.sorted.bam -r ${sample_name}.fastq -g $params.reference --modbam-output ${shard_name}.modifications.sorted.bam -t $task.cpus
    """
}

process calculate_read_modification_frequency {
    cpus params.threads
    memory '32 G'
    time '1d'
    
    publishDir "${sample_name}_results", mode: 'copy'

    input:
        val sample_name
        file "${sample_name}.modifications.sorted.bam"
    output:
        val sample_name, emit: sample_name
        path "${sample_name}.read_modifications.tsv", emit: readmod
    shell:
    """
    $params.mbtools read-frequency ${sample_name}.modifications.sorted.bam > ${sample_name}.read_modifications.tsv
    """
}

process calculate_reference_frequency {
    cpus params.threads
    memory '32 G'
    time '1d'
    
    publishDir "${sample_name}_results", mode: 'copy'

    input:
        val sample_name
        file "${sample_name}.modifications.sorted.bam"
    output:
        val sample_name, emit: sample_name
        path "${sample_name}.reference_modifications.tsv", emit: readmod
    shell:
    """
    $params.mbtools reference-frequency ${sample_name}.modifications.sorted.bam > ${sample_name}.reference_modifications.tsv
    """
}

process calculate_fragmentation {
    cpus params.threads
    memory '32 G'
    time '1d'
    
    publishDir "${sample_name}_results", mode: 'copy'

    input:
        val sample_name
        file "${sample_name}.read_modifications.tsv"
    output:
        file "${sample_name}.fragmentation.ratios.tsv"
    shell:
    """
    python $params.fragmentation -o ${sample_name}.fragmentation.ratios.tsv -s 100 151 -l 151 221 -b 5000000 ${sample_name}.read_modifications.tsv

    """
}
process split_bams{
    cpus params.threads
    memory '32 G'
    time '1d'

    input:
        val sample_name
        file "${sample_name}.sorted.bam"
    output:
        tuple file "shard*.bam", file "shard*.bai"
    shell:
    """
    module load picard
    java -jar \$PICARDJAR SplitSamByNumberOfReads I="${sample_name}.sorted.bam O=$workDir SPLIT_TO_N_READS=$params.segment_nanopolish CREATE_INDEX=true
    """

}
process merge_bams {
    cpus params.threads
    memory '32 G'
    time '1d'
    
    publishDir "${sample_name}_results", mode: 'copy'

    input:
        val sample_name
        tuple bams
    output:
        path "${sample_name}.modifications.sorted.bam", emit: modbam
    shell:
    """
    samtools merge -o ${sample_name}.modifications.sorted.bam" ${bams}
    """
}


workflow pipeline {
    take:
        run_directory
    main:
        
        dir = file("data/*")

        allFiles = dir.list()
        for( def file : allFiles ) {
            println file
        }

        dir.eachFile {
            println "${item.getName()}"
        }

        /*align_output.bam = run_directory + "**.bam"*/
        /*align_output.bai = run_directory + "**.bai"*/
        /*fastq = run_directory + "**.fastq"*/
        /*nanopolish_index = "**.index*".collect(run_directory + it)*/

        /*nanopolish_index.view()*/

        /*if(params.segment_nanopolish > 0) {*/
            /*bams_split = split_bams(align_output.bam)*/
            /*fn =  nanopolish_shard(sample_name,*/
                                        /*it[0],*/
                                        /*merge_fastq_output,*/
                                        /*it[1][1], it[1][0],*/
                                        /*run_directory, index_output)*/
            /*bams_split.out.view()*/
            /*c = Channel.fromFilePairs("${sample_name}_results/shards/shard*.{bai,bam}")*/
            /*modbam_output = merge_bams(basecall_output.sample_name, c.collect(fn))*/


        /*} else {*/
            /*modbam_output = nanopolish_call_methylation(basecall_output.sample_name,*/
                                                     /*merge_fastq_output, align_output.bam,*/
                                                     /*align_output.bai, basecall_output.run_directory, index_output)*/
        /*}*/
        /*read_frequency_output = calculate_read_modification_frequency(modbam_output.sample_name, modbam_output.modbam)*/
        /*referency_frequency_output = calculate_reference_frequency(modbam_output.sample_name, modbam_output.modbam)*/
        /*fragmentation_output = calculate_fragmentation(read_frequency_output.sample_name, read_frequency_output.readmod)*/
}


workflow {

    runs = Channel.fromPath( 'data/*', type: 'dir')

    // determine sample name for each input
    input = runs.map( {[it.getSimpleName, it]} )
    input.view()

    pipeline(runs)

}
