#!/usr/bin/env nextflow

samples = Channel.fromPath('data/*', type: 'dir')

sample_names = Channel.create()

bams = Channel.create()
bais = Channel.create()
sample_files = Channel.create()

samples.map { 
        sample_names << it.getSimpleName()
        it.eachFileRecurse { item -> 
        if (item.getName().startsWith(it.getSimpleName())) {
            sample_files << item
            if (item.getName().endsWith('.bam')){ 
                println "${item.getName()}"
                bams << item
            }
            if (item.getName().endsWith('.bai')){ 
                println "${item.getName()}"
            }
            if (item.getName().endsWith('.fastq.gz')){ 
                println "${item.getName()}"
            }
            if (item.getName().endsWith('.index')){ 
                println "${item.getName()}"
            }
            if (item.getName().endsWith('.index.gz')){ 
                println "${item.getName()}"
            }
            if (item.getName().endsWith('.index.fai')){ 
                println "${item.getName()}"
            }
            if (item.getName().endsWith('.index.readdb')){ 
                println "${item.getName()}"
            }
        } 
    }
}


process split_bams {
    cpus params.threads
    memory '1 G'
    time '1d'

    input:
        val sample_name from sample_names
        file "$sample_name*.bam.bai" from sample_files
        file "${sample_name}*.bam" from sample_files
    output:
        /*tuple file("shard*.bam"), file("shard*.bai")*/
    shell:
    """
    echo $sample_name*.bam
    """
/*    module load picard*/
    /*java -jar \$PICARDJAR SplitSamByNumberOfReads I=${sample_name}.bam O=$workDir SPLIT_TO_N_READS=$params.segment_nanopolish CREATE_INDEX=true*/
    /*"""*/

}

//TODO: Map each file in channel to hash map

//TODO: write process to split bams and create
// set of channels from those bams


/*process nanopolish_shard {*/
    /*cpus params.threads*/
    /*memory '32 G'*/
    /*time '3d'*/
    
    /*publishDir "${sample_name}_results/shards", mode: 'copy'*/

    /*input:*/
        /*val sample_name from samples*/
        /*val shard_name*/
        /*file fastq*/
        /*file bam*/
        /*file bai*/
        /*tuple file(nano_index, nano_index.gz, nano_index.fai, nano_index.readdb) */
    /*output:*/
        /*val sample_name, emit: sample_name*/
        /*path "${shard_name}.modifications.sorted.bam", emit: modbam*/
    /*shell:*/
    /*"""*/
    /*export HDF5_PLUGIN_PATH=/.mounts/labs/simpsonlab/users/jsimpson/code/cfdna_pipeline/etc/*/
    /*$params.nanopolish call-methylation -b ${sample_name}.sorted.bam -r ${sample_name}.fastq -g $params.reference --modbam-output ${shard_name}.modifications.sorted.bam -t $task.cpus*/
    /*"""*/
/*}*/

/*process calculate_read_modification_frequency {*/
    /*cpus params.threads*/
    /*memory '32 G'*/
    /*time '1d'*/
    
    /*publishDir "${sample_name}_results", mode: 'copy'*/

    /*input:*/
        /*val sample_name*/
        /*file "${sample_name}.modifications.sorted.bam"*/
    /*output:*/
        /*val sample_name, emit: sample_name*/
        /*path "${sample_name}.read_modifications.tsv", emit: readmod*/
    /*shell:*/
    /*"""*/
    /*$params.mbtools read-frequency ${sample_name}.modifications.sorted.bam > ${sample_name}.read_modifications.tsv*/
    /*"""*/
/*}*/

/*process calculate_reference_frequency {*/
    /*cpus params.threads*/
    /*memory '32 G'*/
    /*time '1d'*/
    
    /*publishDir "${sample_name}_results", mode: 'copy'*/

    /*input:*/
        /*val sample_name*/
        /*file "${sample_name}.modifications.sorted.bam"*/
    /*output:*/
        /*val sample_name, emit: sample_name*/
        /*path "${sample_name}.reference_modifications.tsv", emit: readmod*/
    /*shell:*/
    /*"""*/
    /*$params.mbtools reference-frequency ${sample_name}.modifications.sorted.bam > ${sample_name}.reference_modifications.tsv*/
    /*"""*/
/*}*/

/*process calculate_fragmentation {*/
    /*cpus params.threads*/
    /*memory '32 G'*/
    /*time '1d'*/
    
    /*publishDir "${sample_name}_results", mode: 'copy'*/

    /*input:*/
        /*val sample_name*/
        /*file "${sample_name}.read_modifications.tsv"*/
    /*output:*/
        /*file "${sample_name}.fragmentation.ratios.tsv"*/
    /*shell:*/
    /*"""*/
    /*python $params.fragmentation -o ${sample_name}.fragmentation.ratios.tsv -s 100 151 -l 151 221 -b 5000000 ${sample_name}.read_modifications.tsv*/

    /*"""*/
/*}*/
/*process merge_bams {*/
    /*cpus params.threads*/
    /*memory '32 G'*/
    /*time '1d'*/
    
    /*publishDir "${sample_name}_results", mode: 'copy'*/

    /*input:*/
        /*val sample_name*/
        /*tuple bams*/
    /*output:*/
        /*path "${sample_name}.modifications.sorted.bam", emit: modbam*/
    /*shell:*/
    /*"""*/
    /*samtools merge -o ${sample_name}.modifications.sorted.bam" ${bams}*/
    /*"""*/
/*}*/


