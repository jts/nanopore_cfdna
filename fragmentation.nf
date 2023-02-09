
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
