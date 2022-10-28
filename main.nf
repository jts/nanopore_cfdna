#!/usr/bin/env nextflow
// Structure based on https://github.com/epi2me-labs/wf-template/blob/master/main.nf

nextflow.enable.dsl = 2

process simulate {
    memory '32 G'
    time '1d'

    publishDir "methylomes", mode: 'copy'

    input:
        each coverage
        each p01
        each p11
        each lung_proportion
        each atlas
        each repeat
    output:
        val coverage
        val p01
        val p11
        val lung_proportion
        val atlas
        val repeat
        path "${coverage}-${p01}-${p11}-${lung_proportion}-${atlas}-methylome.tsv"
    shell:
    if ( atlas == "bermanAtlas")
        """
        nanomix_em simulate --atlas ${projectDir}/atlases/${atlas}.tsv --p01 $p01 --p11 $p11 --coverage $coverage --lung_proportion $lung_proportion\
                --region_size 1 \
                --name "${coverage}-${p01}-${p11}-${lung_proportion}-${atlas}-methylome.tsv"
        """
    else
        """
        nanomix_em simulate --atlas ${projectDir}/atlases/${atlas}.tsv --p01 $p01 --p11 $p11 --coverage $coverage --lung_proportion $lung_proportion\
                --region_size 5 \
                --name "${coverage}-${p01}-${p11}-${lung_proportion}-${atlas}-methylome.tsv"
        """
}
process deconvolute {
    cpus params.threads
    memory '64 G'
    time '1d'

    publishDir "deconvolution", mode: 'copy'
    
    input:
        val coverage
        val p01
        val p11
        val lung_proportion
        val atlas
        val repeat
        path "${coverage}-${p01}-${p11}-${lung_proportion}-${atlas}-methylome.tsv"
        each model
        each init_sigma
    output:
        tuple val(coverage),
        val(p01),
        val(p11),
        val(lung_proportion),
        val(atlas),
        val(init_sigma),
        (stdout)
    script:
    if ( model == 'mmse' )
        """
        head -n1 ${coverage}-${p01}-${p11}-${lung_proportion}-${atlas}-methylome.tsv > header.tsv
        tail -n+2 ${coverage}-${p01}-${p11}-${lung_proportion}-${atlas}-methylome.tsv | bedtools merge -c 4,5 -o sum -i - | cat header.tsv - | cut -f1-5 > aggregated_methylome.tsv
        nanomix.py --model $init_sigma --atlas ${projectDir}/atlases/${atlas}.tsv aggregated_methylome.tsv --p11 $p11 --p01 $p01 > init_sigma.tsv

        nanomix_em deconvolute --atlas ${projectDir}/atlases/${atlas}.tsv  \
                --methylome ${coverage}-${p01}-${p11}-${lung_proportion}-${atlas}-methylome.tsv \
                --p01 $p01 --p11 $p11 --max_iter 2000 --stop_thresh 0.00002 --min_proportion 0.02 --sigma init_sigma.tsv \
                > output.tsv
                grep lung -i output.tsv| grep alveolar -i -v| cut -f2 | head -n1
        """
                /*grep log-likelihood output.tsv | cut -f2*/
    else
        """
        head -n1 ${coverage}-${p01}-${p11}-${lung_proportion}-${atlas}-methylome.tsv > header.tsv
        tail -n+2 ${coverage}-${p01}-${p11}-${lung_proportion}-${atlas}-methylome.tsv | bedtools merge -c 4,5 -o sum -i - | cat header.tsv - | cut -f1-5 > aggregated_methylome.tsv
        nanomix.py --model $model --atlas ${projectDir}/atlases/${atlas}.tsv  \
                aggregated_methylome.tsv\
                --p11 $p11 --p01 $p01 --random_inits \
                > output.tsv
                grep lung -i output.tsv| grep alveolar -i -v| cut -f2 |  head -n1
        """
}
process init_sigma {
    input:
        val coverage
        val p01
        val p11
        val lung_proportion
        val atlas
        val repeat
        path "${coverage}-${p01}-${p11}-${lung_proportion}-${atlas}-methylome.tsv"
        each init
    output:
        path "init_sigma.tsv"
    script:
        """
        head -n1 ${coverage}-${p01}-${p11}-${lung_proportion}-${atlas}-methylome.tsv > header.tsv
        tail -n+2 ${coverage}-${p01}-${p11}-${lung_proportion}-${atlas}-methylome.tsv | bedtools merge -c 4,5 -o sum -i - | cat header.tsv - | cut -f1-5 > aggregated_methylome.tsv
        nanomix.py --model $init --atlas ${projectDir}/atlases/${atlas}.tsv aggregated_methylome.tsv --p11 $p11 --p01 $p01 > init_sigma.tsv
        """
}
process evaluate_deconv {
    cpus params.threads
    memory '64 G'
    time '1d'

    publishDir "class_probs/${coverage}-${p01}-${p11}-${lung_proportion}-${atlas}-methylome/exp${repeat}", mode: 'copy'
    
    input:
        val coverage
        val p01
        val p11
        val lung_proportion
        val atlas
        val repeat
        path "${coverage}-${p01}-${p11}-${lung_proportion}-${atlas}-methylome.tsv"
        each init_sigma
    output:
        path "sigma.tsv"
        path "class_probs.tsv"
        path "incorrect_assignments.tsv"
    shell:
        """
        nanomix_em evaluate --atlas ${projectDir}/atlases/${atlas}.tsv  \
                --methylome ${coverage}-${p01}-${p11}-${lung_proportion}-${atlas}-methylome.tsv \
                --p01 $p01 --p11 $p11 --stop_thresh 0.000005 --min_proportion 0.02 --sigma ${init_sigma}
        """
}
process plot_deconvolution_loss {
    publishDir "plots", mode: 'copy'

    input:
        file lung_proportion  
    output:
        path "*deconvolution_loss.png"
    script:
    """
    ${projectDir}/scripts/plot_deconv_loss.r ${lung_proportion}
    """
}
workflow {
    coverage = Channel.from(0.1, 0.5, 1.0, 5.0, 10.0)
    p01 = Channel.from(0..5).map {it/100}
    p11 = Channel.from(100..95).map {it/100}
    repeats = Channel.from(0..4)
    true_lung_proportion = Channel.from(0.3)
    atlas = Channel.from("bermanAtlas")
    model = Channel.from("mmse")
    init = Channel.from("llse", "nnls", "null")
    lung_proportion = file("lung_proportion.tsv")
    lung_proportion.text = "coverage\tp01\tp11\ttrue_lung_proportion\tatlas\tinit_sigma\tlung_proportion\n"

    methylome = simulate(coverage, p01, p11, true_lung_proportion, atlas, repeats)
    /*evaluate_deconv(methylome, init_sigma)*/
    mixture = deconvolute(methylome, model, init)
    mixture.subscribe { assert it.size() == 7
                        lung_proportion.append(it.join("\t")) }
}
