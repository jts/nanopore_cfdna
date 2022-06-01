#! /.mounts/labs/simpsonlab/sw/miniconda3/envs/cfdna/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
library(karyoploteR)
library(glue)


chromosomes <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY')


unpack_sample <- function(s){
    s = tail(strsplit(s, "/")[[1]], n=1)
    sample_name = sub('.fragmentation.ratios.tsv', '', s)
    sample_name
}
i = 1
for (arg in args) {
    sample = unpack_sample(arg)
    pdf <- glue('{sample}.fragmentome.pdf') 
    pdf(pdf, width=14, height=6, pointsize=10)
    pp <- getDefaultPlotParams(plot.type=4)
    pp$data1min <- -5
    pp$data1max <- 5
    kp <- plotKaryotype(plot.type=4,main=glue('Ratio of Short (100-150bp) to Long (151-220bp) Fragments (normalized)\n{sample}'), plot.params = pp)
    kpDataBackground(kp, data.panel = 1)
    df <- read.table(arg, header=TRUE, sep='\t')
    for (chromosome in chromosomes){
        i <- which(df[, 1] == chromosome)
        x <- df[i,2]
        y <- df[i,4]
        kpLines(kp, chr=chromosome, x=x, y=y, data.panel=1)
    }
    kpAxis(kp, data.panel= 1, ymin=-5, ymax=5)
    dev.off()
}




