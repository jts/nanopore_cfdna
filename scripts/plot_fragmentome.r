#! /.mounts/labs/simpsonlab/sw/miniconda3/envs/cfdna/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
library(karyoploteR)


chromosomes <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY')

pdf <- '~/jbroadbent/results/cfdna/fragmentome.pdf' 
pdf(pdf, width=18, height=6, pointsize=10)
pp <- getDefaultPlotParams(plot.type=3)
pp$data1min <- -5 
pp$data1max <- 5
pp$data1height <- 10
pp$data2min <- -5 
pp$data2max <- 5
pp$data2height <- 10
kp <- plotKaryotype(plot.type=3, main='Fragmentome Plot', plot.params = pp)
kpDataBackground(kp, data.panel = 1)
kpDataBackground(kp, data.panel = 2)

frag_ratio_normal <- args[2]
df <- read.table(frag_ratio_normal, header=TRUE, sep='\t')
j <- 0
for (chromosome in chromosomes){
	i <- which(df[, 1] == chromosome)
	x <- df[i,2]
	y <- df[i,4]
	kpLines(kp, chr=chromosome, x=x, y=y, data.panel=1)
	j <- j + 1
	kpAddLabels(kp, labels='control', data.panel = 1)
}

frag_ratio_cancer <- args[1] 
df <- read.table(frag_ratio_cancer, header=TRUE, sep='\t')
j <- 0
for (chromosome in chromosomes){
	i <- which(df[, 1] == chromosome)
	x <- df[i,2]
	y <- - df[i,4]
	kpLines(kp, chr=chromosome, x=x, y=y, data.panel=2)
	j <- j + 1
	kpAddLabels(kp, labels='tumour', data.panel = 2)
}

dev.off()


