#! /.mounts/labs/simpsonlab/sw/miniconda3/envs/cfdna/bin/Rscript
library(ggplot2)
library(glue)
library(RColorBrewer)
args = commandArgs(trailingOnly=TRUE)
sample = strsplit(args[1], '_')[[1]][1]
print(sample)
title = args[1]
#file <- glue('~/jbroadbent/results/cfdna/{sample}.acchist.png')

df <- read.table(args[1], header=TRUE, sep='\t')


#sample_size <- sprintf("N = %d", nrow(df))
#title = glue('Sample {sample}, {sample_size}')
#ggplot(data = df, aes(x = acc))+
    #geom_histogram(binwidth = 1, color = "black", fill = "blue") +
    #xlim(70, 100)+
    #labs(title=title, x='Accuracy')

#ggsave(file, width=12, height=12)

file <- glue('~/jbroadbent/results/cfdna/{sample}.accreadlen.png')
# Set color palette for 2D heatmap
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

c <- cor(df$read_length, df$acc)
title <- sprintf("N = %d r = %.3f", nrow(df), c)

ggplot(df, aes(read_length, acc)) +
	geom_bin2d(bins=50) + scale_fill_gradientn(colors=r, trans="log10") +
	xlab("Read Length (bp)") +
    # reverse y axis
    scale_y_reverse() +
	ylab("Error Rate") +
    xlim(0, 1000) +
    # breaks every 150 bp in the x axis with labels
	theme_bw(base_size=14) +
	ggtitle(glue('Sample {sample}: {title}'))
# save plot
ggsave(file, width=12, height=12)


