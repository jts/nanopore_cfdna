#! /.mounts/labs/simpsonlab/sw/miniconda3/envs/cfdna/bin/Rscript
library(ggplot2)
library(data.table)
library(ggpubr)
args = commandArgs(trailingOnly = TRUE)

df = read.table(args[1], sep='\t', header=TRUE)
df = df[df$coverage== 1.0, ]
df = df[(
        (df$p11 == 1.0) 
        #(df$p11 == 0.97) |
        #(df$p11 == 0.95)
        )&(
        #(df$p01 == 0.05) |
        #(df$p01 == 0.03) |
        (df$p01 == 0.00)
        ),]
atlas = df$atlas[1]
pseudozero = exp(-6)
df$lung_proportion[df$lung_proportion == 0] <- pseudozero
plot = ggplot(data = df, aes(x=true_lung_proportion, y=lung_proportion, color=model)) +
    #geom_smooth(method = "lm", aes(fill=model) ) + 
    #geom_boxplot(width=.5) +
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width=.1, size=1.5) +
    stat_summary(fun = "mean", geom = "point", size=3, fill = 'black') +
    stat_cor(method = "pearson") +
    labs(title = "Methylation Deconvolution of Simulated Monocyte/Lung Mixture",
         subtitle = atlas,
         x = "True Lung Proportion",
         color = "Model",
         y = "True Lung Proportion") +
    theme(axis.title.y = element_text(angle=90),
          plot.title = element_text(size = 20) ) +
    scale_x_log10(limits = c(pseudozero, 0.4) ) +
    scale_y_log10(limits = c(pseudozero, 0.4) ) +
    #scale_x_continuous(breaks = seq(0.0, 0.5, 0.1) ) +
    #scale_y_continuous(breaks = seq(0.0, 0.5, 0.1) ) + 
    geom_abline(aes(intercept = 0, slope=1)) +
    facet_grid(cols = vars(p01), rows = vars(p11))
ggsave("scatter.png", width=14, height=10)
