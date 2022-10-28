#! /.mounts/labs/simpsonlab/sw/miniconda3/envs/cfdna/bin/Rscript
library(ggplot2)
library(data.table)
library(dplyr)
args = commandArgs(trailingOnly = TRUE)
df = read.table(args[1], sep='\t', header=TRUE)
df <- df[df$true_lung_proportion == 0.3, ]
df <- df %>%
    group_by(coverage, p01, p11, init_sigma, true_lung_proportion, atlas) %>%
    summarise(mean_lung_proportion = mean(lung_proportion))
atlas = df$atlas[1]
plot = ggplot(df, aes(p01, p11) ) +
    geom_tile(aes(fill=mean_lung_proportion), color="grey50") +
    labs(title = "Methylation Deconvolution of Simulated 70:30 Monocyte/Lung Mixture",
         subtitle = atlas,
         x = "p01 (miscalled methylation)",
         y = "p11 (correctly called methylation)",
         fill = "Lung Proportion") +
    theme(axis.text.x = element_text(angle=90),
          axis.title.y = element_text(angle=90),
          plot.title = element_text(size = 20),
          strip.text.y = element_text(size=15),
          strip.text.x = element_text(size=15)) +
    scale_x_continuous(breaks = seq(0.0, 0.1, 0.01), expand=c(0,0)) +
    scale_y_continuous(breaks = seq(0.9, 1.0, 0.01), expand=c(0,0)) + 
    scale_fill_gradient2(low="red", high="blue", mid="white", midpoint=0.3)+
    facet_grid(rows = vars(init_sigma), cols = vars(coverage)) 
ggsave("likelihood_raster.png", width=16, height=10)
