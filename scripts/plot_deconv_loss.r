#! /.mounts/labs/simpsonlab/sw/miniconda3/envs/cfdna/bin/Rscript
library(ggplot2)
library(data.table)
args = commandArgs(trailingOnly = TRUE)

df = read.table(args[1], sep='\t', header=TRUE)
df <- df[df$true_lung_proportion == 0.3, ]
atlas = df$atlas[1]

#p11_vec <- unique(df$p01)
#for (p11 in p11_vec) {
    #plot = ggplot(data = df[(df$p01==p11),], aes(x=p11, y=lung_proportion, color=model)) +
        #geom_smooth(method = "lm") + 
        #geom_point() +
        #labs(title = "Methylation Deconvolution of Simulated Lung/Monocyte Mixture",
             #x = "p11",
             #y = "Predicted Lung Proportion") +
        #ylim(0.0,0.7) +
        #geom_hline(aes(yintercept = true_lung_proportion)) +
        #facet_wrap(vars(coverage))
    #ggsave(paste(p11, "deconvolution_loss.png", sep="-"), width=14, height=10)
#}


plot = ggplot(data = df[(df$p11==1.0),], aes(x=p01, y=lung_proportion, color=model)) +
    geom_smooth(method = "lm") + 
    geom_point() +
    labs(title = "Methylation deconvolution of 30%lung/70%monocyte in-silico mixture",
         subtitle = atlas,
         x = "Nanopore Miscall Rate (p01)",
         y = "Predicted Lung Proportion") +
    ylim(0.0,0.5) +
    geom_hline(aes(yintercept = true_lung_proportion)) +
    facet_wrap(vars(coverage))
ggsave("p01-deconv_loss.png", width=14, height=10)
plot = ggplot(data = df[(df$p01==0.0),], aes(x=p11, y=lung_proportion, color=model)) +
    geom_smooth(method = "lm") + 
    geom_point() +
    labs(title = "Methylation deconvolution of 30%lung/70%monocyte in-silico mixture",
         x = "Nanopore Miscall Rate (p11)",
         y = "Predicted Lung Proportion") +
    ylim(0.0,0.5) +
    geom_hline(aes(yintercept = true_lung_proportion)) +
    facet_wrap(vars(coverage))
ggsave("p11-deconv_loss.png", width=14, height=10)
