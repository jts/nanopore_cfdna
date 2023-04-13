#! /.mounts/labs/simpsonlab/sw/miniconda3/envs/cfdna/bin/Rscript
library(ggplot2)
library(data.table)
library(ggpubr)
args = commandArgs(trailingOnly = TRUE)

#keep = c("lung", "lung_alveolar", "lung_endothelial")
keep = c("HCT116")


# create an empty data frame with four columns
df_main = data.frame(matrix(ncol = 4, nrow = 0))
for (f in args) {
    df = read.table(f, sep='\t', header = TRUE)
    df = df[df$cell_type %in% keep,]
    predicted_lung_proportion = sum(as.numeric(df$proportion))
    # seperate f by "/" and take the last element
    # which is the sample name
    sample_name = unlist(strsplit(f, "/"))[length(unlist(strsplit(f, "/")))]
    # seperate sample name by "." and take the first element
    # which is the true lung percentage
    true_lung_percentage_seed = unlist(strsplit(sample_name, "\\."))[1]
    seed = unlist(strsplit(true_lung_percentage_seed, "_"))[2]
    true_lung_percentage = unlist(strsplit(true_lung_percentage_seed, "_"))[1]
    if (true_lung_percentage == "PBMC") {
        true_lung_percentage = 0
    } else if (true_lung_percentage == "HSAEC" || true_lung_percentage == "HCT116") {
        true_lung_percentage = 100
    } else {
        true_lung_percentage = as.numeric(true_lung_percentage)
    }
    true_lung_proportion = true_lung_percentage/100
    atlas = unlist(strsplit(sample_name, "\\."))[2]
    model = unlist(strsplit(sample_name, "\\."))[3]
    # write a row in df_main
    df_main = rbind(df_main, c(atlas, model, seed, true_lung_proportion, predicted_lung_proportion))
}
colnames(df_main) = c("atlas", "model", "seed", "true_lung_proportion", "predicted_lung_proportion")
# change the type of the columns "true_lung_proportion" and "predicted_lung_proportion" to numeric
df_main$true_lung_proportion = as.numeric(df_main$true_lung_proportion)
df_main$predicted_lung_proportion = as.numeric(df_main$predicted_lung_proportion)
# set values below 0.001 to 0
df_main$predicted_lung_proportion[df_main$predicted_lung_proportion < 0.005] = 0.005
df_main$true_lung_proportion[df_main$true_lung_proportion < 0.005] = 0.005
# remove rows with 0s
#df_main = df_main[df_main$predicted_lung_proportion != 0,]
#df_main = df_main[df_main$true_lung_proportion != 0,]

print(df_main)
# print the rows where model is llse and true_lung_proportion is 0.01
#print(df_main[df_main$model == "llse" & df_main$true_lung_proportion == 0.01,])
#print(df_main[df_main$model == "mmse" & df_main$true_lung_proportion == 0.0,])

labels = c("0.00", "0.01", "0.02", "0.04", "0.08", "0.16", "0.32", "1.00")
breaks = c(0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 1.00)

plot = ggplot(data = df_main, aes(x=true_lung_proportion, y=predicted_lung_proportion, color=model)) +
    #geom_smooth(method = "lm", aes(fill=model) ) + 
    #geom_boxplot(width=.5) +
    stat_summary(fun.data = "mean_cl_boot", fun.args = list(conf.int = 0.9), geom = "errorbar", width=0.1, size=1) +
    stat_summary(fun = "mean", geom = "point", size=3, fill = 'black') +
    stat_cor(method = "pearson") +
    labs(title = "Methylation Deconvolution of Simulated Monocyte/Lung Mixture",
         x = "True Lung Proportion",
         color = "Model",
         y = "Predicted Lung Proportion") +
    theme(axis.title.y = element_text(angle=90, size=6),
          plot.title = element_text(size = 10) ) +
    scale_x_continuous(trans=scales::log_trans(base = 2), breaks = breaks, labels = labels, limits = c(0.005, 1.00)) +
    scale_y_continuous(trans=scales::log_trans(base = 2), breaks = breaks, labels = labels, limits = c(0.005, 1.00)) + 
    geom_abline(aes(intercept = 0, slope=1)) +
    # use bw theme
    theme_bw() +
    facet_grid(cols = vars(atlas))
ggsave("scatter.png", width=12, height=6, dpi=1000, units="in")
