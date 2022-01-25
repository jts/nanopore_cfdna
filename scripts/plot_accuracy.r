#! /.mounts/labs/simpsonlab/sw/miniconda3/envs/cfdna/bin/Rscript
library(ggplot2)
library(data.table)
args = commandArgs(trailingOnly = TRUE)

df = read.csv(args[1], row.names=1)

print(head(df))

remove_X <- function(s) {
	as.numeric(sub('X', '', s))
}
#sample_name = sub('.cpgs_deconv_output.csv', '', args[1])
#df = cor(df)
#df <- data.frame(coverage = rownames(df),
		 #corr = df[,1],
		 #sample = sample_name)
#df = data.frame(apply(df, 2, remove_X))

df = data.frame()

for (arg in args){
    sample_name = tail(strsplit(arg, "/")[[1]], n=1)
	sample_name = sub('.cpgs_deconv_output.csv', '', sample_name)
	df2 = read.table(arg, row.names=1, sep=',')
    # Order by descending coverage
    tdf = transpose(df2)
    df2 = transpose(tdf[order(-tdf$V1),])
    colnames(df2) <- df2[1,]
    df2 <- df2[-1,]
	df2 = cor(df2)
    print(df2)
	df2 <- data.frame(coverage = rownames(df2),
			  corr = df2[,1],
			  sample = sample_name)
	df2[,1:2] = (apply(df2[,1:2], 2, remove_X))
	df = data.frame(rbind(df, df2))
}

print(df)

plot = ggplot(df, aes(x=coverage, y=corr, color=sample)) +
	geom_point() +
	geom_smooth(level=0.3, formula='y ~ log(x)') + 

	labs(title = "Pearson correlation of deconvolution output vector with original sample and downsampled coverage",
	     x = "Coverage",
	     y = "Correlation to Original Sample")

ggsave("coverageVaccuracy.png", width=10, height=8)

