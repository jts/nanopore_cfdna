#! /.mounts/labs/simpsonlab/sw/miniconda3/envs/cfdna/bin/Rscript
library(ggplot2)
library(glue)
library(data.table)
args = commandArgs(trailingOnly = TRUE)
remove_X <- function(s) {
	as.numeric(sub('X', '', s))
}

get_sn <- function(s){
    sample_name = tail(strsplit(s, "/")[[1]], n=1)
    sample_name = sub('.cpgs_deconv_output.csv', '', sample_name)
    sample_name = tail(strsplot(sample_name, "_")[[1]], n=1)
    sample_name
}

correlation <- function(args) {
    df = data.frame()
    for (arg in args){
        sample_name = get_sn(arg)
        df2 = read.table(arg, row.names=1, sep=',')

        # Order by descending coverage
        tdf = transpose(df2)
        df2 = transpose(tdf[order(-tdf$V1),])
        colnames(df2) <- df2[1,]
        df2 <- df2[-1,]
        df2 = cor(df2)
        print(df2)
        df2 <- data.frame(coverage = rownames(df2),
                  acc = df2[,1],
                  sample = sample_name)
        df2[,1:2] = (apply(df2[,1:2], 2, remove_X))
        df = data.frame(rbind(df, df2))
    }
}

chi2 <- function(expected, observed){
    sum = 0
    for (i in 1:length(expected)) {
        sum = sum + (observed[i]-expected[i])**2
    }
    print(sqrt(sum/length(expected)))
    sum
}

chi_squared <-function(args) {
    df = data.frame()
    for (arg in args) {
        sample_name = get_sn(arg)
        df2 = read.table(arg, row.names=1, sep=',')

        # Order by descending coverage
        tdf = transpose(df2)
        df2 = transpose(tdf[order(-tdf$V1),])
        colnames(df2) <- df2[1,]
        df2 <- df2[-1,]

        expected = df2[,1]
        df_chi = c()
        print(df2)
        for (j in 1:length(df2)){
            df_chi[j] = chi2(expected, df2[,j])
        }
        print(df_chi)

        df2 <- data.frame(coverage = colnames(df2),
                  acc = df_chi,
                  sample = sample_name)
        df2[,1:2] = (apply(df2[,1:2], 2, remove_X))
        df = data.frame(rbind(df, df2))
    }
    df
}

df = chi_squared(args)
print(df)
sample_name = get_sn(args[1])
plot = ggplot(df, aes(x=coverage, y=acc, color=sample)) +
	geom_point() +
	geom_smooth(se=FALSE) + 

	labs(title = "Standard Deviation of deconvolution output vector with original sample and downsampled coverage",
	     x = "Coverage",
	     y = "Average Error to Original Sample")

write.table(df, file=glue("{sample_name}.coverageVerror.tsv"),
            sep='\t')

ggsave(glue("{sample_name}.coverageVerror.png"), width=14, height=10)

