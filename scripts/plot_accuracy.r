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
    exp_name = tail(strsplit(sample_name, "_")[[1]], n=1)
    sample_name = sub(glue('_{exp_name}'), '', sample_name)
    sample_name
}

get_exp <- function(s){
    exp_name = tail(strsplit(s, "/")[[1]], n=1)
    exp_name = sub('.cpgs_deconv_output.csv', '', exp_name)
    exp_name = tail(strsplit(exp_name, "_")[[1]], n=1)
    exp_name
}

correlation <- function(args) {
    df = data.frame()
    for (arg in args){
        sample_name = get_sn(arg)
        df2 = read.table(arg, row.names=1, sep=',')
        exp_name = get_exp(arg)

        # Order by descending coverage
        tdf = transpose(df2)
        df2 = transpose(tdf[order(-tdf$V1),])
        colnames(df2) <- df2[1,]
        df2 <- df2[-1,]
        df2 = cor(df2)
        print(df2)
        df2 <- data.frame(coverage = rownames(df2),
                  acc = df2[,1],
                  exp = exp_name,
                  sample = sample_name)
        df2[,1:2] = (apply(df2[,1:2], 2, remove_X))
        df = data.frame(rbind(df, df2))
    }
}

stderr_one <- function(expected, observed){
    sum = 0
    for (i in 1:length(expected)) {
        sum = sum + (observed[i]-expected[i])**2
    }
    sqrt((sum)/length(expected))
}

stderr <-function(args) {
    df = data.frame()
    for (arg in args) {
        sample_name = get_sn(arg)
        exp_name = get_exp(arg)
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
            df_chi[j] = stderr_one(expected, df2[,j])
        }
        print(df_chi)

        df2 <- data.frame(coverage = colnames(df2),
                  acc = df_chi,
                  exp = exp_name,
                  sample = sample_name)
        df2[,1:2] = (apply(df2[,1:2], 2, remove_X))
        df = data.frame(rbind(df, df2))
    }
    df
}

for (t in 1:ceiling(length(args)/10) ){

    samples = args[((t-1)*10):(min(length(args), t*10))]

    print(samples)
    df = stderr(samples)
    print(df)
    plot = ggplot(df, aes(coverage, acc, fill=exp)) +
        geom_col(position="dodge") + 
        labs(title = "Standard Deviation of deconvolution output vector with original sample and downsampled coverage",
             x = "Coverage",
             y = "Average Standard Deviation to Original Sample") +
        xlim(0,12)+
        facet_grid(rows=vars(sample))

    write.table(df, file=glue("coverageVstd.{t}.tsv"),
                sep='\t')

    ggsave(glue("coverageVstd.{t}.png"), width=14, height=10)



    #df = correlation(samples)
    #print(df)
    #plot = ggplot(df, aes(coverage, acc, fill=sample)) +
        #geom_col(position="dodge")

        #labs(title = "Correlation of deconvolution output vector with original sample and downsampled coverage",
             #x = "Coverage",
             #y = "Average Correlation to Original Sample") +
        #facet_grid(rows=vars(sample))


    #write.table(df, file=glue("coverageVcorr.{t}.tsv"),
                #sep='\t')

    #ggsave(glue("coverageVcorr.{t}.png"), width=14, height=10)
}
