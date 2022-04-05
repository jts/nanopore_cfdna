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
    sample_name = sub('.deconv_output.csv', '', sample_name)
    exp_name = tail(strsplit(sample_name, "-")[[1]], n=1)
    sample_name = sub(glue('-{exp_name}'), '', sample_name)
    sample_name
}

get_exp <- function(s){
    exp_name = tail(strsplit(s, "/")[[1]], n=1)
    exp_name = sub('.deconv_output.csv', '', exp_name)
    exp_name = tail(strsplit(exp_name, "-")[[1]], n=1)
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
        df2 <- data.frame(coverage = rownames(df2),
                  acc = df2[,1],
                  exp = exp_name,
                  sample = sample_name)
        df2[,1:2] = (apply(df2[,1:2], 2, remove_X))
        print(df2)
        df = data.frame(rbind(df, df2))
    }
    df
}

stderr_one <- function(expected, observed){
    sum = 0
    for (i in 1:length(expected)) {
        sum = sum + (observed[i]-expected[i])**2
    }
    sqrt((sum))
}

stderr <-function(args) {
    df = data.frame()
    for (arg in args) {
        sample_name = get_sn(arg)
        exp_name = get_exp(arg)
        df2 = read.table(arg, sep='\t')

        # Order by descending coverage
        print(df2)
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

    samples = args[(((t-1)*10)+1):(min(length(args), t*10))]

    

    print(samples)
    df = stderr(samples)
    print(df)
    plot = ggplot(df, aes(coverage, acc, fill=exp)) +
        geom_col(position="dodge") + 
        labs(title = "Methylation Deconvolution Accuracy",
             x = "Coverage",
             y = "Euclidean Distance to Original Sample") +
        xlim(0,10)+
        scale_x_continuous(breaks=c(0.1,0.5,1,2,3,4,5,6,7,8,9,10))+
        ylim(0,1) +
        facet_grid(rows=vars(sample))

    write.table(df, file=glue("coverageVdistance.{t}.tsv"),
                sep='\t')

    ggsave(glue("coverageVdistance.{t}.png"), width=14, height=10)
}

# Plot Aggregate
df_all = stderr(args)
df = aggregate(df_all,
               by = list(df_all$coverage, df_all$exp),
               FUN = mean)
print(df_all)
plot = ggplot(df_all, aes(coverage, acc, color=exp)) +
    geom_point() + 
    geom_smooth(level=0.99)+
    labs(title = "Methylation Deconvolution Accuracy",
         x = "Coverage",
         y = "Euclidean Distance to Original Sample") +
    xlim(0,10) +
    scale_x_continuous(breaks=c(0.1,0.5,1,2,3,4,5,6,7,8,9,10))+
    ylim(0, 1)

write.table(df, file=glue("coverageVdistance.agg.tsv"),
            sep='\t')

ggsave(glue("coverageVdistance.agg.png"), width=14, height=10)
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
