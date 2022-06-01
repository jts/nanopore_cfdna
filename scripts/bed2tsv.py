#! /usr/bin/env python
import argparse
import pandas as pd
import pyranges as pr

def decrease_index(df):
    df.Start -= 1
    df.End -= 1
    return df

def increase_index(df):
    df.Start += 1
    df.End += 1
    return df

def reindex_n(df):
    df.Start = df.Start_n
    df.End = df.End_n
    df.Strand = '+'

    return df

def combine(df_sample, df_atlas):
    """
    Get the rows of sites present in the atlas
    reconcile 0-based indexing on the + and - strands

    :param df_atlas: pd.DataFrame of CpG sites with the 'start' index indicating the C
    :param df_sample: pd.DataFrame of CpG sites with the 'end' index indicating the C on the + strand and the 'start' index indicating the C on the - strand
    """
    # compute total calls and modified calls method from Berman et al. 2022
    df_sample.total_reads = round(df_sample.Score * df_sample.BlockCount/1000)
    df_sample.methylated_reads = df_sample.total_reads*df_sample.BlockSizes /100

    # Intersect
    # Assume atlas loci are from the positive strand
    df_atlas.Strand = '+'
    df_atlas = df_atlas.drop_duplicate_positions()
    negative_strand = df_sample.intersect(df_atlas, strandedness='opposite')
    positive_strand = df_sample.apply(increase_index).intersect(df_atlas, strandedness='same')

    # Join and compute total calls, modified calls and modification freq
    df_tsv = positive_strand.join(negative_strand, how='outer', strandedness='opposite', suffix='_n')
    df_tsv = pr.PyRanges(df_tsv.df.replace(-1, 0))
    df_tsv.total_calls = df_tsv.total_reads + df_tsv.total_reads_n
    df_tsv.modified_calls = df_tsv.methylated_reads + df_tsv.methylated_reads_n
    df_tsv.modification_frequency = df_tsv.modified_calls/df_tsv.total_calls

    # Reindex and select columns
    df_tsv = pr.concat([df_tsv[df_tsv.Start > 0], df_tsv[df_tsv.Start <= 0].apply(reindex_n)]).sort(["Chromosome", "Start", "End"])
    df_tsv = df_tsv[['total_calls', 'modified_calls', 'modification_frequency']]

    return df_tsv

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('-atlas', help='bed file of atlas regions')
    parser.add_argument('-o')
    args = parser.parse_args()
    sample_file = args.input
    atlas = args.atlas
    df_sample, df_atlas = pr.read_bed(sample_file), pr.read_bed(atlas)

    df_tsv = combine(df_sample, df_atlas)
    print(df_tsv)

    cols = {'Chromosome':'chr', "Start":'start', 'End':'end'}
    df = df_tsv.df.rename(columns=cols)
    df.to_csv(args.o,
            sep='\t',
            index=False,
            header=True)

if __name__ == '__main__':
    main()
