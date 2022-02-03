#! /usr/bin/env python

"""
Generate cpg values for all cpg's in the reference_atlas.csv from nanopolish
methylation calls
"""

import sys
import os
import argparse
import pandas as pd
import numpy as np
from time import perf_counter
from tqdm import tqdm
import re

script_dir = os.path.realpath(os.path.dirname(__file__))

CPG_LOCI = os.path.join(script_dir, 'cpg_loci.csv')

def slice_cpg(df, c, start, end):
    """
    get cpg's on chromosome c from start to end
    """
    return df.loc[c].loc[start:end].acc

def fill_NN(df_refmod, df_cpgloci):
    """
    do an outer join on df_refmod and df_cpgloci
    fill NA values
    """

    df = df_cpgloci.join(df_refmod, how='outer')

    df = df['freq']

    result = {}

    for c in {c for c,_ in df.index if '_' not in c}:
        temp = df[c].interpolate('nearest')
        df[c] = [i for i in temp]


    return df

def join_df(df_cpgloci, refmods, outfile,
            fill = True,
            offset=0,
            gDNA = False):
    """
    perform join operation on all reference modification files
    """

    df = df_cpgloci

    columns = ['chromosome', 'position', 'modified_frequency']
    column_map = {'chromosome':'chr',
                'position': 'pos',
                'modified_frequency': 'freq'}
    if gDNA:
        columns = ['chromosome', 'start', 'methylated_frequency']
        column_map = {'chromosome':'chr',
                'start': 'pos',
                'methylated_frequency': 'freq'}
    if fill:
        f = open(f'{outfile}.log', 'w')
        f.write('sample_name\tmissing_cpgs\tfilled_cpgs\n')

    for refmod in refmods:
        sample_name = re.search('C\d*.*\dC', refmod)[0][1:-1]
        if gDNA: sample_name = refmod.split('/')[1]
        print(f'Loading {sample_name}...')
        try:
            df_refmod = pd.read_csv(refmod, sep='\t')#, usecols=columns)
        except pd.errors.EmptyDataError:
            continue
        print(df_refmod.head())
        df_refmod.rename(columns=column_map,
                            inplace=True)
        df_refmod['pos'] = df_refmod['pos'] + offset
        df_refmod = df_refmod.set_index(['chr', 'pos'])
        if fill:
            na_count = df.join(df_refmod, how='left').freq.isna().sum()
            print(f'{na_count} NA CpGs found.\
                    Filling NA by Nearest Neighbour...')
            df_refmod = fill_NN(df_refmod, df_cpgloci)
        else:
            df_refmod = df_cpgloci.join(df_refmod, how='inner')['freq']

        print('Joining dataframes...')
        df = df.join(df_refmod, how='left')
        if fill:
            na_filled = na_count - df.freq.isna().sum()
            print(f'{na_filled} CpGs filled')
            f.write(f'{sample_name}\t{na_count}\t{na_filled}\n')

        df.rename(columns={'freq':sample_name},
                            inplace=True)

    if fill: f.close()

    return df.drop_duplicates()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('refmods', nargs='+',
                        help='reference_modifications.tsv file')
    parser.add_argument('-o', help='output file')
    parser.add_argument('--gDNA', action='store_true')
    parser.add_argument('--fill', action='store_true')
    parser.add_argument('--offset', type=int, default=0,
                        help= '+1, -1 offset assistance')
    args = parser.parse_args()

    # Get location of cpg sites
    print('Loading CPG LOCI...')
    df_cpgloci = pd.read_csv(CPG_LOCI)
    df_cpgloci = df_cpgloci.set_index(['chr', 'pos'])
    df_cpgloci = df_cpgloci.sort_index().drop_duplicates()

    # Get cpg modification frequency
    count = perf_counter()
    df_result = join_df(df_cpgloci, args.refmods, args.o,
                        offset = args.offset,
                        gDNA = args.gDNA,
                        fill = args.fill)
    print(f'get_cpgfreq executed in {round(perf_counter()-count, 5)} seconds')
    print(df_result.describe())

    if args.o:
        df_result.to_csv(args.o, index=False)
    else:
        print(df_result)



if __name__ == '__main__':
    main()
