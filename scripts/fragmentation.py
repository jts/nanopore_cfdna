#! /usr/bin/env python

import csv
from collections import defaultdict
import argparse
import numpy as np

chromosomes = {'chr' + str(num) for num in range(1,23)} | {'chrX', 'chrY'}

def intFromString(s):

    return int("".join(filter(str.isdigit, s)))

def roundToBase(x, base):
    """
    Round a number down to the nearest bin size
    """
    return base * (x//base)

def parse_reads(f):
    """
    Get alignment lengths from nanopore read.tsv output
    """

    dtypes = [('chromosome', 'U5'), ('start', int),\
                ('end', int), ('length', float), ('mapq', int)]
    reads = np.loadtxt(f,delimiter='\t',
                        dtype=dtypes,
                        skiprows=1,
                        usecols=(1, 2, 3, 4, 6), #MapQ is col 6
                        )
    return reads[[x in chromosomes for x in reads['chromosome']]]


def get_bins(f):
    """
    Get 5Mb bins from the tsv file in the fragmentation paper
    """



def bin_frag_ratios(reads,
                    bin_size = 5e+6,
                    short_range = (100, 150),
                    long_range = (151, 220)
                    ):
    """
    bin together align_length into 5Mb Bins
    compute the frag_ratios

    :input reads: (numpy.ndarray)
    :param bin_size:

    :return: numpy array of tuples
            (chromosome, start, end, length, frag_ratio, mean read_length, n)

    """
    dtypes = [('chromosome', 'U5'), ('start', int),\
                ('end', int), ('length', float), ('mapq', int)]

    frag_ratios = []
    chromosomes = {x[0] for x in reads}

    reads_dict = {chromosome:[] for chromosome in chromosomes}
    for x in reads:
        reads_dict[x[0]].append(x)

    for chromosome, reads_chr in reads_dict.items():

        # Find max read length
        reads_chr = np.array(reads_chr, dtype=dtypes)
        max_end = max(reads_chr['end'])

        # Put reads into bins
        bin_dict = {Bin:[] \
                for Bin in range(0, int(max_end + bin_size), int(bin_size))}
        for read in reads_chr:
            bin_dict[roundToBase(read['start'], bin_size)].append(read['length'])

        # Compute frag ratio
        for Bin, binned_reads in bin_dict.items():
            n = len(binned_reads)
            if n == 0: continue
            mean_length = np.mean(binned_reads)
            frag_ratio = compute_frag_ratio(binned_reads, short_range, long_range)

            # Append to numpy array
            frag_ratios.append( (chromosome, Bin, Bin + bin_size,\
                                frag_ratio, mean_length, n) )

    dtypes = [('chromosome', 'U5'), ('start', int), ('end', int),\
                ('frag_ratio', float), ('mean_length', float), ('n', int)]

    return np.array(frag_ratios, dtype=dtypes)

def compute_frag_ratio(reads,
                        short_range = (100, 150),
                        long_range = (151, 220)
                        ):
    """
    Given a list of reads, calculate ratio of short to long reads

    param short_range: tuple
    param long_range: tuple
    return ratio: float
    """

    try:
        r = len([x for x in reads if x in range(*short_range)]) /\
            len([x for x in reads if x in range(*long_range)])
    except ZeroDivisionError:
        return np.nan

    return r

def standardize(X):

    from sklearn import preprocessing

    X['frag_ratio'] = preprocessing.scale(X['frag_ratio'])

    return X

def mapq_filter(X, t):
    """
    filter the numpy array for reads with mapq > threshold t
    """

    return X[X['mapq'] >= t]

def coverage_filter(X, t, bin_size):
    """
    filter frag_ratio bins for those with atleast tx coverage
    """
    return X[X['mean_length']*X['n'] > t*bin_size]


def main():

    p = argparse.ArgumentParser()

    p.add_argument("input_file", help="tsv containing nanopore reads")
    p.add_argument('-o', '--output_file', help='filename with .bed extension',
                    default='output.tsv')
    p.add_argument('-s', '--short_range', type=int, nargs='+', default = [100, 151],
                    help='short read length range')
    p.add_argument('-l', '--long_range', type=int, nargs='+', default = [151, 220],
                    help='short read length range')
    p.add_argument('-b', '--bin_size', type=int, default = int(1e+7))
    p.add_argument('-tq', '--mapq_threshold', type=int, default = 40)
    p.add_argument('-tc', '--coverage_threshold', type=int, default = 0)

    args = p.parse_args()


    big_test = '/.mounts/labs/simpsonlab/projects/cfdna/jbroadbent_dev5/210726_AIX_0360_LSK109_run1_results/210726_AIX_0360_LSK109_run1.read_modifications.tsv'
    small_test = '../../projects/cfdna/dev5/analysis/210726_AIX_0360_LSK109_run1.read_testset.tsv'

    # args.input_file = big_test
    # args.output_file= 'junk.tsv'

    # Parse reads from TSV file
    reads = parse_reads(args.input_file)

    reads = mapq_filter(reads, args.mapq_threshold)

    # Compute fragmentation ratios in bins
    X = bin_frag_ratios(reads,
                                    args.bin_size,
                                    args.short_range,
                                    args.long_range)
    # Filter bins for coverage
    # X = coverage_filter(X, args.coverage_threshold, args.bin_size)

    # Standardize across the genome
    X = standardize(X)

    # Sort
    X = np.sort(X, order=['chromosome', 'start'])

    # Save output to tsv
    np.savetxt(args.output_file, X,
                delimiter='\t',
                header='chromosome\tstart\tend\tfragmentation ratio\tmean read length\tnumber of reads',
                fmt='%s')

if __name__ == '__main__':
    main()
