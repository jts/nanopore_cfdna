import csv
import pandas as pd
import pyranges as pr

def write_bed(in_f, out_f):
    strand_dict = {'F':'+', 'R':'-'}
    with open(in_f, 'r') as csvfile:
        reader = csv.reader(csvfile)
        for probe, build, chr, start, strand in reader:
            if build != '37':
                print('probe of not 37 found')
                continue
            print('chr{}\t{}\t{}\t{}\t{}\t{}'.format(chr, start, int(start)+1, probe, 0, strand_dict[strand]))

def make_dic(f):
    d = {}
    with open(f, 'r') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        head = True
        for chr, start, end, name, score, strand in reader:
            if head:
                head = False
                continue
            d[name] = (chr, start, end)

    return d

def f(row, cpgToLocus):
    try:
        row.Chromosome = cpgToLocus[row.name][0]
        row.Start = int(cpgToLocus[row.name][1])
        return row
    except KeyError:
        return None

def main():
    cpgToLocus = make_dic('full_cpg_loci_hg38.bed')

    atlas = pd.read_csv('CpGs.100bp-block.1000.csv', index_col=0)
    atlas['Chromosome'] = None
    atlas['Start'] = None
    atlas = atlas.apply(f, cpgToLocus=cpgToLocus, axis=1)
    atlas = atlas.dropna()
    atlas['Start'] = atlas['Start'].astype(int)
    atlas['End'] = atlas['Start'] + 1
    gr = pr.PyRanges(atlas).sort()
    breakpoint()
    atlas = gr.df.set_index(['Chromosome', 'Start', 'End'])
    atlas.to_csv('BermanAtlas.tsv', sep='\t', index=True)

if __name__ == '__main__':
    main()
