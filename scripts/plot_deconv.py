#! /usr/bin/env python
import numpy as np
import pandas as pd
from scipy import optimize
import argparse
import os.path as op
import sys
from multiprocessing import Pool
import math
import matplotlib.pylab as plt
import matplotlib.cm
import matplotlib.colors

ATLAS_FILE = '/u/jbroadbent/jbroadbent/code/cfdna/meth_atlas/reference_atlas.csv'
OUT_PATH = 'deconv_plot.png'

# Plotting parameters:
NR_CHRS_XTICKS = 30         # number of characters to be printed of the xticks
FIG_SIZE = (15, 7)          # figure size
COLOR_MAP = 'tab20'         # color map. See https://matplotlib.org/users/colormaps.html
#COLOR_MAP = 'Vega10'
# tissues with less than OTHERS_THRESH contribution will be clustered to 'other' (black):
OTHERS_THRESH = 0.05


####################################
#       Plotting methods           #
####################################
def hide_non_blood(df):
    blood_cells = ['tcell', 'erythroblast', 'nkcell', 'bcell', 'progenitor', 'hsc', 'monocyte', 'macrophage', 'eosinophil', 'neutrophil']
    selection = [name not in blood_cells for name in df.index]
    others = df[selection].sum()
    df = df.drop(df.index[tuple([selection])])
    df = df.append(others.rename('other'))
    return df



def hide_small_tissues(df):
    """
    tissues with very small contribution are grouped to the 'other' category.
    :return: The DataFrame with the new category ('other'),
             where the low-contribution tissues are set to 0.
    """
    others = df[df < OTHERS_THRESH].sum()
    df[df < OTHERS_THRESH] = 0.0
    df = df.append(others.rename('other'))
    return df


def gen_bars_colors_hatches(nr_tissues):
    """
    Generate combinations of colors and hatches for the tissues bars
    Every tissue will get a tuple of (color, hatch)
    the last tuple is for the 'other' category, and is always black with no hatch.
    :return: a list of tuples, with length == nr_tissues
    """
    matplotlib.rcParams['hatch.linewidth'] = 0.3
    hatches = [None, 'xxx', '...', 'O', '++'][:nr_tissues // 7]

    nr_colors = int(math.ceil(nr_tissues / len(hatches)) + 1)
    print(nr_colors, hatches)

    # generate bars colors:
    cmap = matplotlib.cm.get_cmap(COLOR_MAP)
    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=float(nr_colors))
    colors = [cmap(norm(k)) for k in range(nr_colors)]
    # for c in range(nr_colors): print(norm(c), cmap(norm(c)), cmap(c))

    def get_i_bar_tuple(i):
        color_ind = i % nr_colors
        hatch_ind = int(i // math.ceil(nr_tissues / len(hatches)))
        return colors[color_ind], hatches[hatch_ind]

    colors_hatches_list = [get_i_bar_tuple(i) for i in range(nr_tissues-1)]
    return colors_hatches_list + [((0, 0, 0, 1), None)]

def g(x):
    if x == 'ct':
        return 0
    else:
        return 1/float(x)


def f(x):
    if x == 'ct':
        return 0
    elif 'relapse' in x:
        x_new = "".join(x.split('_')[:2] + ['_'])
        return sum([ord(x0) for x0 in x_new]) + 0.0001
    elif 'BCR' in x:
        return sum([ord(x0) for x0 in x]) + 500
    elif 'add' in x:
        return 1
    else:
        return sum([ord(x0) for x0 in x])

def plot_res(df, outpath, show=False):

    # df = hide_small_tissues(df)
    # df = hide_non_blood(df)
    nr_tissues, nr_samples = df.shape
    print(df.shape)

    # generate bars colors and hatches:
    colors_hatches = gen_bars_colors_hatches(nr_tissues)

    plt.figure(figsize=FIG_SIZE)
    r = [i for i in range(nr_samples)]
    bottom = np.zeros(nr_samples)
    for c in colors_hatches: print(c)
    for i in range(nr_tissues):
        plt.bar(r, list(df.iloc[i, :]),
                edgecolor='white',
                width=0.85,
                label=df.index[i],
                bottom=bottom,
                color=colors_hatches[i][0],
                hatch=colors_hatches[i][1])
        bottom += np.array(df.iloc[i, :])

    # Custom x axis
    plt.xticks(r, [w[:NR_CHRS_XTICKS] for w in df.columns], rotation='vertical', fontsize=9)
    plt.xlabel("sample")
    plt.xlim(-.6, nr_samples - .4)

    plt.ylabel("Mixture Proportion")
    # Add a legend and a title
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
    plt.title('Deconvolution Results')

    # adjust layout, save and show
    plt.tight_layout(rect=[0, 0, .83, 1])
    plt.savefig(outpath)
    if show:
        plt.show()



####################################
#            main                  #
####################################


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('input')

    parser.add_argument('-name', default=OUT_PATH, help='Output directory')

    args = parser.parse_args()

    df = pd.read_csv(args.input, sep='\t', index_col='ct')
    # df = df.reindex(sorted(df.columns), axis=1)
    # df = df.sort_values('10', ascending=False)

    plot_res(df, args.name)


if __name__ == '__main__':
    main()
