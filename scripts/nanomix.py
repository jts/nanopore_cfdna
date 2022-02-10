#! /usr/bin/env python

import argparse
import numpy
import sys
import csv
import os

from scipy.stats import binom
from scipy.optimize import minimize, nnls, Bounds

script_dir = os.path.dirname(__file__)
ATLAS = os.path.join(script_dir, '..', 'atlases', 'meth_atlas.csv')

class ReferenceAtlas:
    def __init__(self, filename):
        self.cpg_ids = list()
        self.v = dict()
        self.K = None

        with open(filename) as csvfile:
            reader = csv.DictReader(csvfile)
            data = list()
            for row in reader:
                cpg_id = row['CpGs']
                self.cpg_ids.append(cpg_id)
                cell_types = list(row.keys())[1:]
                self.K = len(cell_types)
                r = list()
                for k in cell_types:
                    if k not in self.v:
                        self.v[k] = list()
                    self.v[k].append(float(row[k]))
                    r.append(float(row[k]))
                data.append(r)

        self.A = numpy.array(data).reshape((self.get_num_cpgs(), self.get_num_cell_types()))

    def get_x(self, sigma):
        x = numpy.matmul(self.A, sigma)
        return x

    def get_num_cpgs(self):
        return len(self.cpg_ids)

    def get_cell_types(self):
        return list(self.v.keys())

    def get_num_cell_types(self):
        return len(self.v.keys())

class Sample:
    def __init__(self, name, x_hat, m, t):
        self.name = name
        self.x_hat = x_hat
        self.m = m
        self.t = t

# Binomial model with sequencing errors, when epsilon = 0
# this is the same as the perfect data model
def log_likelihood_sequencing_with_errors(atlas, sigma, sample, epsilon):
    sigma_t = sigma.reshape( (atlas.K, 1) )

    # the solver we use can try values that are outside
    # the constraints we impose, we need to clip here to prevent
    # things from blowing up
    x = numpy.clip(numpy.ravel(atlas.get_x(sigma_t)), 0, 1.0)
    p = x * (1 - epsilon) + (1 - x) * epsilon
    b =  binom.logpmf(sample.m, sample.t, p)

    #print("SigmaT", sigma_t)
    #print("SigmaSum", numpy.sum(sigma_t))
    #print("m", sample.m)
    #print("t", sample.t)
    #print("x", x)
    #print("B", b)
    #print("Sum", numpy.sum(b))
    return numpy.sum(b)

def eq_constraint(x):
    return 1 - numpy.sum(x)

#
# Model wrappers
#
def fit_llse(atlas, sample, epsilon):
    sigma_0 = numpy.array([ [ 1.0 / atlas.K ] * atlas.K ])
    f = lambda x: -1 * log_likelihood_sequencing_with_errors(atlas, x, sample, epsilon)

    bnds = [ (0.0, 1.0) ] * atlas.K
    cons = ({'type': 'eq', 'fun': eq_constraint})
    res = minimize(f, sigma_0, method='SLSQP', options={'maxiter': 10, 'disp':False}, bounds=bnds, constraints=cons)
    return res.x

def fit_nnls(atlas, sample):

    # add sum=1 constraint
    t = numpy.array([1.0] * atlas.K).reshape( (1, K) )
    A = numpy.append(atlas.A, t, axis=0)
    b = numpy.append(sample.x_hat, [1.0], axis=0)
    res = nnls(A, b)
    return res[0]

def fit_nnls_constrained(atlas, sample):
    sigma_0 = numpy.array([ [ 1.0 / atlas.K ] * atlas.K ])
    f = lambda x: numpy.linalg.norm(atlas.A.dot(x) - sample.x_hat)
    bnds = [ (0.0, 1.0) ] * atlas.K
    cons = ({'type': 'eq', 'fun': eq_constraint})
    res = minimize(f, sigma_0, method='SLSQP', options={'maxiter': 10, 'disp':False}, bounds=bnds, constraints=cons)
    return res.x

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--atlas', type=str, default=ATLAS)
    parser.add_argument('--input', required=True, type=str)
    parser.add_argument('--model', default='llse', type=str, help='deconvolution model options: [nnml, llse]')
    args = parser.parse_args()
    atlas = ReferenceAtlas(args.atlas)

    coverage = 10
    epsilon = 0.05

    data = dict()
    with open(args.input) as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            sample_names = list(row.keys())[1:]
            for s in sample_names:
                if s not in data:
                    data[s] = dict()

                if row[s] != "":
                    data[s][row['acc']] = float(row[s])
                else:
                    data[s][row['acc']] = 0.0

    # convert to Samples and run
    Y = []
    for sn in sample_names:
        # get cpg frequenices in order of atlas
        xhat = [data[sn].get(cpg, 0.0) for cpg in atlas.cpg_ids]

        # fill in default coverage values (for now)
        t = [ coverage ] * len(xhat)
        m = [ int(t[i] * xhat[i]) for i in range(0, len(t)) ]

        s = Sample(sn, numpy.array(xhat), numpy.array(m), numpy.array(t))
        if args.model == 'nnls':
            Y.append(fit_nnls_constrained(atlas, s))
        else:
            Y.append(fit_llse(atlas, s, epsilon))
    # output
    print(f"ct,{','.join([sn for sn in sample_names])}")
    for (idx, ct) in enumerate(atlas.get_cell_types()):
        print(f"{ct},{','.join([str(y[idx]) for y in Y])}")

if __name__ == "__main__":
    main()
