#!/usr/bin/python

"""
Compute Spearman's correlation of netprophet scores and scertf alignment scores.
"""

import sys
import os
import argparse
import numpy
import matplotlib.pyplot as plt
import compt_pw_corr
import math
import glob
import os.path

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Compute rank order correlations of inferred and known motifs")
    parser.add_argument('-d', '--dict_conv', dest='dict_conv', type=str, default='resources/np_scertf_names.txt')
    parser.add_argument('-s', '--dir_st', dest='dir_st', type=str)
    parser.add_argument('-n0', '--dir_np0', dest='dir_np0', type=str)
    parser.add_argument('-n1', '--dir_np1', dest='dir_np1', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)

    parsed.dir_st = check_dir(parsed.dir_st)
    parsed.dir_np0 = check_dir(parsed.dir_np0)
    parsed.dir_np1 = check_dir(parsed.dir_np1)

    # get both gene and systematic names of the given tf name
    fire_names = []
    fns = glob.glob(parsed.dir_np0 + "*")
    for fn in fns:
        fire_names.append(os.path.basename(fn))

    names_dict = {}
    lines = open(parsed.dict_conv, 'r').readlines()
    for line in lines:
        if line.split()[1] in fire_names:
            names_dict[line.split()[0]] = line.split()[1]

    # compute correlations 
    corrs_out0 = compt_corrs(parsed.dir_st, parsed.dir_np0, names_dict)
    corrs_out1 = compt_corrs(parsed.dir_st, parsed.dir_np1, names_dict)

    # show bar plot
    bins = numpy.linspace(-0.25, 0.25, 50)
    labels = ["fire_pwm_7mers", "fire_pwm_8mers"]
    plt.hist(corrs_out0, bins, alpha=0.5, label=labels[0])
    plt.hist(corrs_out1, bins, alpha=0.5, label=labels[1])
    plt.legend(loc="upper right")
    plt.xlabel('Bins')
    plt.ylabel('spearman correlation counts')
    plt.show()

def check_dir(dirname):
    if not dirname.endswith('/'):
        dirname += '/'
    return dirname

def parse_file(filename):
    rank = {}
    lines = open(filename, 'r').readlines()
    for line in lines:
        linesplit = line.split()
        if len(linesplit) != 0:
            rank[linesplit[0]] = linesplit[1]
    return rank

def compt_corrs(dir_st, dir_np, names_dict):
    corrs_out = []
    for gene_name, sys_name in names_dict.iteritems():
        # parse scertf and netprophet rankings
        st_rank = parse_file(dir_st + gene_name)
        np_rank = parse_file(dir_np + sys_name)
        temp_corr = compt_pw_corr.compute_corr(st_rank, np_rank)
        if not math.isnan(temp_corr):
            corrs_out.append(temp_corr)
    return corrs_out

if __name__ == "__main__":
    main(sys.argv)
