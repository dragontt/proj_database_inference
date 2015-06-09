#!/usr/bin/python

"""
Compute Spearman's rank order correlation of netprophet scores and bart+np scores.
"""

import numpy
import matplotlib.pyplot as plt
from scipy.stats import rankdata
from scipy.stats import spearmanr

def main():
    # fn_np = "/Users/KANG/proj_database_inference/resources/np_network_orig/np_combined_model.adjmtr"
    # fn_ba = "/Users/KANG/proj_database_inference/resources/np_network_bart_holstege_350k/holstegeBartNp.adjmtr"
    # fn_ba2 = "/Users/KANG/proj_database_inference/resources/np_network_bart_holstege_100k/holstegeBartNp.adjmtr"
    # fn_ba3 = "/Users/KANG/proj_database_inference/resources/np_network_bart_holstege_50k/holstegeBartNp.adjmtr"

    # corrs1 = compt_corrs(fn_np, fn_ba)
    # corrs2 = compt_corrs(fn_np, fn_ba2)
    # corrs3 = compt_corrs(fn_np, fn_ba3)

    # bins = numpy.linspace(-1, 1, 50)
    # # plt.hist(corrs1, bins, alpha=0.5)
    # labels = ["np vs bart_350k", "np vs bart_100k", "np vs bart_50k"]
    # plt.hist(corrs1, bins, alpha=0.35, label=labels[0])
    # plt.hist(corrs2, bins, alpha=0.35, label=labels[1])
    # plt.hist(corrs3, bins, alpha=0.35, label=labels[2])
    # plt.legend(loc="upper right")
    # plt.xlabel('bins')
    # plt.ylabel('spearman correlation counts')
    # plt.show()

    fn_np = "/Users/KANG/proj_database_inference/temp/huNp.adjmtr"
    fn_ba = "/Users/KANG/proj_database_inference/temp/holstegeBartNp.adjmtr"
    
    mtr_np = numpy.loadtxt(fn_np)
    mtr_ba = numpy.loadtxt(fn_ba)
    corrs1 = compt_corrs(mtr_np, mtr_ba, 'exclude_zero')

    fn_np = "/Users/KANG/proj_database_inference/temp/huNp.tsv"
    fn_ba = "/Users/KANG/proj_database_inference/temp/holstegeBartNp.tsv"

    mtr_np = parse_tsv(fn_np)
    mtr_ba = parse_tsv(fn_ba)
    corrs2 = compt_corrs(mtr_np, mtr_ba, 'include_zero')

    bins = numpy.linspace(-1, 1, 50)
    labels = ["350k rank cutoff", "no rank cutoff"]
    plt.hist(corrs1, bins=bins, alpha=0.5, label=labels[0])
    plt.hist(corrs2, bins=bins, alpha=0.5, label=labels[1])
    plt.legend(loc="upper right")
    plt.xlabel('rank order correlation')
    plt.ylabel('counts')
    plt.title("huNp vs holstegeBartNp")
    plt.show()

def compt_corrs(mtr_np, mtr_ba, method):
    corrs = numpy.zeros(len(mtr_np))
    total_intersect = 0

    for i in range(len(mtr_np)):
        if method == 'exclude_zero':
            # find the non-zero indices
            ind_np = numpy.nonzero(mtr_np[i] != 0)
            ind_ba = numpy.nonzero(mtr_ba[i] != 0)
            if len(ind_np[0]) <= 1 or len(ind_ba[0]) <= 1:
                # if the tf binds to no target
                corrs[i] = -1
            else:
                # clean up the index data structure
                ind_np = numpy.squeeze(numpy.asarray(ind_np))
                ind_ba = numpy.squeeze(numpy.asarray(ind_ba))
                if ind_np.shape == (0,) or ind_ba.shape == (0,):
                    corrs[i] = -1
                else:
                    # find the union of 2 arrays
                    ind_new = numpy.union1d(ind_np, ind_ba) 
                    total_intersect += len(numpy.intersect1d(ind_np, ind_ba))
                    arr_np = mtr_np[i][ind_new]
                    arr_ba = mtr_ba[i][ind_new]
                    # rank arrays and compute rank order correlation
                    rank_np = rankdata(arr_np)
                    rank_ba = rankdata(arr_ba)
                    corrs[i] = spearmanr(rank_np, rank_ba)[0]

        elif method == 'include_zero':
            # rank arrays and compute rank order correlation
            total_intersect += len(mtr_np[i])
            rank_np = rankdata(mtr_np[i])
            rank_ba = rankdata(mtr_ba[i])
            corrs[i] = spearmanr(rank_np, rank_ba)[0]

    print "Total intersect tf-target edges " + str(total_intersect)
    return corrs

def parse_tsv(fn):
    lines = open(fn, "r").readlines()

    targets = lines[0].split()
    tfs = [None] * (len(lines)-1)
    
    adjmtr = numpy.ndarray([len(tfs), len(targets)])

    for i in range(1,len(lines)):
        temp = lines[i].split()
        tfs[i-1] = temp[0]
        for j in range(1,len(temp)):
            adjmtr[i-1, j-1] = float(temp[j])

    return adjmtr

if __name__ == "__main__":
    main()
