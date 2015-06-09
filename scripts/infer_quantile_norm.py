#!/usr/bin/python

"""
Infer a motif by quantile normalization of two sets of scores, which are the DNA binding domain 
percent identity and the correlation of database motif-target scoring rankings with netprophet 
score rankings respectively.

Methods: 
standard - standard quantile normalization
prenorm - prenormalize the datasets before standard quantile normalization
bidirect - map distribution to each other, average scores of each normalized pair, and average
           the rankings of the two normalized pairs 
"""

import sys
import argparse
import glob
import os.path
import scipy.stats
import numpy
import quantile_norm as qn
import bidirect_norm as bn
import matplotlib.pyplot as plt
import time

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Infer the most informative motif")
    parser.add_argument('-m', '-method', dest='method', type=str)
    parser.add_argument('-a', '-dir_dbd_aaid', dest='dir_dbd_aaid', type=str)
    parser.add_argument('-n', '-dir_corr_np', dest='dir_corr_np', type=str)
    parser.add_argument('-o', '-dir_output', dest='dir_output', type=str)
    parser.add_argument('-d', '-dict_conv', dest='dict_conv', type=str, default='resources/np_scertf_names.txt')
    parser.add_argument('-c', '-dir_cisbp_rank', dest='dir_cisbp_rank', type=str)
    parser.add_argument('-s', '-dir_scertf_rank', dest='dir_scertf_rank', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)
    parsed.dir_dbd_aaid = check_dir(parsed.dir_dbd_aaid)
    parsed.dir_corr_np = check_dir(parsed.dir_corr_np)
    parsed.dir_output = check_dir(parsed.dir_output)
    parsed.dir_cisbp_rank = check_dir(parsed.dir_cisbp_rank)
    parsed.dir_scertf_rank = check_dir(parsed.dir_scertf_rank)

    """ Infer motifs by quantile normalization of the two sets of scores """
    sys.stdout.write("\rProcessing inference ... ")

    # get motifs names from both dbd and netprophet
    names_aaid = get_names(parsed.dir_dbd_aaid)
    names_np = get_names(parsed.dir_corr_np)
    names_common = list(set(names_aaid) & set(names_np))
    data_out = [None] * len(names_common)

    for i in range(len(names_common)):
        query = names_common[i]
        tic = time.clock()
        # parse aaids and np correlations of the query
        [motifs_dbd_all, aaids_dbd_all] = parse_scores(parsed.dir_dbd_aaid + query + ".aaid")            
        [motifs_np_all, corrs_np_all] = parse_scores(parsed.dir_corr_np + query)

        # sort the scores with the same motif order
        dict_dbd_all = {}
        for j in range(len(motifs_dbd_all)):
            temp_motif = motifs_dbd_all[j].split(':')[0]
            if temp_motif not in dict_dbd_all:
                dict_dbd_all[temp_motif] = aaids_dbd_all[j]
            else:
                if aaids_dbd_all[j] > dict_dbd_all[temp_motif]:
                    dict_dbd_all[temp_motif] = aaids_dbd_all[j]
        aaids_dbd_all = [None] * len(motifs_np_all)
        for j in range(len(motifs_np_all)):
            aaids_dbd_all[j] = dict_dbd_all[motifs_np_all[j]] if motifs_np_all[j] in dict_dbd_all else float(0)

        scores_orig = numpy.array([aaids_dbd_all, corrs_np_all])

        # standard quantile normalization of the socres
        if parsed.method == "standard":
            scores_norm = qn.quantile_norm(scores_orig)
            scores_avg = (scores_norm[0,:] + scores_norm[1,:]) / 2
            rank_scores_avg = scipy.stats.rankdata(scores_avg)

        # quantile normalize the socres with prenormalization
        elif parsed.method == "prenorm":  
            for j in range(2):
                scores_orig[j] = scores_orig[j] / float(numpy.max(scores_orig[j]) - float(numpy.min(scores_orig[j])))
            scores_norm = qn.quantile_norm(scores_orig)
            scores_avg = (scores_norm[0,:] + scores_norm[1,:]) / 2
            rank_scores_avg = scipy.stats.rankdata(scores_avg) 

        # bidirection normalization
        elif parsed.method == "bidirect":
            rank_scores_avg = bn.bidirect_norm(scores_orig)

        # find the motif of the highest ranking 
        index = numpy.argmax(rank_scores_avg)
        motif_max = motifs_np_all[index]
        aaid_max = aaids_dbd_all[index]/100
        corr_np_max = corrs_np_all[index]
        data_out[i] = [query, motif_max, aaid_max, corr_np_max]

        toc = time.clock()
        print i, query, toc-tic
        
    sys.stdout.write("Done\n")
        
    """ Evaluate inference """
    sys.stdout.write("\rEvaluate inference ... ")

    # get dict of name conversion
    dict_names = {}
    lines = open(parsed.dict_conv, "r").readlines()
    for line in lines:
        dict_names[line.split()[1]] = line.split()[0]

    # compute pairwise correlation of query and inferred motif
    for i in range(len(data_out)):
        # only evaluate the tfs available in scertf 
        if data_out[i][1] == None:
            data_out[i].append(-1)
        elif data_out[i][0] not in dict_names.keys():
            data_out[i].append(-2)
        else:
            dict_ranklist_inferred = {}
            lines = open(parsed.dir_cisbp_rank + data_out[i][1], "r").readlines()
            for line in lines:
                dict_ranklist_inferred[line.split()[0]] = float(line.split()[1])
            names_query_ordered = []
            ranklist_query = []
            ranklist_inferred = []
            lines = open(parsed.dir_scertf_rank + dict_names[data_out[i][0]], "r").readlines()
            for line in lines:
                names_query_ordered.append(line.split()[0])
                ranklist_query.append(line.split()[1])
                ranklist_inferred.append(dict_ranklist_inferred[line.split()[0]])
            corr_query_inferred = scipy.stats.spearmanr(ranklist_query, ranklist_inferred)[0]

            if numpy.isnan(corr_query_inferred):
                data_out[i].append(0)
            else:
                data_out[i].append(corr_query_inferred)

    sys.stdout.write("Done\n")

    """ Present and write results """
    # write data
    writer = open(parsed.dir_output + "quantile_normalization_" + parsed.method + ".txt", "w")
    writer.write("#query_motif\tinferred_motif\tdbd_aaid\tnp_cisbp_corr\tscertf_cisbp_corr\n")
    for datum in data_out:
        if datum[1] == None:
            writer.write("%s\tNone\tNone\tNone\tNone\n" % datum[0])
        elif datum[4] == -2:
            writer.write("%s\t%s\t%.5f\t%.5f\tNone(-2)\n" % (datum[0], datum[1], datum[2], datum[3]))
        else:
            writer.write("%s\t%s\t%.5f\t%.5f\t%.5f\n" % (datum[0], datum[1], datum[2], datum[3], datum[4]))
    writer.close()

    # plot histograms
    data_out = numpy.array(data_out)
    data_arr = numpy.transpose(data_out[:, 4]).astype(float)
    data_arr = numpy.delete(data_arr, numpy.where(data_arr == -2))
    plt.figure(num=None, figsize=(10,6), dpi=80)
    plt.hist(data_arr, bins=25)
    plt.title("Motif Inference by Quantile Normalization _" + parsed.method)
    plt.savefig(parsed.dir_output + "quantile_normalization_" + parsed.method + ".png")

def check_dir(file_dir):
    if not file_dir.endswith('/'):
        file_dir += '/'
    return file_dir

def get_names(file_dir):
    names = []
    filenames = glob.glob(file_dir + "/*")
    for filename in filenames:
        filename = os.path.basename(filename).split('.')[0]
        if not filename.startswith('_'):
            names.append(filename)
    return names

def parse_scores(file_dir):
    lines = open(file_dir, "r").readlines()
    motifs_all = [None] * len(lines)
    corrs_all = [None] * len(lines)
    for j in range(len(lines)):
        line = lines[j].split()
        motifs_all[j] = line[0]
        corrs_all[j] = 0 if line[1] == "nan" else float(line[1])
    data = [motifs_all, corrs_all]
    return data

if __name__ == "__main__":
    main(sys.argv)
