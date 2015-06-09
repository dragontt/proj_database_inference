#!/usr/bin/python

"""
Infer a motif by chance, evaluated by scertf correlation.
"""

import sys
import argparse
import glob
import os.path
import scipy.stats
import numpy
import random
import matplotlib.pyplot as plt

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Infer a motif by chance")
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

    """ Infer motifs by chance """
    sys.stdout.write("\rProcessing inference ... ")

    # get motifs names from both dbd and netprophet
    names_aaid = get_names(parsed.dir_dbd_aaid)
    names_np = get_names(parsed.dir_corr_np)
    names_common = list(set(names_aaid) & set(names_np))
    data_out = [None] * len(names_common)

    [motifs_cisbp, foo] = parse_scores(parsed.dir_corr_np + names_common[0])

    for i in range(len(names_common)):
        query = names_common[i]

        # pick a random motif for the query
        index_random = random.randint(0, len(motifs_cisbp))
        motif_random = motifs_cisbp[index_random]
        data_out[i] = [query, motif_random]
        
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
    writer = open(parsed.dir_output + "by_chance.txt", "w")
    writer.write("#query_motif\tinferred_motif\tscertf_cisbp_corr\n")
    for datum in data_out:
        if datum[1] == None:
            writer.write("%s\tNone\tNone(-1)\n" % datum[0])
        elif datum[2] == -2:
            writer.write("%s\t%s\tNone(-2)\n" % (datum[0], datum[1]))
        else:
            writer.write("%s\t%s\t%.5f\n" % (datum[0], datum[1], datum[2]))
    writer.close()

    # plot histograms
    data_out = numpy.array(data_out)
    data_arr = numpy.transpose(data_out[:, 2]).astype(float)
    data_arr = numpy.delete(data_arr, numpy.where(data_arr == -2))
    plt.figure(num=None, figsize=(10,6), dpi=80)
    plt.hist(data_arr, bins=25)
    plt.title("Motif Inference by Chance")
    plt.savefig(parsed.dir_output + "by_chance.png")

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
