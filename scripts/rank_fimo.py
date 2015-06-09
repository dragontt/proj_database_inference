#!/usr/bin/python

import sys
import argparse
import os.path
import glob
from scipy.stats import rankdata
import numpy
import operator

global num_targets
num_targets = 6220

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Sort rankings of fimo outputs")
    parser.add_argument('-i', '--input_dir', dest='input_dir', type=str)
    parser.add_argument('-t', '--target_names', dest='target_names', type=str)
    parser.add_argument('-o', '--output_dir', dest='output_dir', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)

    if not parsed.input_dir.endswith('/'):
        parsed.input_dir += '/'
    if not parsed.output_dir.endswith('/'):
        parsed.output_dir += '/'

    # get target names
    target_names = [0]*num_targets
    lines = open(parsed.target_names, 'r').readlines()
    for i, line in enumerate(lines):
        if not line.startswith('\n'):
            target_names[i] = line.strip()

    # get tf names
    tf_names = []
    filenames = glob.glob(parsed.input_dir + "*.summary")
    for fn in filenames:
        temp_name = os.path.basename(fn).split('.')
        tf_names.append(temp_name[0])

    # sort fimo rankings for each tf
    for i in range(len(tf_names)):
        print 'Processing %s' % tf_names[i]
        lines = open(filenames[i], 'r').readlines()
        dict_rank = {}
        for line in lines:
            temp_name = line.split()[1].strip()
            temp_score = line.split()[3].strip()
            # only parse the aligned target in netprophet
            if temp_name in target_names:   
                dict_rank[temp_name] = abs(float(temp_score))
        # set target score not aligned but in netprophet to 0
        target_name_diff = set(target_names) - set(dict_rank.keys())
        for item in target_name_diff:
            dict_rank[item] = 0
        sorted_rank = compute_rank(dict_rank)

        # write output file
        out_file = parsed.output_dir + tf_names[i]
        writer = open(out_file, 'w')
        for i in range(len(sorted_rank)):
            writer.write("%s\t%.8f\n" % (sorted_rank[i][0], sorted_rank[i][1]))
        writer.close()

def compute_rank(dict_scores):
    # sort np scores
    sorted_scores = sorted(dict_scores.items(), key=operator.itemgetter(1))
    sorted_scores.reverse()
    rank_out = [0]*num_targets
    rank_out[0] = [sorted_scores[0][0], 1]
    flag = False
    for i in range(1, num_targets):
        rank_out[i] = [sorted_scores[i][0], i+1]
        # deal with tie rankings
        if sorted_scores[i][1] == sorted_scores[i-1][1]:
            if i != num_targets-1:
                if not flag:
                    temp_index = [i-1, i]
                    temp_rank = [rank_out[i-1][1], rank_out[i][1]]
                    flag = True
                else:
                    temp_index.append(i)
                    temp_rank.append(rank_out[i][1])
            else:
                if flag:
                    temp_index.append(i)
                    temp_rank.append(rank_out[i][1])
                    temp_mean_rank = numpy.mean(temp_rank)
                    for j in range(len(temp_index)):
                        rank_out[temp_index[j]][1] = temp_mean_rank
        else:
            if flag:
                temp_mean_rank = numpy.mean(temp_rank)
                for j in range(len(temp_index)):
                    rank_out[temp_index[j]][1] = temp_mean_rank
                flag = False
    return rank_out

if __name__ == "__main__":
    main(sys.argv)
