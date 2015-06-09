#!/usr/bin/python

"""
Query a list of scertf tfs, combine the scertf pwm fimo scan scores.
"""
# Local process: pwd
# /Users/KANG/proj_motifcomparison/scritpts
# python combine_scertf_fimo_score.py -t ../resources/np_network_orig/tf.orfs -g ../resources/np_network_orig/orfs -f ../output/scertf_fimo/ -o ../output/combined_output_orig/scertf_pwm/binding.adjmtr

import sys
import argparse
import os.path
import numpy
import subprocess

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Combine scertf pwm fimo scan scores.")
    parser.add_argument('-t', '--fn_tfs', dest='fn_tfs', type=str)
    parser.add_argument('-g', '--fn_targets', dest='fn_targets', type=str)
    parser.add_argument('-f', '--dir_fimo', dest='dir_fimo', type=str)    
    parser.add_argument('-o', '--fn_adjmtr', dest='fn_adjmtr', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)

    # check directory
    if not parsed.dir_fimo.endswith("/"):
        parsed.dir_fimo += "/"

    # get the lists of tfs, targets
    tfs = get_list(parsed.fn_tfs)
    targets = get_list(parsed.fn_targets)
    adjmtr = numpy.zeros([len(tfs), len(targets)])

    # build the adjmtr
    for i in range(len(tfs)):
        # get fimo scores for the scertf pwm
        fn_pwm = parsed.dir_fimo + tfs[i] + ".summary"
        if os.path.isfile(fn_pwm):
            dict_scores = get_fimo_scores(fn_pwm)
            for j in range(len(targets)):
                t = targets[j]
                adjmtr[i, j] = dict_scores[t] if t in dict_scores else 0
        
    # write adjmtr file
    numpy.savetxt(parsed.fn_adjmtr, adjmtr, fmt='%0.15f', delimiter='\t')
    # write_adjmtr(adjmtr, parsed.fn_adjmtr)

def get_list(fn):
    lines = open(fn, "r").readlines()
    l = [None] * len(lines)
    for i in range(len(lines)):
        l[i] = lines[i].split()[0]
    return l

def get_fimo_scores(fn):
    lines = open(fn, "r").readlines()
    d = {}
    for line in lines:
        temp_name = line.split()[1].strip()
        temp_score = float(line.split()[3].strip())
        d[temp_name] = temp_score
    return d

def write_adjmtr(adjmtr, fn):
    writer = open(fn, "w")
    for i in range(len(adjmtr)):
        for j in range(len(adjmtr[i])):
            if adjmtr[i,j] == 0:
                writer.write("%d\t" % adjmtr[i,j])
            else:
                writer.write("%0.15f\t" % adjmtr[i,j])
        writer.write("\n")
    writer.close()

if __name__ == "__main__":
    main(sys.argv)
