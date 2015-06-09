#!/usr/bin/python

"""
Query a list of database inference results, combine the inferred pwm fimo scan scores,
which were preprocessed in the motif inference step.
"""

import sys
import argparse
import glob
import os.path
import numpy

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Combine inferred pwm fimo scan scores.")
    parser.add_argument('-f', '--dir_fimo', dest='dir_fimo', type=str)   
    parser.add_argument('-t', '--fn_tfs', dest='fn_tfs', type=str)
    parser.add_argument('-g', '--fn_targets', dest='fn_targets', type=str) 
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
    fns = glob.glob(parsed.dir_fimo + "*.summary")
    for fn in fns:
        tf = os.path.basename(fn).split(".")[0].strip()
        index = tfs.index(tf)
        dict_scores = get_fimo_scores(fn)
        for j in range(len(targets)):
            t = targets[j]
            adjmtr[index, j] = dict_scores[t] if t in dict_scores else 0
        
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
                writer.write("0.\t")
            else:
                writer.write("%0.15f\t" % adjmtr[i,j])
        writer.write("\n")
    writer.close()

if __name__ == "__main__":
    main(sys.argv)
