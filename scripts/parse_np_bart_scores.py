#!/usr/bin/env python

"""
Given a netprophet bart network adjacency matrix, output an adjmtr with top informative edges 
(consistent to original netprophet edge count) and output a file for each regulator with the 
scores for that regulator's targets.
"""

import sys
import os
import argparse
import numpy
from scipy.stats import rankdata

def parse_args(argv):
    ''' A method for taking in command line arguments and specifying
    help strings. '''
    parser = argparse.ArgumentParser(description="Parse a NetProphet adjacency matrix")
    parser.add_argument('-n', '--name', dest='name', type=str)
    parser.add_argument('-i', '--dir_input', dest='dir_input', type=str)
    parser.add_argument('-o', '--dir_output', dest='dir_output', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)

    parsed.dir_input = check_dir(parsed.dir_input)
    parsed.dir_output = check_dir(parsed.dir_output)

    # edge count: 50k, 100k, 350k (original NetProphet)
    np_edge_count = 100000

    # zero out low rank edges
    lines = open(parsed.dir_input + parsed.name + ".tsv", "r").readlines()

    targets = lines[0].split()
    tfs = [None] * (len(lines)-1)
    
    adjmtr = numpy.ndarray([len(tfs), len(targets)])

    for i in range(1,len(lines)):
        temp = lines[i].split()
        tfs[i-1] = temp[0]
        for j in range(1,len(temp)):
            adjmtr[i-1, j-1] = float(temp[j])

    rankmtr = adjmtr.shape[0]*adjmtr.shape[1] + 1 - rankdata(adjmtr)
    indices_low_rank = rankmtr > np_edge_count
    indices_low_rank = numpy.reshape(indices_low_rank, (len(tfs), len(targets)))
    adjmtr[indices_low_rank] = 0

    # write data filtered adjmtr
    writer = open(parsed.dir_output + parsed.name + ".adjmtr", "w")
    for i in range(adjmtr.shape[0]):
        for j in range(adjmtr.shape[1]):
            if adjmtr[i, j] == 0:
                writer.write("0\t")
            else:
                writer.write("%0.2f\t" % adjmtr[i, j])
        writer.write("\n")
    writer.close()

    # write individual tf to target score files
    for i in range(adjmtr.shape[0]):
        writer = open(parsed.dir_output + tfs[i], "w")
        writer.write("#target\tscore\n")
        for j in range(adjmtr.shape[1]):
            if adjmtr[i, j] == 0:
                writer.write("%s\t0\n" % targets[j])
            else:
                writer.write("%s\t%0.2f\n" % (targets[j], adjmtr[i, j]))
        writer.close()

def check_dir(fd):
    if not fd.endswith('/'):
        fd += '/'
    return fd

if __name__ == "__main__":
    main(sys.argv)
