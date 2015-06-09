#!/usr/bin/python

import sys
import argparse
import glob
import os.path
from scipy.stats import rankdata
import numpy

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Compute rankings of netprophet scores in range [0,1]")
    parser.add_argument('-i', '--input_dir', dest='input_dir', type=str, default='')
    parser.add_argument('-o', '--output_dir', dest='output_dir', type=str)
    parser.add_argument('-m', '--method', dest='method', type=str, default='use_abs_nozero')
    # optional method: use_abs, use_sign, use_abs_nozero, use_sign_nozero
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)

    parsed.input_dir = check_dir(parsed.input_dir)
    parsed.output_dir = check_dir(parsed.output_dir)

    # get tf names in both netprophet and scertf
    fns = glob.glob(parsed.input_dir + "*")
    for fn in fns:
        tf = os.path.basename(fn)
        print tf

        # get scores
        lines = open(fn, "r").readlines()
        targets = []
        scores = []
        for i in range(1,len(lines)):
            target = lines[i].split()[0]
            score = float(lines[i].split()[1])
            if parsed.method == "use_abs":
                targets.append(target)
                scores.append(abs(score))
            elif parsed.method == "use_sign":
                targets.append(target)
                scores.append(score)
            elif parsed.method == "use_abs_nozero":
                if score != 0:
                    targets.append(target)
                    scores.append(abs(score))
            elif parsed.method == "use_sign_nozero":
                if score != 0:
                    targets.append(target)
                    scores.append(score)

        # rank the scores
        rankings = rankdata(scores)
        rankings = len(scores) + 1 - rankings
        sorted_indices = numpy.argsort(rankings)

        # write data
        writer = open(parsed.output_dir + tf, "w")
        for i in range(len(sorted_indices)):
            writer.write("%s\t%.2f\n" % (targets[sorted_indices[i]], rankings[sorted_indices[i]]))
        writer.close()

def check_dir(directory):
    if not directory.endswith("/"):
        directory += "/"
    return directory

if __name__ == "__main__":
    main(sys.argv)
