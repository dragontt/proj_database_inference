#!/usr/bin/python

"""
Quantitatively evaluate the overall results of the inferred motif evaluation.
"""

import sys
import argparse
import glob
import os.path
import numpy
import matplotlib.pyplot as plt

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Evaluate the overall motif inference")
    parser.add_argument('-i1', '-dir_input1', dest='dir_input1', type=str)
    parser.add_argument('-i2', '-dir_input2', dest='dir_input2', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)

    # scertf correlation cutoff
    cutoff = 0.5

    # compute median of scertf correlation
    out1 = parse_all_data(check_dir(parsed.dir_input1), cutoff)
    out2 = parse_all_data(check_dir(parsed.dir_input2), cutoff)

    # plot histogram
    labels = ["np_orig", "np_bart_holstege"]
    # plt.hist(corrs_out0, bins, alpha=0.5, label=labels[0])
    # plt.hist(corrs_out1, bins, alpha=0.5, label=labels[1])
    plt.figure(num=None, figsize=(10,10), dpi=80)
    plt.subplot(2,1,1)
    # plt.bar(indices, medians)
    plt.bar(out1[0], out1[2], alpha=0.5, color="blue", label=labels[0])
    plt.bar(out2[0], out2[2], alpha=0.5, color="green", label=labels[1])
    plt.xticks([])
    plt.title("Median of the Scertf Correlation by Methods")
    plt.legend(loc="upper left")
    plt.subplot(2,1,2)
    # plt.bar(indices, count_corr50)
    # plt.xticks(indices, methods, rotation=90)
    plt.bar(out1[0], out1[3], alpha=0.5, color="blue", label=labels[0])
    plt.bar(out2[0], out2[3], alpha=0.5, color="green", label=labels[1])
    plt.xticks(out1[0], out1[1], rotation=90)
    plt.title("Count of the Scertf Correlation over " + str(cutoff) + " by Methods")
    plt.show()
    # plt.savefig(parsed.dir_output + "_all_resutls.png")

def check_dir(file_dir):
    if not file_dir.endswith('/'):
        file_dir += '/'
    return file_dir

def parse_all_data(directory, cutoff):
    fns = glob.glob(directory + "*.txt")
    indices = numpy.arange(len(fns))
    methods = [None] * len(fns)
    medians = [None] * len(fns)
    count_corr50 = [None] * len(fns)
    for i in range(len(fns)):
        methods[i] = os.path.basename(fns[i]).split('.')[0]
        [medians[i], count_corr50[i]] = parse_data(fns[i], cutoff)
    out = [indices, methods, medians, count_corr50]
    return out

def parse_data(fn, cutoff):
    lines = open(fn, "r").readlines()
    data = [None] * (len(lines)-1)
    count = 0
    for i in range(1,len(lines)):
        linesplit = lines[i].split()
        datum = linesplit[len(linesplit)-1]
        if datum.startswith("None"):
            data[i-1] = -1 
        else:
            data[i-1] = float(datum)
            if float(datum) >= cutoff:
                count += 1
    med = numpy.median(numpy.array(data))
    return [med, count]

if __name__ == "__main__":
    main(sys.argv)
