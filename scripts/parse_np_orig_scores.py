#!/usr/bin/env python

"""
Given a netprophet adjacency list, outputs a file for each
regulator with the scores for that regulator's targets.
"""

import sys
import os
import argparse
import csv

def parse_args(argv):
    ''' A method for taking in command line arguments and specifying
    help strings. '''
    parser = argparse.ArgumentParser(description="Used to parse a netprophet adjacency list to generate a target score files for each TF.")
    # examples with default values and types, and short and long args
    parser.add_argument('netprophet_adjlst', metavar='netprophet_adjlst',
                        help='the .adjlst output of netprophet')
    parser.add_argument('-o', '--output_directory', dest='output_directory',
                        type=str, default='./processed_hits/')
    parser.add_argument('-t', '--targets', dest='targets',
                        type=str, default=None)
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    """ Prints the input string + newline to standard out. """
    sys.stderr.write(st + "\n")

def main(argv):
    """ The main module should take in inputs from the command line,
    carry out the 'CLI' functionality of the script, then write the
    results to stdout. """
    parsed = parse_args(argv)

    assert parsed.targets != None, 'need to specify a list of target genes!'

    # use a dictionary to hold all the target scores for each
    # target
    target_scores = {}

    # make a list of target genes
    target_genes = set([line.strip() for line in open(parsed.targets)])
    targets_unscored = {}
    all_tfs = []

    reader = open(parsed.netprophet_adjlst)
    for line in csv.DictReader(reader, delimiter='\t'):
        # read the items we care about into variables
        TF_name = line['REGULATOR']
        target_name = line['TARGET']
        score = line['SCORE']
        sign = line['CSIGN']
        # modify score based on sign
        if float(sign) < 0:
            score = str(-1*float(score))

        # check if we have an entry for the TF_name already
        if TF_name not in target_scores:
            target_scores[TF_name] = []
            targets_unscored[TF_name] = set(target_genes)
            all_tfs.append(TF_name)
        # keep track of the target_name, score pair for this TF
        target_scores[TF_name].append((target_name, score))
        targets_unscored[TF_name].remove(target_name)
    # add zero scores for each TF:
    for tf in target_scores:
        for unscored in targets_unscored[tf]:
            target_scores[tf].append((unscored, '0'))

    # check if the output directory exists, and make it if it doesn't
    if not parsed.output_directory.endswith('/'):
        parsed.output_directory += '/'
    if not os.path.exists(parsed.output_directory):
        os.makedirs(parsed.output_directory)

    # loop through our collected scores, outputting files for each TF
    # containing all the scores
    for tf in target_scores:
        file_name = os.path.join(parsed.output_directory,
                                 tf)
        output_text = "target\tscore"
        for item in target_scores[tf]:
            output_text += "\n" + item[0] + "\t" + item[1]
        writer = open(file_name, 'w')
        writer.write(output_text)
        writer.close()


if __name__ == "__main__":
    main(sys.argv)
