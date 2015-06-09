#!/usr/bin/python

"""
Parse all percent identity matrix files, and output a ordered list of percent identities
of the query motif DBD to all database motif DBDs.
"""

import sys
import os
import argparse
import glob

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Parse and combine percent identities.")
    parser.add_argument('-i', '-dir_input', dest='dir_input', type=str)
    parser.add_argument('-o', '-fn_output', dest='fn_output', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed

def main(argv):
    parsed = parse_args(argv)

    if not parsed.dir_input.endswith('/'):
        parsed.dir_input += '/'

    # parse all percent identity matrix files
    filenames = glob.glob(parsed.dir_input + '*.pim')
    data = [None]*len(filenames)
    for i, filename in enumerate(filenames):
        lines = open(filename, "r").readlines()
        data[i] = [lines[2].split()[0], float(lines[2].split()[1])]

    # sort percent identity
    data = sorted(data, key=lambda a:a[1])
    data.reverse()

    # write file
    writer = open(parsed.fn_output, "w")
    for datum in data:
        writer.write("%s\t%0.5f\n" % (datum[0], datum[1]))
    writer.close()

if __name__ == "__main__":
    main(sys.argv)
