#!/usr/bin/python

"""
Generate fasta files of the query motif's DBD to all database DBDs individually.
Each file has 2 redudant query DBD sequence, and 1 database DBD sequence.
"""

import sys
import os
import argparse

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Generate paired dbd")
    parser.add_argument('-m', '-query_motif', dest='query_motif', type=str)
    parser.add_argument('-f1', '-dbd_fasta1', dest='dbd_fasta1', type=str)
    parser.add_argument('-f2', '-dbd_fasta2', dest='dbd_fasta2', type=str)
    parser.add_argument('-o', '-dir_output', dest='dir_output', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed

def main(argv):
    parsed = parse_args(argv)

    if not parsed.dir_output.endswith('/'):
        parsed.dir_output += '/'

    # get the dbd sequence of query motif from dbd_fasta1 file
    lines = open(parsed.dbd_fasta1, "r").readlines()
    for i in range(len(lines)/2):
        temp_motif = lines[i*2].split('>')[1].split(':')[0]
        if temp_motif == parsed.query_motif:
            query_header = lines[i*2].strip()
            query_dbd = lines[i*2+1].strip()
            break
        if i == len(lines):
            print "No query motif found in fasta file."
            sys.exit()

    # get all dbd sequences from the database in dbd_fasta2 file
    lines = open(parsed.dbd_fasta2, "r").readlines()
    num_dbds = len(lines)/2
    all_headers = [None]*num_dbds
    all_dbds = [None]*num_dbds
    for i in range(num_dbds):
        all_headers[i] = lines[i*2].strip()
        all_dbds[i] = lines[i*2+1].strip()

    # create temp file of paired query dbd and datatbase dbd in fasta format
    for i in range(num_dbds):
        temp_name = all_headers[i].split('>')[1].split(':')
        writer = open(parsed.dir_output + parsed.query_motif + '_' + \
            temp_name[0] + '_' + temp_name[1] + '.fasta', "w")
        writer.write("%s\n%s\n%s\n%s\n>void\n%s\n" % \
            (query_header, query_dbd, all_headers[i], all_dbds[i], query_dbd))
        writer.close()

if __name__ == "__main__":
    main(sys.argv)
