#!/usr/bin/python

import sys
import argparse

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Convert scertf's gene names to systematic names")
    parser.add_argument('input_fn', metavar='input_fn')
    parser.add_argument('-d', '--lookup_dict', dest='lookup_dict', type=str)
    parser.add_argument('-o', '--output_fn', dest='output_fn', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)

    lines = open(parsed.lookup_dict, 'r').readlines()
    name_dict = {}
    for i, line in enumerate(lines):
        if i > 0:
            linesplit = line.split('\t')
            if len(linesplit[0]) != 0:
                name_dict[linesplit[0]] = linesplit[1].strip()

    writer = open(parsed.output_fn, 'w')

    lines = open(parsed.input_fn, 'r').readlines()
    for line in lines:
        gene_name = line.strip()
        if len(gene_name) != 0:
            if gene_name.startswith('Y') and (gene_name.endswith('W') or gene_name.endswith('C')):
                writer.write('%s\t%s\n' % (gene_name, gene_name))
            else:
                writer.write('%s\t%s\n' % (gene_name, name_dict[gene_name]))
    writer.close()

if __name__ == "__main__":
    main(sys.argv)