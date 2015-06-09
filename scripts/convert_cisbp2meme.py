#!/usr/bin/python

import sys
import os
import argparse

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Convert inferred motif pwms in CIS_BP database to meme format")
    parser.add_argument('cisbp_dir', metavar='cisbp_dir', help='Directory of inferred motifs')
    parser.add_argument('-o', '--output_dir', dest='output_dir', type=str, default='../resources/motif_meme')
    parser.add_argument('-t', '--tf_names', dest='tf_names', type=str, default='../resources/tf_names.txt')
    parser.add_argument('-b', '--background_frequency', dest='background_frequency', \
        type=str, default='A 0.25 C 0.25 G 0.25 T 0.25')
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)

    # write a list tf names
    writer0 = open(parsed.tf_names, "w")

    for input_fn in os.listdir(parsed.cisbp_dir):
        # get motif names
        if not input_fn.startswith("."):
            motif_name = os.path.splitext(input_fn)[0]
            writer0.write("%s\n" % motif_name)

            # parse motif file
            lines = open("%s/%s" % (parsed.cisbp_dir, input_fn),"r").readlines()
            motif_freq = []
            for i, line in enumerate(lines):
                if i > 0:
                    motif_freq.append([float(line.split()[1]), float(line.split()[2]),
                     float(line.split()[3]), float(line.split()[4])])

            # write motif to output file
            writer = open(parsed.output_dir + "/" + motif_name, "w")
            writer.write("MEME version 4.9.1\n\nALPHABET= ACGT\n\nstrands: + -\n")
            # writer.write("\nBackground letter frequencies\n%s\n" % parsed.background_frequency)
            writer.write("\nMOTIF %s\n" % motif_name)
            writer.write("letter-probability matrix: alength= 4 w= %d nsites= 20 E= 0\n" % len(motif_freq))
            for i in range(0, len(motif_freq)):
                for j in range(0, 4):
                    writer.write("%.3f\t" % motif_freq[i][j])
                writer.write("\n")
            writer.close()
    writer0.close()

def isfloat(var):
    try: 
        float(var)
        return True
    except ValueError:
        return False

if __name__ == "__main__":
    main(sys.argv)
