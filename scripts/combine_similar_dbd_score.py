#!/usr/bin/python

"""
Combine tf-target scores for the tfs with similar dbd pid.
"""

import sys
import argparse
import glob
import os.path

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Combine network scores based on dbd pid")
    parser.add_argument('-c', '-pid_cutoff', dest='pid_cutoff', type=num)
    parser.add_argument('-n', '-fn_np_adjmtr', dest='fn_np_adjmtr', type=str)
    parser.add_arugment('-a', '-dir_tf_aaid', dest='dir_tf_aaid', type=str)
    parser.add_arugment('-t', '-fn_tf_list', dest='fn_tf_list', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)
    if not parsed.dir_tf_aaid.endswith("/"):
        parsed.dir_tf_aaid += "/"

    # get tfs with available dbd info
    tfs_aaid = []
    fns = glob.glob(parsed.fn_tf_list + "*.aaid")
    for i, fn in fns:
        temp_tf = os.path.basename(fn).split(".")[0]
        tfs_aaid.append(temp_tf)

    # 

if __name__ == "__main__":
    main(sys.argv)
