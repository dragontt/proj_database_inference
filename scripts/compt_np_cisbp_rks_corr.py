#!/usr/bin/python

"""
Compute Spearman's rank order correlation of the NetProphet score rankings and the 
database (CIS-BP) motif-target scan score rankings.
"""

import sys
import argparse
import os
import glob
import os.path
import compt_pw_corr
import time

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Compute rank order correlations of np scores and datrabase motif-target score rankings")
    parser.add_argument('-n', '--dir_np', dest='dir_np', type=str, help='NetProphet rankings')
    parser.add_argument('-c', '--dir_cisbp', dest='dir_cisbp', type=str, help="Database motif-target score rankings")
    parser.add_argument('-o', '--dir_output', dest='dir_output', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)

    parsed.dir_np = check_dir(parsed.dir_np)
    parsed.dir_cisbp = check_dir(parsed.dir_cisbp)
    parsed.dir_output = check_dir(parsed.dir_output)

    # get the list of query motif names and their ranked lists
    query_tf_fns = glob.glob(parsed.dir_np + "/*")
    query_tfs = []
    query_tfs_skip = []
    ranks_query = []
    for i in range(len(query_tf_fns)):
        temp_tf = os.path.basename(query_tf_fns[i])
        if os.stat(query_tf_fns[i]).st_size == 0:
            query_tfs_skip.append(temp_tf)
        else:
    	   query_tfs.append(temp_tf)
    	   ranks_query.append(parse_rank(query_tf_fns[i]))
    print "TF count:", len(ranks_query), "Skip TF:", query_tfs_skip

    # get the list of database motif names and their ranked lists
    db_motif_fns = glob.glob(parsed.dir_cisbp + "/*")
    db_motifs = [None] * len(db_motif_fns)
    ranks_db = [None] * len(db_motifs)
    for i in range(len(db_motif_fns)):
    	db_motifs[i] = os.path.basename(db_motif_fns[i])
    	ranks_db[i] = parse_rank(db_motif_fns[i])
    db_motifs = [x for x in db_motifs if x is not None]
    ranks_db = [x for x in ranks_db if x is not None]
    print "Database motif count:", len(ranks_db)

	# compute ranked list correlations
    for i in range(len(query_tfs)):
    	tic = time.clock()
    	writer = open(parsed.dir_output + query_tfs[i], "w")
    	for j in range(len(db_motifs)):
    		corr = compt_pw_corr.compute_corr(ranks_db[j], ranks_query[i])
    		writer.write("%s\t%.5f\n" % (db_motifs[j], corr))
    	writer.close()
    	toc = time.clock()
    	print query_tfs[i], i, "/", len(query_tfs)-1, "Time elapsed:", toc-tic, "sec"

def check_dir(dirname):
    if not dirname.endswith('/'):
        dirname += '/'
    return dirname

def parse_rank(filename):
    rank = {}
    lines = open(filename, 'r').readlines()
    for line in lines:
        linesplit = line.split()
        rank[linesplit[0]] = linesplit[1]
    return rank

if __name__ == "__main__":
    main(sys.argv)
