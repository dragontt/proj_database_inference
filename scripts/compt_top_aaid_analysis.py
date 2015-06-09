#!/usr/bin/python

"""
Analyze the quality of the inferred motif with the highest DNA binding domain sequence percent identity values 
in terms of its scertf alignment correlation ranking and distance. 
"""

import sys
import os
import argparse
import glob
import numpy
import scipy.stats
import matplotlib.pyplot as plt

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Analyze the quality of the inferred motif with the highest percent identity.")
    parser.add_argument('dir_inputs', metavar='dir_inputs')
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)
    parsed.dir_inputs = check_dir(parsed.dir_inputs)

    # get genetic and systematic names of tfs 
    tf_gen_names = []
    tf_sys_names = []
    for (dirpath, dirnames, filenames) in os.walk(parsed.dir_inputs):
        for filename in filenames:
            if filename.endswith('.txt') and not filename.startswith('_'):
                tf_gen_names.append(filename.split('_')[0])
                tf_sys_names.append(filename.split('_')[1])

    # initialize saving data
    out_vals = []
    out_pcts = []

    writer_val = open(parsed.dir_inputs + '_analysis_val.txt', 'w')
    writer_pct = open(parsed.dir_inputs + '_analysis_pct.txt', 'w')
    header = ['avg_top_st_corr', 'avg_top_st_rank', 'avg_aaid_dist_max_st', 'max_aaid_st_corr', 'max_aaid_st_rank', 'max_aaid_dist_max_st']
    writer_val.write('#tf_name\t%s\t%s\t%s\t%s\t%s\t%s\n' % (header[0], header[1], header[2], header[3], header[4], header[5]))
    writer_pct.write('#tf_name\t%s\t%s\t%s\t%s\t%s\t%s\n' % (header[0], header[1], header[2], header[3], header[4], header[5]))

    for i in range(len(tf_gen_names)):
        # get aaid and scertf correlations of all tfs
        gen_name = tf_gen_names[i]
        sys_name = tf_sys_names[i]
        filename = glob.glob(parsed.dir_inputs + gen_name + '*.txt')
        lines = open(filename[0], 'r').readlines()
        num_motifs = len(lines) - 1
        aaid = []
        st_corr = []
        for line in lines:
            if not line.startswith('#'):
                aaid.append(float(line.split()[0]))
                st_corr.append(float(line.split()[1]))
        st_rank = num_motifs + 1 - scipy.stats.rankdata(st_corr)
        
        # get the most correlated motif with scertf
        max_st_index = numpy.argmax(st_corr)
        max_aaid = aaid[max_st_index]
        max_st_corr = st_corr[max_st_index]

        # analyze the motifs with the highest percent identity value
        top_aaid_val = max(aaid)
        top_aaid_indices = [i for i, val in enumerate(aaid) if val == top_aaid_val]
        top_aaid = []
        top_st_corr = []
        top_st_rank = []
        for index in top_aaid_indices:
            top_aaid.append(aaid[index])
            top_st_corr.append(st_corr[index])
            top_st_rank.append(st_rank[index])
        
        avg_top_st_corr = numpy.mean(top_st_corr)
        avg_top_st_rank = numpy.mean(top_st_rank)
        pct_avg_top_st_corr = avg_top_st_corr/max_st_corr
        pct_avg_top_st_rank = avg_top_st_rank/num_motifs

        avg_aaid_dist_max_st = numpy.mean([max_aaid - top_aaid[j] for j in range(len(top_aaid))])
        if max_aaid == 0:
            pct_avg_aaid_dist_max_st = 0
        else:
            pct_avg_aaid_dist_max_st = avg_aaid_dist_max_st/max_aaid

        # analyze the motif with the best scertf correlation among the motifs with the highest percent idenetity
        max_aaid_index = top_aaid_indices[0]
        for i in range(2, len(top_aaid_indices)):
            if st_corr[max_aaid_index] < st_corr[top_aaid_indices[i]]:
                max_aaid_index = top_aaid_indices[i]

        max_aaid_st_corr = st_corr[max_aaid_index]
        max_aaid_st_rank = st_rank[max_aaid_index]
        pct_max_aaid_st_corr = max_aaid_st_corr/max_st_corr
        pct_max_aaid_st_rank = max_aaid_st_rank/num_motifs

        max_aaid_dist_max_st = aaid[max_aaid_index] - max_aaid
        if max_aaid == 0:
            pct_max_aaid_dist_max_st = 0
        else:
            pct_max_aaid_dist_max_st = max_aaid_dist_max_st/max_aaid

        # save data
        temp_vals = [avg_top_st_corr, avg_top_st_rank, avg_aaid_dist_max_st, max_aaid_st_corr, max_aaid_st_rank, max_aaid_dist_max_st]
        temp_pcts = [pct_avg_top_st_corr, pct_avg_top_st_rank, pct_avg_aaid_dist_max_st, pct_max_aaid_st_corr, pct_max_aaid_st_rank, pct_max_aaid_dist_max_st]
        out_vals.append(temp_vals)
        out_pcts.append(temp_pcts)
        writer_val.write('%s\t%.5f\t%.0f\t%.5f\t%.5f\t%.0f\t%.5f\n' % \
            (sys_name, temp_vals[0], temp_vals[1], temp_vals[2], temp_vals[3], temp_vals[4], temp_vals[5]))
        writer_pct.write('%s\t%.5f\t%.0f\t%.5f\t%.5f\t%.0f\t%.5f\n' % \
            (sys_name, temp_pcts[0], temp_pcts[1], temp_pcts[2], temp_pcts[3], temp_pcts[4], temp_pcts[5]))

    # format ouput values and percentages
    out_vals = numpy.array(out_vals)
    out_pcts = numpy.array(out_pcts)

    # plot histograms
    num_bins = 20
    plt.figure(num=None, figsize=(18, 10), dpi=80)
    for i in range(6):
        plt.subplot(2, 3, i+1)
        plt.hist(out_vals[:,i], bins=num_bins, color='blue')
        plt.title('%s' % header[i])
    plt.savefig(parsed.dir_inputs + '_analysis_vals_hist.png')

    plt.figure(num=None, figsize=(18, 10), dpi=80)
    for i in range(6):
        plt.subplot(2, 3, i+1)
        plt.hist(out_pcts[:,i], bins=num_bins, color='blue')
        plt.title('%s' % header[i])
    plt.savefig(parsed.dir_inputs + '_analysis_pcts_hist.png')

def check_dir(dirname):
    if not dirname.endswith('/'):
        dirname += '/'
    return dirname

if __name__ == "__main__":
    main(sys.argv)
