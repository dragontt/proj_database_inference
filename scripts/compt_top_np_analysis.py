#!/usr/bin/python

"""
Analyze the quality of the inferred motifs in the top 5 percent correlations with netprophet score rankings,
in terms of their average ranking in correlations with scertf alignments and average distance (in the space of 
netprophet correlation) to the motif with the highest correlation to the scertf alignment. And analyze the
highest correlated inferred motif in netprophet in scertf alignment correlation ranking and distance. 
"""

import sys
import os
import argparse
import glob
import numpy
import scipy.stats
import matplotlib.pyplot as plt

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Analyze the quality of the inferred motif with top correlation with netprophet score rankings.")
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
    header = ['avg_top_st_corr', 'avg_top_st_rank', 'avg_np_dist_max_st', 'max_np_st_corr', 'max_np_st_rank', 'max_np_dist_max_st']
    writer_val.write('#tf_name\t%s\t%s\t%s\t%s\t%s\t%s\n' % (header[0], header[1], header[2], header[3], header[4], header[5]))
    writer_pct.write('#tf_name\t%s\t%s\t%s\t%s\t%s\t%s\n' % (header[0], header[1], header[2], header[3], header[4], header[5]))

    for i in range(len(tf_gen_names)):
        # get np and scertf correlations of all tfs
        gen_name = tf_gen_names[i]
        sys_name = tf_sys_names[i]
        filename = glob.glob(parsed.dir_inputs + gen_name + '*.txt')
        lines = open(filename[0], 'r').readlines()
        num_motifs = len(lines) - 1
        np_corr = []
        st_corr = []
        for line in lines:
            if not line.startswith('#'):
                np_corr.append(float(line.split()[0]))
                st_corr.append(float(line.split()[1]))
        st_rank = num_motifs + 1 - scipy.stats.rankdata(st_corr)
        
        # get the most correlated motif with scertf
        max_st_index = numpy.argmax(st_corr)
        max_np_corr = np_corr[max_st_index]
        max_st_corr = st_corr[max_st_index]

        # analyze top 5 percent motifs correlated with netprophet
        top_np_indices = numpy.argsort(numpy.array(np_corr))[::-1][0:int(num_motifs*0.05)+1]
        top_np_corr = []
        top_st_corr = []
        top_st_rank = []
        for index in top_np_indices:
            top_np_corr.append(np_corr[index])
            top_st_corr.append(st_corr[index])
            top_st_rank.append(st_rank[index])
        
        avg_top_st_corr = numpy.mean(top_st_corr)
        avg_top_st_rank = numpy.mean(top_st_rank)
        pct_avg_top_st_corr = avg_top_st_corr/max_st_corr
        pct_avg_top_st_rank = avg_top_st_rank/num_motifs

        avg_np_dist_max_st = numpy.mean([max_np_corr - top_np_corr[j] for j in range(len(top_np_corr))])
        pct_avg_np_dist_max_st = avg_np_dist_max_st/max_np_corr

        # analyze the most correlated motif with netprophet
        max_np_index = top_np_indices[0]

        max_np_st_corr = st_corr[max_np_index]
        max_np_st_rank = st_rank[max_np_index]
        pct_max_np_st_corr = max_np_st_corr/max_st_corr
        pct_max_np_st_rank = max_np_st_rank/num_motifs

        max_np_dist_max_st = np_corr[max_np_index] - max_np_corr
        pct_max_np_dist_max_st = max_np_dist_max_st/max_np_corr

        # save data
        temp_vals = [avg_top_st_corr, avg_top_st_rank, avg_np_dist_max_st, max_np_st_corr, max_np_st_rank, max_np_dist_max_st]
        temp_pcts = [pct_avg_top_st_corr, pct_avg_top_st_rank, pct_avg_np_dist_max_st, pct_max_np_st_corr, pct_max_np_st_rank, pct_max_np_dist_max_st]
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
