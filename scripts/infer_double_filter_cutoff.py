#!/usr/bin/python

"""
Infer a motif by first filtering the motifs with DNA binding domain percent identity higher than 
or equal to 50%, and picking the most informative motif that has the highest correlation with 
netprophet score rankings from the filtered motifs.
"""

import sys
import argparse
import glob
import os.path
import scipy.stats
import numpy
import matplotlib.pyplot as plt

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Infer the most informative motif")
    parser.add_argument('-a', '-dir_dbd_aaid', dest='dir_dbd_aaid', type=str)
    parser.add_argument('-n', '-dir_corr_np', dest='dir_corr_np', type=str)
    parser.add_argument('-o', '-dir_output', dest='dir_output', type=str)
    parser.add_argument('-d', '-dict_conv', dest='dict_conv', type=str, default='resources/np_scertf_names.txt')
    parser.add_argument('-c', '-dir_cisbp_rank', dest='dir_cisbp_rank', type=str)
    parser.add_argument('-s', '-dir_scertf_rank', dest='dir_scertf_rank', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed

def errprint(st):
    sys.stderr.write(st + "\n")

def main(argv):
    parsed = parse_args(argv)
    parsed.dir_dbd_aaid = check_dir(parsed.dir_dbd_aaid)
    parsed.dir_corr_np = check_dir(parsed.dir_corr_np)
    parsed.dir_output = check_dir(parsed.dir_output)
    parsed.dir_cisbp_rank = check_dir(parsed.dir_cisbp_rank)
    parsed.dir_scertf_rank = check_dir(parsed.dir_scertf_rank)

    # cutoff of aaid
    aaid_cutoffs = [50, 40, 30]
    data_arr = [None] * len(aaid_cutoffs)

    for k in range(len(aaid_cutoffs)):

        """ Infer motifs by double filtering """
        aaid_cutoff = aaid_cutoffs[k]
        
        sys.stdout.write("\rProcessing inference ... ")

        # get motifs names from both dbd and netprophet
        names_aaid = get_names(parsed.dir_dbd_aaid)
        names_np = get_names(parsed.dir_corr_np)
        names_common = list(set(names_aaid) & set(names_np))
        count_not_inferred = 0
        data_out = [None] * len(names_common)

        for i in range(len(names_common)):
            query = names_common[i]
            # parse and filter aaids of the query
            motifs_filtered = []
            aaid_filtered = []
            lines = open(parsed.dir_dbd_aaid + query + ".aaid", "r").readlines()
            for j in range(len(lines)):
                line = lines[j].split()
                if float(line[1]) >= aaid_cutoff:
                    motifs_filtered.append(line[0].split(':')[0])
                    aaid_filtered.append(float(line[1]))
                else:
                    break
            
            if not motifs_filtered:
                count_not_inferred += 1
                data_out[i] = [query, None, -1, -1]
            else:
                # parse np corrs of the query
                lines = open(parsed.dir_corr_np + query, "r").readlines()
                motifs_np_all = [None] * len(lines)
                corrs_np_all = [None] * len(lines)
                indices_np_filtered = []
                for j in range(len(lines)):
                    line = lines[j].split()
                    motifs_np_all[j] = line[0]
                    corrs_np_all[j] = 0 if line[1] == "nan" else float(line[1])
                    if line[0] in motifs_filtered:
                        indices_np_filtered.append(j)

                # temp_rank = scipy.stats.rankdata(corrs_np_all)
                # rank_corrs_np_all = len(corrs_np_all) + 1 - temp_rank

                # find the most informative motif 
                index_np_max = indices_np_filtered[0]
                motif_np_max = motifs_np_all[index_np_max]
                corr_np_max = corrs_np_all[index_np_max]
                for j in range(1,len(indices_np_filtered)):
                    temp_motif = motifs_np_all[indices_np_filtered[j]]
                    temp_corr = corrs_np_all[indices_np_filtered[j]]
                    if corrs_np_all[indices_np_filtered[j]] > corr_np_max:
                        index_np_max = indices_np_filtered[j]
                        motif_np_max = motifs_np_all[index_np_max]
                        corr_np_max =corrs_np_all[index_np_max]
                # rank_corr_np_max = rank_corrs_np_all[index_np_max]
                aaid_max = aaid_filtered[motifs_filtered.index(motif_np_max)]
                data_out[i] = [query, motif_np_max, aaid_max, corr_np_max]
        
        sys.stdout.write("Done\n")
            
        """ Evaluate inference """
        sys.stdout.write("\rEvaluating inference ... ")

        # get dict of name conversion
        dict_names = {}
        lines = open(parsed.dict_conv, "r").readlines()
        for line in lines:
            dict_names[line.split()[1]] = line.split()[0]

        # compute pairwise correlation of query and inferrred motif
        for i in range(len(data_out)):
            # only evaluate the tfs available in scertf 
            if data_out[i][1] == None:
                data_out[i].append(-1)
            elif data_out[i][0] not in dict_names.keys():
                data_out[i].append(-2)
            else:
                dict_ranklist_inferred = {}
                lines = open(parsed.dir_cisbp_rank + data_out[i][1], "r").readlines()
                for line in lines:
                    dict_ranklist_inferred[line.split()[0]] = float(line.split()[1])
                names_query_ordered = []
                ranklist_query = []
                ranklist_inferred = []
                lines = open(parsed.dir_scertf_rank + dict_names[data_out[i][0]], "r").readlines()
                for line in lines:
                    names_query_ordered.append(line.split()[0])
                    ranklist_query.append(line.split()[1])
                    ranklist_inferred.append(dict_ranklist_inferred[line.split()[0]])
                corr_query_inferred = scipy.stats.spearmanr(ranklist_query, ranklist_inferred)[0]

                if numpy.isnan(corr_query_inferred):
                    data_out[i].append(0)
                else:
                    data_out[i].append(corr_query_inferred)

        sys.stdout.write("Done\n")

        """ Present and write results """
        # print messages
        print "DBD percent identity cutoff:", aaid_cutoff
        print "Inferred motif count:", len(names_common)-count_not_inferred,'/', len(names_common)

        # write data
        writer = open(parsed.dir_output + "dbl_flt_pid_cutoff_" + str(aaid_cutoff) + ".txt", "w")
        writer.write("#query_motif\tinferred_motif\tdbd_aaid\tnp_cisbp_corr\tscertf_cisbp_corr\n")
        for datum in data_out:
            if datum[1] == None:
                writer.write("%s\tNone\tNone(-1)\tNone(-1)\tNone(-1)\n" % datum[0])
            elif datum[4] == -2:
                writer.write("%s\t%s\t%.5f\t%.5f\tNone(-2)\n" % (datum[0], datum[1], datum[2], datum[3]))
            else:
                writer.write("%s\t%s\t%.5f\t%.5f\t%.5f\n" % (datum[0], datum[1], datum[2], datum[3], datum[4]))
        writer.close()

        data_out = numpy.array(data_out)
        data_arr[k] = numpy.transpose(data_out[:, 4]).astype(float)
        data_arr[k] = numpy.delete(data_arr[k], numpy.where(data_arr[k] == -2))

    # plot histograms
    plt.figure(num=None, figsize=(10,18), dpi=80)
    for k in range(len(aaid_cutoffs)):
        plt.subplot(len(aaid_cutoffs), 1, k+1)
        plt.hist(data_arr[k], bins=25)
        plt.title("Scertf Evaluation of Inference with DBD PID Cutoff " + str(aaid_cutoffs[k]) + "%")
    # plt.show()
    plt.savefig(parsed.dir_output + "dbl_flt_pid_cutoff.png")

def check_dir(file_dir):
    if not file_dir.endswith('/'):
        file_dir += '/'
    return file_dir

def get_names(file_dir):
    names = []
    filenames = glob.glob(file_dir + "/*")
    for filename in filenames:
        filename = os.path.basename(filename).split('.')[0]
        if not filename.startswith('_'):
            names.append(filename)
    return names

if __name__ == "__main__":
    main(sys.argv)
