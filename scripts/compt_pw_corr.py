#!/usr/bin/python

"""
Compute spearman's rank order correlation, given two dictionaries of rankings. 
infer_rank: rankings of inferred motif to s.cerevisiae promoters alignment score
konwn_rank: rankings of known motif (netprophet or scertf) to s.cerevisiae promoters alignment score
"""

import operator
import scipy.stats

def compute_corr(infer_rank, known_rank):
    # sort rankings of inferred motifs
    list_infer_rank = [0]*len(known_rank)
    list_known_rank = [0]*len(known_rank)
    
    sorted_known_rank = sorted(known_rank.items(), key=operator.itemgetter(1))
    sorted_known_rank.reverse()

    for i in range(len(sorted_known_rank)):
        list_known_rank[i] = sorted_known_rank[i][1]
        if sorted_known_rank[i][0] in infer_rank:
            list_infer_rank[i] = infer_rank[sorted_known_rank[i][0]]
    list_infer_rank = list(scipy.stats.rankdata(list_infer_rank))

    # compute spearman's rank order correlation
    temp = scipy.stats.spearmanr(list_infer_rank, list_known_rank)
    corr = temp[0]
    
    return corr