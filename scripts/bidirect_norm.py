#!/usr/bin/python

""" 
Bidirectional normalization maps distributions and combines rankings.
Specifically, the algorithm maps Distribution A to B (served as the base), then combines the scores of normalized 
Distribution A and base B, and computes a ranked list of the combined scores. Repeat the same procedure to map
Distribution B to A, yielding another ranked list. Finally average the two ranked lists as the normalized output.
"""

from scipy.stats import rankdata
import numpy

def bidirect_norm(data):
	""" data must be in numpy.ndarray type """
	rank1 = unidirect_norm(data[0], data[1])
	rank2 = unidirect_norm(data[1], data[0])
	rank_out = rankdata(rank1 + rank2)
	return rank_out

def unidirect_norm(A, B):
	# map distribution in one direction and rank the combined scores
	N = len(A)
	A_rank = rankdata(A)
	B_rank = rankdata(B)
	B_dict = {}
	for i in range(N):
		B_dict[B_rank[i]] = B[i]
	A_norm = numpy.zeros(N)
	for i in range(N):
		if A_rank[i] not in B_dict.keys():
			my_key = min(B_dict.keys(), key = lambda score:abs(score-A_rank[i]))
		else:
			my_key = A_rank[i]
		A_norm[i] = B_dict[my_key]
	AB = A_norm + B
	AB_rank = rankdata(AB)
	return AB_rank
