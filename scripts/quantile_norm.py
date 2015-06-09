#!/usr/bin/python

""" 
Quantile normalization maps multiple distributions to be identical in their stats properties.
"""

from scipy.stats import rankdata
import numpy

def quantile_norm(data):
	""" data must be in numpy.ndarray type """
	N = len(data)
	M = len(data[0])

	# rank and sort data values
	ranks = [None] * N
	for k in range(N):
		ranks[k] = rankdata(data[k], method="min")
		data[k] = numpy.array(sorted(data[k]))

	# normalize sorted data and map rank to values
	normalized = [None] * M
	for i in range(M):
		normalized[i] = numpy.mean(data[:,i])

	data = numpy.zeros((N,M))
	for k in range(N):
		for i in range(M):
			index = ranks[k][i]-1
			index = int(index)
			data[k][i] = normalized[index]
	
	return data
