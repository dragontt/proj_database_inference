#!/usr/bin/python

import numpy

lines = open('../resources/yeast_dbd_np_scertf_names.txt', 'r').readlines()
common_names = []
for line in lines:
	if len(line) > 0:
		common_names.append(line.split()[0])

lines = open('../resources/yeast.dbd.aligned.pim' ,'r').readlines()
dbd_names = []
scores = []
index_del = []
for line in lines:
	if len(line.strip()) > 0 and not line.startswith('#'):
		temp_line = line.split()
		temp_index = temp_line[0].strip(':')	# the index starts from 1
		temp_name = temp_line[1].split('_')[0]
		temp_scores = temp_line[2:len(temp_line)]
		dbd_names.append(temp_name)
		scores.append(temp_scores)
		if not temp_name in common_names:
			index_del.append(int(temp_index))

for i in range(len(index_del)):
	scores = numpy.delete(scores, (index_del[i]-1-i), axis=0)
	scores = numpy.delete(scores, (index_del[i]-1-i), axis=1)
	dbd_names.remove(dbd_names[index_del[i]-i])

writer = open('../resources/yeast.dbd.aligned.pim_filtered', 'w')
for i in range(len(scores)):
	for j in range(len(scores[i])):
		writer.write('%s\t' % scores[i][j])
	writer.write('\n')
writer.close()

writer = open('../resources/yeast_dbd_names_filtered.txt', 'w')
for i in range(len(dbd_names)):
	writer.write('%s\n' % dbd_names[i])
writer.close()
