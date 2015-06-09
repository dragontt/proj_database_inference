#!/usr/bin/python

import numpy

motifs = []
lines = open("../resources/cisbp_all_motifs_no_dir_yeast.txt", "r").readlines()
for line in lines:
	motifs.append(line.split()[0])

writer = open("../resources/cisbp_all_dbds_no_dir_yeast.fasta", "w")
flag = False
lines = open("../resources/cisbp_all_dbds.fasta", "r").readlines()
for i, line in enumerate(lines):
	if i % 2 == 0:
		temp_header = line
		temp_motif = line.split('>')[1].split(':')[0]
		if temp_motif in motifs:
			flag = True
		else:
			flag = False
	else:
		if flag:
			writer.write("%s%s" % (temp_header, line))
writer.close()
