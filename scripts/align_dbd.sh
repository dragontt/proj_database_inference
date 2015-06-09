#! /bin/bash

# Compute DBD percent identity of one query motif's DBD to all motif DBDs in the database. 
# First compile the fasta files and make a temp file of paired query DBD to individual databse DBD.
# Then compute the percent idenity of those paired DBDs. Lastly, create a ordered list of aaids of 
# the query motif DBD to all database motif DBDs.

# QUERY_MOTIF=$1
# DIR_QUERY_DBDS=$2
# DIR_DATABASE_DBDS=$3

# DIR_PROJ=$HOME/proj_motifcomparison

# # create temp directories
# DIR_TEMP=$DIR_PROJ/output/cisbp_dbd_align/${QUERY_MOTIF}_temp
# mkdir -p $DIR_TEMP

# # create files of paired query DBD to individual datatbase DBD
# python $DIR_PROJ/scripts/generate_paired_dbd.py -m $QUERY_MOTIF -f1 $DIR_QUERY_DBDS -f2 $DIR_DATABASE_DBDS -o $DIR_TEMP

# # align and compute percent identiy of each dbd pair
# for filename_full in $DIR_TEMP/*.fasta; do
# 	filename_partial=$(basename $filename_full)
# 	# echo $filename_partial
# 	/usr/local/bin/clustalo --infile ${filename_full} --outfile $DIR_TEMP/${QUERY_MOTIF}_${filename_partial}.out \
# 	--seqtype protein --distmat-out $DIR_TEMP/${QUERY_MOTIF}_${filename_partial}.pim --full --percent-id --force
# done

# # parse all percent identity matrix files, and output a ordered list of aadis
# python $DIR_PROJ/scripts/parse_paired_aaid.py -i $DIR_TEMP -o $DIR_PROJ/output/cisbp_dbd_align/$QUERY_MOTIF.aaid

# rm -r $DIR_TEMP

QUERY_MOTIF=$1
DIR_QUERY_DBDS=$2

DIR_PROJ=$HOME/proj_motifcomparison

# create temp directories
DIR_TEMP=$DIR_PROJ/output/tf_dbd_align/${QUERY_MOTIF}_temp
mkdir -p $DIR_TEMP

# create files of paired query DBD to individual datatbase DBD
python $DIR_PROJ/scripts/generate_paired_dbd.py -m $QUERY_MOTIF -f1 $DIR_QUERY_DBDS -f2 $DIR_QUERY_DBDS -o $DIR_TEMP

# align and compute percent identiy of each dbd pair
for filename_full in $DIR_TEMP/*.fasta; do
	filename_partial=$(basename $filename_full)
	# echo $filename_partial
	/usr/local/bin/clustalo --infile ${filename_full} --outfile $DIR_TEMP/${QUERY_MOTIF}_${filename_partial}.out \
	--seqtype protein --distmat-out $DIR_TEMP/${QUERY_MOTIF}_${filename_partial}.pim --full --percent-id --force
done

# parse all percent identity matrix files, and output a ordered list of aadis
python $DIR_PROJ/scripts/parse_paired_aaid.py -i $DIR_TEMP -o $DIR_PROJ/output/tf_dbd_align/$QUERY_MOTIF.aaid

rm -r $DIR_TEMP
