#!/usr/bin/python

"""
Query CIS-BP MySQL database for inferred motif PWMs and DNA binding domain sequences.
"""

import sys
import os
import argparse
import mysql.connector
import unicodedata
import numpy

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Query the cisbp mysql database")
    parser.add_argument('-o', '-fn_output', dest='fn_output', type=str)
    parsed = parser.parse_args(argv[1:])
    return parsed

def main(argv):
    parsed = parse_args(argv)

    """ connection to mysql """ 
    cnx = mysql.connector.connect(user="root", password="root", database="cis_bp")
    cursor = cnx.cursor()

    """ query of all motif dbd sequences """ 
    # # query = ("SELECT Motif_ID, MotifFeature_Sequence FROM motif_features WHERE Motif_ID IN \
    # #     (SELECT Motif_ID FROM motif_all WHERE Evidence = 'I' and Species != 'Saccharomyces_cerevisiae')")
    # query = ("SELECT * FROM motif_features")
    # cursor.execute(query)

    # writer = open(parsed.fn_output, "w")
    # for item in cursor:
    #     temp_id = unicodedata.normalize('NFKD', item[1]).encode('ascii', 'ignore') 
    #     temp_from_pos = item[3]
    #     temp_to_pos = item[4] 
    #     temp_seq = unicodedata.normalize('NFKD', item[5]).encode('ascii', 'ignore')  
    #     writer.write(">%s:%d-%d\n%s\n" % (temp_id, temp_from_pos, temp_to_pos, temp_seq))
    # writer.close()  

    """ query of tfs with no direct evidence, associated with inferred motifs """
    # tf_info_dict = {}
    # fam_dict = {}
    # tf_motif_dict = {}

    # query_id = ("SELECT TF_ID, TF_Name, DBID, TF_Species, Family_ID FROM tfs WHERE TF_Status = 'I'")
    # query_dbd = ("SELECT Family_ID, DBDs FROM tf_families WHERE Family_ID IN \
    #     (SELECT Family_ID FROM tfs WHERE TF_Status = 'I')")
    # query_motif = ("SELECT TF_ID, Motif_ID FROM motif_all WHERE TF_ID IN \
    #     (SELECT TF_ID FROM tfs WHERE TF_Status = 'I')")

    # cursor.execute(query_dbd)
    # for item in cursor:
    #     temp_famid = unicodedata.normalize('NFKD', item[0]).encode('ascii', 'ignore') 
    #     temp_dbd = unicodedata.normalize('NFKD', item[1]).encode('ascii', 'ignore') 
    #     fam_dict[temp_famid] = temp_dbd

    # cursor.execute(query_id)
    # for item in cursor:
    #     temp_id = unicodedata.normalize('NFKD', item[0]).encode('ascii', 'ignore')
    #     temp_name = unicodedata.normalize('NFKD', item[1]).encode('ascii', 'ignore') 
    #     temp_dbid = unicodedata.normalize('NFKD', item[2]).encode('ascii', 'ignore') 
    #     temp_spec = unicodedata.normalize('NFKD', item[3]).encode('ascii', 'ignore') 
    #     temp_famid = unicodedata.normalize('NFKD', item[4]).encode('ascii', 'ignore')  
    #     tf_info_dict[temp_id] = [temp_name, temp_dbid, temp_spec, fam_dict[temp_famid]]

    # cursor.execute(query_motif)
    # for item in cursor:
    #     temp_id = unicodedata.normalize('NFKD', item[0]).encode('ascii', 'ignore') 
    #     temp_motif = unicodedata.normalize('NFKD', item[1]).encode('ascii', 'ignore') 
    #     if temp_id not in tf_motif_dict:
    #         temp_list = [temp_motif]
    #     else:
    #         temp_list = tf_motif_dict[temp_id]
    #         temp_list.append(temp_motif)
    #     tf_motif_dict[temp_id] = temp_list

    # writer = open(parsed.fn_output, "w")
    # tf_key_sorted = sorted(tf_info_dict.keys())
    # for tf in tf_key_sorted:
    #     temp_info = tf_info_dict[tf]
    #     temp_motifs = tf_motif_dict[tf]
    #     writer.write(">%s: %s, %s, %s, %s\n" % (tf, temp_info[0], temp_info[1], temp_info[2], temp_info[3]))
    #     for i, temp_motif in enumerate(temp_motifs):
    #         writer.write("%s" % temp_motif)
    #         if i != len(temp_motifs)-1:
    #             writer.write(", ")
    #         else:
    #             writer.write("\n")
    # writer.close()

    """ query all species """
    # query = ("SELECT DISTINCT(Species), Count(*) FROM motif_all WHERE Evidence = 'D' GROUP BY Species")
    # cursor.execute(query)
    # writer = open(parsed.fn_output, "w")
    # for item in cursor:
    #     spec = unicodedata.normalize('NFKD', item[0]).encode('ascii', 'ignore')
    #     count = item[1]
    #     print spec
    #     writer.write("%s\t%d\n" % (spec, count))
    # writer.close()

    """ query all motifs have direct evidence and eukaryote but not yeast """
    # lines = open("/Users/KANG/proj_motifcomparison/resources/cisbp_all_species.txt", "r").readlines()
    # dict_spec = {}
    # for i, line in enumerate(lines):
    #     if i > 0:
    #         temp_spec = line.split()[0]
    #         temp_cell = line.split()[1]
    #         dict_spec[temp_spec] = temp_cell

    # query = ("SELECT Motif_ID, Species FROM motif_all WHERE Evidence = 'D'")
    # cursor.execute(query)
    # motifs = []
    # species = []
    # for item in cursor:
    #     temp_motif = unicodedata.normalize('NFKD', item[0]).encode('ascii', 'ignore')
    #     temp_species = unicodedata.normalize('NFKD', item[1]).encode('ascii', 'ignore')
    #     if temp_species != 'Saccharomyces_cerevisiae' and dict_spec[temp_species] == 'eukary':
    #         motifs.append(temp_motif)
    #         species.append(temp_species)
    # data = numpy.column_stack((motifs, species))
    # data = data[numpy.argsort(data[:,0])]

    # writer = open(parsed.fn_output, "w")
    # for temp_data in data:
    #     writer.write("%s\t%s\n" % (temp_data[0], temp_data[1]))
    # writer.close()

    """ query all motifs without the ones with direct yeast evidence """
    query = ("SELECT DISTINCT(Motif_ID) FROM motifs WHERE TF_ID IN \
        (SELECT TF_ID FROM tfs WHERE TF_Species != 'Saccharomyces_cerevisiae')")
    cursor.execute(query)
    writer = open(parsed.fn_output, "w")
    for item in cursor:
        temp_motif = unicodedata.normalize('NFKD', item[0]).encode('ascii', 'ignore')
        writer.write("%s\n" % temp_motif)
    writer.close()

    """ query the species with more than 5 direct motif evidence """
    # query = ("SELECT Species, Spec_Count FROM \
    #     (SELECT Species, count(*) AS Spec_Count FROM motif_all WHERE Evidence = 'D' GROUP BY Species) \
    #     AS spec_dir_evid WHERE spec_dir_evid.Spec_Count >= 5")
    # cursor.execute(query)
    # writer = open(parsed.fn_output, "w")
    # for item in cursor:
    #     spec = unicodedata.normalize('NFKD', item[0]).encode('ascii', 'ignore')
    #     writer.write("%s\n" % spec)
    # writer.close()

    """ close mysql connection """
    cursor.close()
    cnx.close()

if __name__ == "__main__":
    main(sys.argv)
