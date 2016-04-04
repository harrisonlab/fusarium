#!/usr/bin/python

'''
These commands were used to parse results from Blast searches of FoC pathogen
unique genes against the FoL genome as well as BLast searches of known SIX genes
against the FoL genome. These commands take information on location of blast
hits & suppliment this information with information from other files detailing
genes intersected in the target genome as well as information on the query gene
such as its location in the queries' genome.
'''


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
from sets import Set
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()
ap.add_argument('--blast_csv',required=True,type=str,help='The blast_pipe.sh results file')
ap.add_argument('--FoL_intersected_genes',required=True,type=str,help='A bed file of FoL genes intersecting Blast hits')
# ap.add_argument('--FoC_interescted_reblast',required=True,type=str,help='A bed file of FoC genes intersecting the location of reciprocal blast hits')

conf = ap.parse_args()


with open(conf.blast_csv) as f:
    blast_csv_lines = f.readlines()

with open(conf.FoL_intersected_genes) as f:
    FoL_intersect_lines = f.readlines()

# with open(conf.FoC_interescted_reblast) as f:
#     FoC_reblast_lines = f.readlines()

column_list=[]

#-----------------------------------------------------
# Step 2
# Read Blast csv file and store in a dictionary of hits.
#-----------------------------------------------------
# The csv file may contain columns showing relationship
# between query sequences. The first line of the file
# contains column headers. The number of column headers
# containg 'Grp' can be counted and this number of
# columns skipped when reading lines.
#
# Not all columns will contain blast hits. The number
# of columns in a line must be counted to check if a
# blast hit is present.


blast_id_set = Set([])
i=0
blast_dict = defaultdict(list)
for line in blast_csv_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    if "ID" in split_line[0]:
        for column in split_line:
            if "Grp" in column:
                i+=1
        hit_contig=i+4
        hit_stand=i+9
        hit_start=i+10
        hit_end=i+11
        extract_list=itemgetter(hit_contig, hit_start, hit_end)(split_line)
        continue
    blast_id=split_line[0]
    blast_id_set.add(blast_id)
    column_list = ["", "", "", ""]
    # print(str(len(split_line)) + "-" + str(hit_contig))
    if len(split_line) > hit_contig:
        # print "pass"
        column_list=itemgetter(hit_contig, hit_start, hit_end, hit_stand)(split_line)
    for column in column_list:
        blast_dict[blast_id].append(column)

#-----------------------------------------------------
# Step 3
# Build a dictionary of intersected genes
#-----------------------------------------------------

intersect_dict = defaultdict(list)
for line in FoL_intersect_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    if "gene" in split_line[11]:
        blast_id = split_line[8].strip('";')
        blast_id = blast_id.replace('ID=', '').replace('_BlastHit_1', '')
        column_list=itemgetter(12, 13, 15, 17)(split_line)
        for column in column_list:
            intersect_dict[blast_id].append(column)

#-----------------------------------------------------
# Step 4
# Build a dictionary of FoC genes intersecting the hit
# location of reciprocal blast queries.
#-----------------------------------------------------
#
# FoC_reblast_dict = defaultdict(list)
# for line in FoC_reblast_lines:
#     line = line.rstrip()
#     split_line = line.split("\t")
#     if "transcript" in split_line[11]:
#         query_id = split_line[8].strip('";')
#         query_id = query_id.replace('ID=', '').replace('_BlastHit_1', '').replace('extracted_hit_of_', '')
#         column_list = []
#         # column_list=itemgetter(12, 13, 15, 17)(split_line)
#         intersect_id = split_line[17]
#         # for column in column_list:
#             # FoC_reblast_dict[query_id].append(column)
#             # print column
#         # print (intersect_id + " - " + query_id)
#         if query_id in intersect_id:
#             # print "match"
#             FoC_reblast_dict[query_id]=[intersect_id, "match"]
#         # Bedtools may intersect a number of transcripts with the blast hit.
#         # This keeps the best transcript.
#         elif "match" in FoC_reblast_dict[query_id]:
#             pass
#         else:
#             # print "nomatch"
#             FoC_reblast_dict[query_id]=[intersect_id, "no_match"]
#
#

#-----------------------------------------------------
# Step 5
# Print final table of information on query, blast
# results and genes intersecting blast results
#-----------------------------------------------------

print ("\t".join(["query_id", "hit_FoL_contig", "hit_start", "hit_end", "hit_strand", "FoL_gene_start", "FoL_gene_end", "FoL_strand", "FoL_gene_ID"]))

for blast_id in blast_id_set:
    useful_columns=[blast_id]
    # useful_columns.extend(FoC_genes_dict[blast_id])
    useful_columns.extend(blast_dict[blast_id])
    # if FoC_reblast_dict[blast_id]:
    #     useful_columns.extend(FoC_reblast_dict[blast_id])
    # else:
    #     useful_columns.extend(["", ""])
    if intersect_dict[blast_id]:
        useful_columns.extend(intersect_dict[blast_id])
    else:
        # useful_columns.extend([".", ".", ".", "no gene"])
        useful_columns.extend(["", "", "", ""])
    # print ("\t".join(useful_columns))
    print useful_columns
