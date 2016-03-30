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
ap.add_argument('--FoC_genes_gff',required=True,type=str,help='A gff file of the genes from the FoC')
conf = ap.parse_args()


with open(conf.blast_csv) as f:
    blast_csv_lines = f.readlines()

with open(conf.FoL_intersected_genes) as f:
    FoL_intersect_lines = f.readlines()

with open(conf.FoC_genes_gff) as f:
    FoC_genes_lines = f.readlines()

column_list=[]

#-----------------------------------------------------
# Step 2
# Read Blast csv file and store in a dictionary of hits.
#-----------------------------------------------------

blast_id_set = Set([])
i=0
blast_dict = defaultdict(list)
for line in blast_csv_lines:
    line = line.rstrip()
    split_line = line.split()
    if "ID" in split_line[0]:
        # print line
        for column in split_line:
            if "Grp" in column:
                i+=1
        # print i
        useful_start=i+1
        hit_contig=i+4
        # print hit_contig
        hit_stand=i+9
        hit_start=i+10
        hit_end=i+11
        # extract_list=split_lines[hit_contig, hit_start, hit_end]
        extract_list=itemgetter(hit_contig, hit_start, hit_end)(split_line)
        # print extract_list
        continue
    # print line
    blast_id=split_line[0]
    blast_id_set.add(blast_id)
    if len(split_line) > hit_contig:
        column_list=itemgetter(hit_contig, hit_start, hit_end, hit_stand)(split_line)
    # useful_columns=split_line[useful_start:]
    # print useful_columns
    # blast_dict[blast_id].append(blast_id)
    for column in column_list:
        blast_dict[blast_id].append(column)

#-----------------------------------------------------
# Step 2
# Build a dictionary of intersected genes
#-----------------------------------------------------

# blast_id_set = Set([])

intersect_dict = defaultdict(list)
for line in FoL_intersect_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    # print line
    # print split_line[11]
    if "gene" in split_line[11]:
        # print line
        blast_id = split_line[8].strip('";')
        blast_id = blast_id.replace('ID=', '').replace('_BlastHit_1', '')
        # blast_id_set.add(blast_id)
        # print blast_id
        column_list=itemgetter(12, 13, 15, 17)(split_line)
        for column in column_list:
            intersect_dict[blast_id].append(column)

#-----------------------------------------------------
# Step 3
# Append co-ordinates from the FoC genome, showing
# source of original Blast queries
#-----------------------------------------------------

FoC_genes_dict = defaultdict(list)
for line in FoC_genes_lines:
    if "gff-version" in line:
        continue
    line = line.rstrip()
    split_line = line.split()
    gene_id=split_line[8]
    if gene_id in blast_id_set:
        column_list=itemgetter(0, 3, 4, 6)(split_line)
        for column in column_list:
            FoC_genes_dict[gene_id].append(column)


#---
print ("\t".join(["query_id", "FoC_contig", "FoC_gene_start", "FoC_gene_end", "FoC_gene_strand", "hit_FoL_contig", "hit_start", "hit_end", "hit_strand", "FoL_gene_start", "FoL_gene_end", "FoL_strand", "FoL_gene_ID"]))

for blast_id in blast_id_set:
    useful_columns=[blast_id]
    useful_columns.extend(FoC_genes_dict[blast_id])
    useful_columns.extend(blast_dict[blast_id])
    useful_columns.extend(intersect_dict[blast_id])
    print ("\t".join(useful_columns))
