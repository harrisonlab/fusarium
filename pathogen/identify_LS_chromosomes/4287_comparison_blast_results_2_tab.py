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
ap.add_argument('--FoC_interescted_reblast',required=True,type=str,help='A bed file of FoC genes intersecting the location of reciprocal blast hits')
ap.add_argument('--FoC_SigP',required=True,type=str,help='A file containing a list of signal-peptide containing genes')
ap.add_argument('--FoC_TM_list',required=True,type=str,help='A file containing a list of transmembrane containing genes')
ap.add_argument('--FoC_MIMP_list',required=True,type=str,help='A file containing a list of genes within 2kb of MIMPs')
ap.add_argument('--FoC_effectorP',required=True,type=str,help='A file containing results of effectorP')


# ap.add_argument('--FoC_expression',required=True,type=str,help='A file containing details of gene expression')

conf = ap.parse_args()


with open(conf.blast_csv) as f:
    blast_csv_lines = f.readlines()

with open(conf.FoL_intersected_genes) as f:
    FoL_intersect_lines = f.readlines()

with open(conf.FoC_genes_gff) as f:
    FoC_genes_lines = f.readlines()

with open(conf.FoC_interescted_reblast) as f:
    FoC_reblast_lines = f.readlines()

with open(conf.FoC_SigP) as f:
    FoC_SigP_lines = f.readlines()

with open(conf.FoC_TM_list) as f:
    FoC_tmhmm_lines = f.readlines()

with open(conf.FoC_MIMP_list) as f:
    FoC_mimp_lines = f.readlines()

with open(conf.FoC_effectorP) as f:
    FoC_effectorP_lines = f.readlines()

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
        # extract_list=itemgetter(hit_contig, hit_start, hit_end)(split_line)
        continue
    blast_id=split_line[0]
    blast_id_set.add(blast_id)
    # column_list = ["no hit", ".", ".", "."]
    column_list = ["", "", "", ""]
    if len(split_line) > hit_contig:
        column_list=itemgetter(hit_contig, hit_start, hit_end, hit_stand)(split_line)
    for column in column_list:
        blast_dict[blast_id].append(column)

#-----------------------------------------------------
# Step 3
# Build a dictionary of intersected genes
# The final column of the input file containing
# gff annotation information was split into two columns
#-----------------------------------------------------

intersect_dict = defaultdict(list)
for line in FoL_intersect_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    if "gene" in split_line[11]:
        blast_id = split_line[8].strip('";')
        blast_id = blast_id.replace('ID=', '').replace('_BlastHit_1', '')
        column_list = itemgetter(12, 13, 15)(split_line)
        for column in column_list:
            intersect_dict[blast_id].append(column)

        feature_info_list = "".join(split_line[17]).split(";")
        FoL_gene_ID = ""
        FoL_feature_desc = ""
        for feature in feature_info_list:
            if "gene_id=" in feature:
                FoL_gene_ID = re.sub(r"^.*=", '', feature)
            if "description=" in feature:
                FoL_feature_desc = re.sub(r"^.*=", '', feature)
        intersect_dict[blast_id].extend([FoL_gene_ID, FoL_feature_desc])

#-----------------------------------------------------
# Step 4
# Append co-ordinates from the FoC genome, showing
# source of original Blast queries
#-----------------------------------------------------

FoC_genes_dict = defaultdict(list)
for line in FoC_genes_lines:
    if "gff-version" in line:
        continue
    line = line.rstrip()
    split_line = line.split("\t")
    gene_id=split_line[8]
    column_list = ["", "", "", ""]
    if gene_id in blast_id_set:
        column_list=itemgetter(0, 3, 4, 6)(split_line)
        for column in column_list:
            FoC_genes_dict[gene_id].append(column)

#-----------------------------------------------------
# Step 5
# Build a dictionary of FoC genes intersecting the hit
# location of reciprocal blast queries.
#-----------------------------------------------------

FoC_reblast_dict = defaultdict(list)
for line in FoC_reblast_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    if "transcript" in split_line[11]:
        query_id = split_line[8].strip('";')
        query_id = query_id.replace('ID=', '').replace('_BlastHit_1', '').replace('extracted_hit_of_', '')
        column_list = []
        intersect_id = split_line[17]
        if query_id in intersect_id:
            FoC_reblast_dict[query_id]=[intersect_id, "match"]
        elif "match" in FoC_reblast_dict[query_id]:
            pass
        else:
            FoC_reblast_dict[query_id]=[intersect_id, ""]

#-----------------------------------------------------
# Steps 6 & 7
# Build a dictionary of FoC genes containing Signal
# peptides
#
# Extend this dictionary with information on FoC genes
# containing transmembrane domains.
#
# If a protein contains a signal peptide and does not
# contain a transmembrane domain, it is treated as a
# putative secreted protein.
#-----------------------------------------------------

FoC_secreted_dict = defaultdict(list)
for line in FoC_SigP_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0].replace('>', '')
    FoC_secreted_dict[gene_id] = split_line[1:]

for line in FoC_tmhmm_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0]
    num_helices = split_line[4].replace('PredHel=', '')
    if num_helices == '0':
        FoC_secreted_dict[gene_id].extend(["NO", num_helices])
        if FoC_secreted_dict[gene_id][0] == "YES":
            FoC_secreted_dict[gene_id].append("Yes")
        else:
            FoC_secreted_dict[gene_id].append("")
    else:
        FoC_secreted_dict[gene_id].extend(["YES", num_helices, ""])

#-----------------------------------------------------
# Step 8
# Build a dictionary of FoC genes within 2Kb of MIMPs
#-----------------------------------------------------

FoC_mimp_set =  Set([])

for line in FoC_mimp_lines:
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0]
    FoC_mimp_set.add(gene_id)

#-----------------------------------------------------
# Step 9
# Build a dictionary of FoC genes that have been
# predicted as effectors by effectorP
#-----------------------------------------------------

FoC_effectorP_dict = defaultdict(list)
First = True
for line in FoC_effectorP_lines:
    # print line
    if First == True:
        First = False
        continue
    line = line.rstrip("\n")
    split_line = line.split("\t")
    gene_id = split_line[0]
    # print gene_id
    # FoC_effectorP_dict[gene_id] = [split_line[1], split_line[2]]
    if gene_id in blast_id_set:
        split_line[1] = split_line[1].replace('Effector', 'Yes').replace("Non-effector", "")
        column_list = itemgetter(1, 2)(split_line)
    else:
        column_list = ["", ""]
    FoC_effectorP_dict[gene_id].extend(column_list)


#-----------------------------------------------------
# Step 10
# Print final table of information on query, blast
# results and genes intersecting blast results
#-----------------------------------------------------

print ("\t".join(["query_id", "FoC_contig", "FoC_gene_start", "FoC_gene_end", "FoC_gene_strand", "SigP", "P-value", "cleavage_site", "TransMem_protein", "No._helices", "Secreted", "MIMP_in_2Kb", "effectorP", "P-value", "hit_FoL_contig", "hit_start", "hit_end", "hit_strand", "reblast_hit", "reblast match", "FoL_gene_start", "FoL_gene_end", "FoL_strand", "FoL_gene_ID", "FoL_gene_description"]))

for blast_id in blast_id_set:
    useful_columns=[blast_id]
    useful_columns.extend(FoC_genes_dict[blast_id])
    useful_columns.extend(FoC_secreted_dict[blast_id])
    # useful_columns.extend(FoC_tmhmm_dict[blast_id])
    mimp_col=""
    if blast_id in FoC_mimp_set:
        mimp_col="Yes"
    useful_columns.append(mimp_col)
    useful_columns.extend(FoC_effectorP_dict[blast_id])
    useful_columns.extend(blast_dict[blast_id])
    if FoC_reblast_dict[blast_id]:
        useful_columns.extend(FoC_reblast_dict[blast_id])
    else:
        useful_columns.extend(["", ""])
    if intersect_dict[blast_id]:
        useful_columns.extend(intersect_dict[blast_id])
    else:
        useful_columns.extend(["", "", "", ""])
    print ("\t".join(useful_columns))
