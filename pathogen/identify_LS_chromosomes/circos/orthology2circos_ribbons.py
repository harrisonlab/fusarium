#!/usr/bin/python

'''
This program is used to identify single copy orthologs between two species in
an ortholoMCL orthology analysis. It then takes these single copy orthologs
and extracts gene locations from gff files. It outputs data in a format allowing
synteny plots to be drawn in circos.

'''

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument('--orthology',required=True,type=str,help='The text file containing orthogroups from OrthoMCL')
ap.add_argument('--name1',required=True,type=str,help='The name of organism 1 as was used in orthoMCL analysis')
ap.add_argument('--gff1',required=True,type=str,help='A gff file of genes from organism 1')
ap.add_argument('--name2',required=True,type=str,help='The name of organism 2 as was used in orthoMCL analysis')
ap.add_argument('--gff2',required=True,type=str,help='A gff file of genes from organism 2')

conf = ap.parse_args()

with open(conf.orthology) as f:
    orthology_lines = f.readlines()

with open(conf.gff1) as f:
    gff1_lines = f.readlines()

with open(conf.gff2) as f:
    gff2_lines = f.readlines()

strain1 = conf.name1
strain2 = conf.name2

gff1_dict = defaultdict(list)
gff2_dict = defaultdict(list)

#-----------------------------------------------------
# Step 2
# Store gene locations for organism 1 in a dictionary
#-----------------------------------------------------

for line in gff1_lines:
    if line.startswith('#'):
        pass
    else:
        line = line.rstrip()
        splitline = line.split("\t")
        if splitline[2] == 'mRNA' or splitline[2] == 'transcript':
            Col9 = splitline[8]
            split_col9 = Col9.split(';');
            ID_part = [ attribute for attribute in split_col9 if 'ID' in attribute ]
            ID = ID_part[0]
            ID_split = ID.split('=')
            genename = ID_split[-1]
            genename.replace('transcript:','')
            gff1_dict[genename] = [splitline[0], splitline[3], splitline[4]]

#-----------------------------------------------------
# Step 2
# Store gene locations for organism 2 in a dictionary
#-----------------------------------------------------

for line in gff2_lines:
    if line.startswith('#'):
        pass
    else:
        line = line.rstrip()
        splitline = line.split("\t")
        if splitline[2] == 'mRNA' or splitline[2] == 'transcript':
            Col9 = splitline[8]
            split_col9 = Col9.split(';');
            ID_part = [ attribute for attribute in split_col9 if 'ID' in attribute ]
            ID = ID_part[0]
            ID_split = ID.split('=')
            genename = ID_split[-1]
            genename = genename.replace('transcript:','')
            # print genename
            gff2_dict[genename] = [splitline[0], splitline[3], splitline[4]]

#-----------------------------------------------------
# Step 2
# Identify single copy genes in orthogroups
#-----------------------------------------------------


strain1_OrthoMCL = strain1 + '\|'
strain2_OrthoMCL = strain2 + '\|'
# print strain1_OrthoMCL

for line in orthology_lines:
    line = line.rstrip()
    # split_line = line.split("\t")
    matches1 = re.finditer(strain1_OrthoMCL, line)
    matches2 = re.finditer(strain2_OrthoMCL, line)
    # matches1 = re.finditer(strain1 + '\|', line)
    # matches2 = re.finditer(strain2 + '\|', line)
    # print matches1
    num_hits1 = sum(1 for m in matches1)
    num_hits2 = sum(1 for m in matches2)
    # print num_hits1
    # print num_hits1
    if ( num_hits1 == 1 and num_hits2 == 1 ):
        split_line = line.split(" ")
        orthogroup = split_line[0]
        # print orthogroup
        # print line
        for gene in split_line:
            if strain1 + "|" in gene:
                # print gene
                gene = gene.replace(strain1 + "|",'')
                # strain1_gene_set.add(gene)
                # print gff1_dict[gene]
                outline_part1 = gff1_dict[gene]
                outline_part1[0] = strain1 + "_" + outline_part1[0]
            elif strain2 + "|" in gene:
                # print gene
                gene = gene.replace(strain2 + "|",'').replace("-p1", "")
                # print gene
                # print gff2_dict[gene]
                outline_part2 = gff2_dict[gene]
                # print outline_part2
                outline_part2[0] = strain2 + "_" + outline_part2[0]
                # strain2_gene_set.add(gene)
        # outline =[ strain1 + "_" + outline_part1[0], outline_part1[1], outline_part1[2], strain2 + "_" + outline_part2[0], outline_part2[1], outline_part2[2] ]
        # print "\t".join(outline)
        print "\t".join(outline_part1) + "\t" + "\t".join(outline_part2)
        # gene = g in split_line if strain1 in g
            # print (gene + "\n")
        # print orthogroup
        # single_copy_orthogroups.add(orthogroup)

#-----------------------------------------------------
# Step 3
# Import variables & load input files
#-----------------------------------------------------
