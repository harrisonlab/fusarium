#!/usr/bin/python

'''
This program is used to identify orthologs present on a given chromomse in
organism1, within its own genome and within the genome of organism2. It then
takes these orthologs and extracts gene locations from gff files. It outputs
data in a format allowing synteny plots to be drawn in circos.
'''

#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
from collections import defaultdict
from sets import Set

ap = argparse.ArgumentParser()
ap.add_argument('--chr1',required=True,type=str,nargs='+',help='The contigs of organism2 from which orthogroups should be extracted.')
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

contig_list = conf.chr1

gff1_dict = defaultdict(list)
gff2_dict = defaultdict(list)
specified_genes_dict = defaultdict(list)

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
            genename = genename.replace('transcript:','').replace('T0','')
            gff1_dict[genename] = [splitline[0], splitline[3], splitline[4]]
            specified_contig = False
            # print 'sweet'
            specified_contig = [ x for x in contig_list if x == splitline[0] ]
            if len(specified_contig) >= 1 :
                specified_genes_dict[genename] = [splitline[0], splitline[3], splitline[4]]
                # print genename
                # print 'dude'

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
            genename = genename.replace('transcript:','').replace('T0','')
            gff2_dict[genename] = [splitline[0], splitline[3], splitline[4]]

#-----------------------------------------------------
# Step 2
# Identify single copy genes in orthogroups
#-----------------------------------------------------


strain1_OrthoMCL = strain1 + '\|'
strain2_OrthoMCL = strain2 + '\|'

pairset1 = Set([])
pairset2 = Set([])

for line in orthology_lines:
    line = line.rstrip()
    # print line
    if strain1 + "|" in line:
        # print line
        split_line = line.split(" ")
        orthogroup = split_line.pop(0)
        # print orthogroup
        for gene in split_line:
            # print gene
            if strain1 + "|" in gene:
                gene = gene.replace(strain1 + "|", "")
                # print gene
                if specified_genes_dict[gene]:
                    outline_part1 = specified_genes_dict[gene]
                    outline_part1[0] = strain1 + "_" + outline_part1[0]
                    for partner in split_line:
                        if (( strain1 + "|" in partner
                            or strain2 + "|" in partner)
                            and
                            ( gene not in partner )
                            ):
                            if strain1 + "|" in partner:
                                partner = partner.replace(strain1 + "|", "")
                                pair_list = [gene, partner]
                                pair_list.sort()
                                pair = "-".join(pair_list)
                                if pair in pairset1:
                                    pass
                                else:
                                    pairset1.add(pair)
                                    outline_part2 = gff1_dict[partner]
                                    outline_part2[0] = strain1 + "_" + outline_part2[0]
                                    print "\t".join(outline_part1) + "\t" + "\t".join(outline_part2)
                                    # print pair
                            elif strain2 + "|" in partner:
                                # print "badgers"
                                partner = partner.replace(strain2 + "|", "")
                                pair_list = [gene, partner]
                                pair_list.sort()
                                pair = "-".join(pair_list)
                                if pair in pairset2:
                                    pass
                                else:
                                    pairset2.add(pair)
                                    outline_part2 = gff2_dict[partner]
                                    outline_part2[0] = strain2 + "_" + outline_part2[0]
                                    print "\t".join(outline_part1) + "\t" + "\t".join(outline_part2)
                                    # print pair

                    # print "awesome"

#
#
#     matches1 = re.finditer(strain1_OrthoMCL, line)
#     matches2 = re.finditer(strain2_OrthoMCL, line)
#     num_hits1 = sum(1 for m in matches1)
#     num_hits2 = sum(1 for m in matches2)
#     if ( num_hits1 >= 1 and num_hits2 >= 1 ):
#         split_line = line.split(" ")
#         orthogroup = split_line[0]
#         # print orthogroup
#         # print line
#         for gene in split_line:
#             if strain1 + "|" in gene:
#                 # print gene
#                 gene = gene.replace(strain1 + "|",'')
#                 # strain1_gene_set.add(gene)
#                 # print gff1_dict[gene]
#                 outline_part1 = gff1_dict[gene]
#                 outline_part1[0] = strain1 + "_" + outline_part1[0]
#             elif strain2 + "|" in gene:
#                 # print gene
#                 gene = gene.replace(strain2 + "|",'')
#                 # print gff2_dict[gene]
#                 outline_part2 = gff2_dict[gene]
#                 # print outline_part2
#                 outline_part2[0] = strain2 + "_" + outline_part2[0]
#                 # strain2_gene_set.add(gene)
#         # outline =[ strain1 + "_" + outline_part1[0], outline_part1[1], outline_part1[2], strain2 + "_" + outline_part2[0], outline_part2[1], outline_part2[2] ]
#         # print "\t".join(outline)
#         print "\t".join(outline_part1) + "\t" + "\t".join(outline_part2)
#         # gene = g in split_line if strain1 in g
#             # print (gene + "\n")
#         # print orthogroup
#         # single_copy_orthogroups.add(orthogroup)
#
# #-----------------------------------------------------
# # Step 3
# # Import variables & load input files
# #-----------------------------------------------------
