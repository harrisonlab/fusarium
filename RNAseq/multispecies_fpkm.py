#!/usr/bin/python

'''
/home/armita/git_repos/emr_repos/scripts/fusarium/RNAseq/multispecies_fpkm.py --mummer_tsv analysis/genome_alignment/mummer/F.oxysporum/fo47/fo47_vs_Fus2/fo47_vs_Fus2_coords.tsv --ref_genes gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3 --query_genes assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_transcripts_parsed.gff3 --ref_fpkm_files alignment/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_72hrs_rep1/fpkm/genes.fpkm_tracking alignment/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_72hrs_rep2/fpkm/genes.fpkm_tracking alignment/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_72hrs_rep3/fpkm/genes.fpkm_tracking --query_fpkm_files alignment/F.oxysporum/fo47/FO47_72hrs_rep1/fpkm/genes.fpkm_tracking alignment/F.oxysporum/fo47/FO47_72hrs_rep2/fpkm/genes.fpkm_tracking alignment/F.oxysporum/fo47/FO47_72hrs_rep3/fpkm/genes.fpkm_tracking | less
'''

import sys,argparse
from collections import defaultdict
from sets import Set
import numpy as np

#######################################
#            Import variables         #
#                                     #
#                                     #
#######################################

ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)

ap.add_argument('--mummer_tsv',required=True,type=str,help='mummer alignment output in tsv')
ap.add_argument('--ref_genes',required=True,type=str,help='gff file of gene models from the reference organism')
ap.add_argument('--query_genes',required=True,type=str,help='gff file of gene models from the query organism')
ap.add_argument('--ref_fpkm_files',required=True,nargs= '+', type=str,help='cufflinks fpkm_tracking files giving fpkm by gene')
ap.add_argument('--query_fpkm_files',required=True,nargs= '+', type=str,help='cufflinks fpkm_tracking files giving fpkm by gene')



conf = ap.parse_args() #sys.argv

#
# fpkm_dict = defaultdict(list)
# CAZY_set = Set()
# effP_set = Set()
# transposon_set = Set()
# metabolite_set = Set()
ref_alignment_dict = defaultdict(list)
query_alignment_dict = defaultdict(list)
intervals_by_ref_contig_dict = defaultdict(list)
intervals_by_query_contig_dict = defaultdict(list)
relative_coord_dict = defaultdict(list)
ref_fpkm_dict = defaultdict(list)
query_fpkm_dict = defaultdict(list)

#######################################
#    Import fpkm for gene models      #
#                                     #
#                                     #
#######################################
#

for fpkm_file in conf.ref_fpkm_files:
    with open(fpkm_file) as f:
        for line in f:
            if line.startswith('tracking_id'):
                continue
            line = line.rstrip()
            split_line = line.split()
            gene_id = split_line[0]
            fpkm = split_line[9]
            ref_fpkm_dict[gene_id].append(fpkm)

for fpkm_file in conf.query_fpkm_files:
    with open(fpkm_file) as f:
        for line in f:
            if line.startswith('tracking_id'):
                continue
            line = line.rstrip()
            split_line = line.split()
            gene_id = split_line[0]
            fpkm = split_line[9]
            query_fpkm_dict[gene_id].append(fpkm)


#######################################
#    Import coordinates of Mummer     #
#               alignment             #
#                                     #
#######################################
# Store alignment lines in a dictionary
# ordered by contig name and by the
# start location of the alignment
# This allows easy extraction of relavent
# lines later.

with open(conf.mummer_tsv) as f:
    startreading = False
    for line in f:
        line = line.rstrip()
        if startreading == False and line.startswith('[S1]'):
            startreading = True
        elif startreading == False:
            continue
        else:
            split_line = line.split()
            ref_start = split_line[0]
            query_start = split_line[2]
            ref_contig = split_line[10]
            query_contig = split_line[11]
            key = ":".join([ref_contig, ref_start])
            ref_alignment_dict[key] = line
            key = ":".join([query_contig, query_start])
            query_alignment_dict[key] = line

def sort_alignment_dict_func(key):
    split_key = key.split('_')
    contig = split_key[1]
    start = key.split(':')[1]
    return int(contig), int(start)

keys = ref_alignment_dict.keys()
sorted_keys = sorted(keys, key=lambda x: sort_alignment_dict_func(x))

#######################################
#    Import reference Gff lines       #
#    and sort by contig               #
#                                     #
#######################################
# import reference gff lines and sort
# by contig and start location
# these are stored for later reference

gff_lines = []
with open(conf.ref_genes) as f:
    for line in f:
        line = line.rstrip()
        split_line = line.split()
        if line.startswith('#'):
            continue
        elif 'gene' in split_line[2]:
            gff_lines.append(line)

def sort_gff_func(line):
    split_line = line.split('\t')
    col1 = split_line[0]
    contig = col1.split('_')[1]
    start = split_line[3]
    return int(contig), int(start)

sorted_ref_gff_lines = sorted(gff_lines, key=lambda x: sort_gff_func(x))


#######################################
#    Build a dictionary of alignment  #
#     intervals, stored by contig     #
#######################################
# intervals for alignemnts to
# reference contigs are stored in a
# dictionary to allow easy identification
# of whether a reference gene is
# within an aligned region.

for key in sorted_keys:
    alignment_line = ref_alignment_dict[key]
    split_line = alignment_line.split()
    ref_start = split_line[0]
    ref_end = split_line[1]
    query_start = split_line[2]
    query_end = split_line[3]
    ref_length = split_line[4]
    query_length = split_line[5]
    ref_contig = split_line[10]
    query_contig = split_line[11]
    intervals_by_ref_contig_dict[ref_contig].append([ref_start, ref_end])
    intervals_by_query_contig_dict[query_contig].append([query_start, query_end])
# print(intervals_by_ref_contig_dict)


#
# # print sorted_keys
# for key in sorted_keys:
#     alignment_line = alignment_dict[key]
#     split_line = alignment_line.split()
#     ref_start = split_line[0]
#     ref_end = split_line[1]
#     query_start = split_line[2]
#     query_end = split_line[3]
#     ref_length = split_line[4]
#     query_length = split_line[5]
#     ref_contig = split_line[10]
#     query_contig = split_line[11]
#
#     print "\t".join([ref_contig, ref_start, ref_end, query_contig, query_start, query_end])


#######################################
#    Store genes contained within     #
#           Aligned regions           #
#######################################
# For each gene, the contig it is
# present on is identified. The
# appropriate list of alignment boundaries
# is extracted from a dictionary. If the
# gene's atrt and stop co-ordinates are
# within an aligned region gene is stored
# for later reference. It's relative
# co-ordinates to the start of the alignment
# are used for reference, as these allow
# comparison to genes from the query
# organism.

for ref_gff_line in sorted_ref_gff_lines:
    split_ref_gff_line = ref_gff_line.split()
    contig = split_ref_gff_line[0]
    gene_start = split_ref_gff_line[3]
    gene_end = split_ref_gff_line[4]
    features = split_ref_gff_line[8]
    split_features = features.split(';')
    gene_id = split_features[0]
    gene_id = gene_id.replace('ID=', '')
    # print "\t".join([contig, gene_start, gene_end, gene_id])
    interval_list = intervals_by_ref_contig_dict[contig]
    alignment_start = 'none'
    for (lower, upper) in interval_list:
        if int(lower) <= int(gene_start) <= int(upper) and int(lower) <= int(gene_end) <= int(upper):
            alignment_start = lower
            # print "\t".join([lower, upper])
            # print "\t".join([contig, gene_start, gene_end, gene_id])
            break
    if alignment_start == 'none':
        continue
    key = ":".join([contig, alignment_start])
    alignment_line = ref_alignment_dict[key]
    split_alignment_line = alignment_line.split()
    ref_start = split_alignment_line[0]
    ref_end = split_alignment_line[1]
    query_start = split_alignment_line[2]
    query_end = split_alignment_line[3]
    ref_length = split_alignment_line[4]
    query_length = split_alignment_line[5]
    ref_contig = split_alignment_line[10]
    query_contig = split_alignment_line[11]
    # print "\t".join([ref_contig, ref_start, ref_end, query_contig, query_start, query_end])
    # print "\t".join([contig, gene_start, gene_end, gene_id])
# ---
# Identify relative location of gene to alignment location
# ---
    relative_start = int(gene_start) - int(ref_start)
    relative_end = int(gene_end) - int(ref_start)
    # print "\t".join([contig, str(relative_start), str(relative_end), gene_id])
    # fpkm = 'unknown'
    fpkm = ":".join(ref_fpkm_dict[gene_id])
    key = "_".join([contig, str(relative_start), str(relative_end)])
    relative_coord_dict[key] = [contig, str(relative_start), str(relative_end), gene_id, fpkm]

# print relative_coord_dict




#######################################
#    Import query Gff lines       #
#    and sort by contig               #
#                                     #
#######################################
# import query gff lines and sort
# by contig and start location
# these are stored for later reference

query_gff_lines = []
with open(conf.query_genes) as f:
    for line in f:
        line = line.rstrip()
        split_line = line.split()
        if line.startswith('#'):
            continue
        elif 'mRNA' in split_line[2]:
        # elif 'gene' in split_line[2]:
            query_gff_lines.append(line)

def sort_query_gff_func(line):
    split_line = line.split('\t')
    col1 = split_line[0]
    contig = col1.split('.')[1]
    start = split_line[3]
    return int(contig), int(start)
    # split_line = line.split('\t')
    # col1 = split_line[0]
    # contig = col1.split('_')[1]
    # start = split_line[3]
    # return int(contig), int(start)

sorted_query_gff_lines = sorted(query_gff_lines, key=lambda x: sort_query_gff_func(x))



#######################################
#    Store query genes contained within     #
#           Aligned regions           #
#######################################
# For each query gene, the contig it is
# present on is identified. The
# appropriate list of alignment boundaries
# is extracted from a dictionary. If the
# gene's atrt and stop co-ordinates are
# within an aligned region gene is stored
# for later reference. It's relative
# co-ordinates to the start of the alignment
# are used for reference, as these allow
# comparison to genes from the query
# organism.

for gff_line in sorted_query_gff_lines:
    split_gff_line = gff_line.split()
    contig = split_gff_line[0]
    gene_start = split_gff_line[3]
    gene_end = split_gff_line[4]
    features = split_gff_line[8]
    split_features = features.split(';')
    gene_id = split_features[0]
    gene_id = gene_id.replace('ID=', '')
    # fo_47 gene_id's from the gff file contain transcript numbers
    gene_id = gene_id[:-2]
    # print gene_id
    # print "\t".join([contig, gene_start, gene_end, gene_id])
    interval_list = intervals_by_query_contig_dict[contig]
    alignment_start = 'none'
    for (lower, upper) in interval_list:
        if int(lower) <= int(gene_start) <= int(upper) and int(lower) <= int(gene_end) <= int(upper):
            alignment_start = lower
            # print "\t".join([lower, upper])
            # print "\t".join([contig, gene_start, gene_end, gene_id])
            break
    if alignment_start == 'none':
        continue
    key = ":".join([contig, alignment_start])
    alignment_line = query_alignment_dict[key]
    split_alignment_line = alignment_line.split()
    ref_start = split_alignment_line[0]
    ref_end = split_alignment_line[1]
    query_start = split_alignment_line[2]
    query_end = split_alignment_line[3]
    ref_length = split_alignment_line[4]
    query_length = split_alignment_line[5]
    ref_contig = split_alignment_line[10]
    query_contig = split_alignment_line[11]
    # print "\t".join([ref_contig, ref_start, ref_end, query_contig, query_start, query_end])
    # print "\t".join([contig, gene_start, gene_end, gene_id])
# ---
# Identify relative location of gene to alignment location
# ---
    relative_start = int(gene_start) - int(query_start)
    relative_end = int(gene_end) - int(query_start)
    # print "\t".join([contig, str(relative_start), str(relative_end), gene_id])
    # fpkm = 'unknown'
    fpkm = ":".join(query_fpkm_dict[gene_id])
    key = "_".join([ref_contig, str(relative_start), str(relative_end)])
    # relative_coord_dict[key] = [gene_id, fpkm]
    if relative_coord_dict[key]:
        # print "\t".join([ref_contig, ref_start, ref_end, query_contig, query_start, query_end])
        print "\t".join([contig, str(relative_start), str(relative_end), gene_id, fpkm]) + "\t" + "\t".join(relative_coord_dict[key])
