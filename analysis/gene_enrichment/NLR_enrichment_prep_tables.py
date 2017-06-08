#!/usr/bin/python

'''
This program extracts Interproscan domains, superfamily domains and Pfam
domains associated with NLR proteins for genes and builds Fishcers
contingeny tables showing numbers of annotated genes on specified contigs
in comparison to other specified contigs
'''


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
import numpy
from sets import Set
# from Bio import SeqIO
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()
ap.add_argument('--interpro',required=True,type=str,help='Interproscan annotations in tsv')
ap.add_argument('--set1_name',required=True,type=str,help='Name for contigs in set 1')
ap.add_argument('--contig_set1',required=True,nargs='+',type=str,help='List of contigs in set 1')
ap.add_argument('--set2_name',required=True,type=str,help='Name for contigs in set 2')
ap.add_argument('--contig_set2',required=True,nargs='+',type=str,help='List of contigs in set 2')
ap.add_argument('--outdir',required=True,type=str,help='Output directory for results')

conf = ap.parse_args()

with open(conf.interpro) as f:
    interpro_lines = f.readlines()


set1_name = conf.set1_name
set1_contigs = set(conf.contig_set1)
set2_name = conf.set2_name
set2_contigs = set(conf.contig_set2)

print set1_name
print set1_contigs
print set2_name
print set2_contigs

# gff1_dict = defaultdict(list)
# specified_genes_dict = defaultdict(list)

#-----------------------------------------------------
# Step 2
#
#-----------------------------------------------------

pfam_domain_dict = {
# N-terminal domains:
'PF00168' : "C2",
'PF17109' : "Goodbye",
'PF14479': "Helo",
'PF17111' : "Helo-like",
'PF06985' : "HET",
'PF01734' : "Patatin",
'PF00082' : "Peptidase_S8",
'PF11558' : "PFD",
'PF00069' : "PKinase",
'PF01048' : "PNP_UDP",
'PF04607' : "RelA_Spot",
'PF17107' : 'SesA',
'PF17046' : 'sesB-like',
# Nucleotide binding domains:
'PF00931' : 'NB-ARC',
'PF05729' : 'NACHT',
'SSF52540' : 'P-loop_NTPase', # P-loop NTPase superfamily (CL0023)
'IPR027417' : 'P-loop_NTPase',
'PF13401' : 'AAA',
# C-terminal repeat domains:
"SSF48452" : 'TPR', # TPR superfamily (CL0020)
"PF00515" : 'TPR', # TPR1
"PF07719" : 'TPR', # TPR2
"PF07720" : 'TPR', # TPR3
"PF07721" : 'TPR', # TPR4
"PF12688" : 'TPR', # TPR5
"PF13174" : 'TPR', # TPR6
"PF13176" : 'TPR', # TPR7
"PF13181" : 'TPR', # TPR8
'PF13371' : 'TPR', # TPR9
'PF13374' : 'TPR', # TPR10
"PF13414" : 'TPR', # TPR11
"PF13424" : 'TPR', # TPR12
# "PF07719" : 'TPR', # Can't find entry for TPR_13
"PF13428" : 'TPR', # TPR14
"PF13429" : 'TPR', # TPR15
"PF13432" : 'TPR', # TPR16
'PF13431' : 'TPR', # TPR17
"PF13512" : 'TPR', # TPR18
'PF14559' : 'TPR', # TPR19
"PF14561" : 'TPR', # TPR20
"PF09976" : 'TPR', # TPR21
"PS50005" : 'TPR', # Prosite profile for TPR
'IPR013026' : 'TPR', # Interproscan annotation for TPR
'SM00386' : 'half_TPR',
'SSF48403' : 'Ankyrin', # Ankyrin superfamily (CL0465)
'PF12796' : 'Ankyrin',
'SSF50978' : 'WD40-like', # WD40-like superfamily (no code)
'PF00400' : 'WD40'
}

N_terminal_list = set(['PF00168', 'PF17109', 'PF14479', 'PF17111', 'PF06985', 'PF01734', 'PF00082', 'PF11558', 'PF00069', 'PF01048', 'PF04607', 'PF17046', 'PF17107'])
NBD_list = set(['SSF52540', 'IPR027417', 'PF13401', 'PF00931', 'PF05729']) # The P-loop containing protein is tested first so that NACHT or NBARC annotations will replace it
C_terminal_list = set([
 'SSF48452', 'SSF48403', 'SSF50978',
 'SM00386', "PF00515", "PF07719", "PF07720",
 "PF07721", "PF12688", "PF13174", "PF13176", "PF13181", 'PF13371',
 'PF13374', "PF13414", "PF13424", "PF13428", "PF13429", "PF13432",
 'PF13431', "PF13512", 'PF14559', "PF14561", "PF09976", "PS50005",
 'IPR013026', 'PF12796',
 'PF00400'
 ])

set1_count = 0
set2_count = 0
set1_dict = defaultdict(int)
set2_dict = defaultdict(int)
NLR_structure_dict1 = defaultdict(int)
NLR_structure_dict2 = defaultdict(int)
seen_NLR_set = set()

for line in interpro_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    gene_id = split_line[0]
    m = re.findall("IPR......|PF.....|SSF.....", line)
    # print line
    # print m
    if m:
        NLR_set = set(m)
    else:
        # NLR_set = set(["no_annotation"])
        continue
    # print NLR_set
    NLR_hits = [ x for x in NLR_set if x in pfam_domain_dict]
    # if NLR_hits == []:
    #     continue

    NLR_hit_set = set([])
    for term in NLR_hits:
        NLR_hit_set.add(pfam_domain_dict[term])

    # print NLR_hits
    if split_line[1] in set1_contigs:
        # print NLR_hit_set
        for NLR in NLR_hit_set:
            set1_dict[NLR] += 1
            seen_NLR_set.add(NLR)
        set1_count += 1
        NLR_structure = "-".join(NLR_hit_set)
        NLR_structure_dict1[NLR_structure] += 1
        # [seen_NLR_set.add(x) for x in NLR_hit_set]
    elif split_line[1] in set2_contigs:
        for NLR in NLR_hit_set:
            set2_dict[NLR] += 1
            seen_NLR_set.add(NLR)
        set2_count += 1
        NLR_structure = "-".join(NLR_hit_set)
        NLR_structure_dict2[NLR_structure] += 1

        # [seen_NLR_set.add(x) for x in NLR_hit_set]

print set1_count
print set2_count

for NLR in seen_NLR_set:
    set1_NLR = set1_dict[NLR]
    set2_NLR = set2_dict[NLR]
    set1_others = set1_count - set1_NLR
    set2_others = set2_count - set2_NLR
    outline1 = "\t".join([str(NLR), str(set1_NLR), str(set2_NLR)]) + "\n"
    outline2 = "\t".join(["Other genes", str(set1_others), str(set2_others)]) + "\n"

    outfile = "".join([conf.outdir, "/", NLR, "_fischertable.txt"])
    # print outfile
    o = open(outfile, 'w')
    o.write("".join([outline1, outline2]))
    o.close()

for key in NLR_structure_dict1.keys():
    print key + "\t" + str(NLR_structure_dict1[key])
for key in NLR_structure_dict2.keys():
    print key + "\t" + str(NLR_structure_dict2[key])
