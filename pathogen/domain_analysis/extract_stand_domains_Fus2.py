#!/usr/bin/python

'''
This program is used to extract information on NLR genes containing STAND domains
from summary tables of Fusarium genes. It outputs an oragsnims genes containing
NLR genes (identified by N-terminal STAND domain, a nucleotide binding domain
that is either NACHT or NB-ARC and a C-terminal domain that is an Ankyrin repeat,
a WD domain or a TPR domain)
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
ap.add_argument('--summary_table',required=True,type=str,help='The summary table from Fo_build_annot_table.py')


conf = ap.parse_args()


with open(conf.summary_table) as f:
    summary_table_lines = f.readlines()

pfam_domain_dict = {
# N-terminal domains:
'PF00168' : "C2",
'PF17109' : "Goodbye",
'PF14479': "Helo",
'PF17111' : "Helo-like",
'PF06985' : "HET",
'PF01734' : "Patatin",
'PF00082' : "Peptidase S8",
'PF11558' : "PFD",
'PF00069' : "PKinase",
'PF01048' : "PND UDP",
'PF04607' : "RelA Spot",
'PF17107' : 'SesA',
'PF17046' : 'sesB-like',
# Nucleotide binding domains:
'PF00931' : 'NB-ARC',
'PF05729' : 'NACHT',
'SSF52540' : 'P-loop NTPase', # P-loop NTPase superfamily (CL0023)
'IPR027417' : 'P-loop NTPase',
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
'SM00386' : 'half TPR',
'SSF48403' : 'Ankyrin', # Ankyrin superfamily (CL0465)
'PF12796' : 'Ankyrin',
'SSF50978' : 'WD40-like', # WD40-like superfamily (no code)
'PF00400' : 'WD40'
}

#-----------------------------------------------------
# Step 2
# Extract interproscan information on NLR candidate
# proteins
#-----------------------------------------------------

NLR_prot_dict = defaultdict(list)
# path_contig_dict = defaultdict(list)
# nonpath_contig_dict = defaultdict(list)
Fus2_contig_dict = defaultdict(list)
i = 0
for line in summary_table_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    gene_id = split_line[0]
    Fus2_contig = split_line[1]
    if len(split_line) >= 32:
        interpro_info = split_line[31]
        # print interpro_info
        orthogroup = split_line[15]
        orthogroup_content = split_line[17]
        Fol_contig = split_line[21]

        N_terminal_list = ['PF00168', 'PF17109', 'PF14479', 'PF17111', 'PF06985', 'PF01734', 'PF00082', 'PF11558', 'PF00069', 'PF01048', 'PF04607', 'PF17046', 'PF17107']
        NBD_list = ['SSF52540', 'IPR027417', 'PF13401', 'PF00931', 'PF05729'] # The P-loop containing protein is tested first so that NACHT or NBARC annotations will replace it
        C_terminal_list = [
         'SSF48452', 'SSF48403', 'SSF50978',
         'SM00386', "PF00515", "PF07719", "PF07720",
         "PF07721", "PF12688", "PF13174", "PF13176", "PF13181", 'PF13371',
         'PF13374', "PF13414", "PF13424", "PF13428", "PF13429", "PF13432",
         'PF13431', "PF13512", 'PF14559', "PF14561", "PF09976", "PS50005",
         'IPR013026', 'PF12796',
         'PF00400'
         ]
        N_terminal = 'NA'
        NBD = 'NA'
        C_terminal = 'NA'
        for x in N_terminal_list:
            if x in interpro_info:
                domain = pfam_domain_dict[x]
                N_terminal = domain

        for x in NBD_list:
            if x in interpro_info:
                domain = pfam_domain_dict[x]
                NBD = domain

        for x in C_terminal_list:
            if x in interpro_info:
                domain = pfam_domain_dict[x]
                C_terminal = domain

        # if any( y != 'NA' for y in [N_terminal, NBD, C_terminal] ):
        if all( y != 'NA' for y in [N_terminal, NBD] ):
        # if all( y != 'NA' for y in [N_terminal, NBD, C_terminal] ):
            combination = N_terminal + "---" + NBD + "---" + C_terminal

            # store_line = gene_id + "\t" + orthogroup + "\t" + Fol_contig + "\t" + interpro_info
            store_line = gene_id + "\t" + orthogroup + "\t" + orthogroup_content
            NLR_prot_dict[combination].append(store_line)
            i += 1
            Fus2_contig_dict[Fus2_contig].append((combination + "\t" + store_line))
            # path_contig_list = [3, 6, 14, 15]
            # nonpath_contig_list = [1, 2, 4, 5, 7, 8, 9, 10, 11, 12, 13]
            # if any( str(z) == Fol_contig for z in path_contig_list ):
            #     path_contig_dict[Fol_contig].append((combination + "\t" + store_line))
            # elif any( str(z) == Fol_contig for z in nonpath_contig_list ):
            #     nonpath_contig_dict[Fol_contig].append((combination + "\t" + store_line))
            # else:
            #     nonpath_contig_dict['16'].append((combination + "\t" + store_line))

#-----------------------------------------------------
# Step 3
# Print genes containing NLR domains
#-----------------------------------------------------

for k in sorted(NLR_prot_dict, key=lambda k: len(NLR_prot_dict[k]), reverse=True):
    num_genes = len(NLR_prot_dict[k])
    print k + "\t" + str(num_genes)

print( "\n\nthe total number of NLR genes was:\t" + str(i) )

print("\n\nthe following genes were found in FoC contigs:")

# for k in sorted(Fus2_contig_dict, key=lambda k: ):
for k in sorted(Fus2_contig_dict.keys()):
    print("\nFus2 contig " + str(k) + " (" + str(len(Fus2_contig_dict[k])) + " genes)")
    print("\n".join(Fus2_contig_dict[k]))


# print("\n\nthe following genes were found in FoL LS contigs:")

# for k in sorted(path_contig_dict, key=int):
#     print("\nFoL contig " + str(k) + " (" + str(len(path_contig_dict[k])) + " genes)")
#     print("\n".join(path_contig_dict[k]))


# print("\n\nthe following genes were found in FoL core contigs:")
#
# for k in sorted(nonpath_contig_dict, key=int):
#     if k == '16':
#         print("\nFoL contig unhit (" + str(len(nonpath_contig_dict[k])) + " genes)")
#     else:
#         print("\nFoL contig " + str(k) + " (" + str(len(nonpath_contig_dict[k])) + " genes)")
#     print("\n".join(nonpath_contig_dict[k]))
