#!/usr/bin/python

'''
This program is used to extract information on NLR genes containing STAND domains
from summary tables of FoL genes. It outputs an oragsnims genes containing
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
import numpy
from sets import Set
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()
ap.add_argument('--summary_table',required=True,type=str,help='The summary table from Fo)build_annot_table.py')
ap.add_argument('--N_term_out',required=True,type=str,help='The output file detailing N terminal domains by chromosome')
ap.add_argument('--NBD_out',required=True,type=str,help='The output file detailing NBD domains by chromosome')
ap.add_argument('--C_term_out',required=True,type=str,help='The output file detailing C terminal domains by chromosome')

conf = ap.parse_args()

file_N_term_out = open(conf.N_term_out, 'w')
file_NBD_out = open(conf.NBD_out, 'w')
file_C_term_out = open(conf.C_term_out, 'w')

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
contig_dict = defaultdict(list)
gene_count_dict = defaultdict(int)
i = 0
for line in summary_table_lines:
    line = line.rstrip()
    split_line = line.split("\t")
    gene_id = split_line[0]
    if len(split_line) >= 34:
        interpro_info = split_line[33]
        orthogroup = split_line[15]
        orthogroup_content = split_line[17]
        Fol_contig = split_line[1]
        gene_start = split_line[2]
        gene_end = split_line[3]

        FoL_Chr1_LS_start = int(5963944)
        FoL_Chr2_LS_start = int(4903160)
        if str(1) == Fol_contig:
            if (int(gene_end) <= FoL_Chr1_LS_start):
                Fol_contig = "1_core"
            else:
                Fol_contig = "1_LS"
        elif str(2) == Fol_contig:
            if (int(gene_end) <= FoL_Chr2_LS_start):
                Fol_contig = "2_core"
            else:
                Fol_contig = "2_LS"
        gene_count_dict[Fol_contig] += 1

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
        combination = N_terminal + "---" + NBD + "---" + C_terminal
        store_line = gene_id + "\t" + orthogroup + "\t" + orthogroup_content

        if any( y != 'NA' for y in [N_terminal, NBD, C_terminal] ):
            NLR_prot_dict[combination].append(store_line)
            i += 1
            contig_dict[Fol_contig].append((combination + "\t" + store_line))

#-----------------------------------------------------
# Step 3
# Print genes containing NLR domains
#-----------------------------------------------------

for k in sorted(NLR_prot_dict, key=lambda k: len(NLR_prot_dict[k]), reverse=True):
    num_genes = len(NLR_prot_dict[k])
    print k + "\t" + str(num_genes)

print( "\n\nthe total number of NLR genes was:\t" + str(i) )

#-----------------------------------------------------
# Step 4
# Print genes by domain and by chromosome
#-----------------------------------------------------



print("\n\nthe following genes were found in FoL regions:")

for FoL_contig in sorted(contig_dict):
    N_terminal_dict = defaultdict(list)
    NBD_dict = defaultdict(list)
    C_terminal_dict = defaultdict(list)
    genes_in_contig = contig_dict[FoL_contig]
    total_genes = gene_count_dict[FoL_contig]
    print("\n\nFoL contig " + str(FoL_contig) + " (" + str(len(genes_in_contig)) + "/ " + str(total_genes) +" genes)")
    file_N_term_out.write("\n\nFoL contig " + str(FoL_contig) + " (" + str(len(genes_in_contig)) + "/ " + str(total_genes) +" genes)" + "\n")
    file_NBD_out.write("\n\nFoL contig " + str(FoL_contig) + " (" + str(len(genes_in_contig)) + "/ " + str(total_genes) +" genes)" + "\n")
    file_C_term_out.write("\n\nFoL contig " + str(FoL_contig) + " (" + str(len(genes_in_contig)) + "/ " + str(total_genes) +" genes)" + "\n")
    for line in genes_in_contig:
        split_line = line.split("\t")
        combination = split_line[0].split("---")
        N_terminal = combination[0]
        NBD = combination[1]
        C_terminal = combination[2]
        if N_terminal != 'NA':
            N_terminal_dict[N_terminal].append(line)
        if NBD != 'NA':
            NBD_dict[NBD].append(line)
        if C_terminal != 'NA':
            C_terminal_dict[C_terminal].append(line)
    print("N-terminal")
    # for keys in sorted(N_terminal_dict, key=lambda keys: len(N_terminal_dict[keys]), reverse=True):
    for keys in ["PKinase", "HET", "Peptidase S8", "PND UDP", "Helo-like", "Helo", "C2", "Patatin", "Goodbye", "PFD", "RelA Spot", 'SesA', 'sesB-like']:
        num_genes = len(N_terminal_dict[keys])
        frac_genes = numpy.divide(float(num_genes), float(total_genes)) * 100
        # print ("\t" + str(keys) + "\t" + str(num_genes) + "\t" + str(frac_genes))
        file_N_term_out.write("\t" + FoL_contig + "\t" + str(keys) + "\t" + str(num_genes) + "\t" + str(frac_genes) + "\n")
    print("NBD")
    # for keys in sorted(NBD_dict, key=lambda keys: len(NBD_dict[keys]), reverse=True):
    for keys in ['P-loop NTPase', 'NACHT', 'NB-ARC', 'AAA']:
        num_genes = len(NBD_dict[keys])
        frac_genes = numpy.divide(float(num_genes), float(total_genes)) * 100
        # print ("\t" + str(keys) + "\t" + str(num_genes) + "\t" + str(frac_genes))
        file_NBD_out.write("\t" + FoL_contig + "\t" + str(keys) + "\t" + str(num_genes) + "\t" + str(frac_genes) + "\n")
    print("C-terminal")
    # for keys in sorted(C_terminal_dict, key=lambda keys: len(C_terminal_dict[keys]), reverse=True):
    for keys in ['Ankyrin', 'TPR', 'WD40', 'WD40-like', 'half TPR']:
        num_genes = len(C_terminal_dict[keys])
        frac_genes = numpy.divide(float(num_genes), float(total_genes)) * 100
        # print ("\t" + str(keys) + "\t" + str(num_genes) + "\t" + str(frac_genes))
        file_C_term_out.write("\t" + FoL_contig + "\t" + str(keys) + "\t" + str(num_genes) + "\t" + str(frac_genes) + "\n")

    del N_terminal_dict
    del NBD
    del C_terminal_dict
