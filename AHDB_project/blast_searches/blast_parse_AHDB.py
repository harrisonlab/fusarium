#!/usr/bin/python

'''
A modification of the blast_parse.py script stored in the pathogen/blast
repository. This program also accepts fasta files of the target genomes and
extracts cds flanking each blast hit of interest.

This summarises results from blast2csv from multiple files, building a matrix
of blast hits for each genome searched against. Percentage identify accross the
blast hit and e-value of the hit need to be over a given threshold.
'''


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import os
import sys
import argparse
import re
import numpy as np
from sets import Set
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

ap = argparse.ArgumentParser()
ap.add_argument('--blast_csv',required=True,nargs='+',type=str,help='Blast2csv output')
ap.add_argument('--headers',required=True,nargs='+',type=str,help='names of target genomes for use as headers in the outpit file. Order must be the same as the Blast2csv output')
ap.add_argument('--genomes',required=True,nargs='+',type=str,help='list of target genome fasta files in the same order as blasthits and headers.')
ap.add_argument('--identity',required=True,type=float,help='Threshold of percentage identity accross the entire query length')
ap.add_argument('--evalue',required=True,type=float,help='Threshold of E-value for any blast hit')
ap.add_argument('--out_prefix',required=True,type=str,help='Prefix for output files.')
conf = ap.parse_args()

csv_files = conf.blast_csv
headers = conf.headers
ref_genomes = conf.genomes
thresh_perc_id = conf.identity
thresh_eval = conf.evalue
prefix = conf.out_prefix

# print len(headers)
# print len(csv_files)
if len(headers) != len(csv_files):
    raise ValueError('different number of header names from blast_csv files')

hits_dict = defaultdict(list)

num = 0
fasta_dict = defaultdict(list)
# out_records = []
for organism, csv_file, genome_f in zip(headers, csv_files, ref_genomes):
    # print "\t".join([organism, csv_file, genome_f])
    num += 1
    # print num
    # header = csv_file.split("/")[-1].replace(".csv", "")
    gene_list = []
    with open(csv_file) as f:
        csv_lines = f.readlines()
    fasta_obj = list(SeqIO.parse(genome_f, "fasta"))
    # fasta_obj_parsed = []
    for record in fasta_obj:
        record_id = record.id.split()[0]
        # record.id = record_id[4:]
        # print record.id
    print organism
    record_dict = SeqIO.to_dict(fasta_obj)
    # record_dict = SeqIO.index(genome_f, "fasta")
    for line in csv_lines:
        line = line.rstrip()
        if line.startswith('ID'):
            headers_line = line
        else:
            split_line = line.split("\t")
            query_ID = split_line[0]
            query_seq = split_line[1]
            query_lgth = int(split_line[2])
            num_hits = int(split_line[3])
            over_threshold = 0
            gene_list.append(query_ID)
            if num_hits > 0:
                hits_cols = split_line[4:]
                for i in range(num_hits):
                    hit_contig = hits_cols.pop(0)
                    hit_eval = float(hits_cols.pop(0))
                    hit_lgth = int(hits_cols.pop(0))
                    hit_perc_lgth = float(hits_cols.pop(0))
                    hit_total_perc_id = float(hits_cols.pop(0))
                    hit_strand = hits_cols.pop(0)
                    hit_start = int(hits_cols.pop(0))
                    hit_end = int(hits_cols.pop(0))
                    hit_seq = hits_cols.pop(0)
                    perc_identity_accross_hit = np.divide(hit_total_perc_id, hit_perc_lgth)
                    if perc_identity_accross_hit >= thresh_perc_id and hit_eval <= thresh_eval:
                        over_threshold += 1
                        # print hit_contig
                        # print record_dict.keys()
                        # print str(thresh_perc_id) + "\t" + str(hit_perc_id)
                        # print str(thresh_eval) + "\t" + str(hit_eval)

                        hit_contig_rec = record_dict[hit_contig]
                        # if hit_contig not in record_dict.keys():
                        #     continue
                        hit_contig_rec = record_dict[hit_contig]
                        # hit_contig_seq = fasta_obj[hit_contig]
                        contig_name = hit_contig_rec.name
                        hit_contig_seq = hit_contig_rec.seq
                        start_end_list = [hit_start, hit_end]
                        x = min(start_end_list) -500
                        y = max(start_end_list) +500
                        if x < 0:
                            x = 0
                        if y > len(hit_contig_seq):
                            y = len(hit_contig_seq)
                        hit_header = "_".join([organism, str(contig_name),str(x),str(y)])
                        hit_seq_expanded = str(hit_contig_seq[x:y])
                        # print hit_seq_expanded
                        hit_record = "\n".join([">" + hit_header, hit_seq_expanded])
                        # hit_record = SeqRecord(hit_seq_expanded, id=hit_header, name='', description='')
                        # print hit_record.format("fasta")
                        # print hit_record
                        # out_records.append(hit_record)
                        fasta_dict[query_ID].append(hit_record)
                        # fasta_dict[query_ID].append("".join([">", hit_header, "\n", hit_seq])
            hits_dict[query_ID].append(str(over_threshold))


# print num
outlines = ["Query ID's\t" + "\t".join(headers)]
# print outline
# outline = "Query ID's\t" + "\t".join(csv_files)
# print outline
for key in gene_list:
    outlines.append(key + "\t" + "\t".join(hits_dict[key]))
    # print outline
    # print len(hits_dict[key])
out_name = prefix + ".csv"
f = open(out_name, 'w')
f.write("\n".join(outlines))
f.close()

# l = fasta_dict.keys()
# print l

for key in gene_list:
    # print key
    if key in fasta_dict:
        # print "badgers"
        out_name = prefix+ "_" + key.replace('|', '_')  + "_hits.fa"
        f = open(out_name, 'w')
        f.write("\n".join(fasta_dict[key]))
        f.close()


        # print "\n".join(fasta_dict[query_ID])

        # for seq_record in fasta_dict[query_ID]:
        #     print seq_record
#             # print seq_record.format("fasta")
#         SeqIO.write(fasta_dict[query_ID], "tmp.fa", "fasta")
# SeqIO.write(out_records, "tmp.fa", "fasta")
