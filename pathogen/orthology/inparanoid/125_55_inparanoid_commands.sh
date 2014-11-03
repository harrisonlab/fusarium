#!/bin/bash
# #----------------------------------------------
# #	OrthoMCL
# #				Commands: 
# #----------------------------------------------
# 
# 
# 
# cd /home/groups/harrisonlab/project_files/fusarium
# mkdir test
# cp /home/armita/prog/orthomcl/orthomclSoftware-v2.0.9/doc/OrthoMCLEngine/Main/orthomcl.config.template .
# 
# 
# vi orthomcl.config.template
# 
# orthomclInstallSchema orthomcl.config.template
# 
# ## this config assumes a mysql database named 'orthomcl'.  adjust according
# ## to your situation.
# # dbVendor=mysql
# # dbConnectString=dbi:mysql:orthomcl
# # dbLogin=armita1
# # dbPassword=strawberry2018
# # similarSequencesTable=SimilarSequences
# # orthologTable=Ortholog
# # inParalogTable=InParalog
# # coOrthologTable=CoOrtholog
# # interTaxonMatchView=InterTaxonMatch
# # percentMatchCutoff=50
# # evalueExponentCutoff=-5
# # oracleIndexTblSpc=NONE
# 
# orthomclAdjustFasta 125 ../gene_pred/augustus/F.oxysporum_fsp_cepae/125/125_aug_out.aa '|'
# 
# orthomclAdjustFasta 55 ../gene_pred/augustus/F.oxysporum_fsp_cepae/55/55_aug_out.aa '|'
# 
# makeblastdb -in 125.fasta -dbtype prot -out 125
# 
# makeblastdb -in 55.fasta -dbtype prot -out 55
# 
# blastp -db 125 -query 55.fasta -outfmt 8 > 55_vs_125.txt
# 
# blastp -db 55 -query 125.fasta -outfmt 8 > 125_vs_55.txt
# 
# blastp -db 125 -query 125.fasta -outfmt 8 > 125_vs_125.txt
# 
# blastp -db 55 -query 55.fasta -outfmt 8 > 55_vs_55.txt
# 
# orthomclBlastParser 125_vs_55.txt 55.fasta > 55_sim_seqs.txt
# 
# orthomclBlastParser 55_vs_125.txt 125.fasta > 125_sim_seqs.txt
# 


#----------------------------------------------
#	INPARANOID
#			commands:
#----------------------------------------------
# cd /home/groups/harrisonlab/project_files/fusarium/inparanoid_test
# 
# cp /home/armita/prog/inparanoid_4.1/ .
# 
# cd inparanoid_4.1/
# 
# mv 125.fasta 125_fasta
# mv 55.fasta 55_fasta
# 
# ./inparanoid.pl 125_fasta 55_fasta
# grep '>' 125.fasta | sed 's/>//' > 125_seqs.txt
# grep '>' 55.fasta | sed 's/>//' > 55_seqs.txt
# cut -f3 table.125_fasta-55_fasta | sed 's/ 125/\n125/g' | cut -d' ' -f1 | tail -n+2 | cat - 125_seqs.txt | sort | uniq -u | wc -l
# #	125 contains 707 non-othologous genes
# cut -f3 table.125_fasta-55_fasta | sed 's/ 125/\n125/g' | cut -d' ' -f1 | tail -n+2 | cat - 125_seqs.txt | sort | uniq -u > 125_uniq_vs_55.txt
# cat 125_uniq_vs_55.txt | cut -d '|' -f2 | xargs -I{} grep -w {} ../../gene_pred/augustus/F.oxysporum_fsp_cepae/125/125_aug.gff | sed 's/AUGUSTUS/125_uniq_vs_55/' > 125_uniq_vs_55.gff
# # upon loading into Geneious many of these 707 transcripts were found to be incomplete proteins on the ends of contigs
# /home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/filter_gff_StartStop.pl 125_uniq_vs_55.gff > 125_uniq_vs_55_filtered.gff
# # filtering the .gff file for complete genes (that had both a start and a stop codon annotations) resulted in 407 unique genes.
# 
# cut -f4 table.125_fasta-55_fasta | sed 's/ 55/\n55/g' | cut -d' ' -f1 | tail -n+2 | cat - 55_seqs.txt | sort | uniq -u | wc -l
# #	55 conatins 433 non-orthologous genes
# cut -f4 table.125_fasta-55_fasta | sed 's/ 55/\n55/g' | cut -d' ' -f1 | tail -n+2 | cat - 55_seqs.txt | sort | uniq -u > 55_uniq_vs_125.txt
# cat 55_uniq_vs_125.txt | cut -d '|' -f2 | xargs -I{} grep -w {} ../../gene_pred/augustus/F.oxysporum_fsp_cepae/55/55_aug.gff | sed 's/AUGUSTUS/55_uniq_vs_125/' > 55_uniq_vs_125.gff
# # upon loading into Geneious many of these 433 transcripts were found to be incomplete proteins on the ends of contigs
# /home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/filter_gff_StartStop.pl 55_uniq_vs_125.gff > 55_uniq_vs_125_filtered.gff
# # filtering the .gff file for complete genes (that had both a start and a stop codon annotations) resulted in 325 unique genes.


set -- 125 55 A23 A28 D2 Fus2 HB17 PG
for a; do 
	shift
	for b; do
	printf "%s - %s\n" "$a" "$b"
	qsub /home/armita/git_repos/emr_repos/tools/pathogen/orthology/inparanoid/sub_inparanoid.sh gene_pred/augustus/F.oxysporum_fsp_cepae/"$a"/"$a"_aug_out.aa gene_pred/augustus/F.oxysporum_fsp_cepae/"$b"/"$b"_aug_out.aa gene_pred/augustus/F.oxysporum_fsp_cepae/"$a"/"$a"_aug.gff gene_pred/augustus/F.oxysporum_fsp_cepae/"$b"/"$b"_aug.gff
	done 
done





#----------------------------------------------
#	Inparanoid 
#				Adapted for cluster
#----------------------------------------------

makeblastdb -in 125_fasta2 -dbtype prot -out 125

makeblastdb -in 55_fasta2 -dbtype prot -out 55

echo '' > 55_vs_125.txt
echo '' > 125_vs_55.txt
echo '' > 125_vs_125.txt
echo '' > 55_vs_55.txt

blastp -db 125 -query 55_fasta2 -outfmt 5 -num_threads 16 -comp_based_stats 3 | ./blast_parser.pl 40 >> 55_vs_125.xml

blastp -db 55 -query 125_fasta2 -outfmt 5 -num_threads 16 -comp_based_stats 3 | ./blast_parser.pl 40 >> 125_vs_55.xml

blastp -db 125 -query 125_fasta2 -outfmt 5 -num_threads 16 -comp_based_stats 3 | ./blast_parser.pl 40 >> 125_vs_125.xml

blastp -db 55 -query 55_fasta2 -outfmt 5 -num_threads 16 -comp_based_stats 3 | ./blast_parser.pl 40 >> 55_vs_55.xml

cp 55_vs_125.xml 55_vs_125
cp 55_vs_55.xml 55_vs_55
cp 125_vs_55.xml 125_vs_55
cp 125_vs_125.xml 125_vs_125
