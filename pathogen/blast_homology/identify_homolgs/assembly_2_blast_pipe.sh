cp /home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_1/125_S3_L001_R1_001.fastq.gz 125/F/.
cp /home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_1/125_S3_L001_R2_001.fastq.gz 125/R/.
cp /home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_2/55_S3_L001_R1_001.fastq.gz 55/F/.
cp /home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_2/55_S3_L001_R2_001.fastq.gz 55/R/.
cp /home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_2/A23_S1_L001_R1_001.fastq.gz A23/F/.
cp /home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_2/A23_S1_L001_R2_001.fastq.gz A23/R/.
cp /home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_2/A28_S2_L001_R1_001.fastq.gz A28/F/.
cp /home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_2/A28_S2_L001_R2_001.fastq.gz A28/R/.
cp /home/groups/harrisonlab/raw_data/raw_seq/fusarium/warwick_seqs/D2/s_5_1_sequence.txt D2/F/.
cp /home/groups/harrisonlab/raw_data/raw_seq/fusarium/warwick_seqs/D2/s_5_2_sequence.txt D2/R/.
cp /home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_1/FUS2_S2_L001_R1_001.fastq.gz Fus2/F/.
cp /home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_1/FUS2_S2_L001_R2_001.fastq.gz Fus2/R/.
cp /home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_2/HB17_S4_L001_R1_001.fastq.gz HB17/F/.
cp /home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_2/HB17_S4_L001_R2_001.fastq.gz HB17/R/.
cp /home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_1/PG_S1_L001_R1_001.fastq.gz PG/F/.
cp /home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_1/PG_S1_L001_R2_001.fastq.gz PG/R/.

gunzip raw_dna/paired/F.oxysporum_fsp_cepae/*/*/*.gz

#PROGRAM="/home/armita/git_repos/emr_repos/tools/seq_tools/initial_pipe/assembly_pipe.sh"
#SCRIPT_DIR="/home/armita/git_repos/emr_repos/tools/seq_tools/initial_pipe"
#F_READ="raw_dna/paired/F.oxysporum_fsp_cepae/125/F/125_S3_L001_R1_001.fastq"
#R_READ="raw_dna/paired/F.oxysporum_fsp_cepae/125/R/125_S3_L001_R2_001.fastq"
#qsub $PROGRAM $F_READ $R_READ 61 600 $SCRIPT_DIR

cp raw_dna/paired/F.oxysporum_fsp_cepae/D2/F/s_5_1_sequence.txt raw_dna/paired/F.oxysporum_fsp_cepae/D2/F/s_5_1_sequence.fastq
cp raw_dna/paired/F.oxysporum_fsp_cepae/D2/R/s_5_2_sequence.txt raw_dna/paired/F.oxysporum_fsp_cepae/D2/R/s_5_2_sequence.fastq
gzip raw_dna/paired/F.oxysporum_fsp_cepae/D2/F/s_5_1_sequence.txt 
gzip raw_dna/paired/F.oxysporum_fsp_cepae/D2/R/s_5_2_sequence.txt
mv raw_dna/paired/F.oxysporum_fsp_cepae/D2/F/s_5_1_sequence.txt.gz raw_dna/paired/F.oxysporum_fsp_cepae/D2/F/s_5_1_sequence.fastq.gz
mv raw_dna/paired/F.oxysporum_fsp_cepae/D2/R/s_5_2_sequence.txt.gz raw_dna/paired/F.oxysporum_fsp_cepae/D2/R/s_5_2_sequence.fastq.gz

qsub /home/armita/git_repos/emr_repos/tools/seq_tools/initial_pipe/assembly_pipe.sh raw_dna/paired/F.oxysporum_fsp_cepae/125/F/125_S3_L001_R1_001.fastq raw_dna/paired/F.oxysporum_fsp_cepae/125/R/125_S3_L001_R2_001.fastq 61 600 /home/armita/git_repos/emr_repos/tools/seq_tools/initial_pipe
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/initial_pipe/assembly_pipe.sh raw_dna/paired/F.oxysporum_fsp_cepae/55/F/55_S3_L001_R1_001.fastq raw_dna/paired/F.oxysporum_fsp_cepae/55/R/55_S3_L001_R2_001.fastq 61 600 /home/armita/git_repos/emr_repos/tools/seq_tools/initial_pipe
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/initial_pipe/assembly_pipe.sh raw_dna/paired/F.oxysporum_fsp_cepae/A23/F/A23_S1_L001_R1_001.fastq raw_dna/paired/F.oxysporum_fsp_cepae/A23/R/A23_S1_L001_R2_001.fastq 61 600 /home/armita/git_repos/emr_repos/tools/seq_tools/initial_pipe
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/initial_pipe/assembly_pipe.sh raw_dna/paired/F.oxysporum_fsp_cepae/A28/F/A28_S2_L001_R1_001.fastq raw_dna/paired/F.oxysporum_fsp_cepae/A28/R/A28_S2_L001_R2_001.fastq 61 300 /home/armita/git_repos/emr_repos/tools/seq_tools/initial_pipe
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/initial_pipe/assembly_pipe.sh raw_dna/paired/F.oxysporum_fsp_cepae/D2/F/s_5_1_sequence.fastq raw_dna/paired/F.oxysporum_fsp_cepae/D2/R/s_5_2_sequence.fastq 61 300 /home/armita/git_repos/emr_repos/tools/seq_tools/initial_pipe
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/initial_pipe/assembly_pipe.sh raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/F/FUS2_S2_L001_R1_001.fastq raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/R/FUS2_S2_L001_R2_001.fastq 61 300 /home/armita/git_repos/emr_repos/tools/seq_tools/initial_pipe
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/initial_pipe/assembly_pipe.sh raw_dna/paired/F.oxysporum_fsp_cepae/HB17/F/HB17_S4_L001_R1_001.fastq raw_dna/paired/F.oxysporum_fsp_cepae/HB17/R/HB17_S4_L001_R2_001.fastq 61 300 /home/armita/git_repos/emr_repos/tools/seq_tools/initial_pipe
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/initial_pipe/assembly_pipe.sh raw_dna/paired/F.oxysporum_fsp_cepae/PG/F/PG_S1_L001_R1_001.fastq raw_dna/paired/F.oxysporum_fsp_cepae/PG/R/PG_S1_L001_R2_001.fastq 61 300 /home/armita/git_repos/emr_repos/tools/seq_tools/initial_pipe

grep -A1 '>' analysis/blast_homology/six_genes/six-appended.fa | cut -d '|' -f5 > analysis/blast_homology/six_genes/six-appended_parsed.fa
vi analysis/blast_homology/six_genes/six-appended_parsed_tmp.fa # Removed initial ' ' from line and added '>' 
sed 's/\s/_/g' analysis/blast_homology/six_genes/six-appended_parsed_tmp.fa > analysis/blast_homology/six_genes/six-appended_parsed.fa

qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/PHIbase/PHI_36_accessions.fa protein assembly/velvet/F.oxysporum_fsp_cepae/125/125_assembly.41/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/PHIbase/PHI_36_accessions.fa protein assembly/velvet/F.oxysporum_fsp_cepae/55/55_assembly.41/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/PHIbase/PHI_36_accessions.fa protein assembly/velvet/F.oxysporum_fsp_cepae/A23/A23_assembly.41/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/PHIbase/PHI_36_accessions.fa protein assembly/velvet/F.oxysporum_fsp_cepae/A28/A28_assembly.41/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/PHIbase/PHI_36_accessions.fa protein assembly/velvet/F.oxysporum_fsp_cepae/D2/D2_assembly.41/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/PHIbase/PHI_36_accessions.fa protein assembly/velvet/F.oxysporum_fsp_cepae/Fus2/Fus2_assembly.41/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/PHIbase/PHI_36_accessions.fa protein assembly/velvet/F.oxysporum_fsp_cepae/HB17/HB17_assembly.41/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/PHIbase/PHI_36_accessions.fa protein assembly/velvet/F.oxysporum_fsp_cepae/PG/PG_assembly.41/sorted_contigs.fa

qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/six_genes/six-appended_parsed.fa dna assembly/velvet/F.oxysporum_fsp_cepae/125/125_assembly.41/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/six_genes/six-appended_parsed.fa dna assembly/velvet/F.oxysporum_fsp_cepae/55/55_assembly.41/sorted_contigs.fa 
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/six_genes/six-appended_parsed.fa dna assembly/velvet/F.oxysporum_fsp_cepae/A23/A23_assembly.41/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/six_genes/six-appended_parsed.fa dna assembly/velvet/F.oxysporum_fsp_cepae/A28/A28_assembly.41/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/six_genes/six-appended_parsed.fa dna assembly/velvet/F.oxysporum_fsp_cepae/D2/D2_assembly.41/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/six_genes/six-appended_parsed.fa dna assembly/velvet/F.oxysporum_fsp_cepae/Fus2/Fus2_assembly.41/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/six_genes/six-appended_parsed.fa dna assembly/velvet/F.oxysporum_fsp_cepae/HB17/HB17_assembly.41/sorted_contigs.fa
qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/six_genes/six-appended_parsed.fa dna assembly/velvet/F.oxysporum_fsp_cepae/PG/PG_assembly.41/sorted_contigs.fa

/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/blast_homology/identify_differentials/homolog_grp_srt.pl analysis/blast_homology/F.oxysporum_fsp_cepae/125/125_PHI_36_accessions.fa_homologs.csv analysis/blast_homology/F.oxysporum_fsp_cepae/55/55_PHI_36_accessions.fa_homologs.csv analysis/blast_homology/F.oxysporum_fsp_cepae/A23/A23_PHI_36_accessions.fa_homologs.csv analysis/blast_homology/F.oxysporum_fsp_cepae/A28/A28_PHI_36_accessions.fa_homologs.csv analysis/blast_homology/F.oxysporum_fsp_cepae/D2/D2_PHI_36_accessions.fa_homologs.csv analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2/Fus2_PHI_36_accessions.fa_homologs.csv analysis/blast_homology/F.oxysporum_fsp_cepae/HB17/HB17_PHI_36_accessions.fa_homologs.csv analysis/blast_homology/F.oxysporum_fsp_cepae/PG/PG_PHI_36_accessions.fa_homologs.csv 
