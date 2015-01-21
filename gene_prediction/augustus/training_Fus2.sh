#Fusarium evidence-based gene prediction

# download and organise RNAseq reads
scp -r andrewarmitage@nero.wsbc.warwick.ac.uk:/home/shared/Fusarium/*.zip timecourse/experiment1/.
mv timecourse/experiment1/Data/Intensities/BaseCalls/1_S1_L001_R1_001.fastq.gz timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/0/F/.
mv timecourse/experiment1/Data/Intensities/BaseCalls/1_S1_L001_R2_001.fastq.gz timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/0/R/.
mv timecourse/experiment1/Data/Intensities/BaseCalls/4hr-root_S7_L001_R1_001.fastq.gz timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/4/F/.
mv timecourse/experiment1/Data/Intensities/BaseCalls/4hr-root_S7_L001_R2_001.fastq.gz timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/4/R/.
mv timecourse/experiment1/Data/Intensities/BaseCalls/8hr-root_S8_L001_R1_001.fastq.gz timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/8/F/.
mv timecourse/experiment1/Data/Intensities/BaseCalls/8hr-root_S8_L001_R2_001.fastq.gz timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/8/R/.
mv timecourse/experiment1/Data/Intensities/BaseCalls/3_S2_L001_R1_001.fastq.gz timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/16/F/.
mv timecourse/experiment1/Data/Intensities/BaseCalls/3_S2_L001_R2_001.fastq.gz timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/16/R/.
mv timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/4/F/4_S3_L001_R1_001.fastq.gz timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/24.1/F/.
mv timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/4/R/4_S3_L001_R2_001.fastq.gz timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/24.1/R/.
mv timecourse/experiment1/Data/Intensities/BaseCalls/36hr-root_S10_L001_R1_001.fastq.gz timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/36/F/.
mv timecourse/experiment1/Data/Intensities/BaseCalls/36hr-root_S10_L001_R2_001.fastq.gz timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/36/R/.
mv timecourse/experiment1/Data/Intensities/BaseCalls/6_S4_L001_R1_001.fastq.gz timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/48/F/.
mv timecourse/experiment1/Data/Intensities/BaseCalls/6_S4_L001_R2_001.fastq.gz timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/48/R/.
mv timecourse/experiment1/Data/Intensities/BaseCalls/7_S5_L001_R1_001.fastq.gz timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/72/F/.
mv timecourse/experiment1/Data/Intensities/BaseCalls/7_S5_L001_R2_001.fastq.gz timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/72/R/.
mv timecourse/experiment1/Data/Intensities/BaseCalls/8_S6_L001_R1_001.fastq.gz timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/96/F/.
mv timecourse/experiment1/Data/Intensities/BaseCalls/8_S6_L001_R2_001.fastq.gz timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/96/R/.

# perform qc of RNAseq timecourse data
for FILE_PATH in $(ls -d timecourse/experiment1/F.oxysporum_fsp_cepae/Fus2/*); do 
	echo $FILE_PATH 
	F_IN=$(ls $FILE_PATH/F/*.gz)
	R_IN=$(ls $FILE_PATH/R/*.gz)
	ILLUMINA_ADAPTERS=/home/armita/git_repos/emr_repos/tools/seq_tools/ilumina_full_adapters.fa
	qsub /home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc/rna_qc_fastq-mcf.sh $F_IN $R_IN $ILLUMINA_ADAPTERS RNA
done

# Train a gene model to Fus2
# 1. Extract fusarium reads from RNAseq experiment
# Login to worker node and perform alignment of raw reads to the Fus2 genome

cd /home/groups/harrisonlab/project_files/fusarium
mkdir alignments
cd alignments/
bowtie2-build ../repeat_masked/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_49/Fus2_combined_49_contigs_unmasked.fa Fus2_bowtie_index
bowtie2 -x Fus2_bowtie_index -1 ../qc_rna/paired/Fus2/96/F/96_qc_F.fastq -2 ../qc_rna/paired/Fus2/96/R/96_qc_R.fastq -S 96hr_Fus2.sam -p 8 --al aligned_unpaired.fastq --al-conc aligned_paried.fastq
# 3695824 reads; of these:
#   3695824 (100.00%) were paired; of these:
#     3258388 (88.16%) aligned concordantly 0 times
#     416333 (11.26%) aligned concordantly exactly 1 time
#     21103 (0.57%) aligned concordantly >1 times
#     ----
#     3258388 pairs aligned concordantly 0 times; of these:
#       38039 (1.17%) aligned discordantly 1 time
#     ----
#     3220349 pairs aligned 0 times concordantly or discordantly; of these:
#       6440698 mates make up the pairs; of these:
#         6413654 (99.58%) aligned 0 times
#         23246 (0.36%) aligned exactly 1 time
#         3798 (0.06%) aligned >1 times
# 13.23% overall alignment rate
mkdir -p F.oxysporum_fsp_cepae/Fus2/96
mv * F.oxysporum_fsp_cepae/Fus2/96/.
