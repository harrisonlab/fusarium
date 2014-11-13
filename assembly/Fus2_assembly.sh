#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 12
#$ -l virtual_free=6G


# fastq-mcf /home/armita/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/F/s_6_1_sequence.txt raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/R/s_6_2_sequence.txt -o qc_dna/paired/F.oxysporum_fsp_cepae/Fus2/Fus2_70bp_F.fastq -o qc_dna/paired/F.oxysporum_fsp_cepae/Fus2/Fus2_70bp_R.fastq -C 1000000 -u -k 20 -t 0.01 -q 30 -p 5

CUR_DIR=$PWD
FORWARD_FILE1=$1
REVERSE_FILE1=$2
FORWARD_FILE2=$3
REVERSE_FILE2=$4
WORK_DIR=/tmp/Fus2_combined

mkdir -p $WORK_DIR
velveth $WORK_DIR 41,61,2 -fastq -shortPaired -separate $FORWARD_FILE1 $REVERSE_FILE1 -shortPaired2 -separate $FORWARD_FILE2 $REVERSE_FILE2
#velveth $WORK_DIR 61,101,2 -fastq -shortPaired -separate $FORWARD_FILE1 $REVERSE_FILE1 -shortPaired2 -separate $FORWARD_FILE2 $REVERSE_FILE2
for hash_directory in $(ls -d $WORK_DIR*); do
	cd $hash_directory
	velvetg . -exp_cov 133 -cov_cutoff 26 -ins_length1 300 -ins_length2 600 -min_contig_lgth 500
	cd ../
done
cp -r $WORK_DIR* $CUR_DIR/assembly/velvet/F.oxysporum_fsp_cepae/Fus2/.

rm -r $TMPDIR

exit

