#Fus2_Assembly
==========

This document details the commands used to assemble and annotate the Fus2 genome.

Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/fusarium

The following is a summary of the work presented in this Readme.

The following processes were applied to Alternaria genomes prior to analysis:
Data qc
Genome assembly
Repeatmasking
Gene prediction
Functional annotation


#Building of directory structure
```shell
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_3/
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/F/
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/R/
```

Sequence data was moved into the appropriate directories

```shell
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	cp $RawDatDir/fusarium/HAPI_seq_1/FUS2_S2_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/F/.
	cp $RawDatDir/fusarium/HAPI_seq_1/FUS2_S2_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/R/.
	cp $RawDatDir/raw_data/raw_seq/fusarium/warwick_seqs/fus2/s_6_1_sequence.txt $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/F/.
	cp $RawDatDir/raw_data/raw_seq/fusarium/warwick_seqs/fus2/s_6_2_sequence.txt $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/R/.
```

One set of sequence reads were in .txt format and unzipped. These were renamed and zipped.
```shell
	mv $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/F/s_6_1_sequence.txt $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/F/s_6_1_sequence.fastq
	mv $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/R/s_6_2_sequence.txt $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/R/s_6_2_sequence.fastq
	gzip $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/F/s_6_1_sequence.fastq
	gzip $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/R/s_6_2_sequence.fast
```

#Data qc

programs:
  fastqc
  fastq-mcf
  kmc

Data quality was visualised using fastqc:
```shell
	for RawData in $(ls raw_dna/paired/*/*/Fus2/*.fastq.gz); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```
