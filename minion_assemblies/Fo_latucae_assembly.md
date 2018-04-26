# F. oxysporum latucae assemblies
==========

This document details the commands used to assemble and annotate the genomes.

Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/fusarium

The following is a summary of the work presented in this Readme.

The following processes were applied to Fusarium genomes prior to analysis:
Data qc
Genome assembly
Repeatmasking
Gene prediction
Functional annotation


# 0. Building of directory structure

### Minion Data

```bash
	RawDatDir=/data/seq_data/minion/2018/20180222_AJ516/AJ516
  RawDatDir=/data/seq_data/minion/2018/20180418_AJ16/AJ16
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	mkdir -p $ProjectDir/raw_dna/minion/F.oxysporum_fsp_lactucae/AJ516

  RawDatDir=/data/seq_data/minion/2018/20180222_AJ520/AJ520
	# The scond run was output to a subfolder of the 20180222_AJ520 folder
  # RawDatDir=
  ProjectDir=/home/groups/harrisonlab/project_files/fusarium
  mkdir -p $ProjectDir/raw_dna/minion/F.oxysporum_fsp_lactucae/AJ520
```


Data was basecalled again using Albacore 2.02 on the minion server:

```bash
screen -a
ssh nanopore@nanopore

# upgrade albacore
wget https://mirror.oxfordnanoportal.com/software/analysis/ont_albacore-2.2.7-cp34-cp34m-manylinux1_x86_64.whl
~/.local/bin/read_fast5_basecaller.py --version
pip3 install --user ont_albacore-2.1.10-cp34-cp34m-manylinux1_x86_64.whl --upgrade
~/.local/bin/read_fast5_basecaller.py --version

mkdir FoLatucae_26-04-18
cd FoLatucae_26-04-18

# Oxford nanopore AJ516 Run 1
Organism=F.oxysporum_fsp_lactucae
Strain=AJ516
Date=22-02-18
FlowCell="FLO-MIN106"
Kit="SQK-RBK001"
RawDatDir=/data/seq_data/minion/2018/20180222_AJ516/AJ516/GA10000/reads
OutDir=~/FoLatucae_26-04-18/$Organism/$Strain/$Date

mkdir -p $OutDir
cd $OutDir
~/.local/bin/read_fast5_basecaller.py \
  --flowcell $FlowCell \
  --kit $Kit \
  --input $RawDatDir \
  --recursive \
  --worker_threads 24 \
  --save_path albacore_v2.2.7 \
  --output_format fastq,fast5 \
  --reads_per_fastq_batch 4000

	# Oxford nanopore AJ516 Run 2
	Organism=F.oxysporum_fsp_lactucae
	Strain=AJ516
	Date=18-04-18
	FlowCell="FLO-MIN106"
	Kit="SQK-RBK001"
	RawDatDir=/data/seq_data/minion/2018/20180418_AJ16/AJ16/GA10000/reads
	OutDir=/data/scratch/nanopore_tmp_data/FoLatucae_26-04-18/$Organism/$Strain/$Date
	mkdir -p $OutDir

	mkdir -p $OutDir
	cd $OutDir
	~/.local/bin/read_fast5_basecaller.py \
	  --flowcell $FlowCell \
	  --kit $Kit \
	  --input $RawDatDir \
	  --recursive \
	  --worker_threads 24 \
	  --save_path albacore_v2.2.7 \
	  --output_format fastq,fast5 \
	  --reads_per_fastq_batch 4000

		# Oxford nanopore AJ520 Run 1
		Organism=F.oxysporum_fsp_lactucae
		Strain=AJ520
		Date=22-02-18
		FlowCell="FLO-MIN106"
		Kit="SQK-RBK001"
		RawDatDir=/data/seq_data/minion/2018/20180222_AJ520/AJ520/GA20000/reads
		OutDir=/data/scratch/nanopore_tmp_data/FoLatucae_26-04-18/$Organism/$Strain/$Date
		mkdir -p $OutDir

		mkdir -p $OutDir
		cd $OutDir
		~/.local/bin/read_fast5_basecaller.py \
		  --flowcell $FlowCell \
		  --kit $Kit \
		  --input $RawDatDir \
		  --recursive \
		  --worker_threads 24 \
		  --save_path albacore_v2.2.7 \
		  --output_format fastq,fast5 \
		  --reads_per_fastq_batch 4000

			# Oxford nanopore AJ516 Run 2
			Organism=F.oxysporum_fsp_lactucae
			Strain=AJ520
			Date=18-04-18
			FlowCell="FLO-MIN106"
			Kit="SQK-RBK001"
			# The second run of AJ520 reads have not yet been copied over.
			# RawDatDir=/data/seq_data/minion/2018/20180222_AJ520/AJ520/GA30000/reads
			OutDir=/data/scratch/nanopore_tmp_data/FoLatucae_26-04-18/$Organism/$Strain/$Date
			mkdir -p $OutDir

			mkdir -p $OutDir
			cd $OutDir
			~/.local/bin/read_fast5_basecaller.py \
			  --flowcell $FlowCell \
			  --kit $Kit \
			  --input $RawDatDir \
			  --recursive \
			  --worker_threads 24 \
			  --save_path albacore_v2.2.7 \
			  --output_format fastq,fast5 \
			  --reads_per_fastq_batch 4000

	#
  # cat Alt_albacore_v2.10_demultiplexed/workspace/pass/barcode01/*.fastq | gzip -cf > Alt_albacore_v2.10_barcode01.fastq.gz
  # cat Alt_albacore_v2.10_demultiplexed/workspace/pass/barcode02/*.fastq | gzip -cf > Alt_albacore_v2.10_barcode02.fastq.gz
  # OutDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.2.7
  # mkdir -p $OutDir
  # chmod +w $OutDir
  # cp Alt_albacore_v2.10_barcode0*.fastq.gz $OutDir/.
  # chmod +rw $OutDir/Alt_albacore_v2.10_barcode0*.fastq.gz
  # scp Alt_albacore_v2.10_barcode*.fastq.gz armita@192.168.1.200:$OutDir/.
  tar -cz -f Alt_albacore_v2.10_demultiplexed.tar.gz Alt_albacore_v2.10_demultiplexed
  OutDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.2.7
  mv Alt_albacore_v2.10_demultiplexed.tar.gz $OutDir/.
  chmod +rw $OutDir/Alt_albacore_v2.10_demultiplexed.tar.gz
```

Sequence data was moved into the appropriate directories

```bash
	RawDatDir=/home/miseq_data/minion/2017/FON-63_2017-05-22		
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	cp $RawDatDir/all_reads_albacore1.1.1.fastq.gz $ProjectDir/raw_dna/minion/F.oxysporum_fsp_narcissi/FON_63/.
```

### MiSeq data

```bash
  RawDatDir=/data/seq_data/miseq/2018/RAW/180213_M04465_0067_000000000-BJ4DR/Data/Intensities/BaseCalls
  ProjectDir=/home/groups/harrisonlab/project_files/fusarium
  OutDir=$ProjectDir/raw_dna/paired/F.oxysporum_fsp_lactucae/AJ516
  mkdir -p $OutDir/F
  mkdir -p $OutDir/R
  cd $OutDir/F
  cp -s $RawDatDir/AJ516_S1_L001_R1_001.fastq.gz .
  cd $OutDir/R
  cp -s $RawDatDir/AJ516_S1_L001_R2_001.fastq.gz .
  OutDir=$ProjectDir/raw_dna/paired/F.oxysporum_fsp_lactucae/AJ520
  mkdir -p $OutDir/F
  mkdir -p $OutDir/R
  cd $OutDir/F
  cp -s $RawDatDir/AJ520_S2_L001_R1_001.fastq.gz .
  cd $OutDir/R
  cp -s $RawDatDir/AJ520_S2_L001_R2_001.fastq.gz .
  cd $ProjectDir
```


#### QC of MiSeq data

programs:
  fastqc
  fastq-mcf
  kmc

Data quality was visualised using fastqc:
```bash
	for RawData in $(ls raw_dna/paired/*/*/*/*.fastq.gz | grep 'lactucae'); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

```bash
	for StrainPath in $(ls -d raw_dna/paired/*/* | grep 'lactucae'); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
		IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
		ReadsF=$(ls $StrainPath/F/*.fastq*)
		ReadsR=$(ls $StrainPath/R/*.fastq*)
		echo $ReadsF
		echo $ReadsR
		qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	done
```

# Identify sequencing coverage
<!--
For Minion data:
```bash
	for RawData in $(ls qc_dna/minion/*/*/*q.gz | grep 'lactucae'); do
		echo $RawData;
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
		GenomeSz=60
		OutDir=$(dirname $RawData)
		qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
	done
```

```bash
  for StrainDir in $(ls -d qc_dna/minion/*/* | grep 'Stocks4'); do
    Strain=$(basename $StrainDir)
    printf "$Strain\t"
    for File in $(ls $StrainDir/*.txt); do
      echo $(basename $File);
      cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
    done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
  done
```
MinION coverage was:
```
Stocks4	65.05
``` -->

For Miseq data:
```bash
	for RawData in $(ls qc_dna/paired/*/*/*/*q.gz | grep 'lactucae'); do
		echo $RawData;
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
		qsub $ProgDir/run_fastqc.sh $RawData;
		GenomeSz=60
		OutDir=$(dirname $RawData)
		qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
	done
```

```bash
	for StrainDir in $(ls -d qc_dna/paired/*/* | grep 'lactucae'); do
		Strain=$(basename $StrainDir)
		printf "$Strain\t"
		for File in $(ls $StrainDir/*/*.txt); do
			echo $(basename $File);
			cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
		done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
	done
```

Miseq coverage was:
```
Stocks4	58.66
```
