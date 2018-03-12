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
<!--
```bash
	RawDatDir=/home/miseq_data/minion/2017/FON-63_2017-05-22
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	mkdir -p $ProjectDir/raw_dna/minion/F.oxysporum_fsp_narcissi/FON_63
```

Sequence data was moved into the appropriate directories

```bash
	RawDatDir=/home/miseq_data/minion/2017/FON-63_2017-05-22		
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	cp $RawDatDir/all_reads_albacore1.1.1.fastq.gz $ProjectDir/raw_dna/minion/F.oxysporum_fsp_narcissi/FON_63/.
``` -->

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
