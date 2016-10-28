#HAPI
==========

Scripts used for the analysis of Fusarium genomes
Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/fusarium

The following is a summary of the work presented in this Readme.

The following processes were applied to Fusarium genomes prior to analysis:
Data qc
Genome assembly
Repeatmasking
Gene prediction
Functional annotation

Analyses performed on these genomes involved BLAST searching for:


ls contigs were identified using:
Alignment of raw reads to assembled genomes
Assembly of remaining reads


#Building of directory structure
```bash
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_3/
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	mkdir -p $ProjectDir/raw_dna/paired/F.proliferatum/A8/F
	mkdir -p $ProjectDir/raw_dna/paired/F.proliferatum/A8/R
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_narcissi/N139/F
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_narcissi/N139/R
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG3/F
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG3/R
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG18/F
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG18/R
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160115_M01678_0041_000000000-AEMMK
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/FOP1/F
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/FOP1/R
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/FOP5/F
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/FOP5/R
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/L5/F
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/L5/R
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/HB6/F
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/HB6/R
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160119_M016878_0042-AGDN7
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/A1-2/F
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/A1-2/R
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/A13/F
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/A13/R
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/CB3/F
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/CB3/R
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/HB6/F
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/HB6/R
	mkdir -p $ProjectDir/raw_dna/paired/F.avenaceum/PG8/F
	mkdir -p $ProjectDir/raw_dna/paired/F.avenaceum/PG8/R
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160401_M004465_0007-AGKF2
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/FOP2/F
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/FOP2/R
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160725_M04465_0020-AP1N0
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/AF064/F
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/AF064/R
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/CH5-2/F
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/CH5-2/R
```

Sequence data was moved into the appropriate directories

```bash
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_3
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	cp $RawDatDir/Fproliferatum_S1_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.proliferatum/A8/F/.
	cp $RawDatDir/Fproliferatum_S1_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.proliferatum/A8/R/.
	cp $RawDatDir/FoxysporumN139_S2_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_narcissi/N139/F/.
	cp $RawDatDir/FoxysporumN139_S2_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_narcissi/N139/R/.
	cp $RawDatDir/FoxysporumPG3_S3_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG3/F/.
	cp $RawDatDir/FoxysporumPG3_S3_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG3/R/.
	cp $RawDatDir/FoxysporumPG18_S4_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG18/F/.
	cp $RawDatDir/FoxysporumPG18_S4_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG18/R/.
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160115_M01678_0041_000000000-AEMMK
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	cp $RawDatDir/FOP1_S1_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/FOP1/F/.
	cp $RawDatDir/FOP1_S1_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/FOP1/R/.
	cp $RawDatDir/FOP5_S2_L001_R1_001.fastq.gz  $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/FOP5/F/.
	cp $RawDatDir/FOP5_S2_L001_R2_001.fastq.gz  $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/FOP5/R/.
	cp $RawDatDir/L5_S3_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/L5/F/.
	cp $RawDatDir/L5_S3_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/L5/R/.
	cp $RawDatDir/HB6_S4_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/HB6/F/.
	cp $RawDatDir/HB6_S4_L001_R2_001.fastq.gz  $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/HB6/R/.
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160119_M016878_0042-AGDN7
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	cp $RawDatDir/A12_S1_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/A1-2/F/.
	cp $RawDatDir/A12_S1_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/A1-2/R/.
	cp $RawDatDir/A13_S2_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/A13/F/.
	cp $RawDatDir/A13_S2_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/A13/R/.
	cp $RawDatDir/CB3_S3_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/CB3/F/.
	cp $RawDatDir/CB3_S3_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/CB3/R/.
	cp $RawDatDir/HB6_S5_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/HB6/F/.
	cp $RawDatDir/HB6_S5_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/HB6/R/.
	cp $RawDatDir/PG8_S4_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.avenaceum/PG8/F/.
	cp $RawDatDir/PG8_S4_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.avenaceum/PG8/R/.
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160401_M004465_0007-AGKF2
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	cp $RawDatDir/FOP2_S1_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/FOP2/F/.
	cp $RawDatDir/FOP2_S1_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/FOP2/R/.
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/raw_reads/160725_M04465_0020-AP1N0
	cp $RawDatDir/AF064_S3_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/AF064/F/.
	cp $RawDatDir/AF064_S3_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/AF064/R/.
	cp $RawDatDir/CH5-2_S2_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/CH5-2/F/.
	cp $RawDatDir/CH5-2_S2_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/CH5-2/R/.
```

This process was repeated for RNAseq data:

```bash
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/fusarium/rna_seq
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDB/F
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDB/R
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_CzapekDox/F
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_CzapekDox/R
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_GlucosePeptone/F
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_GlucosePeptone/R
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDA/F
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDA/R

	cp $RawDatDir/4_S1_L001_R1_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDB/F/.
	cp $RawDatDir/4_S1_L001_R2_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDB/R/.
	cp $RawDatDir/6_S2_L001_R1_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_CzapekDox/F/.
	cp $RawDatDir/6_S2_L001_R2_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_CzapekDox/R/.
	cp $RawDatDir/7_S3_L001_R1_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_GlucosePeptone/F/.
	cp $RawDatDir/7_S3_L001_R2_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_GlucosePeptone/R/.
	cp $RawDatDir/9_S4_L001_R1_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDA/F/.
	cp $RawDatDir/9_S4_L001_R2_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDA/R/.
```



#Data qc

programs:
  fastqc
  fastq-mcf
  kmc

Data quality was visualised using fastqc:
```bash
	for RawData in $(ls raw_dna/paired/*/*/*/*.fastq.gz); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

Trimming was performed on data to trim adapters from
sequences and remove poor quality data. This was done with fastq-mcf

Firstly, those strains with more than one run were identified:

```bash
	for Strain in $(ls -d raw_dna/paired/*/*); do
	NumReads=$(ls $Strain/F/*.gz | wc -l);
		if [ $NumReads -gt 1 ]; then
			echo "$Strain";
			echo "$NumReads";
		fi;
	done
```

```
	raw_dna/paired/F.oxysporum_fsp_cepae/Fus2
	2
	raw_dna/paired/F.oxysporum_fsp_cepae/HB6
	2
```

Trimming was first performed on all strains that had a single run of data:

```bash
	for StrainPath in $(ls -d raw_dna/paired/*/* | grep -v -e 'Fus2' -e 'HB6'); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
		IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
		ReadsF=$(ls $StrainPath/F/*.fastq*)
		ReadsR=$(ls $StrainPath/R/*.fastq*)
		echo $ReadsF
		echo $ReadsR
		qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	done
```

Trimming was then performed for strains with multiple runs of data

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
	IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
	echo "Fus2"
	StrainPath=raw_dna/paired/F.oxysporum_fsp_cepae/Fus2
	ReadsF=$(ls $StrainPath/F/s_6_1_sequence.fastq.gz)
	ReadsR=$(ls $StrainPath/R/s_6_2_sequence.fastq.gz)
	qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	StrainPath=raw_dna/paired/F.oxysporum_fsp_cepae/Fus2
	ReadsF=$(ls $StrainPath/F/FUS2_S2_L001_R1_001.fastq.gz)
	ReadsR=$(ls $StrainPath/R/FUS2_S2_L001_R2_001.fastq.gz)
	qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	echo "HB6"
	StrainPath=raw_dna/paired/F.oxysporum_fsp_cepae/HB6
	ReadsF=$(ls $StrainPath/F/HB6_S4_L001_R1_001.fastq.gz)
	ReadsR=$(ls $StrainPath/R/HB6_S4_L001_R2_001.fastq.gz)
	qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	StrainPath=raw_dna/paired/F.oxysporum_fsp_cepae/HB6
	ReadsF=$(ls $StrainPath/F/HB6_S5_L001_R1_001.fastq.gz)
	ReadsR=$(ls $StrainPath/R/HB6_S5_L001_R2_001.fastq.gz)
	qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
```


Data quality was visualised once again following trimming:
```bash
	for RawData in $(ls qc_dna/paired/*/*/*/*.fq.gz | grep 'FOP2'); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

kmer counting was performed using kmc
This allowed estimation of sequencing depth and total genome size

This was performed for strains with single runs of data

```bash
	for TrimPath in $(ls -d raw_dna/paired/*/* | grep -v -e 'Fus2' -e 'HB6' | grep 'FOP2'); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		TrimF=$(ls $TrimPath/F/*.fastq*)
		TrimR=$(ls $TrimPath/R/*.fastq*)
		echo $TrimF
		echo $TrimR
		qsub $ProgDir/kmc_kmer_counting.sh $TrimF $TrimR
	done
```

and for strains with muiltiple runs of data:

```bash
	for TrimPath in $(ls -d raw_dna/paired/*/Fus2); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		TrimF1=$(ls $TrimPath/F/s_6_1_sequence.fastq.gz)
		TrimR1=$(ls $TrimPath/R/s_6_2_sequence.fastq.gz)
		echo $TrimF1
		echo $TrimR1
		TrimF2=$(ls $TrimPath/F/FUS2_S2_L001_R1_001.fastq.gz)
		TrimR2=$(ls $TrimPath/R/FUS2_S2_L001_R2_001.fastq.gz)
		echo $TrimF2
		echo $TrimR2
		qsub $ProgDir/kmc_kmer_counting.sh $TrimF1 $TrimR1 $TrimF2 $TrimR2
	done
	for TrimPath in $(ls -d raw_dna/paired/*/HB6); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		TrimF1=$(ls $TrimPath/F/HB6_S4_L001_R1_001.fastq.gz)
		TrimR1=$(ls $TrimPath/R/HB6_S4_L001_R2_001.fastq.gz)
		echo $TrimF1
		echo $TrimR1
		TrimF2=$(ls $TrimPath/F/HB6_S5_L001_R1_001.fastq.gz)
		TrimR2=$(ls $TrimPath/R/HB6_S5_L001_R2_001.fastq.gz)
		echo $TrimF2
		echo $TrimR2
		qsub $ProgDir/kmc_kmer_counting.sh $TrimF1 $TrimR1 $TrimF2 $TrimR2
	done
```

mode kmer abundance prior to error correction was reported using the following
commands:

```bash
for File in $(ls qc_dna/kmc/*/*/*_true_kmer_summary.txt); do
basename $File;
tail -n3 $File | head -n1 ;
done
```

```
	PG8_true_kmer_summary.txt
	The mode kmer abundance is:  53 < fail
	125_true_kmer_summary.txt
	The mode kmer abundance is:  46 < pass
	55_true_kmer_summary.txt
	The mode kmer abundance is:  29 < pass
	A1-2_true_kmer_summary.txt
	The mode kmer abundance is:  16 < pass
	A13_true_kmer_summary.txt
	The mode kmer abundance is:  22 < pass
	A23_true_kmer_summary.txt
	The mode kmer abundance is:  33 < pass
	A28_true_kmer_summary.txt
	The mode kmer abundance is:  36 < pass
	CB3_true_kmer_summary.txt
	The mode kmer abundance is:  21 < pass
	D2_true_kmer_summary.txt
	The mode kmer abundance is:  11 < pass
	Fus2_true_kmer_summary.txt
	The mode kmer abundance is:  109 < 2lib
	HB17_true_kmer_summary.txt
	The mode kmer abundance is:  27 < pass
	HB6_true_kmer_summary.txt
	The mode kmer abundance is:  91 < 2 lib
	PG_true_kmer_summary.txt
	The mode kmer abundance is:  58 < pass
	N139_true_kmer_summary.txt
	The mode kmer abundance is:  26 < pass
	FOP1_true_kmer_summary.txt
	The mode kmer abundance is:  32 <- fail
	FOP5_true_kmer_summary.txt
	The mode kmer abundance is:  62 <- fail
	L5_true_kmer_summary.txt
	The mode kmer abundance is:  35 <- pass
	PG18_true_kmer_summary.txt
	The mode kmer abundance is:  24 <- pass
	PG3_true_kmer_summary.txt
	The mode kmer abundance is:  28 <- pass
	A8_true_kmer_summary.txt
	The mode kmer abundance is:  21 <- pass
```

#Assembly

Assembly was performed with:
* Spades

## Spades Assembly



```bash
	for StrainPath in $(ls -d qc_dna/paired/*/* | grep -v -e 'Fus2' -e 'HB6' | grep 'FOP2'); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
		Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
		Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
		F_Read=$(ls $StrainPath/F/*.fq.gz)
		R_Read=$(ls $StrainPath/R/*.fq.gz)
		OutDir=assembly/spades/$Organism/$Strain
		Jobs=$(qstat | grep 'submit_SPA' | grep 'qw' | wc -l)
		while [ $Jobs -gt 1 ]; do
			sleep 5m
			printf "."
			Jobs=$(qstat | grep 'submit_SPA' | grep 'qw' | wc -l)
		done		
		printf "\n"
		echo $F_Read
		echo $R_Read
		qsub $ProgDir/submit_SPAdes.sh $F_Read $R_Read $OutDir correct 10
	done
```

Assembly for PG8, FOP1 and FOP5 failed due to a lack of memory, as such the assembly was
resubmitted with more RAM.

```bash
	for StrainPath in $(ls -d qc_dna/paired/*/* | grep -e 'FOP5' -e 'PG8' -e 'FOP1'); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
		Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
		Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
		F_Read=$(ls $StrainPath/F/*.fq.gz)
		R_Read=$(ls $StrainPath/R/*.fq.gz)
		OutDir=assembly/spades/$Organism/$Strain
		echo $F_Read
		echo $R_Read
		qsub $ProgDir/submit_SPAdes_HiMem.sh $F_Read $R_Read $OutDir correct 10
	done
```

Assemblies were submitted for genomes with data from multiple sequencing runs:

```bash
for StrainPath in $(ls -d qc_dna/paired/F.*/HB6); do
  echo $StrainPath
    ProgDir=/home/ransoe/git_repos/tools/seq_tools/assemblers/spades/multiple_libraries
    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
    echo $Strain
    echo $Organism
    TrimF1_Read=$(ls $StrainPath/F/HB6_S4_L001_R1_001_trim.fq.gz);
    TrimR1_Read=$(ls $StrainPath/R/HB6_S4_L001_R2_001_trim.fq.gz);
    TrimF2_Read=$(ls $StrainPath/F/HB6_S5_L001_R1_001_trim.fq.gz);
    TrimR2_Read=$(ls $StrainPath/R/HB6_S5_L001_R2_001_trim.fq.gz);
    echo $TrimF1_Read
    echo $TrimR1_Read
    echo $TrimF2_Read
    echo $TrimR2_Read
    OutDir=assembly/spades/$Organism/$Strain
    qsub $ProgDir/subSpades_2lib.sh $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $OutDir correct 10
  done
	for StrainPath in $(ls -d qc_dna/paired/F.*/Fus2); do
	  echo $StrainPath
	    ProgDir=/home/ransoe/git_repos/tools/seq_tools/assemblers/spades/multiple_libraries
	    Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
	    Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
	    echo $Strain
	    echo $Organism
	    TrimF1_Read=$(ls $StrainPath/F/s_6_1_sequence_trim.fq.gz);
	    TrimR1_Read=$(ls $StrainPath/R/s_6_2_sequence_trim.fq.gz);
	    TrimF2_Read=$(ls $StrainPath/F/FUS2_S2_L001_R1_001_trim.fq.gz);
	    TrimR2_Read=$(ls $StrainPath/R/FUS2_S2_L001_R2_001_trim.fq.gz);
	    echo $TrimF1_Read
	    echo $TrimR1_Read
	    echo $TrimF2_Read
	    echo $TrimR2_Read
	    OutDir=assembly/spades/$Organism/$Strain
	    qsub $ProgDir/subSpades_2lib_HiMem.sh $TrimF1_Read $TrimR1_Read $TrimF2_Read $TrimR2_Read $OutDir correct 10
	  done
```

Quast

```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta | grep 'FOP2'); do
	# for Assembly in $(ls assembly/spades/*/*/*/contigs_min_500bp_renamed.fasta | grep 'Fus2_edited_v2'); do
  # for Assembly in $(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
    qsub $ProgDir/sub_quast.sh $Assembly $OutDir
  done
```

<!-- The results of quast were shown using the following commands:

```bash
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs/report.txt); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev);
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev);
    echo;
    echo $Organism;
    echo $Strain;
    cat $Assembly;
  done
```

The output of this analysis is in the assembly/quast_results.txt file of this
git repository. -->


Contigs were renamed in accordance with ncbi recomendations.

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
	for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```


A Bioproject and Biosample was made with NCBI genbank for submission of genomes.
Following the creation of these submissions, the .fasta assembly was uploaded
through the submission portal. A note was provided requesting that the assembly
be run through the contamination screen to aid a more detailed resubmission in
future. The returned FCSreport.txt was downloaded from the NCBI webportal and
used to correct the assembly to NCBI standards.

NCBI reports (FCSreport.txt) were manually downloaded to the following loactions:

```bash
	for Assembly in $(ls assembly/spades/*/*/*/contigs_min_500bp_renamed.fasta | grep -w -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139'); do
		Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
		NCBI_report_dir=genome_submission/$Organism/$Strain/initial_submission
		mkdir -p $NCBI_report_dir
	done
```

These downloaded files were used to correct assemblies:

```bash
for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp_renamed.fasta | grep -w -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
NCBI_report=$(ls genome_submission/$Organism/$Strain/initial_submission/FCSreport.txt)
OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
mkdir -p $OutDir
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file $NCBI_report > $OutDir/log.txt
done
```

Quast was used to collect details on these assemblies

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta | grep -w -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
echo "$Organism - $Strain"
OutDir=assembly/spades/$Organism/$Strain/ncbi_edits
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

<!--
The Fus2 was manually edited as SIX9 was noted to be split over 2 contigs. These
contigs were manually joined in geneious and exported back to the cluster
at the location indicated below. These contigs were renamed.

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
  touch tmp.csv
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs/*_edited_*.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
``` -->




# Repeatmasking

Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
for BestAss in $(ls assembly/spades/*/*/ncbi_edits/contigs_min_500bp_renamed.fasta | grep -w -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139'); do
# for BestAss in $(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs.fasta); do
Strain=$(echo $BestAss | rev | cut -f3 -d '/' | rev)
Organism=$(echo $BestAss | rev | cut -f4 -d '/' | rev)
OutDir=repeat_masked/$Organism/"$Strain"_ncbi/ncbi_submission
qsub $ProgDir/rep_modeling.sh $BestAss $OutDir
qsub $ProgDir/transposonPSI.sh $BestAss $OutDir
done
```

The published non-pathogen genome for isolate FO47 was also repeatmasked as
this isolate was also used in experimental work.

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
	Fo47Ass=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs.fasta
	for BestAss in $(ls $Fo47Ass); do
		qsub $ProgDir/rep_modeling.sh $BestAss
		qsub $ProgDir/transposonPSI.sh $BestAss
	done
```

Fo47 masked contigs were noted to still contain additional text following spaces.
This disrupted downstream programs.
```bash
for File in $(ls repeat_masked/F.oxysporum/fo47/*/fo47_contigs_*masked.fa); do
sed -i 's/ of Fusarium oxysporum Fo47//g' $File
done
```

The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and harmasked files.

```bash
for File in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep 'ncbi' | grep -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/*/*_contigs_hardmasked.fa | grep 'ncbi' | grep -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
# echo "Number of masked bases:"
# cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[N]", ".")}' | cut -f2 -d ' '
done
```

```
repeat_masked/F.oxysporum_fsp_cepae/125_ncbi/ncbi_submission/125_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
3967501
repeat_masked/F.oxysporum_fsp_cepae/A13_ncbi/ncbi_submission/A13_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
4946858
repeat_masked/F.oxysporum_fsp_cepae/A23_ncbi/ncbi_submission/A23_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
3546886
repeat_masked/F.oxysporum_fsp_cepae/A28_ncbi/ncbi_submission/A28_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
4456082
repeat_masked/F.oxysporum_fsp_cepae/CB3_ncbi/ncbi_submission/CB3_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
2927012
repeat_masked/F.oxysporum_fsp_cepae/PG_ncbi/ncbi_submission/PG_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
3038447
repeat_masked/F.oxysporum_fsp_narcissi/N139_ncbi/ncbi_submission/N139_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
5671216
repeat_masked/F.proliferatum/A8_ncbi/ncbi_submission/A8_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
1124505
```

The number of bases masked by transposonPSI and Repeatmasker were summarised
using the following commands:

```bash
for RepDir in $(ls -d repeat_masked/F.*/*/ncbi_submission | grep 'ncbi' | grep -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139' -e 'Fus2_canu_new'); do
Strain=$(echo $RepDir | rev | cut -f2 -d '/' | rev)
Organism=$(echo $RepDir | rev | cut -f3 -d '/' | rev)  
RepMaskGff=$(ls $RepDir/*_contigs_hardmasked.gff)
TransPSIGff=$(ls $RepDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
printf "$Organism\t$Strain\n"
printf "The number of bases masked by RepeatMasker:\t"
sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
printf "The number of bases masked by TransposonPSI:\t"
sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
printf "The total number of masked bases are:\t"
cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
echo
done
```
<!--
```
F.oxysporum_fsp_cepae   A13_ncbi
The number of bases masked by RepeatMasker:     4673077
The number of bases masked by TransposonPSI:    1555142
The total number of masked bases are:   4946858

F.oxysporum_fsp_cepae   A23_ncbi
The number of bases masked by RepeatMasker:     3279757
The number of bases masked by TransposonPSI:    1061656
The total number of masked bases are:   3546886

F.oxysporum_fsp_cepae   A28_ncbi
The number of bases masked by RepeatMasker:     4149965
The number of bases masked by TransposonPSI:    1190026
The total number of masked bases are:   4456082

F.oxysporum_fsp_cepae   CB3_ncbi
The number of bases masked by RepeatMasker:     2686212
The number of bases masked by TransposonPSI:    842157
The total number of masked bases are:   2927012

F.oxysporum_fsp_cepae   PG_ncbi
The number of bases masked by RepeatMasker:     2827891
The number of bases masked by TransposonPSI:    865813
The total number of masked bases are:   3038447

F.oxysporum_fsp_narcissi        N139_ncbi
The number of bases masked by RepeatMasker:     5377773
The number of bases masked by TransposonPSI:    1642358
The total number of masked bases are:   5671216

F.proliferatum  A8_ncbi
The number of bases masked by RepeatMasker:     924627
The number of bases masked by TransposonPSI:    258296
The total number of masked bases are:   1124505


``` -->

# Gene Prediction


Gene prediction followed three steps:
	Pre-gene prediction
		- Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
	Gene model training
		- Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline
	Gene prediction
		- Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.

# Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.
```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
cd /home/groups/harrisonlab/project_files/fusarium
for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa | grep -w -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139'); do
#for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa | grep -w 'fo47'); do
echo $Genome;
qsub $ProgDir/sub_cegma.sh $Genome dna;
done
```

Outputs were summarised using the commands:
```bash
	for File in $(ls gene_pred/cegma/F*/*/*_dna_cegma.completeness_report | grep -w -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139'); do
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev);
		Species=$(echo $File | rev | cut -f3 -d '/' | rev);
		printf "$Species\t$Strain\n";
		cat $File | head -n18 | tail -n+4;printf "\n";
	done > gene_pred/cegma/cegma_results_dna_summary.txt

	less gene_pred/cegma/cegma_results_dna_summary.txt
```


#Gene prediction

Gene prediction was performed for Fusarium genomes. Two gene prediction
approaches were used:

Gene prediction using Braker1
Prediction of all putative ORFs in the genome using the ORF finder (atg.pl)
approach.


## Gene prediction 1 - Braker1 gene model training and prediction

Gene prediction was performed using Braker1.

First, RNAseq data was aligned to Fusarium genomes.
* qc of RNA seq data is detailed below:

```bash
 for Folder in $(ls -d raw_rna/paired/F.oxysporum_fsp_cepae/*); do
	 FolderName=$(echo $Folder | rev | cut -f1 -d '/' | rev);
	 echo $FolderName;
	 ls -lh $Folder/F;
	 ls -lh $Folder/R;
	done
```
This contained the following data:
```
	55_72hrs_rep1
	sample013_1.combined.fastq.gz
	sample013_2.combined.fastq.gz
	55_72hrs_rep2
	sample014_1.combined.fastq.gz
	sample014_2.combined.fastq.gz
	55_72hrs_rep3
	sample015_1.combined.fastq.gz
	sample015_2.combined.fastq.gz
	control_72hrs_rep1
	sample002_1.combined.fastq.gz
	sample002_2.combined.fastq.gz
	control_72hrs_rep2
	sample004_1.combined.fastq.gz
	sample004_2.combined.fastq.gz
	control_72hrs_rep3
	sample005_1.combined.fastq.gz
	sample005_2.combined.fastq.gz
	FO47_72hrs_rep1
	sample006_1.combined.fastq.gz
	sample006_2.combined.fastq.gz
	FO47_72hrs_rep2
	sample007_1.combined.fastq.gz
	sample007_2.combined.fastq.gz
	FO47_72hrs_rep3
	sample012_1.combined.fastq.gz
	sample012_2.combined.fastq.gz
	Fus2_0hrs_prelim
	1_S1_L001_R1_001.fastq.gz
	1_S1_L001_R2_001.fastq.gz
	Fus2_16hrs_prelim
	3_S2_L001_R1_001.fastq.gz
	3_S2_L001_R2_001.fastq.gz
	Fus2_24hrs_prelim_rep1
	4_S3_L001_R1_001.fastq.gz
	4_S3_L001_R2_001.fastq.gz
	Fus2_24hrs_prelim_rep2
	Fus2_36hrs_prelim
	36hr-root_S10_L001_R1_001.fastq.gz
	36hr-root_S10_L001_R2_001.fastq.gz
	Fus2_48hrs_prelim
	6_S4_L001_R1_001.fastq.gz
	6_S4_L001_R2_001.fastq.gz
	Fus2_4hrs_prelim
	4hr-root_S7_L001_R1_001.fastq.gz
	4hr-root_S7_L001_R2_001.fastq.gz
	Fus2_72hrs_prelim
	7_S5_L001_R1_001.fastq.gz
	7_S5_L001_R2_001.fastq.gz
	Fus2_72hrs_rep1
	sample016_1.combined.fastq.gz
	sample016_2.combined.fastq.gz
	Fus2_72hrs_rep2
	sample018_1.combined.fastq.gz
	sample018_2.combined.fastq.gz
	Fus2_72hrs_rep3
	sample019_1.combined.fastq.gz
	sample019_2.combined.fastq.gz
	Fus2_8hrs_prelim
	8hr-root_S8_L001_R1_001.fastq.gz
	8hr-root_S8_L001_R2_001.fastq.gz
	Fus2_96hrs_prelim
	8_S6_L001_R1_001.fastq.gz
	8_S6_L001_R2_001.fastq.gz
	Fus2_CzapekDox
	6_S2_L001_R1_001_fastqc  6_S2_L001_R1_001.fastq.gz
	6_S2_L001_R2_001_fastqc  6_S2_L001_R2_001.fastq.gz
	Fus2_GlucosePeptone
	7_S3_L001_R1_001_fastqc  7_S3_L001_R1_001.fastq.gz
	7_S3_L001_R2_001_fastqc  7_S3_L001_R2_001.fastq.gz
	Fus2_PDA
	9_S4_L001_R1_001_fastqc  9_S4_L001_R1_001.fastq.gz
	9_S4_L001_R2_001_fastqc  9_S4_L001_R2_001.fastq.gz
	Fus2_PDB
	4_S1_L001_R1_001_fastqc  4_S1_L001_R1_001.fastq.gz
	4_S1_L001_R2_001_fastqc  4_S1_L001_R2_001.fastq.gz
```

Perform qc of RNAseq timecourse data
```bash
	for num in $(qstat | grep 'rna_qc' | cut -f1 -d ' '); do
		Treatment=$(head -n1 rna_qc_fastq-mcf.sh.o$num | cut -f10 -d '/')
		for FilePath in $(raw_rna/paired/F.oxysporum_fsp_cepae/55_72hrs_rep1); do
			echo $FilePath
			FileF=$(ls $FilePath/F/*.gz)
			FileR=$(ls $FilePath/R/*.gz)
			IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
			ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
			qsub $ProgDir/rna_qc_fastq-mcf.sh $FileF $FileR $IlluminaAdapters RNA
			sleep 10
		done
	done
```

Data quality was visualised using fastqc:
```bash
	for RawData in $(ls qc_rna/paired/F.oxysporum_fsp_cepae/*/*/*.fq.gz); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

#### Aligning

Insert sizes of the RNA seq library were unknown until a draft alignment could
be made. To do this tophat and cufflinks were run, aligning the reads against a
single genome. The fragment length and stdev were printed to stdout while
cufflinks was running.

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep 'FOP2'); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		for RNADir in $(ls -d qc_rna/paired/F.oxysporum_fsp_cepae/*); do
			Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
			echo "$Timepoint"
			FileF=$(ls $RNADir/F/*_trim.fq.gz)
			FileR=$(ls $RNADir/R/*_trim.fq.gz)
			OutDir=alignment/$Organism/$Strain/$Timepoint
			ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
			qsub $ProgDir/tophat_alignment.sh $Assembly $FileF $FileR $OutDir
		done
	done
```
Alignments were concatenated prior to running cufflinks:
Cufflinks was run to produce the fragment length and stdev statistics:

```bash
	for Assembly in $(ls repeat_masked/*/Fo47/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep 'FOP2'); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
		OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim
		echo "$Organism - $Strain"
		mkdir -p $OutDir
		samtools merge -f $AcceptedHits \
		alignment/$Organism/$Strain/55_72hrs_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/55_72hrs_rep2/accepted_hits.bam \
		alignment/$Organism/$Strain/55_72hrs_rep3/accepted_hits.bam \
		alignment/$Organism/$Strain/FO47_72hrs_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/FO47_72hrs_rep2/accepted_hits.bam \
		alignment/$Organism/$Strain/FO47_72hrs_rep3/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_0hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_16hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_24hrs_prelim_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_36hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_48hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_4hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_rep2/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_rep3/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_8hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_96hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_CzapekDox/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_GlucosePeptone/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_PDA/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_PDB/accepted_hits.bam
		cufflinks -o $OutDir/cufflinks -p 8 --max-intron-length 4000 $AcceptedHits 2>&1 | tee $OutDir/cufflinks/cufflinks.log
	done
```

Output from stdout included:
```
	Processed 22484 loci.                        [*************************] 100%
	Map Properties:
	     Normalized Map Mass: 50507412.55
	     Raw Map Mass: 50507412.55
	     Fragment Length Distribution: Empirical (learned)
	                   Estimated Mean: 181.98
	                Estimated Std Dev: 78.39
	[13:02:48] Assembling transcripts and estimating abundances.
	Processed 22506 loci.                        [*************************] 100%
```

The Estimated Mean: 181.98 allowed calculation of of the mean insert gap to be
-20bp 182-(97*2) where 97 was the mean read length. This was provided to tophat
on a second run (as the -r option) along with the fragment length stdev to
increase the accuracy of mapping.


Then Rnaseq data was aligned to each genome assembly:

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139' | grep 'ncbi'); do
# for Assembly in $(ls assembly/merged_canu_spades/*/Fus2/filtered_contigs/Fus2_contigs_renamed.fasta); do
# for Assembly in $(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs.fasta); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for RNADir in $(ls -d qc_rna/paired/F.oxysporum_fsp_cepae/* | grep -e 'FO47_72hrs_rep2'); do
Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
echo "$Timepoint"
FileF=$(ls $RNADir/F/*_trim.fq.gz)
FileR=$(ls $RNADir/R/*_trim.fq.gz)
OutDir=alignment/$Organism/$Strain/$Timepoint
InsertGap='-20'
InsertStdDev='78'
Jobs=$(qstat | grep 'tophat' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'tophat' | grep 'qw' | wc -l)
done
printf "\n"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/tophat_alignment.sh $Assembly $FileF $FileR $OutDir $InsertGap $InsertStdDev
done
done
```


#### Braker prediction

Before braker predictiction was performed, I double checked that I had the
genemark key in my user area and copied it over from the genemark install
directory:

```bash
	ls ~/.gm_key
	cp /home/armita/prog/genemark/gm_key_64 ~/.gm_key
```

<!-- ```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep 'ncbi' | grep -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139'); do
Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
done
printf "\n"
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
for File in alignment/$Organism/$Strain/55_72hrs_rep1/accepted_hits.bam alignment/$Organism/$Strain/55_72hrs_rep2/accepted_hits.bam alignment/$Organism/$Strain/55_72hrs_rep3/accepted_hits.bam alignment/$Organism/$Strain/FO47_72hrs_rep1/accepted_hits.bam alignment/$Organism/$Strain/FO47_72hrs_rep2/accepted_hits.bam alignment/$Organism/$Strain/FO47_72hrs_rep3/accepted_hits.bam alignment/$Organism/$Strain/Fus2_0hrs_prelim/accepted_hits.bam alignment/$Organism/$Strain/Fus2_16hrs_prelim/accepted_hits.bam alignment/$Organism/$Strain/Fus2_24hrs_prelim_rep1/accepted_hits.bam alignment/$Organism/$Strain/Fus2_36hrs_prelim/accepted_hits.bam alignment/$Organism/$Strain/Fus2_48hrs_prelim/accepted_hits.bam alignment/$Organism/$Strain/Fus2_4hrs_prelim/accepted_hits.bam alignment/$Organism/$Strain/Fus2_72hrs_prelim/accepted_hits.bam alignment/$Organism/$Strain/Fus2_72hrs_rep1/accepted_hits.bam alignment/$Organism/$Strain/Fus2_72hrs_rep2/accepted_hits.bam alignment/$Organism/$Strain/Fus2_72hrs_rep3/accepted_hits.bam alignment/$Organism/$Strain/Fus2_8hrs_prelim/accepted_hits.bam alignment/$Organism/$Strain/Fus2_96hrs_prelim/accepted_hits.bam alignment/$Organism/$Strain/Fus2_CzapekDox/accepted_hits.bam alignment/$Organism/$Strain/Fus2_GlucosePeptone/accepted_hits.bam alignment/$Organism/$Strain/Fus2_PDA/accepted_hits.bam alignment/$Organism/$Strain/Fus2_PDB/accepted_hits.bam; do
ls -lh $File;
done
done
``` -->


```bash
	# for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -v 'HB17' | grep 'Fus2' | grep -e 'Fus2_canu_new' -e 'Fus2_merged' | grep 'cepae' | grep 'Fus2_merged'); do
	# for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep 'proliferatum'); do
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139' | grep 'ncbi'); do
Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
done
printf "\n"
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
mkdir -p alignment/$Organism/$Strain/concatenated
samtools merge -f alignment/$Organism/$Strain/concatenated/concatenated.bam \
alignment/$Organism/$Strain/55_72hrs_rep1/accepted_hits.bam \
alignment/$Organism/$Strain/55_72hrs_rep2/accepted_hits.bam \
alignment/$Organism/$Strain/55_72hrs_rep3/accepted_hits.bam \
alignment/$Organism/$Strain/FO47_72hrs_rep1/accepted_hits.bam \
alignment/$Organism/$Strain/FO47_72hrs_rep2/accepted_hits.bam \
alignment/$Organism/$Strain/FO47_72hrs_rep3/accepted_hits.bam \
alignment/$Organism/$Strain/Fus2_0hrs_prelim/accepted_hits.bam \
alignment/$Organism/$Strain/Fus2_16hrs_prelim/accepted_hits.bam \
alignment/$Organism/$Strain/Fus2_24hrs_prelim_rep1/accepted_hits.bam \
alignment/$Organism/$Strain/Fus2_36hrs_prelim/accepted_hits.bam \
alignment/$Organism/$Strain/Fus2_48hrs_prelim/accepted_hits.bam \
alignment/$Organism/$Strain/Fus2_4hrs_prelim/accepted_hits.bam \
alignment/$Organism/$Strain/Fus2_72hrs_prelim/accepted_hits.bam \
alignment/$Organism/$Strain/Fus2_72hrs_rep1/accepted_hits.bam \
alignment/$Organism/$Strain/Fus2_72hrs_rep2/accepted_hits.bam \
alignment/$Organism/$Strain/Fus2_72hrs_rep3/accepted_hits.bam \
alignment/$Organism/$Strain/Fus2_8hrs_prelim/accepted_hits.bam \
alignment/$Organism/$Strain/Fus2_96hrs_prelim/accepted_hits.bam \
alignment/$Organism/$Strain/Fus2_CzapekDox/accepted_hits.bam \
alignment/$Organism/$Strain/Fus2_GlucosePeptone/accepted_hits.bam \
alignment/$Organism/$Strain/Fus2_PDA/accepted_hits.bam \
alignment/$Organism/$Strain/Fus2_PDB/accepted_hits.bam
OutDir=gene_pred/braker/$Organism/"$Strain"_braker
AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
GeneModelName="$Organism"_"$Strain"_braker_new
rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker_new
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
```

Fasta and gff files were extracted from Braker1 output.

```bash
for File in $(ls gene_pred/braker/F.*/*_braker/*/augustus.gff | grep -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139' | grep -e 'canu_new' -e 'ncbi' | grep -v '125'); do
getAnnoFasta.pl $File
OutDir=$(dirname $File)
echo "##gff-version 3" > $OutDir/augustus_extracted.gff
cat $File | grep -v '#' >> $OutDir/augustus_extracted.gff
done
```

The relationship between gene models and aligned reads was investigated. To do
this aligned reads needed to be sorted and indexed:

Note - IGV was used to view aligned reads against the Fus2 genome on my local
machine.

```bash
	InBam=alignment/F.oxysporum_fsp_cepae/Fus2_edited_v2/concatenated/concatenated.bam
	ViewBam=alignment/F.oxysporum_fsp_cepae/Fus2_edited_v2/concatenated/concatenated_view.bam
	SortBam=alignment/F.oxysporum_fsp_cepae/Fus2_edited_v2/concatenated/concatenated_sorted
	samtools view -b $InBam > $ViewBam
	samtools sort $ViewBam $SortBam
	samtools index $SortBam.bam
```

<!--
Cufflinks was run to compare the predicted genes to assembled transcripts:

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -e 'Fus2_canu_new' -e 'Fus2_merged'); do
		Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
		while [ $Jobs -gt 1 ]; do
		sleep 10
		printf "."
		Jobs=$(qstat | grep 'tophat' | grep -w 'r' | wc -l)
		done
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
		OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
		echo "$Organism - $Strain"
		mkdir -p $OutDir
		samtools merge -f $AcceptedHits \
		alignment/$Organism/$Strain/55_72hrs_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/55_72hrs_rep2/accepted_hits.bam \
		alignment/$Organism/$Strain/55_72hrs_rep3/accepted_hits.bam \
		alignment/$Organism/$Strain/FO47_72hrs_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/FO47_72hrs_rep2/accepted_hits.bam \
		alignment/$Organism/$Strain/FO47_72hrs_rep3/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_0hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_16hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_24hrs_prelim_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_36hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_48hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_4hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_rep2/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_rep3/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_8hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_96hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_CzapekDox/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_GlucosePeptone/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_PDA/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_PDB/accepted_hits.bam
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
		qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir/cuflfinks
	# cufflinks -o $OutDir/cufflinks -p 8 --max-intron-length 4000 $AcceptedHits 2>&1 | tee $OutDir/cufflinks/cufflinks.log
	done
```
 -->

<!--
The number of Fo47 genes was determined for comparison to number predicted by Braker (16269):
```bash
	fo47_transcripts=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_transcripts.gtf
	cat $fo47_transcripts | grep 'gene_id' | cut -f2 -d '"' | sort | uniq | wc -l
	# 18191
```
 -->



## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Fistly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi'| grep -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139' | grep -e 'ncbi'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
mkdir -p $OutDir
AcceptedHits=$(ls alignment/$Organism/$Strain/concatenated/concatenated.bam)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
done
```

Secondly, genes were predicted using CodingQuary:

```bash

for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi' | grep -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139' | grep -e 'ncbi' | grep -e 'A8' -e '125' -e 'CB3'); do
Jobs=$(qstat | grep 'sub_cuffli' | grep 'qw' | wc -l)
while [ $Jobs -ge 1 ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'sub_cuffli' | grep 'qw' | wc -l)
done
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/codingquary/$Organism/$Strain
CufflinksGTF=$(ls gene_pred/cufflinks/$Organism/$Strain/concatenated/transcripts.gtf)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
done
```

Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

```bash
	# for BrakerGff in $(ls gene_pred/braker/F.*/*_braker_new/*/augustus.gff3 | grep -w -e 'Fus2'); do
for BrakerGff in $(ls gene_pred/braker/F.*/*/*/augustus.gff3 | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi'| grep -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139' | grep -e 'ncbi' | grep -v '125'); do
Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker_pacbio//g'| sed 's/_braker//g')
Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
CodingQuaryGff=gene_pred/codingquary/$Organism/$Strain/out/PredictedPass.gff3
PGNGff=gene_pred/codingquary/$Organism/$Strain/out/PGN_predictedPass.gff3
AddDir=gene_pred/codingquary/$Organism/$Strain/additional
FinalDir=gene_pred/final_genes/$Organism/$Strain/final
AddGenesList=$AddDir/additional_genes.txt
AddGenesGff=$AddDir/additional_genes.gff
FinalGff=$AddDir/combined_genes.gff
mkdir -p $AddDir
mkdir -p $FinalDir

bedtools intersect -v -a $CodingQuaryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
$ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary

$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuary.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary
cp $BrakerGff $FinalDir/final_genes_Braker.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker
cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta

GffBraker=$FinalDir/final_genes_CodingQuary.gff3
GffQuary=$FinalDir/final_genes_Braker.gff3
GffAppended=$FinalDir/final_genes_appended.gff3
cat $GffBraker $GffQuary > $GffAppended

done
```

The final number of genes per isolate was observed using:
```bash
for DirPath in $(ls -d gene_pred/final_genes/F.*/*/final | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi'| grep -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139'  | grep -e 'ncbi' | grep -v '125'); do
echo $DirPath;
cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
echo "";
done
```
```
gene_pred/final_genes/F.oxysporum_fsp_cepae/A13_ncbi/final
18016
922
18938

gene_pred/final_genes/F.oxysporum_fsp_cepae/A23_ncbi/final
17311
1098
18409

gene_pred/final_genes/F.oxysporum_fsp_cepae/A28_ncbi/final
17404
1469
18873

gene_pred/final_genes/F.oxysporum_fsp_cepae/CB3_ncbi/final
16857
1318
18175

gene_pred/final_genes/F.oxysporum_fsp_cepae/PG_ncbi/final
16864
1223
18087

gene_pred/final_genes/F.oxysporum_fsp_narcissi/N139_ncbi/final
19075
1772
20847

gene_pred/final_genes/F.proliferatum/A8_ncbi/final
15553
42
15595
```

## Manual editing of gene models

Following annotation of gene models and orthology analysis documented below,
some gene models were noted to require manual editing. This was performed and
gene models were re-annotated. In some cases this required removal of genes
annotated by Braker so that the gene models predicted by coding queary could be
favoured. In other cases, gene models were manually edited based upon results
from other Fusarium spp. genomes:

### Fus2:

Gene models requiring editing
 * FTF1 g16934 - delete gene, replace with CodingQuary model.
 * SIX3 g16847, g16190 - A23 one copy
 * SIX5 g16849 - A23 missing intron
 * SIX7 CUFF_4344_1_33 - Gene model good
 * SIX9 g13574, g16273 - premature stop in A23 (end of contig) & one copy in A23 & 125
 * SIX10 g17153 - Gene model good
 * SIX12 PGN_CUFF_4345_1_2 - checked start of Fus2, FoC sequence is too long but FoL reference matches this length - Bring A23 and 125 in line with Fus2 & FoL
 * No gene model

```bash
	Organism=F.oxysporum_fsp_cepae
	Strain=Fus2_canu_new
	OutDir=gene_pred/final_genes/$Organism/$Strain/edited
	mkdir -p $OutDir

	Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
	SourceDir=gene_pred/codingquary/$Organism/$Strain
	BrakerGff=$(ls $SourceDir/final/final_genes_Braker.gff3)
	CodingQuaryGff=$(ls $SourceDir/out/PredictedPass.gff3)
	PGNGff=$(ls $SourceDir/out/PGN_predictedPass.gff3)

	cat $BrakerGff \
	| grep -v -w 'g16934' \
	| grep -v -w 'g13345' \
	> $OutDir/final_genes_Braker_ed.gff3
	printf \
	"contig_19_pilon\tManual_annotation\tgene\t445791\t445964\t0\t+\t.\tID=man_1;
	contig_19_pilon\tManual_annotation\tmRNA\t445791\t445964\t.\t+\t.\tID=man_1.t1;Parent=man_1
	contig_19_pilon\tManual_annotation\tstart_codon\t445791\t445793\t.\t+\t.\tParent=man_1.t1;
	contig_19_pilon\tManual_annotation\tCDS\t445791\t445964\t.\t+\t0\tID=man_1.t1.CDS1;Parent=man_1.t1
	contig_19_pilon\tManual_annotation\texon\t445791\t445964\t.\t+\t.\tID=man_1.t1.exon1;Parent=man_1.t1;
	contig_19_pilon\tManual_annotation\tstop_codon\t445962\t445964\t.\t+\t.\tParent=man_1.t1;\n" \
	>> $OutDir/final_genes_Braker_ed.gff3

	cat $CodingQuaryGff \
	| grep -v -w -e 'NS.00001' -e 'NS.00002' -e 'NS.00003' -e 'NS.04463' -e 'CUFF.11065.1.82' \
	> $OutDir/PredictedPass_ed.gff3
	cp $PGNGff $OutDir/PGN_predictedPass_ed.gff3

	touch $OutDir/manual_annotations.gff3
	#
	# printf \
	# "contig_10_pilon\tManual_annotation\tgene\t1569\t5055\t0\t+\t.\tID=g13345_ps;
	# contig_10_pilon\tManual_annotation\tmRNA\t1569\t5055\t.\t+\t.\tID=g13345_ps.t1;Parent=g13345_ps
	# contig_10_pilon\tManual_annotation\tstart_codon\t1569\t1571\t.\t+\t.\tParent=g13345_ps.t1;
	# contig_10_pilon\tManual_annotation\tCDS\t1569\t2856\t.\t+\t0\tID=g13345_ps.t1.CDS1;Parent=g13345_ps.t1
	# contig_10_pilon\tManual_annotation\texon\t1569\t2856\t.\t+\t.\tParent=g13345_ps.t1
	# contig_10_pilon\tManual_annotation\tintron\t2857\t2905\t.\t+\t.\tParent=g13345_ps.t1
	# contig_10_pilon\tManual_annotation\tCDS\t2906\t5055\t.\t+\t2\tID=g13345_ps.t1.CDS2;Parent=g13345_ps.t1
	# contig_10_pilon\tManual_annotation\texon\t2906\t5055\t.\t+\t.\tID=g13345_ps.t1.exon1;Parent=g13345_ps.t1;
	# contig_10_pilon\tManual_annotation\tstop_codon\t5053\t5055\t.\t+\t.\tParent=g13345_ps.t1;\n" \
	# > $OutDir/pseudogenes.gff3
```

```bash
	Organism=F.oxysporum_fsp_cepae
	Strain=A23
	OutDir=gene_pred/final_genes/$Organism/$Strain/edited
	mkdir -p $OutDir

	Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
	SourceDir=gene_pred/codingquary/$Organism/$Strain
	BrakerGff=$(ls $SourceDir/final/final_genes_Braker.gff3)
	CodingQuaryGff=$(ls $SourceDir/out/PredictedPass.gff3)
	PGNGff=$(ls $SourceDir/out/PGN_predictedPass.gff3)

	cat $BrakerGff \
	| grep -v -w 'g16381' \
	| sed \
	's/contig_818\tAUGUSTUS\tgene\t655\t1273\t0.76\t+\t.\tID=g16444;/contig_818\tAUGUSTUS\tgene\t655\t1188\t0.76\t+\t.\tID=g16444;/g' \
	| sed \
	's/contig_818\tAUGUSTUS\tmRNA\t655\t1273\t0.76\t+\t.\tID=g16444.t1;Parent=g16444/contig_818\tAUGUSTUS\tmRNA\t655\t1188\t0.76\t+\t.\tID=g16444.t1;Parent=g16444/g' \
	| sed \
	's/contig_818\tAUGUSTUS\tintron\t1002\t1057\t0.99\t+\t.\tParent=g16444.t1;/contig_818\tAUGUSTUS\tintron\t1002\t1050\t0.99\t+\t.\tParent=g16444.t1;/g' \
	| sed \
	's/contig_818\tAUGUSTUS\tCDS\t1058\t1273\t0.77\t+\t0\tID=g16444.t1.CDS4;Parent=g16444.t1/contig_818\tAUGUSTUS\tCDS\t1051\t1188\t0.77\t+\t0\tID=g16444.t1.CDS4;Parent=g16444.t1/g' \
	| sed \
	's/contig_818\tAUGUSTUS\texon\t1058\t1273\t.\t+\t.\tID=g16444.t1.exon4;Parent=g16444.t1;/contig_818\tAUGUSTUS\texon\t1051\t1188\t.\t+\t.\tID=g16444.t1.exon4;Parent=g16444.t1;/g' \
	| sed \
	's/contig_818\tAUGUSTUS\tstop_codon\t1271\t1273\t.\t+\t0\tParent=g16444.t1;/contig_818\tAUGUSTUS\tstop_codon\t1186\t1188\t.\t+\t0\tParent=g16444.t1;/g' \
	> $OutDir/final_genes_Braker_ed.gff3

	cp $CodingQuaryGff $OutDir/PredictedPass_ed.gff3
	cat $PGNGff | grep -v -w 'PGN.CUFF.4205' > $OutDir/PGN_predictedPass_ed.gff3

	printf \
	"contig_1793\tManual_annotation\tgene\t2\t219\t1\t+\t.\tID=man_1;
	contig_1793\tManual_annotation\tmRNA\t2\t219\t1\t+\t.\tID=man_1.t1;Parent=man_1
	contig_1793\tManual_annotation\tCDS\t2\t219\t1.0\t+\t0\tID=man_1.t1.CDS1;Parent=man_1.t1
	contig_1793\tManual_annotation\texon\t2\t219\t.\t+\t.\tID=man_1.t1.exon1;Parent=man_1.t1;
	contig_1793\tManual_annotation\tstop_codon\t217\t219\t.\t+\t.\tParent=man_1.t1;\n" \
	>> $OutDir/final_genes_Braker_ed.gff3
```

```bash
	Organism=F.oxysporum_fsp_cepae
	Strain=125
	OutDir=gene_pred/final_genes/$Organism/$Strain/edited
	mkdir -p $OutDir

	Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
	SourceDir=gene_pred/codingquary/$Organism/$Strain
	BrakerGff=$(ls $SourceDir/final/final_genes_Braker.gff3)
	CodingQuaryGff=$(ls $SourceDir/out/PredictedPass.gff3)
	PGNGff=$(ls $SourceDir/out/PGN_predictedPass.gff3)

	cat $BrakerGff \
	| grep -v -w 'g16131' \
	> $OutDir/final_genes_Braker_ed.gff3

	cat $CodingQuaryGff | grep -v -w 'CUFF_14685' >  $OutDir/PredictedPass_ed.gff3
	cp $PGNGff $OutDir/PGN_predictedPass_ed.gff3

	printf \
	"contig_2087\tManual_annotation\tgene\t261\t524\t1\t+\t.\tID=man_1;
	contig_2087\tManual_annotation\tmRNA\t261\t524\t1\t+\t.\tID=man_1.t1;Parent=man_1
	contig_2087\tManual_annotation\tstart_codon\t261\t263\t.\t+\t.\tParent=man_1.t1;
	contig_2087\tManual_annotation\tCDS\t261\t524\t1.0\t+\t0\tID=man_1.t1.CDS1;Parent=man_1.t1
	contig_2087\tManual_annotation\texon\t261\t524\t.\t+\t.\tID=man_1.t1.exon1;Parent=man_1.t1;\n" \
	>> $OutDir/final_genes_Braker_ed.gff3
	printf \
	"contig_972\tManual_annotation\tgene\t333\t506\t1\t-\t.\tID=man_2;
	contig_972\tManual_annotation\tmRNA\t333\t506\t1\t-\t.\tID=man_2.t1;Parent=man_2
	contig_972\tManual_annotation\tstop_codon\t333\t335\t.\t-\t.\tParent=man_2.t1;
	contig_972\tManual_annotation\tCDS\t333\t506\t1.0\t-\t0\tID=man_2.t1.CDS1;Parent=man_2.t1
	contig_972\tManual_annotation\texon\t333\t506\t.\t-\t.\tID=man_2.t1.exon1;Parent=man_2.t1;
	contig_972\tManual_annotation\tstart_codon\t504\t506\t.\t-\t.\tParent=man_2.t1;\n" \
	>> $OutDir/final_genes_Braker_ed.gff3
```

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi'| grep -v -e 'Fus2_canu_new' -e 'A23' -e '125'); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=gene_pred/final_genes/$Organism/$Strain/edited
mkdir -p $OutDir

SourceDir=gene_pred/codingquary/$Organism/$Strain
BrakerGff=$(ls $SourceDir/final/final_genes_Braker.gff3)
CodingQuaryGff=$(ls $SourceDir/out/PredictedPass.gff3)
PGNGff=$(ls $SourceDir/out/PGN_predictedPass.gff3)

cp $BrakerGff $OutDir/final_genes_Braker_ed.gff3
cp $CodingQuaryGff $OutDir/PredictedPass_ed.gff3
cp $PGNGff $OutDir/PGN_predictedPass_ed.gff3
done
```

```bash
for EditDir in $(ls -d gene_pred/final_genes/*/*/edited | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi'| grep -e '125' -e 'A23' -e 'Fus2_canu_new' ); do
Strain=$(echo $EditDir | rev | cut -d '/' -f2 | rev)
Organism=$(echo $EditDir | rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa)
BrakerGff=$EditDir/final_genes_Braker_ed.gff3
CodingQuaryGff=$EditDir/PredictedPass_ed.gff3
PGNGff=$EditDir/PGN_predictedPass_ed.gff3
# ManGff=$EditDir/manual_annotations.gff3
AddDir=$EditDir/additional
FinalDir=gene_pred/final_genes/$Organism/$Strain/final
AddGenesList=$AddDir/additional_genes.txt
AddGenesGff=$AddDir/additional_genes.gff
# FinalGff=$AddDir/combined_genes.gff
mkdir -p $AddDir
mkdir -p $FinalDir

bedtools intersect -v -a $CodingQuaryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
$ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary

$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuary.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary
cp $BrakerGff $FinalDir/final_genes_Braker.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker
# cp $ManGff $FinalDir/final_genes_manual.gff3
# $ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_manual.gff3 $FinalDir/final_genes_manual
cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta

GffBraker=$FinalDir/final_genes_CodingQuary.gff3
GffQuary=$FinalDir/final_genes_Braker.gff3
# GffManual=$FinalDir/final_genes_manual.gff3
GffAppended=$FinalDir/final_genes_appended.gff3
cat $GffBraker $GffQuary > $GffAppended
# cat $GffBraker $GffQuary $GffManual > $GffAppended

done
```

## ORF finder

The genome was searched in six reading frames for any start codon and following
translated identification of a start codon translating sequence until a stop
codon was found. This is based upon the atg.pl script used in paper describing
the P. infestans genome. Additional functionality was added to this script by
also printing ORFs in .gff format.


```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
	for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa | grep -w -e 'Fus2'); do
		qsub $ProgDir/run_ORF_finder.sh $Genome
	done
```

The Gff files from the the ORF finder are not in true Gff3 format. These were
corrected using the following commands:

```bash
	ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
	for ORF_Gff in $(ls gene_pred/ORF_finder/*/*/*_ORF.gff | grep -v '_F_atg_' | grep -v '_R_atg_'); do
		ORF_Gff_mod=$(echo $ORF_Gff | sed 's/_ORF.gff/_ORF_corrected.gff3/g')
		echo ""
		echo "Correcting the following file:"
		echo $ORF_Gff
		echo "Redirecting to:"
		echo $ORF_Gff_mod
		$ProgDir/gff_corrector.pl $ORF_Gff > $ORF_Gff_mod
	done
```

#Functional annotation

## A) Interproscan

Interproscan was used to give gene models functional annotations.
Annotation was run using the commands below:

Note: This is a long-running script. As such, these commands were run using
'screen' to allow jobs to be submitted and monitored in the background.
This allows the session to be disconnected and reconnected over time.

Screen ouput detailing the progress of submission of interporscan jobs
was redirected to a temporary output file named interproscan_submission.log .

```bash
	screen -a
	cd /home/groups/harrisonlab/project_files/fusarium
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Genes in $(ls gene_pred/final_genes/F.*/*/*/final_genes_combined.pep.fasta | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi'| grep -e 'Fus2_canu_new' -e '125' -e 'A23'); do
	echo $Genes
	$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Proteins in $(ls gene_pred/final_genes/F.*/*/*/final_genes_combined.pep.fasta | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi' | grep -e 'Fus2_canu_new' -e '125' -e 'A23'); do
		Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		echo $Strain
		InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
		$ProgDir/append_interpro.sh $Proteins $InterProRaw
	done
```


## B) SwissProt
<!--
```bash
  screen -a
  qlogin
  ProjDir=/home/groups/harrisonlab/project_files/idris
  cd $ProjDir
  for Proteome in $(ls gene_pred/final_genes/F.*/*/*/final_genes_combined.pep.fasta); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=$ProjDir/gene_pred/swissprot/$Species/$Strain/
    mkdir -p $OutDir
    blastp \
    -db /home/groups/harrisonlab/uniprot/swissprot/uniprot_sprot \
    -query $ProjDir/$Proteome \
    -out $OutDir/swissprot_v2015_10_hits.tbl \
    -evalue 1e-100 \
    -outfmt 6 \
    -num_threads 16 \
    -num_alignments 10
  done
``` -->


```bash
for Proteome in $(ls gene_pred/final_genes/F.*/*/*/final_genes_combined.pep.fasta | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi'| grep -e 'Fus2_canu_new' -e '125' -e 'A23'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
OutDir=gene_pred/swissprot/$Organism/$Strain
SwissDbDir=../../uniprot/swissprot
SwissDbName=uniprot_sprot
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
done
```

```bash
	for SwissTable in $(ls gene_pred/swissprot/*/*/swissprot_v2015_10_hits.tbl | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi'| grep -e 'Fus2_canu_new' -e '125' -e 'A23'); do
		# SwissTable=gene_pred/swissprot/Fus2/swissprot_v2015_10_hits.tbl
		Strain=$(echo $SwissTable | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $SwissTable | rev | cut -f3 -d '/' | rev)
		echo "$Organism - $Strain"
		OutTable=gene_pred/swissprot/$Organism/$Strain/swissprot_vJul2016_tophit_parsed.tbl
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
		$ProgDir/swissprot_parser.py --blast_tbl $SwissTable --blast_db_fasta ../../uniprot/swissprot/uniprot_sprot.fasta > $OutTable
	done
```


## C) CAZY proteins

Carbohydrte active enzymes were idnetified using CAZYfollowing recomendations
at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt :

```bash
for Proteome in $(ls gene_pred/final_genes/F.*/*/*/final_genes_combined.pep.fasta | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
OutDir=gene_pred/CAZY/$Organism/$Strain
mkdir -p $OutDir
Prefix="$Strain"_CAZY
CazyHmm=../../dbCAN/dbCAN-fam-HMMs.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/HMMER
qsub $ProgDir/sub_hmmscan.sh $CazyHmm $Proteome $Prefix $OutDir
done
```

The Hmm parser was used to filter hits by an E-value of E1x10-5 or E 1x10-e3 if they had a hit over a length of X %.

Those proteins with a signal peptide were extracted from the list and gff files
representing these proteins made.

```bash
for File in $(ls gene_pred/CAZY/*/*/*CAZY.out.dm | grep -e 'cepae' -e 'proliferatum' -e 'narcissi' | grep 'Fus2'); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $File)
echo "$Organism - $Strain"
ProgDir=/home/groups/harrisonlab/dbCAN
$ProgDir/hmmscan-parser.sh $OutDir/"$Strain"_CAZY.out.dm > $OutDir/"$Strain"_CAZY.out.dm.ps
CazyHeaders=$(echo $File | sed 's/.out.dm/_headers.txt/g')
cat $OutDir/"$Strain"_CAZY.out.dm.ps | cut -f3 | sort | uniq > $CazyHeaders
echo "number of CAZY genes identified:"
cat $CazyHeaders | wc -l
# Gff=$(ls gene_pred/codingquary/$Organism/$Strain/final/final_genes_appended.gff3)
Gff=$(ls gene_pred/final_genes/$Organism/$Strain/final/final_genes_appended.gff3)
CazyGff=$OutDir/"$Strain"_CAZY.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $CazyHeaders $Gff CAZyme ID > $CazyGff

SecretedProts=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/"$Strain"_final_sp_no_trans_mem.aa)
SecretedHeaders=$(echo $SecretedProts | sed 's/.aa/_headers.txt/g')
cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
CazyGffSecreted=$OutDir/"$Strain"_CAZY_secreted.gff
$ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $CazyGff Secreted_CAZyme ID > $CazyGffSecreted
echo "number of Secreted CAZY genes identified:"
cat $CazyGffSecreted | grep -w 'mRNA' | cut -f9 | tr -d 'ID=' | cut -f1 -d ';' > $OutDir/"$Strain"_CAZY_secreted_headers.txt
cat $OutDir/"$Strain"_CAZY_secreted_headers.txt | wc -l
done
```

Note - the CAZY genes identified may need further filtering based on e value and
cuttoff length - see below:

Cols in yourfile.out.dm.ps:
1. Family HMM
2. HMM length
3. Query ID
4. Query length
5. E-value (how similar to the family HMM)
6. HMM start
7. HMM end
8. Query start
9. Query end
10. Coverage

* For fungi, use E-value < 1e-17 and coverage > 0.45

* The best threshold varies for different CAZyme classes (please see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4132414/ for details). Basically to annotate GH proteins, one should use a very relax coverage cutoff or the sensitivity will be low (Supplementary Tables S4 and S9); (ii) to annotate CE families a very stringent E-value cutoff and coverage cutoff should be used; otherwise the precision will be very low due to a very high false positive rate (Supplementary Tables S5 and S10)

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/CAZY
Secreted=gene_pred/final_genes_signalp-4.1/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_final_sp_no_trans_mem_headers.txt
CAZY=gene_pred/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY.out.dm.ps
Gff=gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3
$ProgDir/summarise_CAZY.py --cazy $CAZY --inp_secreted $Secreted --inp_gff $Gff | less -S
```



## D) AntiSMASH

Antismash was run to identify clusters of secondary metabolite genes within
the genome. Antismash was run using the weserver at:
http://antismash.secondarymetabolites.org


<!--
The assembly and Gff annotaiton of gene models was converted into EMBL format prior to submission:

 ```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi'| grep -w 'Fus2_canu_new'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/antismash/$Organism/$Strain
mkdir -p $OutDir
Gff=$(ls gene_pred/final_genes/$Organism/$Strain/final/final_genes_appended.gff3)
seqret -sequence $Assembly -feature -fformat gff -fopenfile $Gff -osformat embl -auto
mv contig_1_pilon.embl $OutDir/"$Strain"_parsed_genome.embl
done
``` -->
<!--
```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi'| grep -w 'Fus2_canu_new'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/antismash/$Organism/$Strain
mkdir -p $OutDir
Gff=$(ls gene_pred/final_genes/$Organism/$Strain/final/final_genes_appended.gff3)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
$ProgDir/antismash_gff2tab.py --gff $Gff --fasta $Assembly --out $OutDir
cp $ProgDir/format_embl.py $OutDir/.
# cp $Assembly $OutDir/assembly.fasta
CurDir=$PWD
cd $OutDir
python format_embl.py
cd $CurDir
rm $OutDir/annotationtable.txt
seqret -sequence $Assembly -feature -fformat gff -fopenfile $Gff -osformat embl -auto
mv contig_1_pilon.embl $OutDir/"$Strain"_parsed_genome.embl
done
```
 -->

Results of web-annotation of gene clusters within the assembly were downloaded to
the following firectories:

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi' | grep -e 'canu_new' -e 'ncbi'); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=analysis/antismash/$Organism/$Strain
mkdir -p $OutDir
done
```

```bash
	for Zip in $(ls analysis/antismash/*/*/*.zip); do
		OutDir=$(dirname $Zip)
		unzip -d $OutDir $Zip
	done
```

```bash
# AntiSmash=analysis/antismash/79c1471f-4a2b-41f7-ba36-18ba94675f59/contig_1_pilon.final.gbk
for AntiSmash in $(ls analysis/antismash/*/*/*/*.final.gbk); do
Organism=$(echo $AntiSmash | rev | cut -f4 -d '/' | rev)
Strain=$(echo $AntiSmash | rev | cut -f3 -d '/' | rev)
OutDir=analysis/antismash/$Organism/$Strain
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
$ProgDir/antismash2gff.py --inp_antismash $AntiSmash > $OutDir/Fus2_secondary_metabolite_regions.gff
printf "Number of clusters detected:\t"
cat $OutDir/Fus2_secondary_metabolite_regions.gff | grep 'antismash_cluster' | wc -l
# GeneGff=gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3
GeneGff=gene_pred/final_genes/$Organism/$Strain/final/final_genes_appended.gff3
bedtools intersect -u -a $GeneGff -b $OutDir/Fus2_secondary_metabolite_regions.gff > $OutDir/metabolite_cluster_genes.gff
cat $OutDir/metabolite_cluster_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' > $OutDir/metabolite_cluster_gene_headers.txt
printf "Number of predicted genes in clusters:\t"
cat $OutDir/metabolite_cluster_gene_headers.txt | wc -l
done
```

#Genomic analysis

## Comparison to FoL 4287

BLast searches were used to identify which genes had homologs on which
chromosomes of the Fusarium lycopersici genome.

```bash
FoLGenomeFa=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum_chromosome_and_additional_contigs.fa
for Proteome in $(ls gene_pred/final_genes/F.*/*/*/final_genes_combined.pep.fasta | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi'| grep -e 'Fus2_canu_new' -e '125' -e 'A23'); do
# for Proteome in $(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_proteins.fasta); do
# for Proteome in $(ls gene_pred/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522_parsed.fa); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=analysis/blast_homology/$Organism/$Strain
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Proteome protein $FoLGenomeFa $OutDir
done
```

Convert top blast hits into gff annotations

```bash
for BlastHitsCsv in $(ls analysis/blast_homology/F.*/*/4287_chromosomal_final_genes_combined.pep.fasta_hits.csv | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi' | grep -e 'Fus2_canu_new' -e '125' -e 'A23'); do
Organism=$(echo $BlastHitsCsv | rev | cut -f3 -d '/' | rev)
Strain=$(echo $BlastHitsCsv | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
HitsGff=$(echo $BlastHitsCsv | sed  's/.csv/.gff/g')
Column2="$Strain"_gene_homolog
NumHits=1
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHitsCsv > $HitsGff
done
```

#### Intersecting blast hits with genes from FoL

```bash
for HitsGff in $(ls analysis/blast_homology/F.*/*/4287_chromosomal_final_genes_combined.pep.fasta_hits.gff | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi' | grep -e 'Fus2_canu_new' -e '125' -e 'A23'); do
Organism=$(echo $HitsGff | rev | cut -f3 -d '/' | rev)
Strain=$(echo $HitsGff| rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
HitsDir=$(dirname $HitsGff)
FoLGenes=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.chr.gff3
FoLIntersect=$HitsDir/4287_chromosomal_final_genes_combined_intersect.bed
bedtools intersect -wo -a $HitsGff -b $FoLGenes > $FoLIntersect
done
```

## Effector genes

Putative pathogenicity and effector related genes were identified within Braker
gene models using a number of approaches:

 * A) From Augustus gene models - Identifying secreted proteins
 * B) From Augustus gene models - Effector identification using EffectorP
 * D) From ORF fragments - Signal peptide & RxLR motif  
 * E) From ORF fragments - Hmm evidence of WY domains  
 * F) From ORF fragments - Hmm evidence of RxLR effectors  


### A) From Augustus gene models - Identifying secreted proteins

Required programs:
 * SignalP-4.1
 * TMHMM

Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
CurPath=$PWD
for Proteome in $(ls gene_pred/final_genes/F.*/*/*/final_genes_combined.pep.fasta | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi' | grep -e 'Fus2_canu_new' -e '125' -e 'A23'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
SplitDir=gene_pred/final_genes_split/$Organism/$Strain
mkdir -p $SplitDir
BaseName="$Organism""_$Strain"_final_preds
$SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
for File in $(ls $SplitDir/*_final_preds_*); do
Jobs=$(qstat | grep 'pred_sigP' | wc -l)
while [ $Jobs -gt 20 ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'pred_sigP' | wc -l)
done
printf "\n"
echo $File
qsub $ProgDir/pred_sigP.sh $File signalp-4.1
done
done
```


The batch files of predicted secreted proteins needed to be combined into a
single file for each strain. This was done with the following commands:
```bash
	for SplitDir in $(ls -d gene_pred/final_genes_split/*/* | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi' | grep -e 'Fus2_canu_new' -e '125' -e 'A23'); do
	Strain=$(echo $SplitDir | rev |cut -d '/' -f1 | rev)
	Organism=$(echo $SplitDir | rev |cut -d '/' -f2 | rev)
	InStringAA=''
	InStringNeg=''
	InStringTab=''
	InStringTxt=''
	SigpDir=final_genes_signalp-4.1
	for GRP in $(ls -l $SplitDir/*_final_preds_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do  
	InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.aa";  
	InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp_neg.aa";  
	InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.tab";
	InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.txt";  
	done
	cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.aa
	cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_neg_sp.aa
	tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.tab
	cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.txt
	done
```

Some proteins that are incorporated into the cell membrane require secretion.
Therefore proteins with a transmembrane domain are not likely to represent
cytoplasmic or apoplastic effectors.

Proteins containing a transmembrane domain were identified:

```bash
for Proteome in $(ls gene_pred/final_genes/F.*/*/*/final_genes_combined.pep.fasta | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi' | grep -e 'Fus2_canu_new' -e '125' -e 'A23'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
qsub $ProgDir/submit_TMHMM.sh $Proteome
done
```

Those proteins with transmembrane domains were removed from lists of Signal
peptide containing proteins

```bash
	for File in $(ls gene_pred/trans_mem/*/*/*_TM_genes_neg.txt | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi' | grep -e 'Fus2_canu_new' -e '125' -e 'A23'); do
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
		echo "$Organism - $Strain"
		TmHeaders=$(echo "$File" | sed 's/neg.txt/neg_headers.txt/g')
		cat $File | cut -f1 > $TmHeaders
		SigP=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp.aa)
		OutDir=$(dirname $SigP)
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
		$ProgDir/extract_from_fasta.py --fasta $SigP --headers $TmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem.aa
		cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | wc -l
	done
```


### B) From Augustus gene models - Effector identification using EffectorP

Required programs:
 * EffectorP.py

```bash
for Proteome in $(ls gene_pred/final_genes/F.*/*/*/final_genes_combined.pep.fasta | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi' | grep -e 'Fus2_canu_new' -e '125' -e 'A23'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
BaseName="$Organism"_"$Strain"_EffectorP
OutDir=analysis/effectorP/$Organism/$Strain
ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation/fungal_effectors
qsub $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
done
```

Those genes that were predicted as secreted and tested positive by effectorP
were identified:

```bash
for File in $(ls analysis/effectorP/*/*/*_EffectorP.txt | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi' | grep -e 'Fus2_canu_new' -e '125' -e 'A23'); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Headers=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_headers.txt/g')
cat $File | grep 'Effector' | cut -f1 > $Headers
Secretome=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp_no_trans_mem.aa)
OutFile=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.aa/g')
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Secretome --headers $Headers > $OutFile
OutFileHeaders=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted_headers.txt/g')
cat $OutFile | grep '>' | tr -d '>' > $OutFileHeaders
cat $OutFileHeaders | wc -l
Gff=$(ls gene_pred/final_genes/$Organism/$Strain/*/final_genes_appended.gff3)
EffectorP_Gff=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.gff/g')
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $OutFileHeaders $Gff effectorP ID > $EffectorP_Gff
done
```

### C) Identification of MIMP-flanking genes

```bash
for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi'); do
Organism=$(echo "$Genome" | rev | cut -d '/' -f4 | rev)
Strain=$(echo "$Genome" | rev | cut -d '/' -f3 | rev)
BrakerGff=$(ls gene_pred/final_genes/$Organism/"$Strain"/final/final_genes_CodingQuary.gff3)
QuaryGff=$(ls gene_pred/final_genes/$Organism/"$Strain"/final/final_genes_Braker.gff3)
OutDir=analysis/mimps/$Organism/$Strain
mkdir -p "$OutDir"
echo "$Organism - $Strain"
ProgDir="/home/armita/git_repos/emr_repos/tools/pathogen/mimp_finder"
$ProgDir/mimp_finder.pl $Genome $OutDir/"$Strain"_mimps.fa $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps.log
$ProgDir/gffexpander.pl +- 2000 $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps_exp.gff
echo "The number of mimps identified:"
cat $OutDir/"$Strain"_mimps.fa | grep '>' | wc -l
bedtools intersect -u -a $BrakerGff -b $OutDir/"$Strain"_mimps_exp.gff > $OutDir/"$Strain"_genes_in_2kb_mimp.gff
bedtools intersect -u -a $QuaryGff -b $OutDir/"$Strain"_mimps_exp.gff >> $OutDir/"$Strain"_genes_in_2kb_mimp.gff
echo "The following transcripts intersect mimps:"
MimpGenesTxt=$OutDir/"$Strain"_genes_in_2kb_mimp.txt
cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' | sort | uniq > $MimpGenesTxt
cat $MimpGenesTxt | wc -l
echo ""
done
```

```
	F.oxysporum_fsp_cepae - Fus2
	The number of mimps identified:
	153
	The following transcripts intersect mimps:
	157
```

```bash
	cat gene_pred/final_genes_signalp-4.1/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_final_sp.aa | grep '>' | sed 's/>//g' | cut -f1 -d '-' | sed -r 's/\s//g' > gene_pred/final_genes_signalp-4.1/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_final_sp_headers.txt
	cat analysis/effectorP/F.oxysporum_fsp_cepae/Fus2/F.oxysporum_fsp_cepae_Fus2_EffectorP.txt | grep 'Effector' | cut -f1 | sed 's/|//g' > analysis/effectorP/F.oxysporum_fsp_cepae/Fus2/F.oxysporum_fsp_cepae_Fus2_EffectorP_headers.txt
	cat gene_pred/trans_mem/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_TM_genes_pos.txt | cut -f1 > gene_pred/trans_mem/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_TM_genes_pos_headers.txt
```


Those genes that were predicted as secreted and within 2Kb of a MIMP
were identified:

```bash
	for File in $(ls analysis/mimps/*/*/*_genes_in_2kb_mimp.gff | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi' | grep -w -v -e 'Fus2_edited_v1' -e 'Fus2_edited_v2' -e 'Fus2_merged' -e 'Fus2'); do
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
		echo "$Organism - $Strain"
		Secretome=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp_no_trans_mem.aa)
		OutFile=$(echo "$File" | sed 's/.gff/_secreted.gff/g')
		SecretedHeaders=$(echo "$Secretome" | sed 's/.aa/_headers.txt/g')
		cat $Secretome | grep '>' | tr -d '>' > $SecretedHeaders
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
		$ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $File secreted_mimp ID > $OutFile
		cat $OutFile | grep -w 'mRNA' | wc -l
	done
```


# 4. Genomic analysis
<!--
## 4.1 Identifcation of protospacers

To facilitate CriprCas editing of the Fusarium oxysporum genome target sites
known as protospacers were identified.

This was done using a published program OPTIMuS as well as a parser script. The
commands to do this were as follows:

```bash
	for GeneSeq in $(ls gene_pred/augustus/F.*/*/*_augustus_preds.codingseq | grep -v 'old'); do
	  Organism=$(echo $GeneSeq | rev | cut -f3 -d '/' | rev)
	  Strain=$(echo $GeneSeq | rev | cut -f2 -d '/' | rev)
		echo $Organism
		echo $Strain
	  ProgDir=~/git_repos/emr_repos/scripts/fusarium_venenatum/OPTIMus
	  OutDir=analysis/protospacers/$Organism/$Strain
	  mkdir -p $OutDir
	  $ProgDir/OPTIMuS_EMR.pl $GeneSeq "threshold" 1 > $OutDir/"$Strain"_protospacer_sites.txt
	  $ProgDir/Optimus2csv.py --inp $OutDir/"$Strain"_protospacer_sites.txt  --out $OutDir/"$Strain"_protospacer_by_gene.csv
	done
```

Unrelated to this project, protospacers were identified in some fasta files for Fiona
```bash
for GeneSeq in $(ls analysis/protospacers/for_fiona/*.fasta | grep -v 'parsed'); do
echo $GeneSeq
OutDir=$(dirname $GeneSeq)
OutName=$(basename $GeneSeq | sed 's/.fasta//g')
# cat $GeneSeq | sed 's/-//g' | tr '///'> $OutDir/"$OutName"_parsed.fasta
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
$ProgDir/fasta_replace_ambiguous_bp.py --inp $GeneSeq > $OutDir/"$OutName"_parsed.fasta
ProgDir=~/git_repos/emr_repos/scripts/fusarium_venenatum/OPTIMus
$ProgDir/OPTIMuS_EMR.pl $OutDir/"$OutName"_parsed.fasta "threshold" 1 > $OutDir/"$OutName"_protospacer_sites.txt
$ProgDir/Optimus2csv.py --inp $OutDir/"$OutName"_protospacer_sites.txt  --out $OutDir/"$OutName"_protospacer_by_gene.csv
done
```
-->

## 4.1 Chracterisation of LS regions

See pathogen/identifying_LS_chromosomes/characterising_LS_contigs.md

## 4.2 Orthology

Orthomcl was used to identify orthologous groups between Fusarium spp. genomes

Genomes were grouped by subspecies and orthology was determined within each
subspecies group. Orthology was also determined between subspecies groups.

| Pathogenic | non-pathogenic | Intermediate |
| ---------- | -------------- | -------------|
| 125        | A28            | 55           |
| A23        | D2             |              |
| Fus2       | PG             |              |



<!--
### 4.2.a) Orthology between pathogenic isolates

The Commands used to run this analysis are shown in
pathogen/orthology/F.oxysporum_fsp.cepae_pathogenic_orthology.md


### 4.2.b) Orthology between non-pathogenic isolates

The Commands used to run this analysis are shown in
pathogen/orthology/F.oxysporum_fsp.cepae_non-pathogenic_orthology.md
 -->

### 4.2.c) Orthology between pathogenic and non-pathogenic isolates

The Commands used to run this analysis are shown in
pathogen/orthology/F.oxysporum_fsp.cepae_pathogen_vs_non-pathogen_orthology.md


<!-- ### 4.2.d) Orthology between all isolates

The Commands used to run this analysis are shown in
pathogen/orthology/F.oxysporum_fsp.cepae_isolates.md -->


## 5.  Searches

## 5.1.A) Identifying SIX genes

Protein sequence of previously characterised SIX genes used to BLAST against
assemblies.

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -w -v 'Fus2_canu_new' | grep 'ncbi'); do
		echo $Assembly
		Query=analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX.fa
		qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
	done
```

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
for BlastHits in $(ls analysis/blast_homology/*/*/*Fo_path_genes_CRX.fa_homologs.csv | grep -v -w 'Fus2_canu_new' | grep 'ncbi'); do
Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
HitsGff=$(echo $BlastHits | sed  's/.csv/.gff/g')
Column2=SIX_homolog
NumHits=12
$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
done
```

	The blast hits were summarised in a single table for all the genomes. The top
	identity of the top blast hit in relation to the enquire query sequence was
	presented for each blast hit.

```bash
OutFile=analysis/blast_homology/Fo_path_genes_CRX_summary.tab
echo "Organism" > tmp2.tab
cat analysis/blast_homology/F.proliferatum/A8/A8_Fo_path_genes_CRX.fa_homologs.csv | cut -f1 >> tmp2.tab
for BlastHits in $(ls analysis/blast_homology/*/*/*Fo_path_genes_CRX.fa_homologs.csv | grep -v -w -e 'Fus2' -e 'Fus2_edited_v1' -e 'Fus2_edited_v2' -e 'Fus2_canu' -e 'Fus2_3' -e '4287' | grep -e 'Fus2_canu_new' -e 'ncbi'); do
Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
echo "$Organism" > tmp.tab
echo "$Strain" >> tmp.tab
cat $BlastHits | cut -f10,19 | tail -n +2 | sed -r "s/\t/,/g" >> tmp.tab
paste tmp2.tab tmp.tab > $OutFile
cp $OutFile tmp2.tab
done
rm tmp.tab
rm tmp2.tab
```

```bash
for HitsGff in $(ls analysis/blast_homology/*/*/*Fo_path_genes_CRX.fa_homologs.gff | grep -v 'trinity' | grep -e 'Fus2_canu_new' -e 'ncbi'); do
Strain=$(echo $HitsGff | rev | cut -f2 -d '/' | rev)
Organism=$(echo $HitsGff | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
GffAppended=$(ls gene_pred/final/$Organism/"$Strain"*/final/final_genes_appended.gff3)
OutDir=$(dirname $HitsGff)
SixIntersect=$OutDir/"$Strain"_Fo_path_genes_CRX.fa_hit_genes.bed
bedtools intersect -wao -a $HitsGff -b $GffAppended > $SixIntersect
bedtools intersect -wao -a $HitsGff -b $GffAppended | cut -f9,18 | grep -v 'Parent'
echo ""
done | tr -d ';' | tr -d '"' | sed 's/ID=//g'> analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX_hit_genes_summary.tab
```

```
F.oxysporum_fsp_cepae - Fus2_canu_new
C5_Seq45_BlastHit_1	g13425
CRX1_Seq49_BlastHit_1	g15739
CRX1_Seq49_BlastHit_2	g11056
CRX1_Seq49_BlastHit_3	g59
CRX1_Seq49_BlastHit_4	g16971
CRX2_Seq50_BlastHit_1	g11056
CRX2_Seq50_BlastHit_2	g15739
CRX2_Seq50_BlastHit_3	g59
MIMP_39_ex.19GJ26_B10-1.seq_BlastHit_1	g15724
MIMP_39_ex.19GJ26_B10-1.seq_BlastHit_2	g17045
MIMP_39_ex.19GJ26_B10-1.seq_BlastHit_3	CUFF_4289_1_31
MIMP_39_ex.19GJ26_B10-1.seq_BlastHit_4	CUFF_327_1_65
MIMP_39_ex.19GJ26_B10-1.seq_BlastHit_5	g12255
MIMP_g29_ex.19GJ26_H12-1.seq_BlastHit_1	g16905
MIMP_g16_ex.19GJ26_B07-1.seq_BlastHit_1	g16179
MIMP_g41_ex.45DG31-1.seq_BlastHit_1	g17144
MIMP_g41_ex.45DG31-1.seq_BlastHit_2	NS_02766
Fusarium_oxysporum_f._sp._lycopersici_isolate_FOL-MM10_secreted_in_xylem_3_(SIX3)_gene,_complete_cds_BlastHit_1	g16847
Fusarium_oxysporum_f._sp._lycopersici_isolate_FOL-MM10_secreted_in_xylem_3_(SIX3)_gene,_complete_cds_BlastHit_2	g16190
Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_5_(SIX5)_gene,_partial_cds_BlastHit_1	g16849
Fusarium_oxysporum_f._sp._lycopersici_secreted_in_xylem_Six7_(SIX7)_mRNA,_complete_cds_BlastHit_1	CUFF_4344_1_33
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six9_(SIX9)_mRNA,_complete_cds_BlastHit_1	g16273
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six9_(SIX9)_mRNA,_complete_cds_BlastHit_2	g13574
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six10_(SIX10)_mRNA,_complete_cds_BlastHit_1	g17153
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six12_(SIX12)_mRNA,_complete_cds_BlastHit_1	PGN_CUFF_4345_1_2
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six14_(SIX14)_mRNA,_complete_cds_BlastHit_1	man_1
```

The orthogroups that these genes belonged to were identified:

```bash
for Gene in $(cat analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX_hit_genes_summary.tab | cut -f2 | tail -n+3 | grep -v -w '.'); do
printf "$Gene - ";
Orthogroup=$(ls analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL_publication/FoC_vs_Fo_vs_FoL_publication_orthogroups.txt);
cat $Orthogroup | grep "Fus2|$Gene.t" | cut -f1 -d ':';
done | cut -f3 -d ' ' | sort | uniq > tmp.txt
```


## 5.1.B) Identifying FTF genes

Previously published FTF genes from Sanchez et al 2016 were blasted against
Fusarium genomes.

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -v -e 'HB17' | grep -w 'Fus2_canu_new'); do
		Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
		Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
		echo $Assembly
		Query=analysis/blast_homology/Fo_path_genes/FTF_cds_Sanchez_et_al_2016.fasta
		OutDir=analysis/FTF/$Organism/$Strain
		mkdir -p $OutDir
		ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
		qsub $ProgDir/run_blast2csv.sh $Query dna $Assembly $OutDir
	done
```

BLAST hits were converted to Gff annotations and intersected with gene models:

```bash
for BlastHits in $(ls analysis/FTF/*/*/*_FTF_cds_Sanchez_et_al_2016.fasta_hits.csv | grep -w 'Fus2_canu_new'); do
Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
OutDir=analysis/FTF/$Organism/$Strain
HitsGff=$(echo $BlastHits | sed  's/.csv/.gff/g')
Column2=FTF_homolog
NumHits=1
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff

GffAppended=$(ls gene_pred/final_genes/$Organism/$Strain/final/final_genes_appended.gff3)
bedtools intersect -wao -a $HitsGff -b $GffAppended > $OutDir/"$Strain"_FTF_hits_intersected.bed
done
```

Fus2 genes g16859 and g10474 were identified as FTF1 and FTF2 homologs.

```bash
cat $OutDir/"$Strain"_FTF_hits_intersected.bed | grep -w 'gene' | cut -f18 | cut -f1 -d ';' | sort | uniq
```

The orthogorups that these genes belonged to was identified:

```bash
	cat analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL_publication/FoC_vs_Fo_vs_FoL_publication_orthogroups.txt | grep -e 'Fus2|g16859.t' -e 'Fus2|g10474.t'
```

```
orthogroup506: A13|g9378.t1 4287|FOXG_09390 A28|g1083.t1 fo47|FOZG_12687T0 PG|g9446.t1 CB3|g6561.t1 Fus2|g10474.t1 A23|g3047.t1 fo47|FOZG_18011T0 Fus2|g13345.t1 PG|g14951.t1 125|g2615.t1 125|g16343.t1 Fus2|g16859.t1 Fus2|NS_02884.t1 A23|g16443.t1 4287|FOXG_14257 4287|FOXG_15059 4287|FOXG_14422 A13|g15293.t1 4287|FOXG_17458 4287|FOXG_14000 Fus2|g13345.t2 4287|FOXG_17084 4287|FOXG_12589 4287|FOXG_12539 4287|FOXG_17123 4287|FOXG_16414
```

This identified orthogroup441 as containing FTF1 & FTF2.

for Strain in Fus2 A13 4287 A28 fo47 PG

<!--
FTF genes were identified by identifying orthogroups containing FoL FTF genes
identified in Sanchez et al 2016.

This led to orthogroups_652 being identified as the FTF  
 -->


## 5.1.C) Identifying Mitochondrial genes in assemblies

 Previously published Fusarium spp. Mitcochondial genes from Al-Reedy 2012 were blasted against
 assembled genomes.

 ```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -v -e 'HB17' | grep -w 'Fus2_canu_new'); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo $Assembly
Query=analysis/blast_homology/Mt_genes/F.spp._mt_prots_Al-Reedy_et_al._2012.fasta
OutDir=analysis/Mt_genes/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query protein $Assembly $OutDir
done
 ```

 BLAST hits were converted to Gff annotations and intersected with gene models:

 ```bash
	 for BlastHits in $(ls analysis/Mt_genes/*/*/*_F.spp._mt_prots_Al-Reedy_et_al._2012.fasta_hits.csv | grep -v -e 'HB17' | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi'); do
		 Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
		 Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
		 echo "$Organism - $Strain"
		 OutDir=analysis/Mt_genes/$Organism/$Strain
		 HitsGff=$(echo $BlastHits | sed  's/.csv/.gff/g')
		 Column2=Mt_gene_homolog
		 NumHits=1
		 ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
		 $ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff

		 GffAppended=$(ls gene_pred/final_genes/$Organism/$Strain/final/final_genes_appended.gff3)
		 bedtools intersect -wao -a $HitsGff -b $GffAppended > $OutDir/"$Strain"_Mt_hits_intersected.bed
	 done
 ```


## 5.2 Identifying PHIbase homologs

The PHIbase database was searched against the assembled genomes using tBLASTx.

```bash
	for Assembly in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
		qsub $ProgDir/blast_pipe.sh analysis/blast_homology/PHIbase/PHI_36_accessions.fa protein $Assembly
	done
```

following blasting PHIbase to the genome, the hits were filtered by effect on
virulence.

First the a tab seperated file was made in the clusters core directory containing
PHIbase. These commands were run as part of previous projects but have been
included here for completeness.
<!--
```bash
	PhibaseDir=/home/groups/harrisonlab/phibase/v3.8
	printf "header\n" > $PhibaseDir/PHI_headers.csv
	cat $PhibaseDir/PHI_accessions.fa | grep '>' | cut -f1 | sed 's/>//g' | sed 's/\r//g' >> $PhibaseDir/PHI_headers.csv
	printf "effect\n" > .$PhibaseDir/PHI_virulence.csv
	cat $PhibaseDir/PHI_accessions.fa | grep '>' | cut -f1 | sed 's/>//g' | rev | cut -f1 -d '|' | rev  >> $PhibaseDir/PHI_virulence.csv
```


```bash
	PhibaseDir=/home/groups/harrisonlab/phibase/v3.8
	PhibaseHeaders=$PhibaseDir/PHI_headers.csv
	PhibaseVirulence=$PhibaseDir/PHI_virulence.csv
	for BlastCSV in $(ls analysis/blast_homology/F*/*/*_PHI_36_accessions.fa_homologs.csv); do
		Strain=$(echo $BlastCSV | rev | cut -f2 -d'/' | rev)
		echo "$Strain"
		OutDir=$(dirname $BlastCSV)
		paste -d '\t' $PhibaseHeaders $PhibaseVirulence $BlastCSV | cut -f-3,1185- > $OutDir/"$Strain"_PHIbase_virulence.csv
		cat $OutDir/"$Strain"_PHIbase_virulence.csv | grep 'NODE_' | cut -f2 | sort | uniq -c | tee $OutDir/"$Strain"_PHIbase_virulence.txt
	done
```
results were:

```
	125
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      1 Increased virulence (Hypervirulence)
	     15 Lethal
	     66 Loss of pathogenicity
	     14 Mixed outcome
	      2 reduced virulence
	    131 Reduced virulence
	     87 Unaffected pathogenicity
	55
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      1 Increased virulence (Hypervirulence)
	     14 Lethal
	     58 Loss of pathogenicity
	     12 Mixed outcome
	      2 reduced virulence
	    121 Reduced virulence
	     83 Unaffected pathogenicity
	A23
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      1 Increased virulence (Hypervirulence)
	     14 Lethal
	     60 Loss of pathogenicity
	     12 Mixed outcome
	      2 reduced virulence
	    123 Reduced virulence
	     79 Unaffected pathogenicity
	A28
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     15 Lethal
	     58 Loss of pathogenicity
	     13 Mixed outcome
	      2 reduced virulence
	    125 Reduced virulence
	     79 Unaffected pathogenicity
	D2
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     16 Lethal
	     62 Loss of pathogenicity
	     12 Mixed outcome
	      2 reduced virulence
	    123 Reduced virulence
	     81 Unaffected pathogenicity
	Fus2
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      1 Increased virulence (Hypervirulence)
	     15 Lethal
	     59 Loss of pathogenicity
	     12 Mixed outcome
	      2 reduced virulence
	    120 Reduced virulence
	     79 Unaffected pathogenicity
	HB17
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      1 Increased virulence (Hypervirulence)
	     14 Lethal
	     60 Loss of pathogenicity
	     13 Mixed outcome
	      2 reduced virulence
	    123 Reduced virulence
	     84 Unaffected pathogenicity
	PG
	      3 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     15 Lethal
	     61 Loss of pathogenicity
	     13 Mixed outcome
	      2 reduced virulence
	    126 Reduced virulence
	     78 Unaffected pathogenicity
	N139
	      4 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     16 Lethal
	     72 Loss of pathogenicity
	     13 Mixed outcome
	      2 reduced virulence
	    149 Reduced virulence
	     93 Unaffected pathogenicity
	PG18
	      4 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     17 Lethal
	     72 Loss of pathogenicity
	     14 Mixed outcome
	      2 reduced virulence
	    153 Reduced virulence
	     97 Unaffected pathogenicity
	PG3
	      4 Chemistry target
	      1 Effector (plant avirulence determinant)
	      2 Increased virulence (Hypervirulence)
	     17 Lethal
	     65 Loss of pathogenicity
	     13 Mixed outcome
	      2 reduced virulence
	    139 Reduced virulence
	     93 Unaffected pathogenicity
	A8
	      4 Chemistry target
	      1 Effector (plant avirulence determinant)
	      4 Increased virulence (Hypervirulence)
	     15 Lethal
	     71 Loss of pathogenicity
	     14 Mixed outcome
	      2 reduced virulence
	    152 Reduced virulence
	     82 Unaffected pathogenicity

```
-->

<!-- The analysis was also performed by blasting the predicted proteins against the
PHIbase database:

The PHIbase database was searched agasinst the assembled genomes using tBLASTx.

```bash
	for Proteins in $(ls gene_pred/braker/F.oxysporum_fsp_cepae/*/*/*.aa); do
		qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh $Proteins protein ../../phibase/v3.8/PHI_accessions.fa
	done
``` -->

### Expressed Genes

As a preliminary analysis of the RNAseq data, highly expressed genes at 72hrs
post infection were identified in Fus2.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
	# samtools merge -f alignment/$Organism/$Strain/concatenated/Fus2_72hpi.bam alignment/$Organism/$Strain/Fus2_72hrs_rep1/accepted_hits.bam alignment/$Organism/$Strain/Fus2_72hrs_rep2/accepted_hits.bam alignment/$Organism/$Strain/Fus2_72hrs_rep3/accepted_hits.bam
	for AcceptedHits in $(ls alignment/*/*/concatenated/concatenated.bam | grep -v -e 'Fus2_edited_v2'); do
		Strain=$(echo $AcceptedHits | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $AcceptedHits | rev | cut -f4 -d '/' | rev)
		OutDir=$(dirname $AcceptedHits)
		echo "$Organism - $Strain"
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
		qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
	done
```

```bash
	for Transcripts in $(ls alignment/F.*/*/concatenated/transcripts.gtf); do
		Strain=$(echo $Transcripts | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Transcripts | rev | cut -f4 -d '/' | rev)
		echo "$Organism - $Strain"
		GeneGff=gene_pred/final_genes/$Organism/$Strain/final/final_genes_appended.gff3
		ExpressedGenes=alignment/$Organism/$Strain/concatenated/expressed_genes.bed
		bedtools intersect -wao -a $Transcripts -b $GeneGff > $ExpressedGenes
	done
```


<!--
```bash
	Organism=F.oxysporum_fsp_cepae
	Strain=Fus2
	samtools merge -f alignment/$Organism/$Strain/concatenated/Fus2_72hpi.bam alignment/$Organism/$Strain/Fus2_72hrs_rep1/accepted_hits.bam alignment/$Organism/$Strain/Fus2_72hrs_rep2/accepted_hits.bam alignment/$Organism/$Strain/Fus2_72hrs_rep3/accepted_hits.bam

	cufflinks -o timecourse/2016_genes/Fus2/72hrs/cufflinks -p 16 --max-intron-length 4000 alignment/$Organism/$Strain/concatenated/Fus2_72hpi.bam

	Transcripts=timecourse/2016_genes/Fus2/72hrs/cufflinks/transcripts.gtf
	GeneGff=gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_edited_v2/final/final_genes_appended.gff3
	ExpressedGenes=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_genes.gff
	bedtools intersect -wao -a $Transcripts -b $GeneGff > $ExpressedGenes
```
-->

<!--
It was noted that not all the top expressed transcripts had gene models. The
transcripts without gene models were identified:

```bash
	ExpressedTranscriptsTxt=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_transcripts.txt
	cat $ExpressedGenes | sort -r -n -t'"' -k6 | cut -f2 -d'"' | uniq | head -n20 > $ExpressedTranscriptsTxt
	ExpressedTranscriptsInter=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_transcripts_intersect.gff
	cat $ExpressedGenes | grep -v -P "\ttranscript_id" | sort -r -n -t'"' -k6 | grep -w -f $ExpressedTranscriptsTxt > $ExpressedTranscriptsInter
	ExpressedTranscriptsGff=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_transcripts.gff
	cat $ExpressedTranscriptsInter | cut -f1-9 | uniq> $ExpressedTranscriptsGff
	cat $ExpressedTranscriptsInter | grep -v -w 'exon' | cut -f2,18 -d'"' --output-delimiter " - " | uniq | less
```

The top 20 expressed transcripts are shown below with genes they correspond to.
9 of the 20 top expressed transcripts have not been predicted as genes. This was
determined, in large, to be a result of transposon being expressed. It was also
noted that cufflinks doesn't always predict direction of a transcript and
therefore when features are intersected, they can not be restricted by strand.

```
	CUFF.4396 - 33 bp - A-rich
	CUFF.2444 - 69 bp - T-rich - repmasked
	CUFF.540 - 60 bp - C-Rich - repmasked
	CUFF.11337 - 92 bp - T-rich - repmasked
	CUFF.14001 - g11792 - Six5
	CUFF.5007 - repmasked
	CUFF.4291 - repmasked
	CUFF.12587 - g11716 - Six3
	CUFF.4090 - repmasked
	CUFF.12214 - g4762
	CUFF.8542 - 103 bp - A-rich - repmasked
	CUFF.13998 - g12630
	CUFF.13635 - g10978
	CUFF.4290 - no gene - Six9
	CUFF.2922 - g12344
	CUFF.9876 - g10974
	CUFF.13435 - g10975
	CUFF.9559 - g10973
	CUFF.12530 - g4720
	CUFF.6539 - g10346
```
-->


<!--
```bash
	TopGenes=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_genes_top20.txt
	cat $ExpressedGenes | sort -r -n -t'"' -k6 | grep -w 'gene' | head -n 20 | cut -f18 > $TopGenes
	InterPro=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2/Fus2_interproscan.tsv
	TopFunctions=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_genes_top20.tsv
	cat $InterPro | grep -w -f $TopGenes > $TopFunctions

	IsolateAbrv=FoC_path_vs_non_path
	WorkDir=analysis/orthology/orthomcl/$IsolateAbrv
	Orthogroups=$WorkDir/"$IsolateAbrv"_orthogroups.txt
	PathDir=$WorkDir/path_orthogroups
	echo "The number of pathogen common genes in the top 20 expressed Fus2 genes at 72hr"
	TopGenes=timecourse/2016_genes/Fus2/72hrs/cufflinks/Fus2_expressed_genes_top20.txt
	cat $TopGenes | grep -w -f $PathDir/Fus2_path_orthogroup_genes.txt
	# secreted
	PathOrthogroupsFus2Secreted=$PathDir/Fus2_path_orthogroup_secreted.txt
	cat $PathOrthogroupsFus2Secreted | grep -w -f $TopGenes
	# near mimps
	MimpGenesTxt=analysis/mimps/F.oxysporum_fsp_cepae/Fus2/Fus2_genes_in_2kb_mimp.txt
  cat $MimpGenesTxt | grep -w -f $TopGenes
	# Intersecting SIX genes

``` -->

### 5.4 Enrichment of genes on LS contigs

Fus2 contigs 10, 14, 16, 17, 19, 20 were identified as lineage specific.

As such interproscan annotation were extracted for these contigs and Functional
enrichments were compared to the entire genome.


```bash

OutDir=analysis/enrichment/LS_regions/F.oxysporum_fsp_cepae/Fus2_canu_new
mkdir -p $OutDir
Gff=gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3
GeneList=$OutDir/Fus2_canu_new_LS_genes_interproscan.txt
cat $Gff | grep -e 'contig_10_pilon' -e 'contig_14_pilon' -e 'contig_16_pilon' -e 'contig_17_pilon' -e 'contig_19_pilon' -e 'contig_20_pilon' | grep 'mRNA' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '='  > $GeneList
InterPro=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv
InterPro_LS=$OutDir/Fus2_canu_new_LS_genes_interproscan.tsv

ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan/feature_extraction
$ProgDir/filter_by_gene.py --inp_interpro $InterPro --inp_txt $GeneList > $InterPro_LS
```

Interproscan tsv files for the entire genome and the LS regions were converted
wego format for analysis by WEGO

```bash
	OutDir=analysis/enrichment/LS_regions/F.oxysporum_fsp_cepae/Fus2_canu_new
	InterPro=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_interproscan.tsv
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan/feature_extraction
	$ProgDir/interpro2wego.py --inp_interpro $InterPro > $OutDir/Fus2_canu_new_GO_WEGO.txt

	InterPro_LS=$OutDir/Fus2_canu_new_LS_genes_interproscan.tsv
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan/feature_extraction
	$ProgDir/interpro2wego.py --inp_interpro $InterPro_LS > $OutDir/Fus2_canu_new_LS_GO_WEGO.txt
```

The WEGO wesite was used to generate summaries of genes on core and LS regions.

## 6. Summarising the Fusarium Proteome

```bash

```
