#HAPI
==========

Scripts used for the analysis of Fusarium genomes
Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/fusarium

The following is a summary of the work presented in this Readme.

The following processes were applied to Alternaria genomes prior to analysis:
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
	for RawData in $(ls qc_dna/paired/*/*/*/*.fq.gz); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

kmer counting was performed using kmc
This allowed estimation of sequencing depth and total genome size

This was performed for strains with single runs of data

```bash
	for TrimPath in $(ls -d raw_dna/paired/*/* | grep -v -e 'Fus2' -e 'HB6'); do
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
	for StrainPath in $(ls -d qc_dna/paired/*/* | grep -v -e 'Fus2' -e 'HB6'); do
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
  for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta); do
  # for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta | grep -e 'Fus2' -e 'HB6' -e 'PG8' -e 'FOP1'); do
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
	for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta | grep -e 'Fus2'); do
  # for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta | grep -e 'Fus2'); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
    $ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file tmp.csv
  done
  rm tmp.csv
```




# Repeatmasking

Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
	for BestAss in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp_renamed.fasta | grep -e 'Fus2'); do
		qsub $ProgDir/rep_modeling.sh $BestAss
		qsub $ProgDir/transposonPSI.sh $BestAss
	done
```

The number of bases masked by transposonPSI and Repeatmasker were summarised
using the following commands:

```bash
	for RepDir in $(ls -d repeat_masked/F.*/*/*); do
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

```
	F.avenaceum	PG8
	The number of bases masked by RepeatMasker:	0
	The number of bases masked by TransposonPSI:	0
	The total number of masked bases are:	0

	F.oxysporum_fsp_cepae	125
	The number of bases masked by RepeatMasker:	3707484
	The number of bases masked by TransposonPSI:	1210934
	The total number of masked bases are:	3952973

	F.oxysporum_fsp_cepae	55
	The number of bases masked by RepeatMasker:	3206161
	The number of bases masked by TransposonPSI:	1031019
	The total number of masked bases are:	3466359

	F.oxysporum_fsp_cepae	A1-2
	The number of bases masked by RepeatMasker:	1937168
	The number of bases masked by TransposonPSI:	816426
	The total number of masked bases are:	2177343

	F.oxysporum_fsp_cepae	A13
	The number of bases masked by RepeatMasker:	4414509
	The number of bases masked by TransposonPSI:	1555142
	The total number of masked bases are:	4712084

	F.oxysporum_fsp_cepae	A23
	The number of bases masked by RepeatMasker:	3156362
	The number of bases masked by TransposonPSI:	1061656
	The total number of masked bases are:	3446078

	F.oxysporum_fsp_cepae	A28
	The number of bases masked by RepeatMasker:	4001386
	The number of bases masked by TransposonPSI:	1189369
	The total number of masked bases are:	4304881

	F.oxysporum_fsp_cepae	CB3
	The number of bases masked by RepeatMasker:	2382071
	The number of bases masked by TransposonPSI:	842157
	The total number of masked bases are:	2630520

	F.oxysporum_fsp_cepae	D2
	The number of bases masked by RepeatMasker:	1363000
	The number of bases masked by TransposonPSI:	594012
	The total number of masked bases are:	1632798

	F.oxysporum_fsp_cepae	Fus2
	The number of bases masked by RepeatMasker:	3619961
	The number of bases masked by TransposonPSI:	1280301
	The total number of masked bases are:	3857605

	F.oxysporum_fsp_cepae	HB17
	The number of bases masked by RepeatMasker:	3385838
	The number of bases masked by TransposonPSI:	1077091
	The total number of masked bases are:	3649652

	F.oxysporum_fsp_cepae	HB6
	The number of bases masked by RepeatMasker:	3216000
	The number of bases masked by TransposonPSI:	995170
	The total number of masked bases are:	3455823

	F.oxysporum_fsp_cepae	PG
	The number of bases masked by RepeatMasker:	2769568
	The number of bases masked by TransposonPSI:	865813
	The total number of masked bases are:	3005423

	F.oxysporum_fsp_narcissi	N139
	The number of bases masked by RepeatMasker:	5404179
	The number of bases masked by TransposonPSI:	1655249
	The total number of masked bases are:	5709023

	F.oxysporum_fsp_pisi	FOP1
	The number of bases masked by RepeatMasker:	0
	The number of bases masked by TransposonPSI:	0
	The total number of masked bases are:	0

	F.oxysporum_fsp_pisi	FOP5
	The number of bases masked by RepeatMasker:	3880611
	The number of bases masked by TransposonPSI:	1313995
	The total number of masked bases are:	4193700

	F.oxysporum_fsp_pisi	L5
	The number of bases masked by RepeatMasker:	1287737
	The number of bases masked by TransposonPSI:	417513
	The total number of masked bases are:	1456488

	F.oxysporum_fsp_pisi	PG18
	The number of bases masked by RepeatMasker:	5349661
	The number of bases masked by TransposonPSI:	1627436
	The total number of masked bases are:	5673770

	F.oxysporum_fsp_pisi	PG3
	The number of bases masked by RepeatMasker:	4686428
	The number of bases masked by TransposonPSI:	1663269
	The total number of masked bases are:	5011872

	F.proliferatum	A8
	The number of bases masked by RepeatMasker:	1065627
	The number of bases masked by TransposonPSI:	278366
	The total number of masked bases are:	1266227
```

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
	for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa); do
		echo $Genome;
		qsub $ProgDir/sub_cegma.sh $Genome dna;
	done
```

Outputs were summarised using the commands:
```bash
	for File in $(ls gene_pred/cegma/F*/*/*_dna_cegma.completeness_report); do
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
	 ls $Folder/F;
	 ls $Folder/R;
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
qsub /home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc/rna_qc_fastq-mcf.sh $FileF $FileR $IlluminaAdapters RNA
sleep 10
done
done
```
```
raw_rna/paired/F.oxysporum_fsp_cepae/55_72hrs_rep1
Your job 6436212 ("rna_qc_fastq-mcf.sh") has been submitted
raw_rna/paired/F.oxysporum_fsp_cepae/55_72hrs_rep2
Your job 6436213 ("rna_qc_fastq-mcf.sh") has been submitted
raw_rna/paired/F.oxysporum_fsp_cepae/55_72hrs_rep3
Your job 6436214 ("rna_qc_fastq-mcf.sh") has been submitted
raw_rna/paired/F.oxysporum_fsp_cepae/control_72hrs_rep2
Your job 6436215 ("rna_qc_fastq-mcf.sh") has been submitted
raw_rna/paired/F.oxysporum_fsp_cepae/control_72hrs_rep3
Your job 6436216 ("rna_qc_fastq-mcf.sh") has been submitted
raw_rna/paired/F.oxysporum_fsp_cepae/FO47_72hrs_rep1
Your job 6436217 ("rna_qc_fastq-mcf.sh") has been submitted
raw_rna/paired/F.oxysporum_fsp_cepae/FO47_72hrs_rep2
Your job 6436218 ("rna_qc_fastq-mcf.sh") has been submitted
raw_rna/paired/F.oxysporum_fsp_cepae/FO47_72hrs_rep3
Your job 6436219 ("rna_qc_fastq-mcf.sh") has been submitted
```


#### Aligning

Then Rnaseq data was aligned to each genome assembly:

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -w -e 'Fus2' -e 'A23' -e 'A28' -e 'D2' -e 'PG' -e '125' | grep -v 'Fus2'); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		for RNADir in $(ls -d qc_rna/paired/F.oxysporum_fsp_cepae/* | grep -v -e '55_72hrs' -e 'FO47_72hrs' -e 'control_72hrs_rep2' -e 'control_72hrs_rep3' | grep 'Fus2_72hrs_rep1'); do
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
```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -w -e 'Fus2' -e 'A23' -e 'A28' -e 'D2' -e 'PG' -e '125'); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		mkdir -p merge alignment/$Organism/$Strain/concatenated
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
		OutDir=gene_pred/braker/$Organism/$Strain
		AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
		GeneModelName="$Organism"_"$Strain"_braker
		rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
		qsub $ProgDir/sub_braker.sh $Assembly $OutDir $AcceptedHits $GeneModelName
	done
```
Fasta and gff files were extracted from Braker1 output.

```bash
for File in $(ls gene_pred/braker/F.*/*/*_braker/augustus.gff ); do
getAnnoFasta.pl $File
OutDir=$(dirname $File)
echo "##gff-version 3" > $OutDir/augustus_extracted.gff
cat $File | grep -v '#' >> $OutDir/augustus_extracted.gff
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
for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa); do
qsub $ProgDir/run_ORF_finder.sh $Genome
done
```

The Gff files produced by ORF finder have some formatting errors. The following
script was run to correct these:

```bash
	ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
	for ORF_Gff in $(ls gene_pred/ORF_finder/*/*/*_ORF.gff | grep -v '_F_atg_' | grep -v '_R_atg_'); do
		ORF_Gff_mod=$(echo $ORF_Gff | sed 's/_ORF.gff/_ORF_mod.gff/g')
		echo ""
		echo "Correcting the following file:"
		echo $ORF_Gff
		echo "Redirecting to:"
		echo $ORF_Gff_mod
		$ProgDir/gff_corrector.pl $ORF_Gff > $ORF_Gff_mod
	done
```


#Functional annotation

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
	for Genes in $(ls gene_pred/braker/F.*/*/*_braker/augustus.aa | grep -w -e 'Fus2' -e 'A23' -e 'A28' -e 'D2' -e 'PG' -e '125' | grep -v 'Fus2'| grep -v '125'); do
		echo $Genes
		$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for StrainPath in $(ls -d gene_pred/interproscan/F.*/*); do
		Strain=$(basename $StrainPath)
		Organism=$(echo $StrainPath | rev | cut -d "/" -f2 | rev)
		echo $Strain
		PredGenes=gene_pred/augustus/"$Organism"/"$Strain"/"$Strain"_augustus_preds.aa
		InterProRaw=gene_pred/interproscan/"$Organism"/"$Strain"/raw
		$ProgDir/append_interpro.sh $PredGenes $InterProRaw
	done
```

#Genomic analysis

## RxLR genes

Putative RxLR genes were identified within Augustus gene models using a number
of approaches:

 * A) From Augustus gene models - Signal peptide & RxLR motif  
 * B) From Augustus gene models - Hmm evidence of WY domains  
 * C) From Augustus gene models - Hmm evidence of RxLR effectors  
 * D) From ORF fragments - Signal peptide & RxLR motif  
 * E) From ORF fragments - Hmm evidence of WY domains  
 * F) From ORF fragments - Hmm evidence of RxLR effectors  


### A) From Augustus gene models - Signal peptide & RxLR motif

Required programs:
 * SigP
 * biopython


Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
	SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	CurPath=$PWD
	for Proteome in $(ls gene_pred/braker/F.*/*/*_braker/augustus.aa | grep -w -e 'Fus2' -e 'A23' -e 'A28' -e 'D2' -e 'PG' -e '125' | grep -v -w 'Fus2'); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		SplitDir=gene_pred/augustus_split/$Organism/$Strain
		mkdir -p $SplitDir
		BaseName="$Organism""_$Strain"_augustus_preds
		$SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
		for File in $(ls $SplitDir/*_augustus_preds_*); do
			Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
			while [ $Jobs -gt 1 ]; do
				sleep 10
				printf "."
				Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
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
	for SplitDir in $(ls -d gene_pred/sigP/F.*/*/split); do
		Strain=$(echo $SplitDir | cut -d '/' -f4)
		Organism=$(echo $SplitDir | cut -d '/' -f3)
		InStringAA=''
		InStringNeg=''
		InStringTab=''
		InStringTxt=''
		for GRP in $(ls -l $SplitDir/*_augustus_preds_*_sp.aa | rev | cut -d '_' -f2 | rev | sort -n); do  
		InStringAA="$InStringAA gene_pred/sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_augustus_preds_$GRP""_sp.aa";  
		InStringNeg="$InStringNeg gene_pred/sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_augustus_preds_$GRP""_sp_neg.aa";  
		InStringTab="$InStringTab gene_pred/sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_augustus_preds_$GRP""_sp.tab";
		InStringTxt="$InStringTxt gene_pred/sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_augustus_preds_$GRP""_sp.txt";  
		done
		cat $InStringAA > gene_pred/sigP/$Organism/$Strain/"$Strain"_aug_sp.aa
		cat $InStringNeg > gene_pred/sigP/$Organism/$Strain/"$Strain"_aug_neg_sp.aa
		tail -n +2 -q $InStringTab > gene_pred/sigP/$Organism/$Strain/"$Strain"_aug_sp.tab
		cat $InStringTxt > gene_pred/sigP/$Organism/$Strain/"$Strain"_aug_sp.txt
	done
```
The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identfy RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
	for Secretome in $(ls gene_pred/sigP/F.*/*/*_aug_sp.aa); do
		ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors;
		Strain=$(echo $Secretome | cut -d '/' -f4);
		Organism=$(echo $Secretome | cut -d '/' -f3) ;
		OutDir=analysis/RxLR_effectors/RxLR_EER_regex_finder/"$Organism"/"$Strain";
		mkdir -p $OutDir;
		printf "\nstrain: $Strain\tspecies: $Organism\n";
		printf "the number of SigP gene is:\t";
		cat $Secretome | grep '>' | wc -l;
		printf "the number of SigP-RxLR genes are:\t";
		$ProgDir/RxLR_EER_regex_finder.py $Secretome > $OutDir/"$Strain"_Aug_RxLR_EER_regex.fa;
		cat $OutDir/"$Strain"_Aug_RxLR_EER_regex.fa | grep '>' | cut -f1 | sed 's/>//g' | sed 's/ //g' > $OutDir/"$Strain"_Aug_RxLR_regex.txt
		cat $OutDir/"$Strain"_Aug_RxLR_regex.txt | wc -l
		printf "the number of SigP-RxLR-EER genes are:\t";
		cat $OutDir/"$Strain"_Aug_RxLR_EER_regex.fa | grep '>' | grep 'EER_motif_start' |  cut -f1 | sed 's/>//g' | sed 's/ //g' > $OutDir/"$Strain"_Aug_RxLR_EER_regex.txt
		cat $OutDir/"$Strain"_Aug_RxLR_EER_regex.txt | wc -l
		printf "\n"
		# ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
		# Col2=RxLR_EER_regex_finder.py
		# GeneNames=$OutDir/"$Strain"_Aug_RxLR_regex.txt
		# GeneModels=gene_pred/augustus/"$Organism"/"$Strain"/"$Strain"_augustus_preds.aa
		# $ProgDir/gene_list_to_gff.pl $GeneNames $GeneModels $Col2 Name > $OutDir/"$Strain"_Aug_RxLR_regex.gff3
	done
```

The results were as follows:

strain: 125	species: F.oxysporum_fsp_cepae  
the number of SigP gene is:	1340  
the number of SigP-RxLR genes are:	13  
the number of SigP-RxLR-EER genes are:	1  

strain: 55	species: F.oxysporum_fsp_cepae  
the number of SigP gene is:	1333  
the number of SigP-RxLR genes are:	12  
the number of SigP-RxLR-EER genes are:	1  

strain: A23	species: F.oxysporum_fsp_cepae  
the number of SigP gene is:	1335  
the number of SigP-RxLR genes are:	13  
the number of SigP-RxLR-EER genes are:	2  

strain: A28	species: F.oxysporum_fsp_cepae  
the number of SigP gene is:	1354  
the number of SigP-RxLR genes are:	16  
the number of SigP-RxLR-EER genes are:	1  

strain: D2	species: F.oxysporum_fsp_cepae  
the number of SigP gene is:	1156  
the number of SigP-RxLR genes are:	8  
the number of SigP-RxLR-EER genes are:	0  

strain: Fus2	species: F.oxysporum_fsp_cepae  
the number of SigP gene is:	1346  
the number of SigP-RxLR genes are:	13  
the number of SigP-RxLR-EER genes are:	1  

strain: HB17	species: F.oxysporum_fsp_cepae  
the number of SigP gene is:	1350  
the number of SigP-RxLR genes are:	13  
the number of SigP-RxLR-EER genes are:	1  

strain: PG	species: F.oxysporum_fsp_cepae  
the number of SigP gene is:	1334  
the number of SigP-RxLR genes are:	11  
the number of SigP-RxLR-EER genes are:	0  

strain: N139	species: F.oxysporum_fsp_narcissi  
the number of SigP gene is:	1765  
the number of SigP-RxLR genes are:	22  
the number of SigP-RxLR-EER genes are:	1  

strain: PG18	species: F.oxysporum_fsp_pisi  
the number of SigP gene is:	1807  
the number of SigP-RxLR genes are:	28  
the number of SigP-RxLR-EER genes are:	1  

strain: PG3	species: F.oxysporum_fsp_pisi  
the number of SigP gene is:	1386  
the number of SigP-RxLR genes are:	11  
the number of SigP-RxLR-EER genes are:	0  

strain: A8	species: F.proliferatum  
the number of SigP gene is:	1975  
the number of SigP-RxLR genes are:	40  
the number of SigP-RxLR-EER genes are:	3  



### B) From Augustus gene models - Hmm evidence of WY domains
Hmm models for the WY domain contained in many RxLRs were used to search gene models predicted with Augustus. These were run with the following commands:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
	for Proteome in $(ls gene_pred/augustus/F.*/*/*_augustus_preds.aa | grep -v '_old'); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_Aug_WY_hmmer.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep 'Initial search space'
		cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
		HmmFasta="$Strain"_Aug_WY_hmmer.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
	done
```

Results were as follows:  

F.oxysporum_fsp_cepae 125
Initial search space (Z):              11203  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae 55  
Initial search space (Z):              11047  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae A23  
Initial search space (Z):              11108  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae A28  
Initial search space (Z):              11131  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae D2  
Initial search space (Z):              10632  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae Fus2  
Initial search space (Z):              11104  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae HB17  
Initial search space (Z):              11126  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae PG  
Initial search space (Z):              11011  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.oxysporum_fsp_narcissi N139  
Initial search space (Z):              14757  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.oxysporum_fsp_pisi PG18  
Initial search space (Z):              14859  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.oxysporum_fsp_pisi PG3  
Initial search space (Z):              11517  [actual number of targets]  
Domain search space  (domZ):               0  [number of targets reported over threshold]  
F.proliferatum A8  
Initial search space (Z):              15034  [actual number of targets]  
Domain search space  (domZ):               1  [number of targets reported over threshold]  



### C) From Augustus gene models - Hmm evidence of RxLR effectors
```bash
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm
	for Proteome in $(ls gene_pred/augustus/F.*/*/*_augustus_preds.aa | grep -v '_old'); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_Aug_RxLR_hmmer.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep 'Initial search space'
		cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
		HmmFasta="$Strain"__Aug_RxLR_hmmer.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
	done
```

F.oxysporum_fsp_cepae 125  
Initial search space (Z):              11203  [actual number of targets]  
Domain search space  (domZ):               3  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae 55  
Initial search space (Z):              11047  [actual number of targets]  
Domain search space  (domZ):               3  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae A23  
Initial search space (Z):              11108  [actual number of targets]  
Domain search space  (domZ):               2  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae A28  
Initial search space (Z):              11131  [actual number of targets]  
Domain search space  (domZ):               1  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae D2  
Initial search space (Z):              10632  [actual number of targets]  
Domain search space  (domZ):               1  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae Fus2  
Initial search space (Z):              11104  [actual number of targets]  
Domain search space  (domZ):               2  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae HB17  
Initial search space (Z):              11126  [actual number of targets]  
Domain search space  (domZ):               2  [number of targets reported over threshold]  
F.oxysporum_fsp_cepae PG  
Initial search space (Z):              11011  [actual number of targets]  
Domain search space  (domZ):               2  [number of targets reported over threshold]  
F.oxysporum_fsp_narcissi N139  
Initial search space (Z):              14757  [actual number of targets]  
Domain search space  (domZ):               1  [number of targets reported over threshold]  
F.oxysporum_fsp_pisi PG18  
Initial search space (Z):              14859  [actual number of targets]  
Domain search space  (domZ):               2  [number of targets reported over threshold]  
F.oxysporum_fsp_pisi PG3  
Initial search space (Z):              11517  [actual number of targets]  
Domain search space  (domZ):               2  [number of targets reported over threshold]  
F.proliferatum A8  
Initial search space (Z):              15034  [actual number of targets]  
Domain search space  (domZ):               2  [number of targets reported over threshold]  

### D) From ORF fragments - Signal peptide & RxLR motif

Required programs:
 * SigP
 * biopython


Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
	SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	CurPath=$PWD
	for Proteome in $(ls gene_pred/ORF_finder/F.*/*/*.aa_cat.fa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		SplitDir=gene_pred/ORF_finder_split/$Organism/$Strain
		mkdir -p $SplitDir
		BaseName="$Organism""_$Strain"_ORF_finder_preds
		$SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
		for File in $(ls $SplitDir/$BaseName_*.fa); do
			Jobs=$(qstat | grep 'pred_sigP' | wc -l)
			while [ $Jobs -ge 32 ]; do
				sleep 10
				printf "."
				Jobs=$(qstat | grep 'pred_sigP' | wc -l)
			done
			printf "\n"
			echo $File
			qsub $ProgDir/pred_sigP.sh $File
		done
	done
```

The batch files of predicted secreted proteins needed to be combined into a
single file for each strain. This was done with the following commands:
```bash
	for SplitDir in $(ls -d gene_pred/ORF_finder_sigP/F.*/*/split); do
	Strain=$(echo $SplitDir | cut -d '/' -f4)
	Organism=$(echo $SplitDir | cut -d '/' -f3)
	InStringAA=''
	InStringNeg=''
	InStringTab=''
	InStringTxt=''
	for GRP in $(ls -l $SplitDir/*_ORF_finder_preds_*_sp.aa | rev | cut -d '_' -f2 | rev | sort -n); do  
	InStringAA="$InStringAA gene_pred/ORF_finder_sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_finder_preds_$GRP""_sp.aa";  
	InStringNeg="$InStringNeg gene_pred/ORF_finder_sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_finder_preds_$GRP""_sp_neg.aa";  
	InStringTab="$InStringTab gene_pred/ORF_finder_sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_finder_preds_$GRP""_sp.tab";
	InStringTxt="$InStringTxt gene_pred/ORF_finder_sigP/$Organism/$Strain/split/"$Organism"_"$Strain"_ORF_finder_preds_$GRP""_sp.txt";  
	done
	cat $InStringAA > gene_pred/ORF_finder_sigP/$Organism/$Strain/"$Strain"_ORF_finder_sp.aa
	cat $InStringNeg > gene_pred/ORF_finder_sigP/$Organism/$Strain/"$Strain"_ORF_finder_neg_sp.aa
	tail -n +2 -q $InStringTab > gene_pred/ORF_finder_sigP/$Organism/$Strain/"$Strain"_ORF_finder_sp.tab
	cat $InStringTxt > gene_pred/ORF_finder_sigP/$Organism/$Strain/"$Strain"_ORF_finder_sp.txt
	done
```
The regular expression R.LR.{,40}[ED][ED][KR] has previously been used to identfy RxLR effectors. The addition of an EER motif is significant as it has been shown as required for host uptake of the protein.

The RxLR_EER_regex_finder.py script was used to search for this regular expression and annotate the EER domain where present.

```bash
	for Secretome in $(ls gene_pred/ORF_finder_sigP/F.*/*/*_ORF_finder_sp.aa); do
		ProgDir=~/git_repos/emr_repos/tools/pathogen/RxLR_effectors;
		Strain=$(echo $Secretome | cut -d '/' -f4);
		Organism=$(echo $Secretome | cut -d '/' -f3) ;
		OutDir=analysis/RxLR_effectors/RxLR_EER_regex_finder/"$Organism"/"$Strain";
		mkdir -p $OutDir;
		printf "\nstrain: $Strain\tspecies: $Organism\n";
		printf "the number of SigP gene is:\t";
		cat $Secretome | grep '>' | wc -l;
		printf "the number of SigP-RxLR genes are:\t";
		$ProgDir/RxLR_EER_regex_finder.py $Secretome > $OutDir/"$Strain"_ORF_RxLR_EER_regex.fa;
		cat $OutDir/"$Strain"_ORF_RxLR_EER_regex.fa | grep '>' | cut -f1 | sed 's/>//g' | sed 's/ //g' > $OutDir/"$Strain"_ORF_RxLR_regex.txt
		cat $OutDir/"$Strain"_ORF_RxLR_regex.txt | wc -l
		printf "the number of SigP-RxLR-EER genes are:\t";
		cat $OutDir/"$Strain"_ORF_RxLR_EER_regex.fa | grep '>' | grep 'EER_motif_start' |  cut -f1 | sed 's/>//g' | sed 's/ //g' > $OutDir/"$Strain"_ORF_RxLR_EER_regex.txt
		cat $OutDir/"$Strain"_ORF_RxLR_EER_regex.txt | wc -l
		printf "\n"
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
		Col2=RxLR_EER_regex_finder.py
		#GeneNames=$OutDir/"$Strain"_ORF_RxLR_regex.txt
		#GeneModels=gene_pred/ORF_finder/"$Organism"/"$Strain"/"$Strain"*.aa_cat.fa
		#$ProgDir/gene_list_to_gff.pl $GeneNames $GeneModels $Col2 Name > $OutDir/"$Strain"_ORF_RxLR_regex.gff3
	done
```

The results were as follows:
```
	strain: 125	species: F.oxysporum_fsp_cepae
	the number of SigP gene is:	15312
	the number of SigP-RxLR genes are:	685
	the number of SigP-RxLR-EER genes are:	17


	strain: 55	species: F.oxysporum_fsp_cepae
	the number of SigP gene is:	16347
	the number of SigP-RxLR genes are:	724
	the number of SigP-RxLR-EER genes are:	17


	strain: A23	species: F.oxysporum_fsp_cepae
	the number of SigP gene is:	16402
	the number of SigP-RxLR genes are:	731
	the number of SigP-RxLR-EER genes are:	17


	strain: A28	species: F.oxysporum_fsp_cepae
	the number of SigP gene is:	16617
	the number of SigP-RxLR genes are:	716
	the number of SigP-RxLR-EER genes are:	16


	strain: D2	species: F.oxysporum_fsp_cepae
	the number of SigP gene is:	16169
	the number of SigP-RxLR genes are:	686
	the number of SigP-RxLR-EER genes are:	15


	strain: Fus2	species: F.oxysporum_fsp_cepae
	the number of SigP gene is:	16262
	the number of SigP-RxLR genes are:	722
	the number of SigP-RxLR-EER genes are:	18


	strain: HB17	species: F.oxysporum_fsp_cepae
	the number of SigP gene is:	16326
	the number of SigP-RxLR genes are:	720
	the number of SigP-RxLR-EER genes are:	17


	strain: PG	species: F.oxysporum_fsp_cepae
	the number of SigP gene is:	16390
	the number of SigP-RxLR genes are:	720
	the number of SigP-RxLR-EER genes are:	15


	strain: N139	species: F.oxysporum_fsp_narcissi
	the number of SigP gene is:	22180
	the number of SigP-RxLR genes are:	986
	the number of SigP-RxLR-EER genes are:	22


	strain: PG18	species: F.oxysporum_fsp_pisi
	the number of SigP gene is:	22146
	the number of SigP-RxLR genes are:	1014
	the number of SigP-RxLR-EER genes are:	14


	strain: PG3	species: F.oxysporum_fsp_pisi
	the number of SigP gene is:	17424
	the number of SigP-RxLR genes are:	740
	the number of SigP-RxLR-EER genes are:	11


	strain: A8	species: F.proliferatum
	the number of SigP gene is:	24561
	the number of SigP-RxLR genes are:	1203
	the number of SigP-RxLR-EER genes are:	27
```

### E) From ORF fragments - Hmm evidence of WY domains
```bash
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer/WY_motif.hmm
	for Proteome in $(ls gene_pred/ORF_finder/F.*/*/*.aa_cat.fa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/RxLR_effectors/hmmer_WY/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_ORF_WY_hmmer.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep 'Initial search space'
		cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
		HmmFasta="$Strain"_ORF_WY_hmmer.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
	done
```
Results were as follows:
```
	F.oxysporum_fsp_cepae 125
	Initial search space (Z):             357389  [actual number of targets]
	Domain search space  (domZ):              11  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae 55
	Initial search space (Z):             356765  [actual number of targets]
	Domain search space  (domZ):              11  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae A23
	Initial search space (Z):             357686  [actual number of targets]
	Domain search space  (domZ):              11  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae A28
	Initial search space (Z):             365675  [actual number of targets]
	Domain search space  (domZ):              10  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae D2
	Initial search space (Z):             352139  [actual number of targets]
	Domain search space  (domZ):               9  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae Fus2
	Initial search space (Z):             353490  [actual number of targets]
	Domain search space  (domZ):              11  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae HB17
	Initial search space (Z):             355148  [actual number of targets]
	Domain search space  (domZ):              11  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae PG
	Initial search space (Z):             356817  [actual number of targets]
	Domain search space  (domZ):              15  [number of targets reported over threshold]
	F.oxysporum_fsp_narcissi N139
	Initial search space (Z):             449092  [actual number of targets]
	Domain search space  (domZ):               9  [number of targets reported over threshold]
	F.oxysporum_fsp_pisi PG18
	Initial search space (Z):             445079  [actual number of targets]
	Domain search space  (domZ):              13  [number of targets reported over threshold]
	F.oxysporum_fsp_pisi PG3
	Initial search space (Z):             384712  [actual number of targets]
	Domain search space  (domZ):              10  [number of targets reported over threshold]
	F.proliferatum A8
	Initial search space (Z):             456646  [actual number of targets]
	Domain search space  (domZ):               6  [number of targets reported over threshold]
```

### F) From ORF fragments - Hmm evidence of RxLR effectors
```bash
	ProgDir=/home/armita/git_repos/emr_repos/scripts/phytophthora/pathogen/hmmer
	HmmModel=/home/armita/git_repos/emr_repos/SI_Whisson_et_al_2007/cropped.hmm
	for Proteome in $(ls gene_pred/ORF_finder/F.*/*/*.aa_cat.fa); do
		Strain=$(echo $Proteome | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		OutDir=analysis/RxLR_effectors/hmmer_RxLR/$Organism/$Strain
		mkdir -p $OutDir
		HmmResults="$Strain"_ORF_RxLR_hmmer.txt
		hmmsearch -T 0 $HmmModel $Proteome > $OutDir/$HmmResults
		echo "$Organism $Strain"
		cat $OutDir/$HmmResults | grep 'Initial search space'
		cat $OutDir/$HmmResults | grep 'number of targets reported over threshold'
		HmmFasta="$Strain"_ORF_RxLR_hmmer.fa
		$ProgDir/hmmer2fasta.pl $OutDir/$HmmResults $Proteome > $OutDir/$HmmFasta
	done
```

The results were as follows:

```
	F.oxysporum_fsp_cepae 125
	Initial search space (Z):             357389  [actual number of targets]
	Domain search space  (domZ):              46  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae 55
	Initial search space (Z):             356765  [actual number of targets]
	Domain search space  (domZ):              46  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae A23
	Initial search space (Z):             357686  [actual number of targets]
	Domain search space  (domZ):              46  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae A28
	Initial search space (Z):             365675  [actual number of targets]
	Domain search space  (domZ):              45  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae D2
	Initial search space (Z):             352139  [actual number of targets]
	Domain search space  (domZ):              41  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae Fus2
	Initial search space (Z):             353490  [actual number of targets]
	Domain search space  (domZ):              46  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae HB17
	Initial search space (Z):             355148  [actual number of targets]
	Domain search space  (domZ):              46  [number of targets reported over threshold]
	F.oxysporum_fsp_cepae PG
	Initial search space (Z):             356817  [actual number of targets]
	Domain search space  (domZ):              51  [number of targets reported over threshold]
	F.oxysporum_fsp_narcissi N139
	Initial search space (Z):             449092  [actual number of targets]
	Domain search space  (domZ):              58  [number of targets reported over threshold]
	F.oxysporum_fsp_pisi PG18
	Initial search space (Z):             445079  [actual number of targets]
	Domain search space  (domZ):              48  [number of targets reported over threshold]
	F.oxysporum_fsp_pisi PG3
	Initial search space (Z):             384712  [actual number of targets]
	Domain search space  (domZ):              35  [number of targets reported over threshold]
	F.proliferatum A8
	Initial search space (Z):             456646  [actual number of targets]
	Domain search space  (domZ):              45  [number of targets reported over threshold]
```

### G) Combined evidence file

Evidence of RxLR effectors from different sources was combined into a single gff
feature file.

This was done using the following commands:
<!--
```bash
		ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff
    # Infiles
    GffAug=gene_pred/augustus_unmasked/P.cactorum/414/414_augustus_preds.gtf
    ORF_RxLRs=analysis/rxlr_atg_unmasked/P.cactorum/414/414_ORF_sp_rxlr.gff3
    ORF_WYs=analysis/hmmer/WY/P.cactorum/414/414_ORF_WY_hmmer.gff3
    # Outfiles
    AugDB=414_aug.db
    WyDB=414_WY.db
    WyID=WY_id.txt
    WyDB_mod=414_WY_note.db
    RxlrDB=414_rxlr.db
    RxlrID=rxlr_id.txt
    RxlrDB_mod=414_rxlr_note.db
    Rxlr_Wy_DB=414_rxlr_WY.db

    OrfMerged=414_rxlr_WY_merged.db
    MergedDB=414_Aug_ORF_merged.db
    FinalDB=414_Aug_ORF.db
    FinalGff=414_Aug_ORF.gff
```

Make a db of aug genes and effector ORFs

```bash
$ProgDir/make_gff_database.py --inp $GffAug --db $AugDB
	$ProgDir/make_gff_database.py --inp $ORF_WYs --db $WyDB
	$ProgDir/make_gff_database.py --inp $ORF_RxLRs --db $RxlrDB
``` -->


## Mimps

The presence of Mimp promotors in Fusarium genomes were identified. This was
done in three steps:
 * Position of Mimps were identified in the genome
 * Genes within 1000bp downstream of the mimp were identified from Augustus
predictions
 * ORFs within 1000bp downstream of the mimp were identified from ORF
predictions

### A) Position of Mimps were identified
Position of Mimps
gff predictions for the position of mimps in the genome were identified
Fus2 was performed separately to ensure mimps are predicted from the correct assembly

```bash
	ProgDir="/home/armita/git_repos/emr_repos/tools/pathogen/mimp_finder"
	for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa); do
		Organism=$(echo "$Genome" | rev | cut -d '/' -f4 | rev)
		Strain=$(echo "$Genome" | rev | cut -d '/' -f3 | rev)
		OutDir=analysis/mimps/$Organism/"$Strain"
		mkdir -p "$OutDir"
		"$ProgDir"/mimp_finder.pl "$Genome" "$OutDir"/"$Strain"_mimps.fa "$OutDir"/"$Strain"_mimps.gff3 > "$OutDir"/"$Strain"_mimps.log
	done
```

```bash
ProgDir="/home/armita/git_repos/emr_repos/tools/pathogen/mimp_finder"
# for Genome in $(ls assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_AssemblyScaffolds.fasta); do
for Genome in $(ls assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Ma_et_al_2010/fusarium_oxysporum_f._sp._lycopersici_mitochondrion_2_contigs.fasta); do
Organism=$(echo "$Genome" | rev | cut -d '/' -f4 | rev)
Strain=$(echo "$Genome" | rev | cut -d '/' -f3 | rev)
OutDir=analysis/mimps/$Organism/"$Strain"
mkdir -p "$OutDir"
"$ProgDir"/mimp_finder.pl "$Genome" "$OutDir"/"$Strain"_mimps.fa "$OutDir"/"$Strain"_mimps.gff3 > "$OutDir"/"$Strain"_mimps.log
done
```

### B) Augustus genes flanking Mimps


<!-- ```bash
	ProgDir=~/git_repos/emr_repos/tools/pathogen/mimp_finder
	for Mimps in $(ls -d analysis/mimps/F.*/*/*_mimps.gff3); do
		Organism=$(echo "$Mimps" | rev | cut -d '/' -f3 | rev)
		Strain=$(echo "$Mimps" | rev | cut -f2 -d '/' | rev)
		MimpDir=$(dirname $Mimps)
		echo "$Mimps"
		"$ProgDir"/gffexpander.pl + 1500 "$Mimps" > "$MimpDir"/"$Strain"_mimps_expanded.gff3
		StrainAugModels=$(ls gene_pred/augustus/$Organism/"$Strain"/*_augustus_preds.gtf)
		StrainORFModels=$(ls gene_pred/ORF_finder/$Organism/"$Strain"/*_ORF_mod.gff)
		bedtools intersect -s -a "$StrainAugModels"  -b "$MimpDir"/"$Strain"_mimps_expanded.gff3 > "$MimpDir"/"$Strain"_mimps_intersected_Aug_genes.gff
		cat "$MimpDir"/"$Strain"_mimps_intersected_Aug_genes.gff | grep 'gene' | rev | cut -f1 -d '=' | rev | sort | uniq > "$MimpDir"/"$Strain"_mimps_intersected_Aug_genes_names.txt
		bedtools intersect -s -a "$StrainORFModels"  -b "$MimpDir"/"$Strain"_mimps_expanded.gff3 > "$MimpDir"/"$Strain"_mimps_intersected_ORF_genes.gff
		cat "$MimpDir"/"$Strain"_mimps_intersected_ORF_genes.gff | grep 'gene' | rev | cut -f1 -d '=' | rev | sort | uniq > "$MimpDir"/"$Strain"_mimps_intersected_ORF_genes_names.txt
		printf "Augustus genes intersected:\t"
		cat "$MimpDir"/"$Strain"_mimps_intersected_Aug_genes_names.txt | wc -l
		printf "ORF fragments intersected:\t"
		cat "$MimpDir"/"$Strain"_mimps_intersected_ORF_genes_names.txt | wc -l
		echo ""
	done
```

Results were as follows:
```
	analysis/mimps/F.oxysporum_fsp_cepae/125/125_mimps.gff3
	Augustus genes intersected:	21
	ORF fragments intersected:	285

	analysis/mimps/F.oxysporum_fsp_cepae/55/55_mimps.gff3
	Augustus genes intersected:	15
	ORF fragments intersected:	254

	analysis/mimps/F.oxysporum_fsp_cepae/A23/A23_mimps.gff3
	Augustus genes intersected:	20
	ORF fragments intersected:	307

	analysis/mimps/F.oxysporum_fsp_cepae/A28/A28_mimps.gff3
	Augustus genes intersected:	5
	ORF fragments intersected:	99

	analysis/mimps/F.oxysporum_fsp_cepae/D2/D2_mimps.gff3
	Augustus genes intersected:	0
	ORF fragments intersected:	87

	analysis/mimps/F.oxysporum_fsp_cepae/Fus2/Fus2_mimps.gff3
	Augustus genes intersected:	20
	ORF fragments intersected:	260

	analysis/mimps/F.oxysporum_fsp_cepae/HB17/HB17_mimps.gff3
	Augustus genes intersected:	29
	ORF fragments intersected:	328

	analysis/mimps/F.oxysporum_fsp_cepae/PG/PG_mimps.gff3
	Augustus genes intersected:	5
	ORF fragments intersected:	170

	analysis/mimps/F.oxysporum_fsp_narcissi/N139/N139_mimps.gff3
	Augustus genes intersected:	15
	ORF fragments intersected:	395

	analysis/mimps/F.oxysporum_fsp_pisi/PG18/PG18_mimps.gff3
	Augustus genes intersected:	5
	ORF fragments intersected:	155

	analysis/mimps/F.oxysporum_fsp_pisi/PG3/PG3_mimps.gff3
	Augustus genes intersected:	3
	ORF fragments intersected:	146

	analysis/mimps/F.proliferatum/A8/A8_mimps.gff3
	Augustus genes intersected:	2
	ORF fragments intersected:	20
``` -->

```bash
	ProgDir=~/git_repos/emr_repos/tools/pathogen/mimp_finder
	for Mimps in $(ls -d analysis/mimps/F.*/*/*_mimps.gff3); do
		Organism=$(echo "$Mimps" | rev | cut -d '/' -f3 | rev)
		Strain=$(echo "$Mimps" | rev | cut -f2 -d '/' | rev)
		OutDir=analysis/mimps_+-2000bp/$Organism/$Strain
		mkdir -p $OutDir
		MimpDir=$(dirname $Mimps)
		echo "$Mimps"
		"$ProgDir"/gffexpander.pl +- 2000 "$Mimps" > "$OutDir"/"$Strain"_mimps_2000bp_expanded.gff3
		StrainAugModels=$(ls gene_pred/augustus/$Organism/"$Strain"/*_augustus_preds.gtf)
		StrainORFModels=$(ls gene_pred/ORF_finder/$Organism/"$Strain"/*_ORF_mod.gff)
		bedtools intersect -a "$StrainAugModels"  -b "$OutDir"/"$Strain"_mimps_2000bp_expanded.gff3 > "$OutDir"/"$Strain"_mimps_intersected_Aug_genes.gff
		#bedtools intersect -s -a "$StrainAugModels"  -b "$OutDir"/"$Strain"_mimps_2000bp_expanded.gff3 > "$OutDir"/"$Strain"_mimps_intersected_Aug_genes.gff
		cat "$OutDir"/"$Strain"_mimps_intersected_Aug_genes.gff | grep 'gene' | rev | cut -f1 -d '=' | rev | sort | uniq > "$OutDir"/"$Strain"_mimps_intersected_Aug_genes_names.txt
		bedtools intersect -s -a "$StrainORFModels"  -b "$OutDir"/"$Strain"_mimps_2000bp_expanded.gff3 > "$OutDir"/"$Strain"_mimps_intersected_ORF_genes.gff
		cat "$OutDir"/"$Strain"_mimps_intersected_ORF_genes.gff | grep 'gene' | rev | cut -f1 -d '=' | rev | sort | uniq > "$OutDir"/"$Strain"_mimps_intersected_ORF_genes_names.txt
		printf "Augustus genes intersected:\t"
		cat "$OutDir"/"$Strain"_mimps_intersected_Aug_genes_names.txt | wc -l
		printf "ORF fragments intersected:\t"
		cat "$OutDir"/"$Strain"_mimps_intersected_ORF_genes_names.txt | wc -l
		echo ""
	done
```

```
	analysis/mimps/F.oxysporum_fsp_cepae/125/125_mimps.gff3
	Augustus genes intersected:	62
	ORF fragments intersected:	767

	analysis/mimps/F.oxysporum_fsp_cepae/55/55_mimps.gff3
	Augustus genes intersected:	47
	ORF fragments intersected:	757

	analysis/mimps/F.oxysporum_fsp_cepae/A23/A23_mimps.gff3
	Augustus genes intersected:	64
	ORF fragments intersected:	755

	analysis/mimps/F.oxysporum_fsp_cepae/A28/A28_mimps.gff3
	Augustus genes intersected:	17
	ORF fragments intersected:	393

	analysis/mimps/F.oxysporum_fsp_cepae/D2/D2_mimps.gff3
	Augustus genes intersected:	11
	ORF fragments intersected:	326

	analysis/mimps/F.oxysporum_fsp_cepae/Fus2/Fus2_mimps.gff3
	Augustus genes intersected:	60
	ORF fragments intersected:	654

	analysis/mimps/F.oxysporum_fsp_cepae/HB17/HB17_mimps.gff3
	Augustus genes intersected:	68
	ORF fragments intersected:	830

	analysis/mimps/F.oxysporum_fsp_cepae/PG/PG_mimps.gff3
	Augustus genes intersected:	16
	ORF fragments intersected:	415

	analysis/mimps/F.oxysporum_fsp_narcissi/N139/N139_mimps.gff3
	Augustus genes intersected:	39
	ORF fragments intersected:	951

	analysis/mimps/F.oxysporum_fsp_pisi/PG18/PG18_mimps.gff3
	Augustus genes intersected:	20
	ORF fragments intersected:	470

	analysis/mimps/F.oxysporum_fsp_pisi/PG3/PG3_mimps.gff3
	Augustus genes intersected:	19
	ORF fragments intersected:	425

	analysis/mimps/F.proliferatum/A8/A8_mimps.gff3
	Augustus genes intersected:	5
	ORF fragments intersected:	58
```

### for FoL reference genome

```bash
ProgDir=~/git_repos/emr_repos/tools/pathogen/mimp_finder
for Mimps in $(ls -d analysis/mimps/F.oxysporum_fsp_lycopersici/4287/*_mimps.gff3); do
Organism=$(echo "$Mimps" | rev | cut -d '/' -f3 | rev)
Strain=$(echo "$Mimps" | rev | cut -f2 -d '/' | rev)
OutDir=analysis/mimps_+-2000bp/$Organism/$Strain
mkdir -p $OutDir
MimpDir=$(dirname $Mimps)
echo "$Mimps"
"$ProgDir"/gffexpander.pl +- 2000 "$Mimps" > "$OutDir"/"$Strain"_mimps_2000bp_expanded.gff3
StrainAugModels=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1.filtered_proteins.ExternalModels.gff3
bedtools intersect -a "$StrainAugModels"  -b "$OutDir"/"$Strain"_mimps_2000bp_expanded.gff3 > "$OutDir"/"$Strain"_mimps_intersected_Aug_genes.gff
cat "$OutDir"/"$Strain"_mimps_intersected_Aug_genes.gff | grep 'gene' | rev | cut -f1 -d '=' | rev | sort | uniq > "$OutDir"/"$Strain"_mimps_intersected_Aug_genes_names.txt
bedtools intersect -s -a "$StrainORFModels"  -b "$OutDir"/"$Strain"_mimps_2000bp_expanded.gff3 > 		printf "Augustus genes intersected:\t"
cat "$OutDir"/"$Strain"_mimps_intersected_Aug_genes_names.txt | wc -l
echo ""
done
```

## Merging effector evidence:

```bash
	for Strain in $(ls -d gene_pred/ORF_finder/F.*/* | rev | cut -f1 -d '/' | rev); do
		Queue=$(qstat | grep 'merge_fus' | sed -E 's/ +/ /g' |  cut -f5 -d ' ' | grep -c 'qw' | sed 's/ //g')
		while [[ $Queue -ge 1 ]]; do
			sleep 10
			printf "."
			Queue=$(qstat | grep 'merge_fus' | sed -E 's/ +/ /g' |  cut -f5 -d ' ' | grep -c 'qw')
		done
		# Node1Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace01' | sed 's/ //g')
		Node2Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace02' | sed 's/ //g')
		# Node3Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace03' | sed 's/ //g')
		Node4Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace04' | sed 's/ //g')
		#Node5Jobs=$(qstat | grep 'merge_fus' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace05' | sed 's/ //g')
		Node6Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace06' | sed 's/ //g')
		Node11Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace11' | sed 's/ //g')
		# while [[ $Node1Jobs -ge 1  && $Node2Jobs -ge 1 && $Node3Jobs -ge 1 && $Node4Jobs -ge 1 ]]; do
		while [[ $Node2Jobs -ge 1 && $Node4Jobs -ge 1 && $Node6Jobs -ge 1 && $Node11Jobs -ge 1 ]]; do
			#while [ [ $Node1Jobs -ge 1 ] && [ $Node2Jobs -ge 1 ] && [ $Node3Jobs -ge 1 ] && [ $Node4Jobs -ge 1 ] && [ $Node5Jobs -ge 1 ] ]; do
			sleep 10
			printf "."
			# Node1Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace01' | sed 's/ //g')
			Node2Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace02' | sed 's/ //g')
			# Node3Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace03' | sed 's/ //g')
			Node4Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace04' | sed 's/ //g')
			#Node5Jobs=$(qstat | grep 'merge_fus' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace05' | sed 's/ //g')
			Node6Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace06' | sed 's/ //g')
			Node11Jobs=$(qstat | grep 'merge' | sed -E 's/ +/ /g' | cut -f 8 -d ' ' | sort | grep -c 'blacklace11' | sed 's/ //g')
		done
		ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/merge_gff;
		echo $Strain;
		Orf_Gff=$(ls gene_pred/ORF_finder/*/$Strain/"$Strain"_ORF.gff);
		Aug_Gff=$(ls gene_pred/augustus/*/"$Strain"/"$Strain"_augustus_preds.gtf | grep -v 'old');
		echo $Orf_Gff;
		echo $Aug_Gff;
		# if [[ $Node1Jobs -lt 1 ]]; then
		# 	qsub -l h=blacklace01.blacklace $ProgDir/merge_fus_effectors.sh $Orf_Gff $Aug_Gff;
		# elif [[ $Node2Jobs -lt 1 ]]; then
		if [[ $Node2Jobs -lt 1 ]]; then
			qsub -l h=blacklace02.blacklace $ProgDir/merge_fus_effectors.sh $Orf_Gff $Aug_Gff;
		# elif [[ $Node3Jobs -lt 1 ]]; then
		# 	qsub -l h=blacklace03.blacklace $ProgDir/merge_fus_effectors.sh $Orf_Gff $Aug_Gff;
		elif [[ $Node4Jobs -lt 1 ]]; then
			qsub -l h=blacklace04.blacklace $ProgDir/merge_fus_effectors.sh $Orf_Gff $Aug_Gff;
			# elif [[ $Node5Jobs -lt 1 ]]; then
			# qsub -l h=blacklace05.blacklace $ProgDir/merge_fus_effectors.sh $Orf_Gff $Aug_Gff;
		elif [[ $Node6Jobs -lt 1 ]]; then
			qsub -l h=blacklace06.blacklace $ProgDir/merge_fus_effectors.sh $Orf_Gff $Aug_Gff;
		elif [[ $Node11Jobs -lt 1 ]]; then
			qsub -l h=blacklace11.blacklace $ProgDir/merge_fus_effectors.sh $Orf_Gff $Aug_Gff;
		else
			echo "Error something has jumped the queue"
		fi
	done
```

# 4. Genomic analysis

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


## 4.2 Orthology

Orthomcl was used to identify orthologous groups between Fusarium spp. genomes

Genomes were grouped by subspecies and orthology was determined within each
subspecies group. Orthology was also determined between subspecies groups.

| Pathogenic | non-pathogenic | Intermediate |
| ---------- | -------------- | -------------|
| 125        | A28            | 55           |
| A23        | D2             |              |
| Fus2       | PG             |              |




### 4.2.a) Orthology between pathogenic isolates

The Commands used to run this analysis are shown in
pathogen/orthology/F.oxysporum_fsp.cepae_pathogenic_orthology.md


### 4.2.b) Orthology between non-pathogenic isolates

The Commands used to run this analysis are shown in
pathogen/orthology/F.oxysporum_fsp.cepae_non-pathogenic_orthology.md


### 4.2.c) Orthology between pathogenic and non-pathogenic isolates

The Commands used to run this analysis are shown in
pathogen/orthology/F.oxysporum_fsp.cepae_pathogen_vs_non-pathogen_orthology.md


### 4.2.d) Orthology between all isolates

The Commands used to run this analysis are shown in
pathogen/orthology/F.oxysporum_fsp.cepae_isolates.md


## 5. BLAST Searches

## 5.1 Identifying SIX genes

Protein sequence of previously characterised SIX genes used to BLAST against
assemlies.

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	Query=analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX.fa
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa); do
		echo $Assembly
		qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
	done
```

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	BlastHits=analysis/blast_homology/FeChina/dip_spades/dip_spades_six-appended_parsed.fa_homologs.csv
	HitsGff=analysis/blast_homology/FeChina/dip_spades/dip_spades_six-appended_parsed.fa_homologs.gff
	Column2=SIX_homolog
	NumHits=5
	$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
```

## 5.2 Identifying PHIbase homologs

The PHIbase database was searched agasinst the assembled genomes using tBLASTx.

```bash
	for Assembly in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa); do
		qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh analysis/blast_homology/PHIbase/PHI_36_accessions.fa protein $Assembly
	done
```

following blasting PHIbase to the genome, the hits were filtered by effect on
virulence.

First the a tab seperated file was made in the clusters core directory containing
PHIbase. These commands were run as part of previous projects but have been
included here for completeness.

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

The analysis was also performed by blasting the predicted proteins against the
PHIbase database:

The PHIbase database was searched agasinst the assembled genomes using tBLASTx.

```bash
	for Proteins in $(ls gene_pred/augustus/F.oxysporum_fsp_cepae/*/*_augustus_preds.aa); do
		qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh $Proteins protein ../../phibase/v3.8/PHI_accessions.fa
	done
```
