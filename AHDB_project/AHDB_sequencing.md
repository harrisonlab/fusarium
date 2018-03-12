# AHDB sequencing
==========

This document details the commands used to assemble and annotate genomes
sequenced as part of the AHDB fusarium project. This excludes the commands used
assemble FoM genome for isolate Stocks4 and FoN isoalte FON_63 whose commands
are documented in the minion_assembly directory of this repository.

Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/fusarium

The following is a summary of the work presented in this Readme.

# 0. Building of directory structure

### Minion Data

```bash
  Organism="combined_samples"
  Strain="FON77_FON139"
	RawDatDir=/data/seq_data/minion/2017/20170920_1225_FON77\ -\ FON139
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	mkdir -p $ProjectDir/raw_dna/minion/$Organism/$Strain
```

Sequence data was moved into the appropriate directories

```bash
	RawDatDir=/data/seq_data/minion/2017/20170920_1225_FON77-FON139
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
  nanopolish index -r $RawDatDir/fast5/pass | gzip -cf > $ProjectDir/raw_dna/minion/$Organism/$Strain/"$Strain"_reads.fa.gz
	cp $RawDatDir/all_reads_albacore1.1.1.fastq.gz $ProjectDir/raw_dna/minion/F.oxysporum_fsp_narcissi/FON_63/.
```

### MiSeq data

```bash
  # 160725
  RawDatDir=/home/scratch/groups/harrisonlab/fusarium/160725_M04465_0020_000000000-AP1N0/160725_M04465_0020_000000000-AP1N0/Data/Intensities/BaseCalls
  ProjectDir=/home/groups/harrisonlab/project_files/fusarium
  OutDir=$ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/AF064
  mkdir -p $OutDir/F
  mkdir -p $OutDir/R
  cd $OutDir/F
  cp -s $RawDatDir/AF064_S3_L001_R1_001.fastq.gz .
  cd $OutDir/R
  cp -s $RawDatDir/AF064_S3_L001_R2_001.fastq.gz .
  OutDir=$ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/CH5-2
  mkdir -p $OutDir/F
  mkdir -p $OutDir/R
  cd $OutDir/F
  cp -s $RawDatDir/CH5-2_S2_L001_R1_001.fastq.gz .
  cd $OutDir/R
  cp -s $RawDatDir/CH5-2_S2_L001_R2_001.fastq.gz .

  # 171006
	RawDatDir=/data/seq_data/miseq/2017/RAW/171006_M04465_0050_000000000-B49T7/Data/Intensities/BaseCalls
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	OutDir=$ProjectDir/raw_dna/paired/F.oxysporum_fsp_narcissi/FON129
	mkdir -p $OutDir/F
	mkdir -p $OutDir/R
  cd $OutDir/F
	cp -s $RawDatDir/FON129_S2_L001_R1_001.fastq.gz .
  cd $OutDir/R
	cp -s $RawDatDir/FON129_S2_L001_R2_001.fastq.gz .
  # 170926
  RawDatDir=/data/seq_data/miseq/2017/RAW/170926_M04465_0049_000000000-B49VF/Data/Intensities/BaseCalls
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	OutDir=$ProjectDir/raw_dna/paired/F.oxysporum_fsp_narcissi/FON77
	mkdir -p $OutDir/F
	mkdir -p $OutDir/R
  cd $OutDir/F
	cp -s $RawDatDir/FON77_S1_L001_R1_001.fastq.gz .
  cd $OutDir/R
	cp -s $RawDatDir/FON77_S1_L001_R2_001.fastq.gz .
  OutDir=$ProjectDir/raw_dna/paired/F.oxysporum_fsp_narcissi/FON81
  mkdir -p $OutDir/F
  mkdir -p $OutDir/R
  cd $OutDir/F
  cp -s $RawDatDir/FON81_S2_L001_R1_001.fastq.gz .
  cd $OutDir/R
  cp -s $RawDatDir/FON81_S2_L001_R2_001.fastq.gz .
  OutDir=$ProjectDir/raw_dna/paired/F.oxysporum_fsp_narcissi/FON89
  mkdir -p $OutDir/F
  mkdir -p $OutDir/R
  cd $OutDir/F
  cp -s $RawDatDir/FON89_S3_L001_R1_001.fastq.gz .
  cd $OutDir/R
  cp -s $RawDatDir/FON89_S3_L001_R2_001.fastq.gz .
  OutDir=$ProjectDir/raw_dna/paired/F.oxysporum_fsp_fragariae/Straw465
  mkdir -p $OutDir/F
  mkdir -p $OutDir/R
  cd $OutDir/F
  cp -s $RawDatDir/Straw465_S4_L001_R1_001.fastq.gz .
  cd $OutDir/R
  cp -s $RawDatDir/Straw465_S4_L001_R2_001.fastq.gz .
  RawDatDir=/data/seq_data/miseq/2017/RAW/171016_M04465_0051_000000000-B43DD/Data/Intensities/BaseCalls
  ProjectDir=/home/groups/harrisonlab/project_files/fusarium
  OutDir=$ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/F81
  mkdir -p $OutDir/F
  mkdir -p $OutDir/R
  cd $OutDir/F
  cp -s $RawDatDir/F81_S3_L001_R1_001.fastq.gz .
  cd $OutDir/R
  cp -s $RawDatDir/F81_S3_L001_R2_001.fastq.gz .
  OutDir=$ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/FOP1-EMR
  mkdir -p $OutDir/F
  mkdir -p $OutDir/R
  cd $OutDir/F
  cp -s $RawDatDir/FOP1-EMR_S2_L001_R1_001.fastq.gz .
  cd $OutDir/R
  cp -s $RawDatDir/FOP1-EMR_S2_L001_R2_001.fastq.gz .
  OutDir=$ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/R2
  mkdir -p $OutDir/F
  mkdir -p $OutDir/R
  cd $OutDir/F
  cp -s $RawDatDir/R2_S1_L001_R1_001.fastq.gz .
  cd $OutDir/R
  cp -s $RawDatDir/R2_S1_L001_R2_001.fastq.gz .
  OutDir=$ProjectDir/raw_dna/paired/F.oxysporum/15-074
  mkdir -p $OutDir/F
  mkdir -p $OutDir/R
  cd $OutDir/F
  cp -s $RawDatDir/15-074_S4_L001_R1_001.fastq.gz .
  cd $OutDir/R
  cp -s $RawDatDir/15-074_S4_L001_R2_001.fastq.gz .


  cd $ProjectDir
```

#### QC of MiSeq data

programs:
  fastqc
  fastq-mcf
  kmc

Data quality was visualised using fastqc:
```bash
  for RawData in $(ls raw_dna/paired/*/*/*/*.fastq.gz | grep -w -e 'FON129' -e 'FON77' -e 'FON81' -e 'FON89' -e 'Straw465' -e 'F81' -e 'FOP1-EMR' -e 'R2' -e '15-074' -e 'AF064' -e 'CH5-2'| grep -e 'AF064' -e 'CH5-2'); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    echo $RawData;
    qsub $ProgDir/run_fastqc.sh $RawData
  done
```

```bash
  for StrainPath in $(ls -d raw_dna/paired/*/* | grep -w -e 'FON129' -e 'FON77' -e 'FON81' -e 'FON89' -e 'Straw465' -e 'F81' -e 'FOP1-EMR' -e 'R2' -e '15-074' -e 'AF064' -e 'CH5-2'| grep -e 'AF064' -e 'CH5-2'); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
    IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
    ReadsF=$(ls $StrainPath/F/*.fastq*)
    ReadsR=$(ls $StrainPath/R/*.fastq*)
    echo $ReadsF
    echo $ReadsR
    qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
  done
```

Data quality was visualised once again following trimming:

```bash
  for TrimData in $(ls qc_dna/paired/*/*/*/*.fq.gz | grep -w -e 'FON129' -e 'FON77' -e 'FON81' -e 'FON89' -e 'Straw465' -e 'F81' -e 'FOP1-EMR' -e 'R2' -e '15-074'); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    echo $RawData;
    qsub $ProgDir/run_fastqc.sh $RawData
  done
```

# Assembly


## Identify sequencing coverage
<!--
For Minion data:
```bash
	for RawData in $(ls qc_dna/minion/*/*/*q.gz | grep 'FON_63'); do
		echo $RawData;
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
		GenomeSz=60
		OutDir=$(dirname $RawData)
		qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
	done
```

```bash
  for StrainDir in $(ls -d qc_dna/minion/*/* | grep 'FON_63'); do
    Strain=$(basename $StrainDir)
    printf "$Strain\t"
    for File in $(ls $StrainDir/*.txt); do
      echo $(basename $File);
      cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
    done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
  done
```
MinION coverage
```
FON_63	122.23
```
-->

For Miseq data:
```bash
for RawData in $(ls qc_dna/paired/*/*/*/*.fq.gz | grep -w -e 'FON129' -e 'FON77' -e 'FON81' -e 'FON89' -e 'Straw465' -e 'F81' -e 'FOP1-EMR' -e 'R2' -e '15-074' -e 'AF064' -e 'CH5-2'| grep -e 'AF064' -e 'CH5-2'); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData;
GenomeSz=60
OutDir=$(dirname $RawData)
qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
done
```

```bash
  for StrainDir in $(ls -d qc_dna/paired/*/* | grep -w -e 'FON129' -e 'FON77' -e 'FON81' -e 'FON89' -e 'Straw465' -e 'F81' -e 'FOP1-EMR' -e 'R2' -e '15-074'); do
    Strain=$(basename $StrainDir)
    printf "$Strain\t"
    for File in $(ls qc_dna/paired/*/"$Strain"/*/*.txt); do
      echo $(basename $File);
      cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
    done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
  done
```

```
15-074	76.3
Straw465	52.3
FON129	98.93
FON77	42.41
FON81	57.62
FON89	57.92
F81	48.03
FOP1-EMR	66.85
R2	55.03
```



## MiSeq Assembly



```bash
for StrainPath in $(ls -d qc_dna/paired/*/* |  grep -w -e 'FON129' -e 'FON77' -e 'FON81' -e 'FON89' -e 'Straw465' -e 'F81' -e 'FOP1-EMR' -e 'R2' -e '15-074' -e 'A1-2' -e 'HB6' -e 'D2' -e 'PG8' -e 'L5' -e 'AF064' -e 'CH5-2'| grep -e 'AF064' -e 'CH5-2'); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*.fq.gz)
R_Read=$(ls $StrainPath/R/*.fq.gz)
OutDir=assembly/spades/$Organism/$Strain
# Jobs=$(qstat | grep 'submit_SPA' | grep 'qw' | wc -l)
# while [ $Jobs -gt 0 ]; do
# sleep 1m
# printf "."
# Jobs=$(qstat | grep 'submit_SPA' | grep 'qw' | wc -l)
# done
echo "$Organism - $Strain"		
echo $F_Read
echo $R_Read
qsub $ProgDir/submit_SPAdes.sh $F_Read $R_Read $OutDir correct 10
# qsub $ProgDir/submit_SPAdes_HiMem.sh $F_Read $R_Read $OutDir correct 10
done
```

Assembly stats were collected using quast

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/spades/F.*/*/contigs.fasta | grep -w -e 'FON129' -e 'FON77' -e 'FON81' -e 'FON89' -e 'Straw465' -e 'F81' -e 'FOP1-EMR' -e 'R2' -e '15-074'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```


# Repeat Masking

Repeat masking was performed on the non-hybrid assembly.

```bash
for Assembly in $(ls assembly/spades/F.*/*/filtered_contigs/contigs_min_500bp.fasta | grep -w -e 'FON129' -e 'FON77' -e 'FON81' -e 'FON89' -e 'Straw465' -e 'F81' -e 'FOP1-EMR' -e 'R2' -e '15-074' -e 'A1-2' -e 'HB6' -e 'D2' -e 'PG8' -e 'L5' | grep -e 'A1-2' -e 'HB6'); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=repeat_masked/$Organism/"$Strain"/filtered_contigs
echo $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
done
```

The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and hardmasked files.

```bash

for File in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep 'FON_63'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/*/*_contigs_hardmasked.fa | grep 'FON_63'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```
```
Number of masked bases:
9075386
```
