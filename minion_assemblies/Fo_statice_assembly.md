# Nectria_haematococca
Commands for the analysis of Nectria haematococca (F. solani) ex. pea genomes.

This document details the commands used to assemble and annotate minion sequence data.

Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/fusarium


# 0. Building of directory structure

### Minion Data

```bash
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	mkdir -p $ProjectDir/raw_dna/minion/F.oxysporum_fsp_statice/Stat10
```


Basecalling from the gridion was used:

```bash
Organism=F.oxysporum_fsp_statice
Strain=Stat10
Date=2018-05-15
RawDatDir=/data/seq_data/minion/2018/20180504_Statice10-180501/Statice10-180501/GA10000
OutDir=$ProjectDir/raw_dna/minion/$Organism/$Strain

cat $RawDatDir/*.fastq | gzip -cf > $OutDir/${Strain}_${Date}_albacore_v2.2.7.fastq.gz
```

<!--
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
``` -->


## Assembly

### Removal of adapters

Splitting reads and trimming adapters using porechop
```bash
	for RawReads in $(ls raw_dna/minion/*/*/*.fastq.gz | grep 'statice'); do
    Organism=$(echo $RawReads| rev | cut -f3 -d '/' | rev)
    Strain=$(echo $RawReads | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
  	OutDir=qc_dna/minion/$Organism/$Strain
  	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
  	qsub $ProgDir/sub_porechop.sh $RawReads $OutDir
  done
```
<!--
#### QC of MiSeq data

programs:
  fastqc
  fastq-mcf
  kmc

Data quality was visualised using fastqc:
```bash
for RawData in $(ls raw_dna/paired/*/*/*/*.fastq.gz | grep 'statice'); do
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
``` -->

# Identify sequencing coverage

For Minion data:
```bash
for RawData in $(ls qc_dna/minion/*/*/*q.gz | grep 'statice'); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
GenomeSz=58
OutDir=$(dirname $RawData)
qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
done
```

```bash
  for StrainDir in $(ls -d qc_dna/minion/*/*); do
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

```
<!--
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
-->

### Read correction using Canu

```bash
for TrimReads in $(ls qc_dna/minion/*/*/*q.gz | grep 'statice'); do
Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
OutDir=assembly/canu-1.6/$Organism/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
qsub $ProgDir/sub_canu_correction.sh $TrimReads 58m $Strain $OutDir
done
```


### Assembbly using SMARTdenovo

```bash
for CorrectedReads in $(ls assembly/canu-1.6/*/*/*.trimmedReads.fasta.gz | grep 'statice'); do
Organism=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $CorrectedReads | rev | cut -f2 -d '/' | rev)
Prefix="$Strain"_smartdenovo
OutDir=assembly/SMARTdenovo/$Organism/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/SMARTdenovo
qsub $ProgDir/sub_SMARTdenovo.sh $CorrectedReads $Prefix $OutDir
done
```


Quast and busco were run to assess the effects of racon on assembly quality:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg | grep 'statice'); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
	echo "$Organism - $Strain"
  OutDir=$(dirname $Assembly)
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
	OutDir=gene_pred/busco/$Organism/$Strain/assembly
	BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
	qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```


Error correction using racon:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg | grep 'statice'); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ReadsFq=$(ls qc_dna/minion/*/$Strain/*q.gz)
Iterations=10
OutDir=$(dirname $Assembly)"/racon_$Iterations"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/racon
qsub $ProgDir/sub_racon.sh $Assembly $ReadsFq $Iterations $OutDir
done
```


```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/*.fasta | grep 'round_10' | grep 'statice'); do
OutDir=$(dirname $Assembly)
echo "" > tmp.txt
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/racon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
```

Quast and busco were run to assess the effects of racon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta | grep 'statice'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/F.*/*/racon*/*.fasta | grep 'statice'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```
```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/A*/*/assembly/*/short_summary_*.txt); do
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```
