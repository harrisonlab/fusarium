# F. oxysporum ex. statice
Commands for the analysis of F. oxysporum ex. statice genomes

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

Nanopolish index had problems with either the fastq or fast5 files from the
gridion - as such they were recalled and used for nanopolish later

```bash
	screen -a
	ssh nanopore@nanopore

	# upgrade albacore
	wget https://mirror.oxfordnanoportal.com/software/analysis/ont_albacore-2.2.7-cp34-cp34m-manylinux1_x86_64.whl
	~/.local/bin/read_fast5_basecaller.py --version
	pip3 install --user ont_albacore-2.2.7-cp34-cp34m-manylinux1_x86_64.whl --upgrade
	~/.local/bin/read_fast5_basecaller.py --version

	mkdir FoStatice_06-09-18
	cd FoStatice_06-09-18
	# Oxford nanopore Stocks4 run1
	Organism=F.oxysporum_fsp_statice
	Strain=Stat10
	Date=2018-05-15
	FlowCell="FLO-MIN106"
	Kit="SQK-RBK001"
	RawDatDir=/data/seq_data/minion/2018/20180504_Statice10-180501/Statice10-180501/GA10000
	OutDir=~/FoStatice_06-09-18/$Organism/$Strain/$Date
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

	cd ~/FoStatice_06-09-18
	cat $OutDir/albacore_v2.2.7/workspace/pass/unclassified/*.fastq | gzip -cf > ${Strain}_${Date}_albacore_v2.2.7.fastq.gz

	tar -cz -f ${Strain}_${Date}_albacore_v2.2.7.tar.gz $OutDir

	FinalDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.2.7
	mkdir -p $FinalDir
	mv ${Strain}_${Date}_albacore_v2.2.7.* $FinalDir/.
	chmod +rw -R $FinalDir
```


### MiSeq data

```bash
  RawDatDir=/data/seq_data/miseq/2018/RAW/180716_M04465_0083_000000000-BM8CL/Data/Intensities/BaseCalls
  ProjectDir=/home/groups/harrisonlab/project_files/fusarium
  OutDir=$ProjectDir/raw_dna/paired/F.oxysporum_fsp_statice/Stat10
  mkdir -p $OutDir/F
  mkdir -p $OutDir/R
  cd $OutDir/F
  cp -s $RawDatDir/Stat10_S1_L001_R1_001.fastq.gz .
  cd $OutDir/R
  cp -s $RawDatDir/Stat10_S1_L001_R2_001.fastq.gz .
  cd $ProjectDir
```


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
for StrainPath in $(ls -d raw_dna/paired/*/* | grep 'statice'); do
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

For Minion data:
```bash
for RawData in $(ls qc_dna/minion/*/*/*q.gz | grep 'statice'); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
GenomeSz=70
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

For Miseq data:
```bash
for RawData in $(ls qc_dna/paired/*/*/*/*q.gz | grep 'statice'); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
qsub $ProgDir/run_fastqc.sh $RawData;
GenomeSz=70
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

```


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


## Nanopolish


<!--
```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta | grep 'Stat10'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ReadsFq=$(ls qc_dna/minion/*/$Strain/*q.gz)
OutDir=qc_dna/minion/$Organism/$Strain/appended
mkdir $OutDir
cat $ReadsFq > $OutDir/${Strain}_appended.fq.gz
done
``` -->

For Stat10
```bash
# for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta | grep -e 'Stat10'); do
# Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
# Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
# echo "$Organism - $Strain"
# # Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# # Note - the full path from home must be used
# ReadDir=raw_dna/nanopolish/$Organism/$Strain
# mkdir -p $ReadDir
# ReadsFq=$(ls raw_dna/minion/*/$Strain/*.fastq.gz)
# Fast5Dir=$(ls -d /data/seq_data/minion/2018/20180504_Statice10-180501/Statice10-180501/GA10000/reads)
# nanopolish index -v -d $Fast5Dir $ReadsFq
# done

for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta | grep -e 'Stat10'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Note - the full path from home must be used

CurDir=$PWD
ScratchDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.2.7

ReadDir=/data2/scratch2/armita/FoN/raw_dna/nanopolish/$Organism/$Strain
mkdir -p $ReadDir
cd $ReadDir
# tar -zxvf $ScratchDir/Stat10_2018-05-15_albacore_v2.2.7.tar.gz
cd $CurDir
Fast5Dir=$(ls -d /data2/scratch2/armita/FoN/raw_dna/nanopolish/F.oxysporum_fsp_statice/Stat10/home/nanopore/FoStatice_06-09-18/F.oxysporum_fsp_statice/Stat10/2018-05-15/albacore_v2.2.7/workspace/pass/unclassified/)

ReadsFq=raw_dna/minion/$Organism/$Strain/${Strain}_appended.fastq.gz
cat $Fast5Dir/*.fastq > $ReadsFq


nanopolish index -d $Fast5Dir $ReadsFq

OutDir=$(dirname $Assembly)/nanopolish
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
# submit alignments for nanoppolish
qsub $ProgDir/sub_minimap2_nanopolish.sh $Assembly $ReadsFq $OutDir/nanopolish
done
```


 Split the assembly into 50Kb fragments an submit each to the cluster for
 nanopolish correction

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta | grep -e 'Stat10'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)/nanopolish
ReadsFq=$(ls raw_dna/minion/*/$Strain/*.fastq.gz | grep -v 'Stat10_2018-05-15')
AlignedReads=$(ls $OutDir/nanopolish/reads.sorted.bam)

NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_makerange.py $Assembly --segment-length 50000 > $OutDir/nanopolish_range.txt

Ploidy=1
echo "nanopolish log:" > $OutDir/nanopolish_log.txt
for Region in $(cat $OutDir/nanopolish_range.txt); do
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
done		
printf "\n"
echo $Region
echo $Region >> $OutDir/nanopolish_log.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
qsub $ProgDir/sub_nanopolish_variants.sh $Assembly $ReadsFq $AlignedReads $Ploidy $Region $OutDir/$Region
done
done
```

A subset of nanopolish jobs needed to be resubmitted as the ran out of RAM

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta | grep -e 'Stat10'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)/nanopolish
# ReadsFq=$(ls raw_dna/minion/*/$Strain/*.fastq.gz | grep '2017-12-03')
ReadsFq=$(ls raw_dna/minion/*/$Strain/*.fastq.gz | grep -v 'Stat10_2018-05-15')
AlignedReads=$(ls $OutDir/nanopolish/reads.sorted.bam)

NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
# python $NanoPolishDir/nanopolish_makerange.py $Assembly --segment-length 50000 > $OutDir/nanopolish_range.txt

Ploidy=1
echo "nanopolish log:" > $OutDir/nanopolish_high_mem_log.txt
ls -lh $OutDir/*/*.fa | grep -v ' 0 ' | cut -f8 -d '/' | sed 's/_consensus.fa//g' > $OutDir/files_present.txt
for Region in $(cat $OutDir/nanopolish_range.txt | grep -vwf "$OutDir/files_present.txt"); do
echo $Region
echo $Region >> $OutDir/nanopolish_high_mem_log.txt
Jobs=$(qstat | grep 'sub_nano_h' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_nano_h' | grep 'qw' | wc -l)
done		
printf "\n"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
qsub $ProgDir/sub_nanopolish_variants_high_mem.sh $Assembly $ReadsFq $AlignedReads $Ploidy $Region $OutDir/$Region
done
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta | grep -e 'Stat10'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=assembly/SMARTdenovo/$Organism/$Strain/nanopolish
mkdir -p $OutDir
# cat "" > $OutDir/"$Strain"_nanoplish.fa
NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
InDir=$(dirname $Assembly)
python $NanoPolishDir/nanopolish_merge.py $InDir/nanopolish/*/*.fa > $OutDir/"$Strain"_nanoplish.fa

echo "" > tmp.txt
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $OutDir/"$Strain"_nanoplish.fa --out $OutDir/"$Strain"_nanoplish_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
```

Quast and busco were run to assess the effects of nanopolish on assembly quality:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/nanopolish/*_nanoplish_min_500bp_renamed.fasta | grep -e 'Stat10'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
# Quast
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
# Busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
  for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary_*.txt | grep -e 'Stat10'); do
  Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
  Version=$(echo $File | rev | cut -d '/' -f2 | rev)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Single=$(cat $File | grep "(S)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Version\t$Complete\t$Single\t$Fragmented\t$Missing\t$Total"
  done
```


Remove unpacked fast5 files:

```bash
# FON129
rm -r /home/groups/harrisonlab/project_files/fusarium_ex_narcissus/home
# FON139
rm -r raw_dna/nanopolish/F.oxysporum_fsp_narcissi/FON139/home
```

### Pilon assembly correction

Assemblies were polished using Pilon

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/nanopolish/*_nanoplish_min_500bp_renamed.fasta  | grep -e 'Stat10'); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
IlluminaDir=$(ls -d qc_dna/paired/*/$Strain)
TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
OutDir=$(dirname $Assembly)/../pilon
Iterations=10
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
done
```

Contigs were renamed
```bash
echo "" > tmp.txt
Assembly=$(ls assembly/SMARTdenovo/*/*/pilon/*.fasta | grep 'Stat10' | grep 'pilon_10')
OutDir=$(dirname $Assembly)
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/pilon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
```

Quast and busco were run to assess the effects of pilon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/*.fasta | grep 'Stat10' | grep 'pilon_min_500bp_renamed.fasta'); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```


```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/*.fasta | grep 'Stat10'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```


```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/F*/*/assembly/*/short_summary_*.txt | grep 'Stat10'); do  
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

# Identifying bacterial contamination in the assembly

Submission to NCBI identified presence of a bacterial contaminant, representing contig1:

```bash
SUBID     	BioProject	BioSample	Organism
--------------------------------------------------------
SUB4837029	PRJNA506536	SAMN10462267	Fusarium oxysporum

[] We ran your sequences through our Contamination Screen. The screen found
contigs that need to be trimmed and/or excluded. Please adjust the
sequences appropriately and then resubmit your sequences. After you remove the
contamination, trim any Ns at the ends of the sequence and remove any sequences
that are shorter than 200 nt and not part of a multi-component scaffold.

Note that hits in eukaryotic genomes to mitochondrial sequences can be ignored
when specific criteria are met. Those criteria are explained below.

Note that mismatches between the name of the adaptor/primer identified in the screen
and the sequencing technology used to generate the sequencing data should not be used
to discount the validity of the screen results as the adaptors/primers of many
different sequencing platforms share sequence similarity.



Screened 68 sequences, 70,551,861 bp.
1 sequence with locations to mask/trim

Trim:
Sequence name, length, span(s), apparent source
contig_1	6847288	1111510..1111595,1111965..1112330,4865118..4867429,5917015..5917152,5935479..5935511,5945808..5945919,5954687..5954720	Bordetella avium 197N,Achromobacter xylosoxidans NH44784-1996
```

The bacterial database was used within deconseq to remove any bactrial contigs:


```bash
  for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/pilon_min_500bp_renamed.fasta | grep 'Stat10'); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    for Exclude_db in "achromobacter"; do
      AssemblyDir=$(dirname $Assembly)
      OutDir=$AssemblyDir/../deconseq_$Exclude_db
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
      qsub $ProgDir/sub_deconseq_no_retain.sh $Assembly $Exclude_db $OutDir
    done
  done
```

Results were summarised using the commands:

```bash
for Exclude_db in "achromobacter"; do
echo $Exclude_db
for File in $(ls assembly/*/*/*/*/log.txt | grep "$Exclude_db"); do
Name=$(echo $File | rev | cut -f3 -d '/' | rev);
Good=$(cat $File |cut -f2 | head -n1 | tail -n1);
Bad=$(cat $File |cut -f2 | head -n3 | tail -n1);
printf "$Name\t$Good\t$Bad\n";
done
done
```

```
achromobacter
Stat10	67	1
```


# Identifying Mitochondrial genes in assemblies

Using a blast based approach of Mt genes:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/../deconseq_achromobacter/contigs_min_500bp_filtered_renamed.fasta | grep 'Stat10'); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo $Assembly
Query=$(ls analysis/blast_homology/Mt_genes/F.spp._mt_prots_Al-Reedy_et_al._2012.fasta)
OutDir=analysis/Mt_genes/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query protein $Assembly $OutDir
done
```

Using an exclusion database with deconseq:


```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/../deconseq_achromobacter/contigs_min_500bp_filtered_renamed.fasta | grep 'Stat10'); do
    Strain=$(echo $Assembly | rev | cut -f5 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f6 -d '/' | rev)
    echo "$Organism - $Strain"
    for Exclude_db in "Fo_mtDNA"; do
      AssemblyDir=$(dirname $Assembly)
      OutDir=$AssemblyDir/../deconseq_$Exclude_db
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
      qsub $ProgDir/sub_deconseq_no_retain.sh $Assembly $Exclude_db $OutDir
    done
  done
```

Results were summarised using the commands:

```bash
for Exclude_db in "Fo_mtDNA"; do
echo $Exclude_db
for File in $(ls assembly/*/*/*/*/log.txt | grep "$Exclude_db" | grep 'Stat10'); do
Name=$(echo $File | rev | cut -f3 -d '/' | rev);
Good=$(cat $File |cut -f2 | head -n1 | tail -n1);
Bad=$(cat $File |cut -f2 | head -n3 | tail -n1);
printf "$Name\t$Good\t$Bad\n";
done
done
```

```
Fo_mtDNA
Stat10	66	1
```

Quast was run on the removed mtDNA:

```bash
for Assembly in $(ls assembly/*/*/*/deconseq_Fo_mtDNA/*_cont.fa | grep 'Stat10'); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```




# Repeat Masking

Repeat masking was performed on the non-hybrid assembly.

```bash
for Assembly in $(ls assembly/*/*/*/deconseq_Fo_mtDNA/contigs_min_500bp_filtered_renamed.fasta | grep 'Stat10'); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=repeat_masked/$Organism/"$Strain"/filtered_contigs
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
done
```

The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and hardmasked files.

```bash

for File in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep 'Stat10'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/*/*_contigs_hardmasked.fa  | grep -e 'Stat10'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```


```
repeat_masked/F.oxysporum_fsp_statice/Stat10/filtered_contigs/Stat10_contigs_softmasked_repeatmasker_TPSI_appended.fa
Number of masked bases:
11270789
```

# Gene prediction

## RNA QC

RNAseq data has been previously qc'd as part of the commands contatined in the
Fusarium repository README document.

## Alignment

Then Rnaseq data was aligned to each genome assembly:

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep 'Stat10'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
# echo "$Organism - $Strain"
# echo "$Organism - $Strain"
for RNADir in $(ls -d qc_rna/paired/F.oxysporum_fsp_cepae/* | grep -v 'control'); do
Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
FileF=$(ls $RNADir/F/*_trim.fq.gz)
FileR=$(ls $RNADir/R/*_trim.fq.gz)
# OutDir=alignment/$Organism/$Strain/$Timepoint
# Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
Jobs=$(squeue -n 'slurm_star.sh' -t 'PD' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
# Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
Jobs=$(squeue -n 'slurm_star.sh' -t 'PD' | wc -l)
done
# printf "\n"
# printf "\n"
echo "$Organism - $Strain - $Timepoint"
# echo $FileF
# echo $FileR
# Prefix=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
# OutDir=../../../../../data/scratch/armita/fusarium/alignment/star/$Organism/$Strain/$Timepoint
OutDir=../../../../../data/scratch/armita/fusarium/alignment/star2/$Organism/$Strain/$Timepoint
# ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
ProgDir=/projects/oldhome/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
# qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
sbatch $ProgDir/slurm_star.sh $Assembly $FileF $FileR $OutDir
# echo "$Strain\t$Timepoint" >> alignment.log
done
done
```

Alignments were concatenated prior to gene prediction. This step was done through
a qlogin session on the cluster.

```bash
for StrainDir in $(ls -d ../../../../../data/scratch/armita/fusarium/alignment/star/*/* | grep 'Stat10'); do
BamFiles=$(ls $StrainDir/*/star_aligmentAligned.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
OutDir=$StrainDir/concatenated
mkdir -p $OutDir
samtools merge -f $OutDir/all_reads_concatenated.bam $BamFiles
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

```bash

for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep 'Stat10'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
AcceptedHits=$(ls ../../../../../data*/scratch*/armita/fusarium/alignment/star/$Organism/$Strain/concatenated/all_reads_concatenated.bam)
OutDir=gene_pred/braker/$Organism/"$Strain"_braker
GeneModelName="$Organism"_"$Strain"_braker_new
rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker_new
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
```

Fasta and gff files were extracted from Braker1 output.

```bash
for File in $(ls gene_pred/braker/*/*_braker/*/augustus.gff3 | grep -e 'FON_63' -e 'Stocks4'); do
getAnnoFasta.pl $File
OutDir=$(dirname $File)
echo "##gff-version 3" > $OutDir/augustus_extracted.gff
cat $File | grep -v '#' >> $OutDir/augustus_extracted.gff
done
```


## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Fistly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep 'Stat10'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
mkdir -p $OutDir
AcceptedHits=$(ls ../../../../../data*/scratch*/armita/fusarium/alignment/star/$Organism/$Strain/concatenated/all_reads_concatenated.bam)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
done
```

Secondly, genes were predicted using CodingQuary:

```bash
# Jobs=$(qstat | grep 'sub_cuffli' | wc -l)
# while [ $Jobs -ge 1 ]; do
# sleep 10m
# printf "."
# Jobs=$(qstat | grep 'sub_cuffli' | wc -l)
# done
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w -e 'Stat10'); do
# Jobs=$(qstat | grep 'sub_cuffli' | wc -l)
# while [ $Jobs -ge 1 ]; do
# sleep 10
# printf "."
# Jobs=$(qstat | grep 'sub_cuffli' | wc -l)
# done
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
for BrakerGff in $(ls gene_pred/braker/F.*/*/*/augustus.gff3 | grep -e 'Stat10'); do
Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker//g')
Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs/*_softmasked_repeatmasker_TPSI_appended.fa)
CodingQuaryGff=gene_pred/codingquary/$Organism/$Strain/out/PredictedPass.gff3
PGNGff=gene_pred/codingquary/$Organism/$Strain/out/PGN_predictedPass.gff3
AddDir=gene_pred/codingquary/$Organism/$Strain/additional
FinalDir=gene_pred/final/$Organism/$Strain/final
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
# -
# This section is edited
$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $AddDir/add_genes_CodingQuary_unspliced.gff3
$ProgDir/correct_CodingQuary_splicing.py --inp_gff $AddDir/add_genes_CodingQuary_unspliced.gff3 > $FinalDir/final_genes_CodingQuary.gff3
# -
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary
cp $BrakerGff $FinalDir/final_genes_Braker.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker
cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta


GffBraker=$FinalDir/final_genes_Braker.gff3
GffQuary=$FinalDir/final_genes_CodingQuary.gff3
GffAppended=$FinalDir/final_genes_appended.gff3
cat $GffBraker $GffQuary > $GffAppended

# cat $BrakerGff $AddDir/additional_gene_parsed.gff3 | bedtools sort > $FinalGff
done
```




In preperation for submission to ncbi, gene models were renamed and duplicate gene features were identified and removed.
 * no duplicate genes were identified


The next step had problems with the masked pacbio genome. Bioperl could not read in
the fasta sequences. This was overcome by wrapping the unmasked genome and using this
fasta file.

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa | grep -e 'Stat10'); do
NewName=$(echo $Assembly | sed 's/_unmasked.fa/_unmasked_wrapped.fa/g')
cat $Assembly | fold > $NewName
done
```


Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

```bash
for GffAppended in $(ls gene_pred/final/*/*/final/final_genes_appended.gff3 | grep -e 'Stat10'); do
Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
FinalDir=gene_pred/final/$Organism/$Strain/final
GffFiltered=$FinalDir/filtered_duplicates.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
$ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered
GffRenamed=$FinalDir/final_genes_appended_renamed.gff3
LogFile=$FinalDir/final_genes_appended_renamed.log
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
$ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
rm $GffFiltered
Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_unmasked_wrapped.fa)
$ProgDir/gff2fasta.pl $Assembly $GffRenamed gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed
# The proteins fasta file contains * instead of Xs for stop codons, these should
# be changed
sed -i 's/\*/X/g' gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.pep.fasta
done
```

```bash
for Gff in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.gff3 | grep -e 'Stat10'); do
	Strain=$(echo $Gff | rev | cut -d '/' -f3 | rev)
	Organism=$(echo $Gff | rev | cut -d '/' -f4 | rev)
	echo "$Strain - $Organism"
	cat $Gff | grep -w 'gene' | wc -l
done
```

```
Stat10 - F.oxysporum_fsp_statice
26985
```

### Fungap

Gene prediction using fungap was also assessed:

This was run from the new cluster

```bash
screen -a

mkdir qc_rna/concat
rm qc_rna/concat/72hpi_CzapekDox_GlucosePeptone_PDA_PDB_RNAseq_*.fq.gz
for RNADir in $(ls -d qc_rna/paired/F.oxysporum_fsp_cepae/* | grep -v -e 'control' -e 'rep2' -e 'rep3' -e 'prelim' -e 'FO47' -e '55'); do
Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
FileF=$(ls $RNADir/F/*_trim.fq.gz)
FileR=$(ls $RNADir/R/*_trim.fq.gz)
echo $RNADir
ls $FileF
ls $FileR
cat $FileF >> qc_rna/concat/72hpi_CzapekDox_GlucosePeptone_PDA_PDB_RNAseq_F.fq.gz
cat $FileR >> qc_rna/concat/72hpi_CzapekDox_GlucosePeptone_PDA_PDB_RNAseq_R.fq.gz
done
gunzip qc_rna/concat/72hpi_CzapekDox_GlucosePeptone_PDA_PDB_RNAseq_*.fq.gz

conda activate fungap

for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa | grep -e 'Stat10'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
FileF=$(ls qc_rna/concat/*_F.fq)
FileR=$(ls qc_rna/concat/*_R.fq)
RefProteome=$(ls assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_AnnotatedProteins.fasta)
ReadF=$(ls qc_rna/concat/*_F.fq)
ReadR=$(ls qc_rna/concat/*_R.fq)
BamAlignment=$(ls ../../../../../data*/scratch*/armita/fusarium/alignment/star/$Organism/$Strain/concatenated/all_reads_concatenated.bam)
Prefix=${Organism}_${Strain}
OutDir=gene_pred/fungap/$Organism/${Strain}
ProgDir=/projects/oldhome/armita/git_repos/emr_repos/tools/gene_prediction/fungap
sbatch $ProgDir/slurm_fungap.sh $Assembly $RefProteome $ReadF $ReadR $BamAlignment $Prefix $OutDir
done

screen -a
cd /projects/oldhome/groups/harrisonlab/project_files/fusarium
mkdir tmp
cd tmp
srun --mem=40gb --cpus-per-task=40 --partition=unlimited --time=4-0:0:0 --pty bash
/home/armita/anaconda2/envs/fungap/bin/gff3_merge -g -n -d ../gene_pred/fungap/F.oxysporum_fsp_statice/Stat10/maker_out/concat/maker_run4/assembly.maker.output/assembly_master_datastore_index.log
/home/armita/anaconda2/envs/fungap/bin/fasta_merge -d ../gene_pred/fungap/F.oxysporum_fsp_statice/Stat10/maker_out/concat/maker_run4/assembly.maker.output/assembly_master_datastore_index.log

```


## Assessing the Gene space in predicted transcriptomes:

```bash
for Assembly in $(ls gene_pred/final/*/*/final/final_genes_appended_renamed.gene.fasta  | grep -e 'Stat10'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
# BuscoDB="Fungal"
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/genes
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
	for File in $(ls gene_pred/busco/*/*/genes/*/short_summary_*.txt); do  
		echo $File;
		cat $File | grep -e '(C)' -e 'Total';
	done
```


#Functional annotation

## A) Interproscan

Interproscan was used to give gene models functional annotations.
Annotation was run using the commands below:

Note: This is a long-running script. As such, these commands were run using
'screen' to allow jobs to be submitted and monitored in the background.
This allows the session to be disconnected and reconnected over time.
<!--
Screen ouput detailing the progress of submission of interporscan jobs
was redirected to a temporary output file named interproscan_submission.log .

```bash
	screen -a
	cd /home/groups/harrisonlab/project_files/fusarium
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Genes in $(ls gene_pred/final_genes/F.*/*/*/final_genes_appended_renamed.pep.fasta | grep -w -e 'FON_63' -e 'Stocks4'); do
	echo $Genes
	$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
for Proteins in $(ls gene_pred/final_genes/F.*/*/*/final_genes_appended_renamed.pep.fasta | grep -w -e 'FON_63' -e 'Stocks4'); do
Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
echo $Strain
InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
$ProgDir/append_interpro.sh $Proteins $InterProRaw
done
```


## B) SwissProt



```bash
for Proteome in $(ls gene_pred/final_genes/F.*/*/*/final_genes_appended_renamed.pep.fasta | grep -w -e 'FON_63' -e 'Stocks4'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
OutDir=gene_pred/swissprot/$Organism/$Strain
SwissDbDir=../../uniprot/swissprot
SwissDbName=uniprot_sprot
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
done
```
-->

# C) Identifying secreted proteins

Required programs:
 * SignalP-4.1
 * TMHMM

Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
CurPath=$PWD
for Proteome in $(ls gene_pred/final/F.*/*/*/final_genes_appended_renamed.pep.fasta | grep 'Stat10'); do
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
for SplitDir in $(ls -d gene_pred/final_genes_split/*/* | grep 'Stat10'); do
Strain=$(echo $SplitDir | rev |cut -d '/' -f1 | rev)
Organism=$(echo $SplitDir | rev |cut -d '/' -f2 | rev)
echo "$Organism - $Strain"
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
for Proteome in $(ls gene_pred/final/F.*/*/*/final_genes_appended_renamed.pep.fasta | grep 'Stat10'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
qsub $ProgDir/submit_TMHMM.sh $Proteome
done
```

Those proteins with transmembrane domains were removed from lists of Signal
peptide containing proteins

```bash
for File in $(ls gene_pred/trans_mem/*/*/*_TM_genes_neg.txt | grep 'Stat10'); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
TmHeaders=$(echo "$File" | sed 's/neg.txt/neg_headers.txt/g')
cat $File | cut -f1 > $TmHeaders
SigP=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp.aa)
OutDir=$(dirname $SigP)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $SigP --headers $TmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem.aa
echo "Number of SigP proteins:"
cat $SigP | grep '>' | wc -l
echo "Number without transmembrane domains:"
cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | wc -l
echo "Number of gene models:"
cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | cut -f1 -d't' | sort | uniq |wc -l
done
```

```
F.oxysporum_fsp_statice - Stat10
Number of SigP proteins:
2608
Number without transmembrane domains:
2084
Number of gene models:
2081
```

### C) Identification of MIMP-flanking genes

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa | grep 'Stat10'); do
Organism=$(echo "$Assembly" | rev | cut -d '/' -f4 | rev)
Strain=$(echo "$Assembly" | rev | cut -d '/' -f3 | rev)
GeneGff=$(ls gene_pred/final/$Organism/"$Strain"/final/final_genes_appended_renamed.gff3)
OutDir=analysis/mimps/$Organism/$Strain
mkdir -p "$OutDir"
echo "$Organism - $Strain"
ProgDir="/home/armita/git_repos/emr_repos/tools/pathogen/mimp_finder"
$ProgDir/mimp_finder.pl $Assembly $OutDir/"$Strain"_mimps.fa $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps.log
$ProgDir/gffexpander.pl +- 2000 $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps_exp.gff
echo "The number of mimps identified:"
cat $OutDir/"$Strain"_mimps.fa | grep '>' | wc -l
bedtools intersect -u -a $GeneGff -b $OutDir/"$Strain"_mimps_exp.gff > $OutDir/"$Strain"_genes_in_2kb_mimp.gff
echo "The following transcripts intersect mimps:"
MimpProtsTxt=$OutDir/"$Strain"_prots_in_2kb_mimp.txt
MimpGenesTxt=$OutDir/"$Strain"_genes_in_2kb_mimp.txt
cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' | sort | uniq > $MimpProtsTxt
cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' | cut -f1 -d '.'| sort | uniq > $MimpGenesTxt
cat $MimpProtsTxt | wc -l
cat $MimpGenesTxt | wc -l
echo ""
done
```

```
F.oxysporum_fsp_statice - Stat10
The number of mimps identified:
163
The following transcripts intersect mimps:
128
128
```

Those genes that were predicted as secreted and within 2Kb of a MIMP
were identified:

```bash
for File in $(ls analysis/mimps/*/*/*_genes_in_2kb_mimp.txt | grep 'Stat10'); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev | sed 's/_chromosomal//g')
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProtsFile=$(echo $File | sed 's/genes/prots/g')
Secretome=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp_no_trans_mem.aa)
OutFile=$(echo "$File" | sed 's/.gff/_secreted.gff/g')
SecretedHeaders=$(echo "$Secretome" | sed 's/.aa/_headers.txt/g')
cat $Secretome | grep '>' | tr -d '>' | sed 's/-p.//g' > $SecretedHeaders
SecretedMimps=$(echo "$File" | sed 's/.txt/_secreted_headers.txt/g')
cat $ProtsFile $SecretedHeaders | cut -f1 | sort | uniq -d > $SecretedMimps
cat $SecretedMimps | wc -l
cat $SecretedHeaders | cut -f1 | cut -f1 -d '.' | sort | uniq | grep -f $File > tmp.txt
cat tmp.txt | wc -l
cat $SecretedHeaders | cut -f1 | cut -f1 -d '.' | sort | uniq | grep -f $File > tmp.txt
cat $SecretedHeaders | cut -f1 | sort | uniq | grep -f $File > tmp2.txt
MimpsGff=$(ls analysis/mimps/$Organism/$Strain/*_mimps.gff)
GenesIn2Kb=$(ls analysis/mimps/$Organism/$Strain/"$Strain"_genes_in_2kb_mimp.gff)
SecretedMimpsGff=$(echo $GenesIn2Kb | sed 's/.gff/_secreted.gff/g')

ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl tmp2.txt $GenesIn2Kb secreted_mimp ID > $SecretedMimpsGff
# cat $SecretedMimpsGff | grep -w 'mRNA' | wc -l
# cat $MimpsGff | grep -w -f $SecretedMimps > $SecretedMimpsGff
done
```

```
F.oxysporum_fsp_statice - Stat10
16
16
```

sequences of these genes were extracted:

```bash
for Headers in $(ls analysis/mimps/F.oxysporum_fsp_statice/*/*_genes_in_2kb_mimp_secreted_headers.txt); do
	echo $Headers
	Organism=$(echo "$Headers" | rev | cut -d '/' -f3 | rev)
	Strain=$(echo "$Headers" | rev | cut -d '/' -f2 | rev)
	Genes=$(ls gene_pred/final/$Organism/"$Strain"/final/final_genes_appended_renamed.gene.fasta)
	HeadersMod=$(echo $Headers | sed 's/_headers.txt/_headers_mod.txt/g')
	cat $Headers | sed "s/\.t.//g" > $HeadersMod
	OutFile=$(echo $Headers | sed 's/_headers.txt/.fa/g')
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
	$ProgDir/extract_from_fasta.py --fasta $Genes --headers $HeadersMod > $OutFile
done
```


## C) CAZY proteins

Carbohydrte active enzymes were idnetified using CAZYfollowing recomendations
at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt :

```bash
for Proteome in $(ls gene_pred/final/F.*/*/*/final_genes_appended_renamed.pep.fasta | grep 'Stat10'); do
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

The Hmm parser was used to filter hits by an E-value of E1x10-5 or E1x10-e3 if they had a hit over a length of X %.

Those proteins with a signal peptide were extracted from the list and gff files
representing these proteins made.

```bash
for File in $(ls gene_pred/CAZY/*/*/*CAZY.out.dm | grep 'Stat10'); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $File)
echo "$Organism - $Strain"
ProgDir=/home/groups/harrisonlab/dbCAN
$ProgDir/hmmscan-parser.sh $OutDir/"$Strain"_CAZY.out.dm > $OutDir/"$Strain"_CAZY.out.dm.ps
CazyHeaders=$(echo $File | sed 's/.out.dm/_headers.txt/g')
cat $OutDir/"$Strain"_CAZY.out.dm.ps | cut -f3 | sort | uniq > $CazyHeaders
echo "number of CAZY proteins identified:"
cat $CazyHeaders | wc -l
# Gff=$(ls gene_pred/codingquary/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
CazyGff=$OutDir/"$Strain"_CAZY.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $CazyHeaders $Gff CAZyme ID > $CazyGff
echo "number of CAZY genes identified:"
cat $CazyGff | grep -w 'gene' | wc -l

SecretedProts=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/"$Strain"_final_sp_no_trans_mem.aa)
SecretedHeaders=$(echo $SecretedProts | sed 's/.aa/_headers.txt/g')
cat $SecretedProts | grep '>' | tr -d '>' > $SecretedHeaders
CazyGffSecreted=$OutDir/"$Strain"_CAZY_secreted.gff
$ProgDir/extract_gff_for_sigP_hits.pl $SecretedHeaders $CazyGff Secreted_CAZyme ID > $CazyGffSecreted
echo "number of Secreted CAZY proteins identified:"
cat $CazyGffSecreted | grep -w 'mRNA' | cut -f9 | tr -d 'ID=' | cut -f1 -d ';' > $OutDir/"$Strain"_CAZY_secreted_headers.txt
cat $OutDir/"$Strain"_CAZY_secreted_headers.txt | wc -l
echo "number of Secreted CAZY genes identified:"
cat $CazyGffSecreted | grep -w 'gene' | wc -l
# cat $OutDir/"$Strain"_CAZY_secreted_headers.txt | cut -f1 -d '.' | sort | uniq | wc -l
done
```

```
F.oxysporum_fsp_statice - Stat10
number of CAZY proteins identified:
1110
number of CAZY genes identified:
1110
number of Secreted CAZY proteins identified:
422
number of Secreted CAZY genes identified:
422
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


### Summary of CAZY families by organism


```bash
for CAZY in $(ls gene_pred/CAZY/*/*/*_CAZY.out.dm.ps | grep 'Stat10'); do
  Strain=$(echo $CAZY | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $CAZY | rev | cut -f3 -d '/' | rev)
  OutDir=$(dirname $CAZY)
  echo "$Organism - $Strain"
  Secreted=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp_no_trans_mem_headers.txt)
  Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/CAZY
  $ProgDir/summarise_CAZY.py --cazy $CAZY --inp_secreted $Secreted --inp_gff $Gff --summarise_family --trim_gene_id 2 --kubicek_2014
done | less -S
```

```
F.oxysporum_fsp_statice - Stat10
B-Galactosidases - 3
A-Galactosidases - 3
Polygalacturonase - 12
A-Arabinosidases - 25
Xylanases - 11
Polygalacturonate lyases - 26
B-Glucuronidases - 4
B-Glycosidases - 14
Cellulases - 21
other - 301
Xyloglucanases - 1
```

## D) AntiSMASH

Antismash was run to identify clusters of secondary metabolite genes within
the genome. Antismash was run using the weserver at:
http://antismash.secondarymetabolites.org

<!--
# Identifying genes on LS regions.


The following contigs were identified as LS in FoN and FoM:
```bash

FoN_core="contig_1 contig_2 contig_3 contig_4 contig_5 contig_6 contig_7 contig_8 contig_10 contig_11 contig_12"
Organism="F.oxysporum_fsp_narcissi"
Strain="FON_63"
echo $FoN_core | sed 's/ /\n/g' > tmp.txt
GenesIn2Kb=$(ls analysis/mimps/$Organism/$Strain/"$Strain"_genes_in_2kb_mimp.gff)
LS_MIMP_genes=$(echo $GenesIn2Kb | sed 's/.gff/_LS.gff/g')
cat $GenesIn2Kb | grep -w -v -f tmp.txt > $LS_MIMP_genes
SecretedMimpsGff=$(ls analysis/mimps/$Organism/$Strain/"$Strain"_genes_in_2kb_mimp_secreted.gff)
LS_MIMP_secreted=$(echo $GenesIn2Kb | sed 's/.gff/_LS.gff/g')
LS_MIMP_secretedHeaders=$(echo $GenesIn2Kb | sed 's/.gff/_LS_headers.txt/g')
cat $SecretedMimpsGff | grep -w -v -f tmp.txt > $LS_MIMP_secreted
cat $LS_MIMP_secreted | grep -w 'mRNA' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '='>  $LS_MIMP_secretedHeaders
cat $LS_MIMP_secreted | grep -w 'gene' | wc -l
Genes=$(ls gene_pred/final_genes/F.*/*/*/final_genes_appended_renamed.cdna.fasta | grep $Strain)
LS_MIMP_secretedFasta=$(echo $GenesIn2Kb | sed 's/.gff/_LS.fa/g')
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $LS_MIMP_secretedHeaders > $LS_MIMP_secretedFasta
cat $LS_MIMP_secretedFasta | grep '>' | wc -l

FoM_core="contig_1 contig_2 contig_3 contig_4 contig_5 contig_6 contig_7 contig_8 contig_9 contig_10 contig_12 contig_17"
Organism="F.oxysporum_fsp_mathioli"
Strain="Stocks4"
echo $FoM_core | sed 's/ /\n/g' > tmp.txt
GenesIn2Kb=$(ls analysis/mimps/$Organism/$Strain/"$Strain"_genes_in_2kb_mimp.gff)
LS_MIMP_genes=$(echo $GenesIn2Kb | sed 's/.gff/_LS.gff/g')
cat $GenesIn2Kb | grep -w -v -f tmp.txt > $LS_MIMP_genes
SecretedMimpsGff=$(ls analysis/mimps/$Organism/$Strain/"$Strain"_genes_in_2kb_mimp_secreted.gff)
LS_MIMP_secreted=$(echo $GenesIn2Kb | sed 's/.gff/_LS.gff/g')
LS_MIMP_secretedHeaders=$(echo $GenesIn2Kb | sed 's/.gff/_LS_headers.txt/g')
cat $SecretedMimpsGff | grep -w -v -f tmp.txt > $LS_MIMP_secreted
cat $LS_MIMP_secreted | grep -w 'mRNA' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '='>  $LS_MIMP_secretedHeaders
cat $LS_MIMP_secreted | grep -w 'gene' | wc -l
Genes=$(ls gene_pred/final_genes/F.*/*/*/final_genes_appended_renamed.cdna.fasta | grep $Strain)
LS_MIMP_secretedFasta=$(echo $GenesIn2Kb | sed 's/.gff/_LS.fa/g')
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $LS_MIMP_secretedHeaders > $LS_MIMP_secretedFasta
cat $LS_MIMP_secretedFasta | grep '>' | wc -l

```

 -->















# 5. Genomic analysis


## 5.1  Identifying SIX gene homologs

## 5.1.a) Performing BLAST searches

BLast searches were performed against the genome:

```bash
  for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep 'Stat10'); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    Query=analysis/blast_homology/six_genes/six-appended_parsed.fa
    qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
  done
```

```bash
cat analysis/blast_homology/F.oxysporum_fsp_statice/Stat10/Stat10_six-appended_parsed.fa_homologs.csv | cut -f1,5,6
```

```
ID	No.hits	Hit
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six10_(SIX10)_mRNA,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six11_(SIX11)_mRNA,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six12_(SIX12)_mRNA,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six13_(SIX13)_mRNA,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six14_(SIX14)_mRNA,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_1_(SIX1)_gene,_complete_cds	1	contig_20
Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_2_(SIX2)_gene,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_isolate_FOL-MM10_secreted_in_xylem_3_(SIX3)_gene,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_isolate_IPO3_secreted_in_xylem_3_(SIX3)_gene,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_isolate_14844_secreted_in_xylem_3_(SIX3)_gene,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_3_(SIX3)_gene,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_SIX3_gene_for_Secreted_in_xylem_3_protein	0
Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_4_(SIX4)_gene,_complete_cds	1	contig_25
Fusarium_oxysporum_f._sp._lycopersici_Six5_mRNA,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_5_(SIX5)_gene,_partial_cds	0
Fusarium_oxysporum_f._sp._passiflorae_SIX6_gene,_complete_cds	1	contig_23
Fusarium_oxysporum_f._sp._niveum_SIX6_gene,_complete_cds	1	contig_23
Fusarium_oxysporum_f._sp._vasinfectum_SIX6_gene,_complete_cds	1	contig_23
Fusarium_oxysporum_f._sp._lycopersici_secreted_in_xylem_Six6_(SIX6)_mRNA,_complete_cds	1	contig_23
Fusarium_oxysporum_f._sp._melonis_isolate_NRRL_26406_secreted_in_xylem_6-like_protein_(SIX6)_gene,_complete_cds	1	contig_23
Fusarium_oxysporum_f._sp._radicis-cucumerinum_isolate_Afu-3_secreted_in_xylem_6-like_protein_(SIX6)_gene,_complete_cds	1	contig_23
Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_6_(SIX6)_gene,_complete_cds	1	contig_23
Fusarium_oxysporum_f._sp._lycopersici_secreted_in_xylem_Six7_(SIX7)_mRNA,_complete_cds	0
Fusarium_oxysporum_f._sp._lilii_isolate_NRRL_28395_secreted_in_xylem_7-like_protein_(SIX7)_gene,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_7_(SIX7)_gene,_complete_cds	0
Fusarium_oxysporum_SIX8_gene,_complete_cds	1	contig_25
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six9_(SIX9)_mRNA,_complete_cds	1	contig_24

```

## 5.1.b) Converting BLAST results to gff annotations

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
for BlastHits in $(ls analysis/blast_homology/*/*/*_six-appended_parsed.fa_homologs.csv | grep 'Stat10'); do
Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)  
Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
HitsGff=analysis/blast_homology/$Organism/$Strain/"$Strain"_six-appended_parsed.fa_homologs.gff
Column2=BLAST_hit
NumHits=3
$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
done
```
