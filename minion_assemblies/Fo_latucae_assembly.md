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
pip3 install --user ont_albacore-2.2.7-cp34-cp34m-manylinux1_x86_64.whl --upgrade
~/.local/bin/read_fast5_basecaller.py --version

mkdir FoLatucae_26-04-18
cd FoLatucae_26-04-18

# Oxford nanopore AJ516 Run 1
Organism=F.oxysporum_fsp_lactucae
Strain=AJ516
Date=22-02-18
FlowCell="FLO-MIN106"
Kit="SQK-LSK108"
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

cat $OutDir/albacore_v2.2.7/workspace/pass/*.fastq | gzip -cf > ${Strain}_${Date}_albacore_v2.2.7.fastq.gz
tar -cz -f ${Strain}_${Date}_albacore_v2.2.7.tar.gz $OutDir

FinalDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.2.7
mkdir -p $FinalDir
mv $OutDir/${Strain}_${Date}_albacore_v2.2.7.* $FinalDir/.
chmod +rw -R $FinalDir


	# OutDir=/data/scratch/nanopore_tmp_data/Fusarium/albacore_v2.2.7
  # mkdir -p $OutDir
  # chmod +w $OutDir
  # cp Alt_albacore_v2.10_barcode0*.fastq.gz $OutDir/.
  # chmod +rw $OutDir/Alt_albacore_v2.10_barcode0*.fastq.gz
  # scp Alt_albacore_v2.10_barcode*.fastq.gz armita@192.168.1.200:$OutDir/.
  # tar -cz -f Alt_albacore_v2.10_demultiplexed.tar.gz Alt_albacore_v2.10_demultiplexed
  # OutDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.2.7
  # mv Alt_albacore_v2.10_demultiplexed.tar.gz $OutDir/.
  # chmod +rw $OutDir/Alt_albacore_v2.10_demultiplexed.tar.gz

	# Oxford nanopore AJ516 Run 2
	Organism=F.oxysporum_fsp_lactucae
	Strain=AJ516
	Date=18-04-18
	FlowCell="FLO-MIN106"
	Kit="SQK-LSK108"
	RawDatDir=/data/seq_data/minion/2018/20180418_AJ16/AJ16/GA10000/reads
	OutDir=~/FoLatucae_26-04-18/$Organism/$Strain/$Date
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

	cd ~/FoLatucae_26-04-18
	cat $OutDir/albacore_v2.2.7/workspace/pass/*.fastq | gzip -cf > ${Strain}_${Date}_albacore_v2.2.7.fastq.gz
	tar -cz -f ${Strain}_${Date}_albacore_v2.2.7.tar.gz $OutDir

	FinalDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.2.7
	mkdir -p $FinalDir
	mv ${Strain}_${Date}_albacore_v2.2.7.* $FinalDir/.
	chmod +rw -R $FinalDir

	# Oxford nanopore AJ520 Run 1
	Organism=F.oxysporum_fsp_lactucae
	Strain=AJ520
	Date=22-02-18
	FlowCell="FLO-MIN106"
	Kit="SQK-RBK001"
	RawDatDir=/data/seq_data/minion/2018/20180222_AJ520/AJ520/GA20000/reads
	OutDir=~/FoLatucae_26-04-18/$Organism/$Strain/$Date
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

	cd ~/FoLatucae_26-04-18
	cat $OutDir/albacore_v2.2.7/workspace/pass/unclassified/*.fastq | gzip -cf > ${Strain}_${Date}_albacore_v2.2.7.fastq.gz

	tar -cz -f ${Strain}_${Date}_albacore_v2.2.7.tar.gz $OutDir

	FinalDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.2.7
	mkdir -p $FinalDir
	mv ${Strain}_${Date}_albacore_v2.2.7.* $FinalDir/.
	chmod +rw -R $FinalDir

	# Oxford nanopore AJ520 Run 2
	Organism=F.oxysporum_fsp_lactucae
	Strain=AJ520
	Date=18-04-18
	FlowCell="FLO-MIN106"
	Kit="SQK-LSK108"
	# The second run of AJ520 reads have not yet been copied over.
	RawDatDir=/data/seq_data/minion/2018/20180426_AJ520_GA30000/GA30000/reads
	OutDir=~/FoLatucae_26-04-18/$Organism/$Strain/$Date
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

		cd ~/FoLatucae_26-04-18
		cat $OutDir/albacore_v2.2.7/workspace/pass/*.fastq | gzip -cf > ${Strain}_${Date}_albacore_v2.2.7.fastq.gz
		tar -cz -f ${Strain}_${Date}_albacore_v2.2.7.tar.gz $OutDir

		FinalDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.2.7
		mkdir -p $FinalDir
		mv ${Strain}_${Date}_albacore_v2.2.7.* $FinalDir/.
		chmod +rw -R $FinalDir
```

Sequence data was moved into the appropriate directories

```bash
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	cd $ProjectDir

	Species="F.oxysporum_fsp_lactucae"
	Strain="AJ516"
	OutDir=$ProjectDir/raw_dna/minion/$Species/$Strain
	mkdir -p $OutDir
	RawDatDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.2.7
	cp -s $RawDatDir/AJ516_18-04-18_albacore_v2.2.7.fastq.gz $OutDir/.
	cp -s $RawDatDir/AJ516_22-02-18_albacore_v2.2.7.fastq.gz $OutDir/.

	Species="F.oxysporum_fsp_lactucae"
	Strain="AJ520"
	OutDir=$ProjectDir/raw_dna/minion/$Species/$Strain
	mkdir -p $OutDir
	RawDatDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.2.7
	cp -s $RawDatDir/AJ520_18-04-18_albacore_v2.2.7.fastq.gz $OutDir/.
	cp -s $RawDatDir/AJ520_22-02-18_albacore_v2.2.7.fastq.gz $OutDir/.
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


## Assembly

### Removal of adapters

Splitting reads and trimming adapters using porechop
```bash
	for RawReads in $(ls raw_dna/minion/*/*/*.fastq.gz | grep 'fsp_lactucae' | grep 'AJ520'); do
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

For Minion data:
```bash
for RawData in $(ls qc_dna/minion/*/*/*q.gz | grep 'lactucae' | grep 'AJ520'); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
GenomeSz=60
OutDir=$(dirname $RawData)
qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
done
```

```bash
  for StrainDir in $(ls -d qc_dna/minion/*/* | grep 'lactucae'| grep 'AJ520'); do
    Strain=$(basename $StrainDir)
    printf "$Strain\t"
    for File in $(ls $StrainDir/*cov.txt); do
      echo $(basename $File);
      cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
    done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
  done
```
MinION coverage was:
```
AJ516   57.57
AJ520	81.4
```

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


### Read correction using Canu

```bash
for ReadDir in $(ls -d qc_dna/minion/*/* | grep 'lactucae'); do
	Strain=$(echo $ReadDir | rev | cut -f1 -d '/' | rev)
	cat $ReadDir/*_trim.fastq.gz > $ReadDir/${Strain}_trim_appended.fastq.gz
done


for TrimReads in $(ls qc_dna/minion/*/*/*_trim_appended.fastq.gz | grep 'AJ516'); do
Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
OutDir=assembly/canu-1.6/$Organism/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
qsub $ProgDir/sub_canu_correction.sh $TrimReads 58m $Strain $OutDir
done

# For AJ520
TrimReads1=$(ls qc_dna/minion/*/*/*_trim.fastq.gz| grep 'AJ520' | head -n1 | tail -n1)
TrimReads2=$(ls qc_dna/minion/*/*/*_trim.fastq.gz| grep 'AJ520' | head -n2 | tail -n1)
Organism=$(echo $TrimReads1 | rev | cut -f3 -d '/' | rev)
Strain=$(echo $TrimReads1 | rev | cut -f2 -d '/' | rev)
echo $TrimReads
OutDir=assembly/canu-1.6/$Organism/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
qsub $ProgDir/sub_canu_correction_2lib.sh $TrimReads1 $TrimReads2 58m $Strain $OutDir
```


### Assembly using SMARTdenovo

```bash
for CorrectedReads in $(ls assembly/canu-1.6/*/*/*.trimmedReads.fasta.gz | grep 'F.oxysporum_fsp_lactucae'| grep 'AJ520'); do
Organism=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $CorrectedReads | rev | cut -f2 -d '/' | rev)
Prefix="$Strain"_smartdenovo
OutDir=assembly/SMARTdenovo/$Organism/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/SMARTdenovo
qsub $ProgDir/sub_SMARTdenovo.sh $CorrectedReads $Prefix $OutDir
# $ProgDir/sub_SMARTdenovo.sh $CorrectedReads $Prefix $OutDir
done
```


Quast and busco were run to assess the effects of racon on assembly quality:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg | grep 'F.oxysporum_fsp_lactucae' | grep 'AJ520'); do
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
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg | grep 'F.oxysporum_fsp_lactucae' | grep 'AJ520'); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ReadsFq=$(ls qc_dna/minion/*/$Strain/*_trim_appended.fastq.gz)
Iterations=10
OutDir=$(dirname $Assembly)"/racon_$Iterations"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/racon
qsub $ProgDir/sub_racon.sh $Assembly $ReadsFq $Iterations $OutDir
done
```

<!-- ```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg | grep 'F.oxysporum_fsp_lactucae' | grep 'AJ516'); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ReadsFq=$(ls qc_dna/minion/*/$Strain/*q.gz)
Iterations=10
OutDir=$(dirname $Assembly)"/racon_$Iterations"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/racon
qsub $ProgDir/sub_racon_2libs.sh $Assembly $ReadsFq $Iterations $OutDir
done
``` -->

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/*.fasta | grep 'round_10' | grep 'F.oxysporum_fsp_lactucae' | grep 'AJ520'); do
OutDir=$(dirname $Assembly)
echo "" > tmp.txt
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/racon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
```

Quast and busco were run to assess the effects of racon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta | grep 'F.oxysporum_fsp_lactucae' | grep 'AJ520'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/F.*/*/racon*/*.fasta | grep 'F.oxysporum_fsp_lactucae' | grep 'AJ520'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```
```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary_*.txt  | grep 'F.oxysporum_fsp_lactucae'); do
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

# Nanopolish


For AJ516

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta | grep -e 'AJ516'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Note - the full path from home must be used
ReadDir=raw_dna/nanopolish/$Organism/$Strain
mkdir -p $ReadDir
ReadsFq=$(ls raw_dna/minion/*/$Strain/*.fastq.gz | grep '18-04-18')
# cat $ReadsFq > $ReadDir/${Strain}_appended.fastq.gz
CurDir=$PWD
TmpDir=/data/scratch/armita/FoL/$Organism/$Strain
mkdir -p $TmpDir
ScratchDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.2.7
# tar -zxvf $ScratchDir/AJ516_18-04-18_albacore_v2.2.7.tar.gz -C $TmpDir
# tar -zxvf $ScratchDir/AJ516_22-02-18_albacore_v2.2.7.tar.gz -C $TmpDir
Fast5Dir1=$(ls -d /data/scratch/armita/FoL/F.oxysporum_fsp_lactucae/AJ516/home/nanopore/FoLatucae_26-04-18/F.oxysporum_fsp_lactucae/AJ516/18-04-18/albacore_v2.2.7/workspace/pass)
# Fast5Dir2=$(ls -d /data/scratch/armita/FoL/F.oxysporum_fsp_lactucae/AJ516/home/nanopore/FoLatucae_26-04-18/F.oxysporum_fsp_lactucae/AJ516/22-02-18/albacore_v2.2.7/workspace/pass)
nanopolish index -d $Fast5Dir1 $ReadsFq
done

for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta | grep -e 'AJ516'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Note - the full path from home must be used
# ReadDir=raw_dna/nanopolish/$Organism/$Strain
# mkdir -p $ReadDir
ReadsFq=$(ls raw_dna/minion/*/$Strain/*.fastq.gz | grep '18-04-18')
OutDir=$(dirname $Assembly)/nanopolish
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
# submit alignments for nanoppolish
qsub $ProgDir/sub_minimap2_nanopolish.sh $Assembly $ReadsFq $OutDir/nanopolish
done
```


For AJ520

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta | grep -e 'AJ520'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Note - the full path from home must be used
ReadDir=raw_dna/nanopolish/$Organism/$Strain
mkdir -p $ReadDir
ReadsFq=$(ls raw_dna/minion/*/$Strain/*.fastq.gz | grep '18-04-18')
# cat $ReadsFq > $ReadDir/${Strain}_appended.fastq.gz
CurDir=$PWD
TmpDir=/data/scratch/armita/FoL/$Organism/$Strain
mkdir -p $TmpDir
ScratchDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.2.7
# tar -zxvf $ScratchDir/AJ520_18-04-18_albacore_v2.2.7.tar.gz -C $TmpDir
Fast5Dir1=$(ls -d /data/scratch/armita/FoL/F.oxysporum_fsp_lactucae/AJ520/home/nanopore/FoLatucae_26-04-18/F.oxysporum_fsp_lactucae/AJ520/18-04-18/albacore_v2.2.7/workspace/pass/)
nanopolish index -d $Fast5Dir1 $ReadsFq
done

for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta | grep -e 'AJ520'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Note - the full path from home must be used
# ReadDir=raw_dna/nanopolish/$Organism/$Strain
# mkdir -p $ReadDir
ReadsFq=$(ls raw_dna/minion/*/$Strain/*.fastq.gz | grep '18-04-18')
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
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta | grep -e 'AJ516' -e 'AJ520' | grep -e 'AJ520'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)/nanopolish
# ReadsFq=$(ls raw_dna/minion/*/$Strain/*.fastq.gz | grep '2017-12-03')
ReadsFq=$(ls raw_dna/minion/*/$Strain/*.fastq.gz | grep -v '22-02-18')
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
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta | grep -e 'AJ516' -e 'AJ520' | grep -e 'AJ520'); do
	Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
	Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
	echo "$Organism - $Strain"
	OutDir=$(dirname $Assembly)/nanopolish
	# ReadsFq=$(ls raw_dna/minion/*/$Strain/*.fastq.gz | grep '2017-12-03')
	ReadsFq=$(ls raw_dna/minion/*/$Strain/*.fastq.gz)
	AlignedReads=$(ls $OutDir/nanopolish/reads.sorted.bam)

	NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
	python $NanoPolishDir/nanopolish_makerange.py $Assembly --segment-length 50000 > $OutDir/nanopolish_range.txt

	Ploidy=1
	echo "nanopolish log:" > $OutDir/nanopolish_high_mem_log.txt
	ls -lh $OutDir/*/*.fa | grep -v ' 0 ' | cut -f8 -d '/' | sed 's/_consensus.fa//g' > $OutDir/files_present.txt
	for Region in $(cat $OutDir/nanopolish_range.txt | grep -vwf "$OutDir/files_present.txt" | head -n1); do
	echo $Region
	echo $Region >> $OutDir/nanopolish_high_mem_log.txt
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
	qsub $ProgDir/sub_nanopolish_variants_high_mem.sh $Assembly $ReadsFq $AlignedReads $Ploidy $Region $OutDir/$Region
	done
	done
```

```bash
for File in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_lactucae/AJ520/racon_10/nanopolish/*/*.fa); do
# echo $File;
Segment=$(echo $File | cut -f7 -d '/');
# echo
Hits=$(cat assembly/SMARTdenovo/F.oxysporum_fsp_lactucae/AJ520/racon_10/nanopolish/nanopolish_range.txt | grep -w "$Segment" | wc -l)
# echo "$Segment - $Hits"
if [[ $Hits < 1 ]]; then
	echo "$Segment - $Hits"
	rm -r assembly/SMARTdenovo/F.oxysporum_fsp_lactucae/AJ520/racon_10/nanopolish/$Segment
fi
done | less
```


```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta | grep -e 'AJ516' -e 'AJ520' | grep -e 'AJ520'); do
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
for Assembly in $(ls assembly/SMARTdenovo/*/*/nanopolish/*_nanoplish_min_500bp_renamed.fasta | grep -e 'AJ516' -e 'AJ520' | grep 'AJ520'); do
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
  for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary_*.txt); do
  Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Single=$(cat $File | grep "(S)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Complete\t$Single\t$Fragmented\t$Missing\t$Total"
  done
```


### Pilon assembly correction

Assemblies were polished using Pilon

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/nanopolish/*_nanoplish_min_500bp_renamed.fasta | grep -e 'AJ516' -e 'AJ520' | grep 'AJ520'); do
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
for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/*.fasta | grep 'pilon_10'| grep -e 'AJ516' -e 'AJ520' | grep 'AJ520'); do
OutDir=$(dirname $Assembly)
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/pilon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
```

Quast and busco were run to assess the effects of pilon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/pilon_min_500bp_renamed.fasta | grep -e 'AJ516' -e 'AJ520' | grep 'AJ520'); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```


```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/*.fasta | grep 'pilon_min_500bp_renamed.fasta' | grep -e 'AJ516' -e 'AJ520' | grep 'AJ520'); do
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
for File in $(ls gene_pred/busco/F*/*/assembly/*/short_summary_*.txt | grep -e 'AJ516' -e 'AJ520'); do
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

```
short_summary_AJ516_smartdenovo.dmo.lay.txt	990	6	694	2041	3725
short_summary_AJ516_smartdenovo_racon_round_1.txt	2678	24	540	507	3725
short_summary_AJ516_smartdenovo_racon_round_2.txt	2757	22	519	449	3725
short_summary_AJ516_smartdenovo_racon_round_3.txt	2775	25	508	442	3725
short_summary_AJ516_smartdenovo_racon_round_4.txt	2817	26	483	425	3725
short_summary_AJ516_smartdenovo_racon_round_5.txt	2812	22	499	414	3725
short_summary_AJ516_smartdenovo_racon_round_6.txt	2819	25	502	404	3725
short_summary_AJ516_smartdenovo_racon_round_7.txt	2796	28	513	416	3725
short_summary_AJ516_smartdenovo_racon_round_8.txt	2800	23	501	424	3725
short_summary_AJ516_smartdenovo_racon_round_9.txt	2803	24	502	420	3725
short_summary_AJ516_smartdenovo_racon_round_10.txt	2828	22	487	410	3725
short_summary_racon_min_500bp_renamed.txt	2828	22	487	410	3725
short_summary_AJ516_nanoplish_min_500bp_renamed.txt	3345	36	183	197	3725

short_summary_AJ520_smartdenovo.dmo.lay.txt	62	0	94	3569	3725
short_summary_AJ520_smartdenovo_racon_round_1.txt	1274	5	321	2130	3725
short_summary_AJ520_smartdenovo_racon_round_2.txt	1386	4	285	2054	3725
short_summary_AJ520_smartdenovo_racon_round_3.txt	1404	7	287	2034	3725
short_summary_AJ520_smartdenovo_racon_round_4.txt	1407	6	268	2050	3725
short_summary_AJ520_smartdenovo_racon_round_5.txt	1430	7	259	2036	3725
short_summary_AJ520_smartdenovo_racon_round_6.txt	1423	8	257	2045	3725
short_summary_AJ520_smartdenovo_racon_round_7.txt	1441	5	240	2044	3725
short_summary_AJ520_smartdenovo_racon_round_8.txt	1412	7	273	2040	3725
short_summary_AJ520_smartdenovo_racon_round_9.txt	1421	8	271	2033	3725
short_summary_AJ520_smartdenovo_racon_round_10.txt	1426	7	265	2034	3725
short_summary_racon_min_500bp_renamed.txt	1426	7	265	2034	3725
short_summary_AJ520_nanoplish_min_500bp_renamed.txt	1731	8	97	1897	3725

```



## Identifying Mitochondrial genes in assemblies
<!--
Using a blast based approach of Mt genes:

```bash
for Assembly in $(ls assembly/merged_SMARTdenovo_spades/*/1166/pilon/pilon_min_500bp_renamed.fasta assembly/SMARTdenovo/A.gaisen/650/pilon/pilon_min_500bp_renamed.fasta); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo $Assembly
Query=$(ls ../../../../home/groups/harrisonlab/project_files/alternaria/analysis/blast_homology/Mt_dna/Mt_genes_A.alternata_Liao_2017.fasta)
OutDir=analysis/Mt_genes/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query protein $Assembly $OutDir
done
```

Using an exclusion database with deconseq:


```bash
  for Assembly in $(ls assembly/merged_SMARTdenovo_spades/*/1166/pilon/pilon_min_500bp_renamed.fasta assembly/SMARTdenovo/A.gaisen/650/pilon/pilon_min_500bp_renamed.fasta); do
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    for Exclude_db in "Aalt_mtDNA"; do
      AssemblyDir=$(dirname $Assembly)
      OutDir=$AssemblyDir/../deconseq_$Exclude_db
      ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
      qsub $ProgDir/sub_deconseq_no_retain.sh $Assembly $Exclude_db $OutDir
    done
  done
```

Results were summarised using the commands:

```bash
for Exclude_db in "Aalt_mtDNA"; do
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
Aalt_mtDNA
1166	22	1
650	23	1
650	27	1
```

Quast was run on the removed mtDNA:

```bash
for Assembly in $(ls assembly/*/A*/*/deconseq_Aalt_mtDNA/*_cont.fa); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```
 -->


# Repeat Masking

Repeat masking was performed on the non-hybrid assembly.

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/pilon_min_500bp_renamed.fasta | grep -e 'AJ516' -e 'AJ520' | grep 'AJ520'); do
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

for File in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep -e 'AJ516' -e 'AJ520'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/*/*_contigs_hardmasked.fa  | grep -e 'AJ516' -e 'AJ520'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```


# Gene prediction

## RNA QC

RNAseq data has been previously qc'd as part of the commands contatined in the
Fusarium repository README document.

## Alignment

Then Rnaseq data was aligned to each genome assembly:

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep 'AJ520'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
# echo "$Organism - $Strain"
# echo "$Organism - $Strain"
for RNADir in $(ls -d qc_rna/paired/F.oxysporum_fsp_cepae/* | grep -v 'control'); do
Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
FileF=$(ls $RNADir/F/*_trim.fq.gz)
FileR=$(ls $RNADir/R/*_trim.fq.gz)
# OutDir=alignment/$Organism/$Strain/$Timepoint
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
done
# printf "\n"
# printf "\n"
echo "$Organism - $Strain - $Timepoint"
# echo $FileF
# echo $FileR
# Prefix=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
OutDir=../../../../../data/scratch/armita/fusarium/alignment/star/$Organism/$Strain/$Timepoint
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
echo "$Strain\t$Timepoint" >> alignment.log
done
done
```

Alignments were concatenated prior to gene prediction. This step was done through
a qlogin session on the cluster.

```bash
	for StrainDir in $(ls -d ../../../../../data/scratch/armita/fusarium/alignment/star/*/* | grep 'AJ520'); do
	BamFiles=$(ls $StrainDir/*/star_aligmentAligned.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
	OutDir=$StrainDir/concatenated
	mkdir -p $OutDir
	samtools merge -f $OutDir/all_reads_concatenated.bam $BamFiles
	done
	for StrainDir in $(ls -d ../../../../../data2/scratch2/armita/fusarium/alignment/star/*/* | grep 'AJ516'); do
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

<!-- ```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -e 'AJ516' -e 'AJ520' | grep 'AJ520'); do
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
``` -->


Gene prediction was run on the new cluster from a screen session with a ssh connection to a worker node:

```bash
screen -a
ssh compute03

ProjDir=/oldhpc/home/groups/harrisonlab/project_files/fusarium

for Assembly in $(ls $ProjDir/repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -e 'AJ516' -e 'AJ520'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"


WorkDir=$HOME/tmp/braker_Fus_${Strain}
OutDir=$ProjDir/gene_pred/braker/$Organism/${Strain}_UTR
AcceptedHits=$(ls /oldhpc/data/scratch/armita/fusarium/alignment/star/*/${Strain}/concatenated/all_reads_concatenated.bam)
GeneModelName="$Organism"_"$Strain"_braker
CurDir=$PWD

mkdir -p $WorkDir
cd $WorkDir

cp $Assembly assembly.fa
cp $AcceptedHits alignedRNA.bam

braker.pl \
  --cores 40 \
  --overwrite \
  --fungus \
  --UTR=on \
  --gff3 \
  --softmasking on \
  --species=$GeneModelName \
  --genome="assembly.fa" \
  --bam="alignedRNA.bam"

mkdir -p $OutDir
cp -r braker/* $OutDir/.
done

# rm -r $WorkDir

```

Fasta and gff files were extracted from Braker1 output.

```bash
for File in $(ls gene_pred/braker/*/*_UTR/*/augustus.hints.gff3 | grep -e 'AJ516' -e 'AJ520' ); do
	Strain=$(echo $File | rev | cut -d '/' -f3 | rev | sed 's/_UTR//g')
	Organism=$(echo $File | rev | cut -d '/' -f4 | rev)
	echo "$Organism - $Strain"
	echo "number of genes:"
	cat $File | grep -v '#' | grep -w 'gene' | wc -l
	echo "number of genes with predicted UTRs"
	cat ${File%.gff3}_utr.gff | grep -v '#' | grep -w 'gene' | wc -l

	getAnnoFasta.pl $File
	OutDir=$(dirname $File)
	echo "##gff-version 3" > $OutDir/augustus_extracted.gff
	cat $File | grep -v '#' >> $OutDir/augustus_extracted.gff
done
```

```
	F.oxysporum_fsp_lactucae - AJ516
	number of genes:
	17374
	number of genes with predicted UTRs
	16260
	F.oxysporum_fsp_lactucae - AJ520
	number of genes:
	16603
	number of genes with predicted UTRs
	15735
```


## Supplementing Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Fistly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -e 'AJ516' -e 'AJ520'); do
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
Jobs=$(qstat | grep 'sub_cuffli' | wc -l)
while [ $Jobs -ge 1 ]; do
sleep 10m
printf "."
Jobs=$(qstat | grep 'sub_cuffli' | wc -l)
done
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -e 'AJ516' -e 'AJ520'); do
Jobs=$(qstat | grep 'sub_cuffli' | wc -l)
while [ $Jobs -ge 1 ]; do
sleep 10
printf "."
Jobs=$(qstat | grep 'sub_cuffli' | wc -l)
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
# for BrakerGff in $(ls gene_pred/braker/F.*/*/*/augustus.gff3 | grep -e 'AJ516' -e 'AJ520'); do
# Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker//g')
# Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
for BrakerGff in $(ls gene_pred/braker/F.*/*/*/augustus.hints.gff3 | grep -e 'AJ516' -e 'AJ520' | grep 'AJ516'); do
Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_UTR//g')
Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs/*_softmasked_repeatmasker_TPSI_appended.fa)
CodingQuaryGff=gene_pred/codingquary/$Organism/$Strain/out/PredictedPass.gff3
PGNGff=gene_pred/codingquary/$Organism/$Strain/out/PGN_predictedPass.gff3
AddDir=gene_pred/codingquary/$Organism/$Strain/additional
# FinalDir=gene_pred/final/$Organism/$Strain/final
FinalDir=gene_pred/final/$Organism/$Strain/publication
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

```
F.oxysporum_fsp_lactucae - AJ516
Genome fasta parsed
NS_02920: upstream edited -
NS_05731: upstream edited -
CUFF_11431_1_62: upstream edited -
NS_10610: upstream edited +
NS_12018: upstream edited +
NS_12300: upstream edited -
NS_12505: upstream edited +
NS_13067: upstream edited -
NS_13068: upstream edited +
NS_13138: upstream edited +
CUFF_12370_1_28: upstream edited -
NS_13176: upstream edited +
NS_13177: upstream edited +
NS_13527: upstream edited -
NS_13642: upstream edited +
NS_13713: upstream edited -
NS_13951: upstream edited -
NS_13953: upstream edited -
NS_13954: upstream edited -
NS_14106: upstream edited +
NS_14274: upstream edited -
NS_14345: upstream edited -
NS_14346: upstream edited -
Genome fasta parsed
g1: upstream edited +
g3102: upstream edited -
g6383: upstream edited -
g6431: upstream edited -
g8760: upstream edited +
g11547: upstream edited -
g11548: upstream edited -
g11773: upstream edited +
g11907: upstream edited -
g13216: upstream edited +
g13719: upstream edited -
g14902: upstream edited +
g16581: upstream edited -
g17306: upstream edited +
g17323: upstream edited +

F.oxysporum_fsp_lactucae - AJ520
Genome fasta parsed
NS_00830: upstream edited -
NS_05752: upstream edited +
PGN_16225: upstream edited -
NS_07976: upstream edited -
NS_07977: upstream edited -
NS_08169: upstream edited -
NS_08170: upstream edited -
NS_08171: upstream edited -
NS_09058: upstream edited +
NS_09059: upstream edited +
NS_09232: upstream edited +
NS_09804: upstream edited -
CUFF_11773_1_12: upstream edited -
NS_10138: upstream edited -
NS_10516: upstream edited +
NS_12584: upstream edited +
NS_12585: upstream edited +
NS_10688: upstream edited +
NS_11166: upstream edited +
NS_11195: upstream edited -
NS_11380: upstream edited -
NS_11381: upstream edited -
NS_11859: upstream edited -
NS_11932: upstream edited +
NS_11963: upstream edited -
NS_12399: upstream edited -
NS_12432: upstream edited +
NS_12470: upstream edited +
NS_12520: upstream edited -
NS_12562: upstream edited +
NS_12567: upstream edited +
NS_12568: upstream edited +
NS_12812: upstream edited +
NS_12813: upstream edited +
Genome fasta parsed
g1518: upstream edited -
g4427: upstream edited -
g6830: upstream edited -
g6860: upstream edited -
g9912: upstream edited -
g9913: upstream edited +
g11071: upstream edited -
g12759: upstream edited -
g14278: upstream edited -
g15167: upstream edited -
g15835: upstream edited -
g16280: upstream edited +
g16332: upstream edited +
g16603: upstream edited -
```


In preperation for submission to ncbi, gene models were renamed and duplicate gene features were identified and removed.
 * no duplicate genes were identified


The next step had problems with the masked pacbio genome. Bioperl could not read in
the fasta sequences. This was overcome by wrapping the unmasked genome and using this
fasta file.

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa | grep -e 'AJ516' -e 'AJ520'); do
NewName=$(echo $Assembly | sed 's/_unmasked.fa/_unmasked_wrapped.fa/g')
cat $Assembly | fold > $NewName
done
```


Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

```bash
# for GffAppended in $(ls gene_pred/final/*/*/final/final_genes_appended.gff3 | grep -e 'AJ516' -e 'AJ520'); do
# Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
# Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
# echo "$Organism - $Strain"
# FinalDir=gene_pred/final/$Organism/$Strain/final
for GffAppended in $(ls gene_pred/final/*/*/publication/final_genes_appended.gff3 | grep -e 'AJ516' -e 'AJ520' | grep 'AJ516'); do
Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev | sed 's/_UTR//g')
Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
FinalDir=gene_pred/final/$Organism/$Strain/publication
GffFiltered=$FinalDir/filtered_duplicates.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
$ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered
GffRenamed=$FinalDir/final_genes_appended_renamed.gff3
LogFile=$FinalDir/final_genes_appended_renamed.log
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
$ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
rm $GffFiltered
Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_unmasked_wrapped.fa)
$ProgDir/gff2fasta.pl $Assembly $GffRenamed $FinalDir/final_genes_appended_renamed
# The proteins fasta file contains * instead of Xs for stop codons, these should
# be changed
sed -i 's/\*/X/g' $FinalDir/final_genes_appended_renamed.pep.fasta
done
```

```bash
for Gff in $(ls gene_pred/final/*/*/publication/final_genes_appended_renamed.gff3  | grep -e 'AJ516' -e 'AJ520'); do
	Strain=$(echo $Gff | rev | cut -d '/' -f3 | rev)
	Organism=$(echo $Gff | rev | cut -d '/' -f4 | rev)
	echo "$Strain - $Organism"
	cat $Gff | grep -w 'gene' | wc -l
done
```

<!--
```
AJ516 - F.oxysporum_fsp_lactucae
22217
AJ520 - F.oxysporum_fsp_lactucae
20614
```
-->

```bash
AJ516 - F.oxysporum_fsp_lactucae
22129
AJ520 - F.oxysporum_fsp_lactucae
20639
```


## Assessing the Gene space in predicted transcriptomes:

```bash
for Assembly in $(ls gene_pred/final/*/*/publication/final_genes_appended_renamed.gene.fasta | grep -e 'AJ516' -e 'AJ520'); do
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
for File in $(ls gene_pred/busco/*/*/genes/*/short_summary_*.txt | grep -e 'AJ516' -e 'AJ520'); do  
echo $File;
cat $File | grep -e '(C)' -e 'Total';
done
```
<!--
```
gene_pred/busco/F.oxysporum_fsp_lactucae/AJ516/genes/run_final_genes_appended_renamed.gene/short_summary_final_genes_appended_renamed.gene.txt
	3670	Complete BUSCOs (C)
	3725	Total BUSCO groups searched
gene_pred/busco/F.oxysporum_fsp_lactucae/AJ520/genes/run_final_genes_appended_renamed.gene/short_summary_final_genes_appended_renamed.gene.txt
	3673	Complete BUSCOs (C)
	3725	Total BUSCO groups searched
```
 -->

```
gene_pred/busco/F.oxysporum_fsp_lactucae/AJ516/genes/run_final_genes_appended_renamed.gene/short_summary_final_genes_appended_renamed.gene.txt
	3667	Complete BUSCOs (C)
	3725	Total BUSCO groups searched
gene_pred/busco/F.oxysporum_fsp_lactucae/AJ520/genes/run_final_genes_appended_renamed.gene/short_summary_final_genes_appended_renamed.gene.txt
	3669	Complete BUSCOs (C)
	3725	Total BUSCO groups searched
```

## FunGap

Gene prediction was run on the new cluster from a screen session with a ssh connection to a worker node:

```bash
screen -a
ssh compute03

ProjDir=/oldhpc/home/groups/harrisonlab/project_files/fusarium

cat /oldhpc/home/groups/harrisonlab/project_files/fusarium/qc_rna/paired/F.oxysporum_fsp_cepae/[^control]*[^prelim][^prelim_rep1]/F/*_trim.fq.gz > $HOME/tmp/concat_1.fastq.gz
cat /oldhpc/home/groups/harrisonlab/project_files/fusarium/qc_rna/paired/F.oxysporum_fsp_cepae/[^control]*[^prelim][^prelim_rep1]/R/*_trim.fq.gz > $HOME/tmp/concat_2.fastq.gz

gunzip $HOME/tmp/concat*.fastq.gz

for Assembly in $(ls $ProjDir/repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -e 'AJ516' -e 'AJ520' | grep 'AJ516'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"

WorkDir=$HOME/tmp/fungap_Fus_${Strain}
OutDir=$ProjDir/gene_pred/fungap/$Organism/${Strain}
AcceptedHits=$(ls /oldhpc/data/scratch/armita/fusarium/alignment/star/*/${Strain}/concatenated/all_reads_concatenated.bam)
RefProteome=$(ls /oldhpc/home/groups/harrisonlab/project_files/fusarium/gene_pred/published/F.oxysporum_fsp_lycopersici/4287_chromosomal/4287_published_genes.fa.pep.fasta)

GeneModelName="$Organism"_"$Strain"_fungap
CurDir=$PWD

mkdir -p $WorkDir
cd $WorkDir

cp $Assembly assembly.fa
cp $AcceptedHits alignedRNA.bam

TrimF=$(ls ../concat_1.fastq)
TrimR=$(ls ../concat_2.fastq)

conda activate fungap

fungap.py \
  --num_cores 40 \
  --output_dir out \
  -1 $TrimF \
	-2 $TrimR \
  --genome_assembly assembly.fa \
  --augustus_species $GeneModelName \
  --sister_proteome $RefProteome \
  --trans_bam "alignedRNA.bam"

mkdir -p $OutDir
cp -r braker/* $OutDir/.
done

# rm -r $WorkDir

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
	for Genes in $(ls gene_pred/final/*/*/publication/final_genes_appended_renamed.pep.fasta | grep -e 'AJ516' -e 'AJ520'); do  
	echo $Genes
	$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
for Proteins in $(ls gene_pred/final/*/*/publication/final_genes_appended_renamed.pep.fasta | grep -e 'AJ516' -e 'AJ520'); do
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
for Proteome in $(ls gene_pred/final/*/*/publication/final_genes_appended_renamed.pep.fasta | grep -e 'AJ516' -e 'AJ520'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
OutDir=gene_pred/swissprot/$Organism/$Strain
SwissDbDir=../../uniprot/swissprot
SwissDbName=uniprot_sprot
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
done
```

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
for Proteome in $(ls gene_pred/final/F.*/*/*/final_genes_appended_renamed.pep.fasta | grep -e 'AJ516' -e 'AJ520'); do
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
for SplitDir in $(ls -d gene_pred/final_genes_split/*/* | grep -e 'AJ516' -e 'AJ520'); do
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
for Proteome in $(ls gene_pred/final/F.*/*/*/final_genes_appended_renamed.pep.fasta | grep -e 'AJ516' -e 'AJ520'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
qsub $ProgDir/submit_TMHMM.sh $Proteome
done
```

Those proteins with transmembrane domains were removed from lists of Signal
peptide containing proteins

```bash
for File in $(ls gene_pred/trans_mem/*/*/*_TM_genes_neg.txt | grep -e 'AJ516' -e 'AJ520'); do
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
F.oxysporum_fsp_lactucae - AJ516
Number of SigP proteins:
1926
Number without transmembrane domains:
1596
Number of gene models:
1593
F.oxysporum_fsp_lactucae - AJ520
Number of SigP proteins:
1850
Number without transmembrane domains:
1530
Number of gene models:
1526
```

### C) Identification of MIMP-flanking genes

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa | grep -e 'AJ516' -e 'AJ520'); do
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
F.oxysporum_fsp_lactucae - AJ516
The number of mimps identified:
191
The following transcripts intersect mimps:
162
160

F.oxysporum_fsp_lactucae - AJ520
The number of mimps identified:
220
The following transcripts intersect mimps:
176
174
```

Those genes that were predicted as secreted and within 2Kb of a MIMP
were identified:

```bash
for File in $(ls analysis/mimps/*/*/*_genes_in_2kb_mimp.txt | grep -e 'AJ516' -e 'AJ520'); do
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
F.oxysporum_fsp_lactucae - AJ516
25
25
F.oxysporum_fsp_lactucae - AJ520
28
33
```

sequences of these genes were extracted:

```bash
for Headers in $(ls analysis/mimps/F.oxysporum_fsp_lactucae/*/*_genes_in_2kb_mimp_secreted_headers.txt); do
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
for Proteome in $(ls gene_pred/final/F.*/*/*/final_genes_appended_renamed.pep.fasta | grep -e 'AJ516' -e 'AJ520'); do
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
for File in $(ls gene_pred/CAZY/*/*/*CAZY.out.dm | grep -e 'AJ516' -e 'AJ520'); do
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
F.oxysporum_fsp_lactucae - AJ516
number of CAZY proteins identified:
964
number of CAZY genes identified:
964
number of Secreted CAZY proteins identified:
408
number of Secreted CAZY genes identified:
408
F.oxysporum_fsp_lactucae - AJ520
number of CAZY proteins identified:
947
number of CAZY genes identified:
947
number of Secreted CAZY proteins identified:
404
number of Secreted CAZY genes identified:
404
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
for CAZY in $(ls gene_pred/CAZY/*/*/*_CAZY.out.dm.ps | grep -e 'AJ516' -e 'AJ520'); do
  Strain=$(echo $CAZY | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $CAZY | rev | cut -f3 -d '/' | rev)
  OutDir=$(dirname $CAZY)
  echo "$Organism - $Strain"
  Secreted=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp_no_trans_mem_headers.txt)
  Gff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/CAZY
  $ProgDir/summarise_CAZY.py --cazy $CAZY --inp_secreted $Secreted --inp_gff $Gff --summarise_family --trim_gene_id 2 --kubicek_2014
done
```

```
F.oxysporum_fsp_lactucae - AJ516
B-Galactosidases - 2
A-Galactosidases - 4
Polygalacturonase - 17
A-Arabinosidases - 22
Xylanases - 10
Polygalacturonate lyases - 23
B-Glucuronidases - 4
B-Glycosidases - 10
Cellulases - 19
Xyloglucanases - 1
other - 295
F.oxysporum_fsp_lactucae - AJ520
B-Galactosidases - 2
A-Galactosidases - 4
Polygalacturonase - 17
A-Arabinosidases - 23
Xylanases - 11
Polygalacturonate lyases - 23
B-Glucuronidases - 4
B-Glycosidases - 13
Cellulases - 21
Xyloglucanases - 1
other - 285
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














## 5.1  Identifying SIX gene homologs

## 5.1.a) Performing BLAST searches

BLast searches were performed against the genome:

```bash
  # for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep -e 'AJ516' -e 'AJ520' | grep 'AJ516'); do
  for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/pilon_min_500bp_renamed.fasta | grep -e 'AJ516' -e 'AJ520' | grep 'AJ520'); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    Query=analysis/blast_homology/six_genes/six-appended_parsed.fa
    qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
  done
```


```bash
cat analysis/blast_homology/F.oxysporum_fsp_lactucae/*/*_six-appended_parsed.fa_homologs.csv | cut -f1,5
```

```
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six10_(SIX10)_mRNA,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six11_(SIX11)_mRNA,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six12_(SIX12)_mRNA,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six13_(SIX13)_mRNA,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six14_(SIX14)_mRNA,_complete_cds	1
Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_1_(SIX1)_gene,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_2_(SIX2)_gene,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_isolate_FOL-MM10_secreted_in_xylem_3_(SIX3)_gene,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_isolate_IPO3_secreted_in_xylem_3_(SIX3)_gene,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_isolate_14844_secreted_in_xylem_3_(SIX3)_gene,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_3_(SIX3)_gene,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_SIX3_gene_for_Secreted_in_xylem_3_protein	0
Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_4_(SIX4)_gene,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_Six5_mRNA,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_5_(SIX5)_gene,_partial_cds	0
Fusarium_oxysporum_f._sp._passiflorae_SIX6_gene,_complete_cds	0
Fusarium_oxysporum_f._sp._niveum_SIX6_gene,_complete_cds	0
Fusarium_oxysporum_f._sp._vasinfectum_SIX6_gene,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_secreted_in_xylem_Six6_(SIX6)_mRNA,_complete_cds	0
Fusarium_oxysporum_f._sp._melonis_isolate_NRRL_26406_secreted_in_xylem_6-like_protein_(SIX6)_gene,_complete_cds	0
Fusarium_oxysporum_f._sp._radicis-cucumerinum_isolate_Afu-3_secreted_in_xylem_6-like_protein_(SIX6)_gene,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_6_(SIX6)_gene,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_secreted_in_xylem_Six7_(SIX7)_mRNA,_complete_cds	0
Fusarium_oxysporum_f._sp._lilii_isolate_NRRL_28395_secreted_in_xylem_7-like_protein_(SIX7)_gene,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_7_(SIX7)_gene,_complete_cds	0
Fusarium_oxysporum_SIX8_gene,_complete_cds	1
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six9_(SIX9)_mRNA,_complete_cds	3
```


## 5.1.b) Converting BLAST results to gff annotations

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
for BlastHits in $(ls analysis/blast_homology/*/*/*_six-appended_parsed.fa_homologs.csv | grep -e 'AJ516' -e 'AJ520'); do
Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)  
Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
HitsGff=analysis/blast_homology/$Organism/$Strain/"$Strain"_six-appended_parsed.fa_homologs.gff
Column2=BLAST_hit
NumHits=3
$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
done
```


## Alignment vs minion genomes

Illumina sequence data from each of the isolates was aligned against the minion
assemblies.


### Read alignment

```bash
for Reference in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa | grep -e 'AJ516' -e 'AJ520'); do
RefStrain=$(echo $Reference | rev | cut -f3 -d '/' | rev)
for StrainPath in $(ls -d qc_dna/paired/*/* | grep -e 'AJ516' -e 'AJ520'); do
Jobs=$(qstat | grep 'sub_bowtie' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_bowtie' | grep 'qw'| wc -l)
done
printf "\n"
echo $StrainPath
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
F_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1)
echo $F_Read
echo $R_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_${RefStrain}
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie.sh $Reference $F_Read $R_Read $OutDir $Strain
done
done
```

### Read coverage

Identify read coverage over each bp

```bash
  for Bam in $(ls analysis/genome_alignment/bowtie/*/*/vs_*/*_aligned_sorted.bam | grep 'vs_650'); do
    Target=$(echo $Bam | rev | cut -f2 -d '/' | rev)
    Strain=$(echo $Bam | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Bam | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain - $Target"
    OutDir=$(dirname $Bam)
    samtools depth -aa $Bam > $OutDir/${Organism}_${Strain}_${Target}_depth.tsv
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/coverage_analysis
    $ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_${Target}_depth.tsv > $OutDir/${Organism}_${Strain}_${Target}_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_${Target}_depth_10kb.tsv
  done
  for Target in "vs_650" "vs_1166"; do
    OutDir=analysis/genome_alignment/bowtie/grouped_${Target}
    mkdir -p $OutDir
    cat analysis/genome_alignment/bowtie/*/*/*/*_*_${Target}_depth_10kb.tsv > $OutDir/${Target}_grouped_depth.tsv
  done
```

### Plot read coverage


```R
library(readr)
setwd("~/Downloads/Aalt/coverage2")

appended_df <- read_delim("~/Downloads/Aalt/coverage2/vs_1166_grouped_depth.tsv", "\t", escape_double = FALSE, col_names = FALSE, col_types = cols(X4 = col_factor(levels = c("675", "97.0013", "97.0016", "650", "648", "24350", "1082", "1164", "635", "743", "1166", "1177"))), trim_ws = TRUE)

myFun <- function(x) {
  c(min = min(x), max = max(x),
    mean = mean(x), median = median(x),
    std = sd(x))
}

colnames(appended_df) <- c("contig","position", "depth", "strain")

appended_df$treatment <- paste(appended_df$strain , appended_df$contig)
tapply(appended_df$depth, appended_df$treatment, myFun)

df2 <- cbind(do.call(rbind, tapply(appended_df$depth, appended_df$treatment, myFun)))
write.csv(df2, '1166_contig_coverage.csv')

appended_df$depth <- ifelse(appended_df$depth > 100, 100, appended_df$depth)

# install.packages("ggplot2")
library(ggplot2)
require(scales)

for (i in 1:22){
contig = paste("contig", i, sep = "_")
p0 <- ggplot(data=appended_df[appended_df$contig == contig, ], aes(x=`position`, y=`depth`, group=1)) +
geom_line() +
labs(x = "Position (bp)", y = "Coverage") +
scale_y_continuous(breaks=seq(0,100,25), limits=c(0,100)) +
facet_wrap(~strain, nrow = 12, ncol = 1, strip.position = "left")
outfile = paste("1166_contig", i, "cov.jpg", sep = "_")
ggsave(outfile , plot = p0, device = 'jpg', path = NULL,
scale = 1, width = 500, height = 500, units = 'mm',
dpi = 150, limitsize = TRUE)
}

```
