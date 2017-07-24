# Stocks_Assembly
==========

This document details the commands used to assemble and annotate the Fus2 genome.

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
```bash
	RawDatDir=/home/miseq_data/minion/2017/MINION_20170424_FNFAB42727_MN18323_sequencing_run_Fusarium_oxysporum_Stocks4
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	mkdir -p $ProjectDir/raw_dna/minion/F.oxysporum/Stocks4
```

Sequence data was moved into the appropriate directories

```bash
	RawDatDir=/home/miseq_data/minion/2017/MINION_20170424_FNFAB42727_MN18323_sequencing_run_Fusarium_oxysporum_Stocks4
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	cp $RawDatDir/all_reads.fastq.gz $ProjectDir/raw_dna/minion/F.oxysporum/Stocks4/.
```

## Assembly


### Canu assembly

f. sp. mathioli

```bash
  Reads=$(ls raw_dna/minion/F.oxysporum/Stocks4/all_reads.fastq.gz)
  GenomeSz="53m"
  Strain=$(echo $Reads | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Reads | rev | cut -f3 -d '/' | rev)
  Prefix="$Strain"
  OutDir="assembly/canu-1.5/$Organism/$Strain"_canu
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/canu
  qsub $ProgDir/submit_canu_minion.sh $Reads $GenomeSz $Prefix $OutDir
```
<!--

# Aligning corrected reads to the longest minion read:

```bash
mkdir -p assembly/test/F.oxysporum/Stocks4
# move the longest minion read to a directory, pretending that it is an assembly
LongestReadFq=$(ls assembly/test/F.oxysporum/Stocks4/longest_read.fastq)
ReadDir=$(dirname $LongestReadFq)
printf ">F.oxysporum_Stocks4_minion\n" > $ReadDir/longest_read.fa
cat $LongestReadFq | head -n2 | tail -n1 >> $ReadDir/longest_read.fa

# Reads=raw_dna/pacbio/F.oxysporum_fsp_cepae/Fus2/extracted/concatenated_pacbio.fastq
# OutDir=analysis/genome_alignment/bwa/F.oxysporum_fsp_cepae/Fus2/vs_Stocks4
Reads=$(ls ../../../../harrir/projects/minion/stocks/fos.correctedReads.fasta.gz)
OutDir=analysis/genome_alignment/bwa/F.oxysporum_fsp_mathioli/Stocks4/vs_Stocks4
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/bwa
qsub $ProgDir/sub_bwa_pacbio.sh $ReadDir/longest_read.fa $Reads $OutDir
```

Sequence data for isolates with a data from two sequencing runs was aligned
against the Fus2 genome

```bash
Reference=$(ls assembly/test/F.oxysporum/Stocks4/longest_read.fa)
for StrainPath in $(ls -d qc_dna/paired/F.*/* | grep -e 'Fus2'); do
echo $StrainPath
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
F1_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n1 | tail -n1);
R1_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n1 | tail -n1);
F2_Read=$(ls $StrainPath/F/*_trim.fq.gz | head -n2 | tail -n1);
R2_Read=$(ls $StrainPath/R/*_trim.fq.gz | head -n2 | tail -n1);
echo $F1_Read
echo $R1_Read
echo $F2_Read
echo $R2_Read
OutDir=analysis/genome_alignment/bowtie/$Organism/$Strain/vs_Stocks4_longest_read
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bowtie/sub_bowtie_2lib.sh $Reference $F1_Read $R1_Read $F2_Read $R2_Read $OutDir $Strain
done
``` -->

# Assembly using SMARTdenovo

```bash
qlogin -pe smp 8
mkdir -p /tmp/last
cd /tmp/last


LongestReadFa=$(ls /home/groups/harrisonlab/project_files/fusarium/assembly/test/F.oxysporum/Stocks4/longest_read.fa)
lastdb stocksread_db $LongestReadFa
cp /home/harrir/projects/minion/stocks/fos.correctedReads.fasta.gz .
gunzip fos.correctedReads.fasta.gz
lastal -Q 0 -P 8 stocksread_db fos.correctedReads.fasta > last_out.txt

Fus2Assembly=$(ls /home/groups/harrisonlab/project_files/fusarium/repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_unmasked.fa)
lastdb Fus2_db $Fus2Assembly
lastal -Q 0 -P 8 Fus2_db $LongestReadFa > longest_read_vs_FoC.txt

# SMARTdevovo assembly
mkdir SMARTdenovo
cd SMARTdenovo
smartdenovo.pl -t 8 ../fos.correctedReads.fasta > wtasm.mak

make -f wtasm.mak 2>&1 | tee SMARTdenovo_run_log.txt
OutDir=/home/groups/harrisonlab/project_files/fusarium/assembly/SMARTdenovo/F.oxysporum_f.sp_mathioli/Stocks4
mkdir -p $OutDir
cp wtasm.dmo.lay.utg $OutDir/Stocks4_SMARTdenovo.fasta

rm /tmp/last
logout
```

Quast for the SMARTdenovo assembly:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_f.sp_mathioli/Stocks4/Stocks4_SMARTdenovo.fasta); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

Busco has replaced CEGMA and was run to check gene space in assemblies

```bash
for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_f.sp_mathioli/Stocks4/Stocks4_SMARTdenovo.fasta); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```

```bash
	for File in $(ls gene_pred/busco/F*/*/assembly/*/short_summary_*.txt); do  
		echo $File;
		cat $File | grep -e '(C)' -e 'Total';
	done
```


## Nanopolish scaffolding

```bash
qlogin -pe smp 8
WorkDir=/tmp/nanopolish
mkdir -p $WorkDir
cd $WorkDir

OutDir=assembly/nanopolish/F.oxysporum_f.sp_mathioli/Stocks4
mkdir -p $OutDir
Fast5Dir=$(ls -d /home/miseq_data/minion/2017/MINION_20170424_FNFAB42727_MN18323_sequencing_run_Fusarium_oxysporum_Stocks4/albacore1.1.1)
nanopolish extract -r $Fast5Dir | gzip -cf > stocks4_reads.fa.gz

Assembly=$(ls /home/groups/harrisonlab/project_files/fusarium/assembly/SMARTdenovo/F.oxysporum_f.sp_mathioli/Stocks4/Stocks4_SMARTdenovo.fasta)
cp $Assembly assembly.fa
bwa index assembly.fa
bwa mem -x ont2d -t 8 assembly.fa  stocks4_reads.fa.gz > aligned.sam
samtools view -b -S aligned.sam | samtools sort - -o tmp > reads.sorted.bam
samtools index reads.sorted.bam

nanopolish eventalign --progress --reads stocks4_reads.fa.gz --bam reads.sorted.bam --genome assembly.fa -t 7 --sam | samtools view -Su - | samtools sort - -o tmp > reads.eventalign.sorted.bam
samtools index reads.eventalign.sorted.bam

--sam -r ectocooler_subset.fasta -b ecto_subset.sorted.bam -g ecto_subset.contigs.fasta --models nanopolish_models.fofn | samtools view -b -S - | samtools sort > ecto_subset.eventalign.sorted.bam


NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_makerange.py assembly.fa > nanopolish_range.txt
cat stocks4_reads.fa.gz | gunzip -cf > stocks4_reads.fa

for line in $(cat nanopolish_range.txt); do
	echo $line
	nanopolish variants -t 8 --reads stocks4_reads.fa --bam reads.sorted.bam --genome assembly.fa --ploidy 1 -w $line --consensus="$line"_consensus.fa --fix-homopolymers --min-candidate-frequency 0.1 > "$line"_variants.txt
done
```

```bash
OutDir=assembly/nanopolish/F.oxysporum_f.sp_mathioli/Stocks4
mkdir -p $OutDir
Fast5Dir=$(ls -d /home/miseq_data/minion/2017/MINION_20170424_FNFAB42727_MN18323_sequencing_run_Fusarium_oxysporum_Stocks4/albacore1.1.1)
# nanopolish extract -r $Fast5Dir | gzip -cf > $OutDir/stocks4_reads.fa.gz

Assembly=$(ls /home/groups/harrisonlab/project_files/fusarium/assembly/SMARTdenovo/F.oxysporum_f.sp_mathioli/Stocks4/Stocks4_SMARTdenovo.fasta)

ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
qsub $ProgDir/sub_bwa_nanopolish.sh $Assembly $OutDir/stocks4_reads.fa.gz $OutDir

qsub $ProgDir/sub_nanopolish_variants.sh $Assembly $OutDir/stocks4_reads.fa.gz <output_directory>

samtools view -b -S $OutDir/aligned.sam | samtools sort - -o tmp > $OutDir/reads.sorted.bam
samtools index reads.sorted.bam
NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_makerange.py assembly.fa > nanopolish_range.txt
cat stocks4_reads.fa.gz | gunzip -cf > stocks4_reads.fa

for line in $(cat nanopolish_range.txt); do
	echo $line
	nanopolish variants -t 8 --reads stocks4_reads.fa --bam reads.sorted.bam --genome assembly.fa --ploidy 1 -w $line --consensus="$line"_consensus.fa --fix-homopolymers --min-candidate-frequency 0.1 > "$line"_variants.txt
done
```






# Assembly will with full trimming of reads:

Splitting reads and trimming adapters using porechop
```bash
	RawReads=raw_dna/minion/F.oxysporum/Stocks4/all_reads.fastq.gz
	OutDir=qc_dna/minion/F.oxysporum/Stocks4
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
	qsub $ProgDir/sub_porechop.sh $RawReads $OutDir
```

Read correction using Canu

```bash
for TrimReads in $(ls qc_dna/minion/F.oxysporum/Stocks4/*_trim.fastq.gz); do
Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
OutDir=assembly/canu-1.5/F.oxysporum_fsp_mathioli/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
qsub $ProgDir/sub_canu_correction.sh $TrimReads 60m $Strain $OutDir
done
```

Assembly using Canu

```bash
for CorrectedReads in $(ls assembly/canu-1.5/F.oxysporum_fsp_mathioli/Stocks4/Stocks4.trimmedReads.fasta.gz); do
Organism=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $CorrectedReads | rev | cut -f2 -d '/' | rev)
OutDir=assembly/canu-1.5/F.oxysporum_fsp_mathioli/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
qsub $ProgDir/sub_canu_assembly_only.sh $CorrectedReads 60m $Strain $OutDir
done

```

Assembbly using SMARTdenovo

```bash
for CorrectedReads in $(ls assembly/canu-1.5/F.oxysporum_fsp_mathioli/Stocks4/Stocks4.trimmedReads.fasta.gz); do
Organism=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $CorrectedReads | rev | cut -f2 -d '/' | rev)
Prefix=$Strain
OutDir=assembly/SMARTdenovo/F.oxysporum_fsp_mathioli/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/SMARTdenovo
qsub $ProgDir/sub_SMARTdenovo.sh $CorrectedReads $Prefix $OutDir
done
```

Quast for the SMARTdenovo assembly:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_mathioli/Stocks4/wtasm.dmo.lay.utg); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

Busco has replaced CEGMA and was run to check gene space in assemblies

```bash
for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_mathioli/Stocks4/wtasm.dmo.lay.utg); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```

```bash
# printf "Organism\tStrain\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/F*/*/assembly/*/short_summary_*.txt | grep 'Stocks4'); do  
# echo $File;
# Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
# Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
# printf "$Organism\t$Strain\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```


Error correction using racon:

```bash
Assembly=$(ls assembly/SMARTdenovo/F.oxysporum_fsp_mathioli/Stocks4/wtasm.dmo.lay.utg)
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ReadsFq=qc_dna/minion/F.oxysporum/Stocks4/all_reads_trim.fastq.gz
OutDir=assembly/SMARTdenovo/F.oxysporum_fsp_mathioli/Stocks4/racon2
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/racon
Iterations=10
qsub $ProgDir/sub_racon.sh $Assembly $ReadsFq $Iterations $OutDir
```

Quast and busco were run to assess the effects of racon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_mathioli/Stocks4/racon/*.fasta | grep 'round_4'); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```


```bash
for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_mathioli/Stocks4/racon2/*.fasta | grep "round_.*.fasta"); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```


```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/F*/*/assembly/*/short_summary_*.txt | grep 'Stocks4'); do  
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
Filename	Complete	Duplicated	Fragmented	Missing	Total
short_summary_wtasm.dmo.lay.txt	285	1	377	3063	3725
short_summary_wtasm_racon_round_1.txt	1940	11	873	912	3725
short_summary_wtasm_racon_round_2.txt	2239	13	744	742	3725
short_summary_wtasm_racon_round_3.txt	2308	16	717	700	3725
short_summary_wtasm_racon_round_4.txt	2281	17	750	694	3725
short_summary_wtasm_racon_round_5.txt	2312	17	709	704	3725
short_summary_wtasm_racon_round_6.txt	2286	14	749	690	3725
short_summary_wtasm_racon_round_7.txt	2306	17	741	678	3725
short_summary_wtasm_racon_round_8.txt	2314	12	740	671	3725
short_summary_wtasm_racon_round_9.txt	2314	19	736	675	3725
short_summary_wtasm_racon_round_10.txt	2350	18	730	645	3725

```


# Assembly correction using nanopolish


```bash
Assembly=$(ls assembly/SMARTdenovo/F.oxysporum_fsp_mathioli/Stocks4/racon2/wtasm_racon_round_10.fasta)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
Fast5Dir=$(ls -d /home/miseq_data/minion/2017/MINION_20170424_FNFAB42727_MN18323_sequencing_run_Fusarium_oxysporum_Stocks4/albacore1.1.1)
ReadDir=raw_dna/nanopolish/$Organism/$Strain
if [ -d $ReadDir ]; then
	echo "reads already extracted"
else
	echo "extracting reads"
	mkdir -p $ReadDir
	CurDir=$PWD
	cd $ReadDir
	nanopolish extract -r $Fast5Dir | gzip -cf > "$Strain"_reads.fa.gz
	cd $CurDir
fi

RawReads=$(ls $ReadDir/"$Strain"_reads.fa.gz)
OutDir=$(dirname $Assembly)
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
# submit alignments for nanoppolish
qsub $ProgDir/sub_bwa_nanopolish.sh $Assembly $RawReads $OutDir/nanopolish
```

 Split the assembly into 50Kb fragments and submit each to the cluster for
 nanopolish correction

```bash
Assembly=$(ls assembly/SMARTdenovo/F.oxysporum_fsp_mathioli/Stocks4/racon2/wtasm_racon_round_10.fasta)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)
RawReads=$(ls raw_dna/nanopolish/$Organism/$Strain/"$Strain"_reads.fa.gz)
AlignedReads=$(ls $OutDir/nanopolish/reads.sorted.bam)

NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_makerange.py $Assembly > $OutDir/nanopolish/nanopolish_range.txt

Ploidy=1
echo "nanopolish log:" > nanopolish_log.txt
for Region in $(cat $OutDir/nanopolish/nanopolish_range.txt | tail -n+21); do
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | wc -l)
done		
printf "\n"
echo $Region
echo $Region >> nanopolish_log.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
qsub $ProgDir/sub_nanopolish_variants.sh $Assembly $RawReads $AlignedReads $Ploidy $Region $OutDir/$Region
done
```

```bash
Assembly=$(ls assembly/SMARTdenovo/F.oxysporum_fsp_mathioli/Stocks4/racon2/wtasm_racon_round_10.fasta)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
OutDir=assembly/SMARTdenovo/$Organism/$Strain/nanopolish
mkdir -p $OutDir
# cat "" > $OutDir/"$Strain"_nanoplish.fa
NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_merge.py assembly/SMARTdenovo/$Organism/$Strain/racon2/*/*.fa > $OutDir/"$Strain"_nanoplish.fa

# for File in $(ls assembly/SMARTdenovo/F.*/*/racon2/*/*.fa | grep ":*-"); do
# 	# echo $File
# 	cat $File >> $OutDir/"$Strain"_nanoplish.fa
# done

echo "" > tmp.txt
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $OutDir/"$Strain"_nanoplish.fa --out $OutDir/"$Strain"_nanoplish_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt

```


Quast and busco were run to assess the effects of racon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_mathioli/Stocks4/nanopolish/Stocks4_nanoplish_min_500bp_renamed.fasta); do
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```


```bash
for Assembly in $(ls assembly/SMARTdenovo/F.oxysporum_fsp_mathioli/Stocks4/nanopolish/Stocks4_nanoplish_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco2.sh $Assembly $BuscoDB $OutDir
done
```
