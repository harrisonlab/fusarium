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
```


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
# Reads=$(ls /home/groups/harrisonlab/project_files/fusarium/raw_dna/minion/F.oxysporum/Stocks4/all_reads.fastq.gz)
Fast5Dir=$(ls -d /home/miseq_data/minion/2017/MINION_20170424_FNFAB42727_MN18323_sequencing_run_Fusarium_oxysporum_Stocks4/albacore1.1.1)
nanopolish extract -r $Fast5Dir | gzip -cf > stocks4_reads.fa.gz

Assembly=$(ls /home/groups/harrisonlab/project_files/fusarium/assembly/SMARTdenovo/F.oxysporum_f.sp_mathioli/Stocks4/Stocks4_SMARTdenovo.fasta)
cp $Assembly assembly.fa
bwa index assembly.fa
# cp $Reads reads.fa.gz
bwa mem -x ont2d -t 8 assembly.fa  stocks4_reads.fa.gz > aligned.sam
# samtools sort $OutDir/aligned.bam -f $OutDir/reads.sorted.bam
samtools view -b -S aligned.sam | samtools sort - -o tmp > reads.sorted.bam
# samtools view -u aligned.bam | samtools sort - -f reads.sorted.bam
samtools index reads.sorted.bam
NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_makerange.py assembly.fa > nanopolish_range.txt
# cat reads.fa.gz | gunzip -cf | head -n1 | cut -f1 -d ' ' | gzip -cf reads_2.fa.gz
cat stocks4_reads.fa.gz | gunzip -cf > stocks4_reads.fa

# cat assembly.fa | cut -f1 > assembly2.fa

for line in $(cat nanopolish_range.txt); do
	echo $line
	nanopolish variants -t 8 --reads stocks4_reads.fa --bam reads.sorted.bam --genome assembly.fa --ploidy 1 -w $line --consensus="$line"_consensus.fa --fix-homopolymers --min-candidate-frequency 0.1 > "$line"_variants.txt
done



cat nanopolish_range.txt | nanopolish variants --consensus polished.{1}.fa -w {1} -r stocks4_reads.fa -b reads.sorted.bam -g assembly2.fa -t 8 --min-candidate-frequency 0.1 -o nanopolish_out.txt

python $NanoPolishDir/nanopolish_makerange.py assembly.fa

# python $NanoPolishDir/nanopolish_makerange.py assembly.fa | parallel --results nanopolish.results -P 8 \
#   nanopolish variants --consensus polished.{1}.fa -w {1} -r $TMPDIR/WT_07-03-17_pass.fastq.gz -b $TMPDIR/reads.sorted.bam -g $TMPDIR/contigs_min_500bp.fasta -t 8 --min-candidate-frequency 0.1



TMPDIR=/tmp/nanopolish
mkdir -p $TMPDIR
bwa index $TMPDIR/contigs_min_500bp.fasta
bwa mem -x ont2d -t 8 $TMPDIR/contigs_min_500bp.fasta $TMPDIR/WT_07-03-17_pass.fastq.gz | samtools sort - -o $TMPDIR/reads.sorted.bam

NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_makerange.py $TMPDIR/contigs_min_500bp.fasta | parallel --results nanopolish.results -P 8 \
    nanopolish variants --consensus polished.{1}.fa -w {1} -r $TMPDIR/WT_07-03-17_pass.fastq.gz -b $TMPDIR/reads.sorted.bam -g $TMPDIR/contigs_min_500bp.fasta -t 8 --min-candidate-frequency 0.1
```
