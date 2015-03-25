#!/bin/bash

# -------------------------
#	Perfrom data QC
# -------------------------
# qc the raw reads before performing an alignment
for Strainz in $(ls -d raw_dna/paired/F.oxysporum_fsp_cepae/*); do 
	ProgPath=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
	F_IN=$(ls $Strainz/F/*.fastq.gz)
	R_IN=$(ls $Strainz/R/*.fastq.gz)
	IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
	SeqType=dna
	qsub "$ProgPath"/rna_qc_fastq-mcf.sh "$F_IN" "$R_IN" "$IlluminaAdapters" "$SeqType"
done

mkdir -p repeat_masked_old/F.oxysporum_fsp_cepae/Fus2/Fus2_assembly.81_repmask
mv repeat_masked/F.oxysporum_fsp_cepae/Fus2/Fus2_assembly.81_repmask repeat_masked_old/F.oxysporum_fsp_cepae/Fus2/Fus2_assembly.81_repmask

# -------------------------
#	Align against assembled genomes
# -------------------------
# For each strain perform an alignment of raw reads against the genomes of all strains.
for Pathz in $(ls -d raw_dna/paired/F.oxysporum_fsp_cepae/*); do  
Strain=$(echo Strainz | cut -d '/' -f4)
echo "using reads for $Strain"
ProgPath=/home/armita/git_repos/emr_repos/tools/pathogen/lineage_specific_regions
F_IN=$(ls $Pathz/F/*.fastq.gz)
R_IN=$(ls $Pathz/R/*.fastq.gz)
for Assemblyz in $(ls repeat_masked/F.oxysporum_fsp_cepae/*/*/*_contigs_unmasked.fa); do
basename $Assemblyz
qsub "$ProgPath"/bowtie2_alignment_pipe.sh	$F_IN $R_IN $Assemblyz
done
done

SummaryFile=assembly/ls_contigs/alignment_summaries.txt
printf "" > "$SummaryFile"
for OUTPUT in $(ls bowtie2_alignment_pipe.sh.e*); do 
ID=$(echo $OUTPUT | rev | cut -d 'e' -f1 | rev | less); 
cat bowtie2_alignment_pipe.sh.o"$ID" | grep -E "Trimmed .* reads .*/F/|Output files: " | sed -e 's/.*\/F\///g' | cut -f1 -d ')' | cut -f2 -d '"' >> "$SummaryFile"; 
cat $OUTPUT >> "$SummaryFile"; 
printf "\n" >> "$SummaryFile"; 
done

# -------------------------
#	Perform alignments against
#	publicly available genomes
# -------------------------

Fo4287_Ec=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Ma_et_al_2010/F.oxysporum_fsp.lycopersici_4287.fasta 
Fo4287_Mt=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Ma_et_al_2010/fusarium_oxysporum_f._sp._lycopersici_mitochondrion_2_contigs.fasta
Fo4287_Out=assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Ma_et_al_2010/F.oxysporum_fsp.lycopersici_4287_appended.fasta
cat $Fo4287_Ec $Fo4287_Mt > $Fo4287_Out

Fo47_Ec=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs.fasta
Fo47_Linear=assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs_linear.fasta
cat $Fo47_Ec | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' > $Fo47_Linear

for Pathz in $(ls -d raw_dna/paired/F.oxysporum_fsp_cepae/D2); do  
Strain=$(echo $Pathz | cut -d '/' -f4)
echo "using reads for $Strain"
ProgPath=/home/armita/git_repos/emr_repos/tools/pathogen/lineage_specific_regions
F_IN=$(ls $Pathz/F/*.fastq.gz)
R_IN=$(ls $Pathz/R/*.fastq.gz)
for Assemblyz in "$Fo4287_Out" "$Fo47_Linear"; do
basename "$Assemblyz"
qsub "$ProgPath"/bowtie2_alignment_pipe.sh	$F_IN $R_IN $Assemblyz
done
done

InFile=assembly/ls_contigs/F.oxysporum_fsp_cepae/PG/vs_PG_assembly.81_repmask_contigs/PG_vs_PG_assembly.81_repmask.bam
# bowtie 2 returns aligned reads when it outputs non-aligned reads
# the number of reads unaligned discordantly and concordantly 
# which are also not aligned singly within the genome is the same 
# number as returned when samtools is used to extract all reads with
# the SAM flag 8
# Building a SMALL index
# 6437691 reads; of these:
#   6437691 (100.00%) were paired; of these:
#     1345863 (20.91%) aligned concordantly 0 times
#     4891988 (75.99%) aligned concordantly exactly 1 time
#     199840 (3.10%) aligned concordantly >1 times
#     ----
#     1345863 pairs aligned concordantly 0 times; of these:
#       1091492 (81.10%) aligned discordantly 1 time
#     ----
#     254371 pairs aligned 0 times concordantly or discordantly; of these:
#       508742 mates make up the pairs; of these:
#         246494 (48.45%) aligned 0 times
#         146460 (28.79%) aligned exactly 1 time
#         115788 (22.76%) aligned >1 times
# samtools view -f 8 "$InFile" > tmp.txt
# cat tmp.txt | grep 'M01678' | wc -l
# 246494
# When looking at the sam flags present in this output file
# the FLAG 8 also pull out other flags that are used to show alignments.
# cat tmp.txt | cut -f2 | sort | uniq -c | sort
#   25434 153 - read2 map -  reads1 unmap +
#   25772 137 - read2 map + read1 munmap +
#   37349 73 -	read1 map + read2 unmap
#   40079 89 - read1 map + read2 unmap -
#   58930 141 -	read2 unmap + read1 unmap +
#   58930 77 -  read1 unmap + read2 unmap +
# FLAGS 141 and 77 will extract the desired unmapped reads.

# for Pathz in $(ls assembly/ls_contigs/F.oxysporum_fsp_cepae/D2/vs_*/*_sorted.bam); do  
# OutFileF=$(echo $Pathz | sed 's/.bam/_unaligned_F.fastq/g')
# OutFileR=$(echo $Pathz | sed 's/.bam/_unaligned_R.fastq/g')
# samtools view -f 77 "$Pathz" > $OutFileF
# samtools view -f 141 "$Pathz" > $OutFileR
# done

# -------------------------
#	Extract unaligned reads
# -------------------------


for Pathz in $(ls assembly/ls_contigs/F.oxysporum_fsp_cepae/*/vs_*/*_sorted.bam); do  
OutFileF=$(echo $Pathz | sed 's/.bam/_unaligned_F.txt/g')
OutFileR=$(echo $Pathz | sed 's/.bam/_unaligned_R.txt/g')
samtools view -f 77 "$Pathz" | cut -f1 > $OutFileF
samtools view -f 141 "$Pathz" | cut -f1 > $OutFileR
done

for Pathz in $(ls -d assembly/ls_contigs/F.oxysporum_fsp_cepae/*/vs_Ma*); do
FileF=$(ls $Pathz/*_F.txt)
FileR=$(ls $Pathz/*_R.txt)
sed -i $FileF -e 's/#CGATGT/#CGATGT\/1/g'
sed -i $FileR -e 's/#CGATGT/#CGATGT\/2/g'
done

for Pathz in $(ls assembly/ls_contigs/F.oxysporum_fsp_cepae/*/vs_*/*_sorted.bam); do  
FilterPath="/home/armita/git_repos/emr_repos/tools/pathogen/lineage_specific_regions"
ThisDir=$(dirname $Pathz)
HeadersF=$(echo $Pathz | sed 's/.bam/_unaligned_F.txt/g')
HeadersR=$(echo $Pathz | sed 's/.bam/_unaligned_R.txt/g')
OutFileF=$(echo $Pathz | sed 's/.bam/_unaligned_F.fastq/g')
OutFileR=$(echo $Pathz | sed 's/.bam/_unaligned_R.fastq/g')
Reads1=$(ls $ThisDir/*vs_*_unaligned.1.fastq)
Reads2=$(ls $ThisDir/*vs_*_unaligned.2.fastq)
$FilterPath/fastq_filter.py $HeadersF $Reads1 > $OutFileF
$FilterPath/fastq_filter.py $HeadersR $Reads2 > $OutFileR
done

# -------------------------
#	Assemble unaligned reads
# -------------------------

ProgPath=/home/armita/git_repos/emr_repos/tools/pathogen/lineage_specific_regions
InsLgth=600

# Assemble 125 regions
for Pathz in $(ls -d assembly/ls_contigs/F.oxysporum_fsp_cepae/125/*); do
FileF=$(ls $Pathz/*_F.fastq)
FileR=$(ls $Pathz/*_R.fastq)
HashLgth=71
ExpCov=40
CutCov=20
qsub $ProgPath/assemble_ls_regions.sh $HashLgth $FileF $FileR $ExpCov $CutCov $InsLgth
done

# 55
for Pathz in $(ls -d assembly/ls_contigs/F.oxysporum_fsp_cepae/55/*); do
	FileF=$(ls $Pathz/*_F.fastq)
	FileR=$(ls $Pathz/*_R.fastq)
	HashLgth=61
	ExpCov=29
	CutCov=14
	qsub $ProgPath/assemble_ls_regions.sh $HashLgth $FileF $FileR $ExpCov $CutCov $InsLgth
done

# A23
for Pathz in $(ls -d assembly/ls_contigs/F.oxysporum_fsp_cepae/A23/*); do
	FileF=$(ls $Pathz/*_F.fastq)
	FileR=$(ls $Pathz/*_R.fastq)
	HashLgth=32
	ExpCov=16
	CutCov=20
	qsub $ProgPath/assemble_ls_regions.sh $HashLgth $FileF $FileR $ExpCov $CutCov $InsLgth
done

# A28
for Pathz in $(ls -d assembly/ls_contigs/F.oxysporum_fsp_cepae/A28/*); do
	FileF=$(ls $Pathz/*_F.fastq)
	FileR=$(ls $Pathz/*_R.fastq)
	HashLgth=61
	ExpCov=35
	CutCov=17
	qsub $ProgPath/assemble_ls_regions.sh $HashLgth $FileF $FileR $ExpCov $CutCov $InsLgth
done

# D2
for Pathz in $(ls -d assembly/ls_contigs/F.oxysporum_fsp_cepae/D2/vs_Ma*); do
	FileF=$(ls $Pathz/*_F.fastq)
	FileR=$(ls $Pathz/*_R.fastq)
	HashLgth=41
	ExpCov=10
	CutCov=5
	qsub $ProgPath/assemble_ls_regions.sh $HashLgth $FileF $FileR $ExpCov $CutCov $InsLgth
done

# Fus2
for Pathz in $(ls -d assembly/ls_contigs/F.oxysporum_fsp_cepae/Fus2/*); do
	FileF=$(ls $Pathz/*_F.fastq)
	FileR=$(ls $Pathz/*_R.fastq)
	HashLgth=49
	ExpCov=40
	CutCov=20
	qsub $ProgPath/assemble_ls_regions.sh $HashLgth $FileF $FileR $ExpCov $CutCov $InsLgth
done

# HB17
for Pathz in $(ls -d assembly/ls_contigs/F.oxysporum_fsp_cepae/HB17/*); do
	FileF=$(ls $Pathz/*_F.fastq)
	FileR=$(ls $Pathz/*_R.fastq)
	HashLgth=51
	ExpCov=26
	CutCov=13
	qsub $ProgPath/assemble_ls_regions.sh $HashLgth $FileF $FileR $ExpCov $CutCov $InsLgth
done

# PG
for Pathz in $(ls -d assembly/ls_contigs/F.oxysporum_fsp_cepae/PG/*); do
	FileF=$(ls $Pathz/*_F.fastq)
	FileR=$(ls $Pathz/*_R.fastq)
	HashLgth=81
	ExpCov=44
	CutCov=22
	qsub $ProgPath/assemble_ls_regions.sh $HashLgth $FileF $FileR $ExpCov $CutCov $InsLgth
done

# -------------------------
#	Summarise assemblies
# -------------------------

for filename in $(ls assembly/ls_contigs/F.oxysporum_fsp_cepae/*/*/*_contigs_ls_*.txt); do
	echo $filename | cut -f6 -d '/' ; 
	cat $filename; 
	echo "";
done > assembly/ls_contigs/assembly_summaries.txt


cat assembly/ls_contigs/assembly_summaries.txt | grep -e "Number of bases in contigs" | grep -v 'kb' | rev |  cut -f1 -d ' ' | rev | less

# -------------------------
#	Perform gene prediction
# -------------------------

#cat assembly/ls_contigs/F.oxysporum_fsp_cepae/125/vs_PG_assembly.81_repmask_contigs/125_vs_PG_assembly.81_repmask_sorted_unaligned_F.fastq | cut -f1 > tmp_unalF.fq 

ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/lineage_specific_regions
for Assemblyz in $(ls assembly/ls_contigs/F.oxysporum_fsp_cepae/*/*/*_contigs_ls_*.fa); do
	echo "$Assemblyz"
	qsub "$ProgDir"/augustus_ls_contig.sh fusarium_graminearum $Assemblyz
done

# -------------------------
#	Summarise gene predictions
# -------------------------

for filename in $(ls gene_pred/ls_contigs/F.oxysporum_fsp_cepae/*/*/*_aug_out.aa); do
	echo $filename | cut -f5,6 -d '/' ; 
	cat $filename | grep '>' | wc -l; 
	echo "";
done > gene_pred/ls_contigs/gene_summaries.txt
