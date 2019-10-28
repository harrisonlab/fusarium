# Genome Annotation

This documentation covers gene prediction, functional annotation and effector
identification performed on FoM and FoN genomes as part of the AHDB Fusarium
project.

## Assembly

Genome assembly was performed on MinION and MiSeq reads for FoM and FoN genomes
and commands are documented in the minion_assembly directory within this
repository.

## Repeatmasking

Repeatmasking was performed on MinION and MiSeq reads for FoM and FoN genomes
and commands are documented in the minion_assembly directory within this
repository.

## Gene prediction


# Gene Prediction


Gene prediction followed three steps:
	Gene model training
		- Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline
	Gene prediction
		- Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.


#Gene prediction

Gene prediction was performed for Fusarium genomes. Two gene prediction
approaches were used:

Gene prediction using Braker1

## Gene prediction 1 - Braker1 gene model training and prediction

Gene prediction was performed using Braker1.

First, RNAseq data was aligned to Fusarium genomes.
* RNAseq data from the Fusarium HAPI project was used
* qc of RNA seq data is detailed in the README file of this repository:


#### Aligning

RNAseq data was aligned to the assemblies using STAR to provide evidence for gene models

```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa | grep -w -e 'FON_63' -e 'Stocks4' | grep -w -e 'Stocks4'); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    for RNADir in $(ls -d qc_rna/paired/F.oxysporum_fsp_cepae/*); do
      FileNum=$(ls $RNADir/F/*_trim.fq.gz | wc -l)
      for num in $(seq 1 $FileNum); do
        while [ $Jobs -gt 1 ]; do
          sleep 1m
          printf "."
          Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
        done
        printf "\n"
        FileF=$(ls $RNADir/F/*_trim.fq.gz | head -n $num | tail -n1)
        FileR=$(ls $RNADir/R/*_trim.fq.gz | head -n $num | tail -n1)
        echo $FileF
        echo $FileR
        Prefix=$(echo $FileF | rev | cut -f1 -d '/' | rev | sed "s/_R.*_trim.fq.gz//g")
        Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
        Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
        echo "$Timepoint"
        OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Prefix
        ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
        qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
      done
    done
  done
```

<!-- Alignment stats were collected:

```bash
for File in $(ls alignment/star/*/*/*/*/star_aligmentLog.final.out| grep -w -e 'FON_63' -e 'Stocks4'); do
Sample=$(echo $File | rev | cut -f2 -d '/' | rev);
ReadNumU=$(cat $File | grep 'Uniquely' | grep 'number' | cut -f2);
ReadPercU=$(cat $File | grep 'Uniquely' | grep '%' | cut -f2);
ReadNumM=$(cat $File | grep 'multiple' | grep 'Number' | cut -f2);
ReadPercM=$(cat $File | grep 'multiple' | grep '%' | cut -f2);
echo -e "$Sample""\t""$ReadNumU""\t""$ReadPercU""\t""$ReadNumM""\t""$ReadPercM";  
done
``` -->


Alignments were concatenated prior to gene prediction. This step was done through
a qlogin session on the cluster.

```bash
	for StrainDir in $(ls -d ls alignment/star/*/* | grep -w -e 'FON_63' -e 'Stocks4' | grep -w -e 'Stocks4'); do
	BamFiles=$(ls $StrainDir/*/*/star_aligmentAligned.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
	OutDir=$StrainDir/concatenated
	mkdir -p $OutDir
	samtools merge -f $OutDir/all_reads_concatenated.bam $BamFiles
	done
```


<!--

Alignments were concatenated prior to running cufflinks:
Cufflinks was run to produce the fragment length and stdev statistics:

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa | grep -w -e 'FON_63' -e 'Stocks4'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
AcceptedHits=$(ls alignment/star/$Organism/$Strain/concatenated/all_reads_concatenated.bam)
OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim
echo "$Organism - $Strain"
mkdir -p $OutDir
cufflinks -o $OutDir/cufflinks -p 8 --max-intron-length 4000 $AcceptedHits 2>&1 | tee $OutDir/cufflinks/cufflinks.log
done
``` -->

#### Braker prediction

Before braker predictiction was performed, I double checked that I had the
genemark key in my user area and copied it over from the genemark install
directory:

```bash
	ls ~/.gm_key
	cp /home/armita/prog/genemark/gm_key_64 ~/.gm_key
```

```bash

for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w -e 'FON_63' -e 'Stocks4' | grep -w -e 'Stocks4'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
AcceptedHits=$(ls alignment/star/$Organism/$Strain/concatenated/all_reads_concatenated.bam)
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
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w -e 'FON_63' -e 'Stocks4' | grep -w -e 'Stocks4'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
mkdir -p $OutDir
AcceptedHits=$(ls alignment/star/$Organism/$Strain/concatenated/all_reads_concatenated.bam)
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
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w -e 'FON_63' -e 'Stocks4' | grep -w -e 'Stocks4'); do
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
for BrakerGff in $(ls gene_pred/braker/F.*/*/*/augustus.gff3 | grep -e 'FON_63' -e 'Stocks4'); do
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

```
gene_pred/final_genes/F.oxysporum_fsp_mathioli/Stocks4/final
19185
1197
20382

gene_pred/final_genes/F.oxysporum_fsp_narcissi/FON_63/final
19243
1443
20686
```


In preperation for submission to ncbi, gene models were renamed and duplicate gene features were identified and removed.
 * no duplicate genes were identified


```bash
for GffAppended in $(ls gene_pred/final/*/*/final/final_genes_appended.gff3); do
Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
FinalDir=gene_pred/final_genes/$Organism/$Strain/final
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
$ProgDir/remove_dup_features.py --inp_gff $GffAppended
$ProgDir/remove_dup_features.py --inp_gff $GffAppended | grep -A2 'Duplicate gene found' | tail -n1 | cut -f2 -d'=' > $FinalDir/filter_list.tmp
GffFiltered=$FinalDir/filtered_duplicates.gff
cat $GffAppended | grep -v -w -f $FinalDir/filter_list.tmp > $GffFiltered
rm $FinalDir/filter_list.tmp
GffRenamed=$FinalDir/final_genes_appended_renamed.gff3
LogFile=$FinalDir/final_genes_appended_renamed.log
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
$ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
rm $GffFiltered

Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs/*_softmasked_repeatmasker_TPSI_appended.fa)
$ProgDir/gff2fasta.pl $Assembly $GffRenamed $FinalDir/final_genes_appended_renamed

# The proteins fasta file contains * instead of Xs for stop codons, these should
# be changed
sed -i 's/\*/X/g' $FinalDir/final_genes_appended_renamed.pep.fasta
done
```

```bash
for Gff in $(ls gene_pred/final_genes/*/*/final/final_genes_appended_renamed.gff3 | grep -w -e 'FON_63' -e 'Stocks4'); do
	Strain=$(echo $Gff | rev | cut -d '/' -f3 | rev)
	Organism=$(echo $Gff | rev | cut -d '/' -f4 | rev)
	echo "$Strain - $Organism"
	cat $Gff | grep -w 'gene' | wc -l
done
```

```
Stocks4 - F.oxysporum_fsp_mathioli
20361
FON_63 - F.oxysporum_fsp_narcissi
20545
```


## Assessing the Gene space in predicted transcriptomes:

```bash
for Assembly in $(ls gene_pred/final_genes/*/*/final/final_genes_appended_renamed.gene.fasta | grep -w -e 'FON_63' -e 'Stocks4'); do
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
for Proteome in $(ls gene_pred/final_genes/F.*/*/*/final_genes_appended_renamed.pep.fasta | grep -w -e 'FON_63' -e 'Stocks4' | grep 'Stocks4'); do
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
for Proteome in $(ls gene_pred/final_genes/F.*/*/*/final_genes_appended_renamed.pep.fasta | grep -w -e 'FON_63' -e 'Stocks4' | grep 'Stocks4'); do
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
for SplitDir in $(ls -d gene_pred/final_genes_split/*/* | grep -w -e 'Stocks4' -e 'FON_63' | grep 'Stocks4'); do
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
for Proteome in $(ls gene_pred/final_genes/F.*/*/*/final_genes_appended_renamed.pep.fasta | grep -w -e 'Stocks4' -e 'FON_63'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
qsub $ProgDir/submit_TMHMM.sh $Proteome
done
```

Those proteins with transmembrane domains were removed from lists of Signal
peptide containing proteins

```bash
for File in $(ls gene_pred/trans_mem/*/*/*_TM_genes_neg.txt | grep -w -e 'Stocks4' -e 'FON_63'); do
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
F.oxysporum_fsp_mathioli - Stocks4
Number of SigP proteins:
1737
Number without transmembrane domains:
1436
Number of gene models:
1436
F.oxysporum_fsp_narcissi - FON_63
Number of SigP proteins:
1880
Number without transmembrane domains:
1555
Number of gene models:
1550
```

### C) Identification of MIMP-flanking genes

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa | grep -w -e 'FON_63' -e 'Stocks4'); do
Organism=$(echo "$Assembly" | rev | cut -d '/' -f4 | rev)
Strain=$(echo "$Assembly" | rev | cut -d '/' -f3 | rev)
GeneGff=$(ls gene_pred/final_genes/$Organism/"$Strain"/final/final_genes_appended_renamed.gff3)
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
F.oxysporum_fsp_mathioli - Stocks4
The number of mimps identified:
124
The following transcripts intersect mimps:
110
110

F.oxysporum_fsp_narcissi - FON_63
The number of mimps identified:
241
The following transcripts intersect mimps:
214
212
```

Those genes that were predicted as secreted and within 2Kb of a MIMP
were identified:

```bash
for File in $(ls analysis/mimps/*/*/*_genes_in_2kb_mimp.txt | grep -w -e 'FON_63' -e 'Stocks4'); do
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
F.oxysporum_fsp_mathioli - Stocks4
15
15
F.oxysporum_fsp_narcissi - FON_63
49
49
```


## C) CAZY proteins

Carbohydrte active enzymes were idnetified using CAZYfollowing recomendations
at http://csbl.bmb.uga.edu/dbCAN/download/readme.txt :

```bash
for Proteome in $(ls gene_pred/final_genes/F.*/*/*/final_genes_appended_renamed.pep.fasta | grep -w -e 'FON_63' -e 'Stocks4'); do
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
for File in $(ls gene_pred/CAZY/*/*/*CAZY.out.dm | grep -w -e 'FON_63' -e 'Stocks4'); do
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
Gff=$(ls gene_pred/final_genes/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
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
F.oxysporum_fsp_mathioli - Stocks4
number of CAZY proteins identified:
902
number of CAZY genes identified:
902
number of Secreted CAZY proteins identified:
377
number of Secreted CAZY genes identified:
377
F.oxysporum_fsp_narcissi - FON_63
number of CAZY proteins identified:
949
number of CAZY genes identified:
949
number of Secreted CAZY proteins identified:
398
number of Secreted CAZY genes identified:
398
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
for CAZY in $(ls gene_pred/CAZY/*/*/*_CAZY.out.dm.ps | grep -w -e 'FON_63' -e 'Stocks4'); do
  Strain=$(echo $CAZY | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $CAZY | rev | cut -f3 -d '/' | rev)
  OutDir=$(dirname $CAZY)
  echo "$Organism - $Strain"
  Secreted=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp_no_trans_mem_headers.txt)
  Gff=gene_pred/final_genes/$Organism/$Strain/final/final_genes_appended_renamed.gff3
  ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/CAZY
  $ProgDir/summarise_CAZY.py --cazy $CAZY --inp_secreted $Secreted --inp_gff $Gff --summarise_family --trim_gene_id 2 --kubicek_2014
done | less -S
```

```
F.oxysporum_fsp_mathioli - Stocks4
B-Galactosidases - 2
A-Galactosidases - 4
Polygalacturonase - 13
A-Arabinosidases - 19
Xylanases - 10
Polygalacturonate lyases - 22
B-Glucuronidases - 4
B-Glycosidases - 10
Cellulases - 19
other - 273
Xyloglucanases - 1
F.oxysporum_fsp_narcissi - FON_63
B-Galactosidases - 2
A-Galactosidases - 4
Polygalacturonase - 12
A-Arabinosidases - 25
Xylanases - 9
Polygalacturonate lyases - 28
B-Glucuronidases - 4
B-Glycosidases - 18
Cellulases - 17
other - 277
Xyloglucanases - 1
```

## D) AntiSMASH

Antismash was run to identify clusters of secondary metabolite genes within
the genome. Antismash was run using the weserver at:
http://antismash.secondarymetabolites.org


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

## 5.1  Identifying SIX gene homologs

## 5.1.a) Performing BLAST searches

BLast searches were performed against the genome:

```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa | grep -w -e 'FON_63' -e 'Stocks4' | grep 'FON_63'); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    Query=analysis/blast_homology/six_genes/six-appended_parsed.fa
    qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
  done
```


```bash
cat analysis/blast_homology/F.oxysporum_fsp_mathioli/*/*_six-appended_parsed.fa_homologs.csv | cut -f1,5
cat analysis/blast_homology/*/FON_63/*_six-appended_parsed.fa_homologs.csv | cut -f1,5

```

```
ID	No.hits
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six10_(SIX10)_mRNA,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six11_(SIX11)_mRNA,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six12_(SIX12)_mRNA,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six13_(SIX13)_mRNA,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six14_(SIX14)_mRNA,_complete_cds	0
Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_1_(SIX1)_gene,_complete_cds	3
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
Fusarium_oxysporum_f._sp._lilii_isolate_NRRL_28395_secreted_in_xylem_7-like_protein_(SIX7)_gene,_complete_cds0
Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_7_(SIX7)_gene,_complete_cds	0
Fusarium_oxysporum_SIX8_gene,_complete_cds	1
Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six9_(SIX9)_mRNA,_complete_cds	3
```


## 5.1.b) Converting BLAST results to gff annotations

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
for BlastHits in $(ls analysis/blast_homology/*/*/*_six-appended_parsed.fa_homologs.csv | grep -w -e 'FON_63' -e 'Stocks4' | grep 'FON_63'); do
Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)  
Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
HitsGff=analysis/blast_homology/$Organism/$Strain/"$Strain"_six-appended_parsed.fa_homologs.gff
Column2=BLAST_hit
NumHits=3
$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
done
```


### E) Effector-like proteins

### EffectorP

```bash
for Proteome in $(ls gene_pred/final/F.*/*/*/final_genes_appended_renamed.pep.fasta | grep -w -e 'FON_63' -e 'Stocks4'); do
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
for File in $(ls analysis/effectorP/*/*/*_EffectorP.txt | grep -w -e 'FON_63' -e 'Stocks4'); do
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
Gff=$(ls gene_pred/final/$Organism/$Strain/*/final_genes_appended_renamed.gff3)
EffectorP_Gff=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.gff/g')
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $OutFileHeaders $Gff effectorP ID > $EffectorP_Gff
cat $EffectorP_Gff | grep -w 'gene' | wc -l
done > tmp.txt
```

```
F.oxysporum_fsp_mathioli - Stocks4
332
332
F.oxysporum_fsp_narcissi - FON_63
385
385
```

#### SSCP

Small secreted cysteine rich proteins were identified within secretomes. These
proteins may be identified by EffectorP, but this approach allows direct control
over what constitutes a SSCP.

```bash
for Secretome in $(ls gene_pred/final_genes_signalp-4.1/*/*/*_final_sp_no_trans_mem.aa | grep -e 'AJ516' -e 'AJ520'); do
Strain=$(echo $Secretome| rev | cut -f2 -d '/' | rev)
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=analysis/sscp/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/sscp
$ProgDir/sscp_filter.py --inp_fasta $Secretome --max_length 300 --threshold 3 --out_fasta $OutDir/"$Strain"_sscp_all_results.fa
cat $OutDir/"$Strain"_sscp_all_results.fa | grep 'Yes' > $OutDir/"$Strain"_sscp.fa
printf "number of SSC-rich genes:\t"
cat $OutDir/"$Strain"_sscp.fa | grep '>' | tr -d '>' | cut -f1 -d '.' | sort | uniq | wc -l
printf "Number of effectors predicted by EffectorP:\t"
EffectorP=$(ls analysis/effectorP/$Organism/$Strain/*_EffectorP_secreted_headers.txt)
cat $EffectorP | wc -l
printf "Number of SSCPs predicted by both effectorP and this approach: \t"
cat $OutDir/"$Strain"_sscp.fa | grep '>' | tr -d '>' > $OutDir/"$Strain"_sscp_headers.txt
cat $OutDir/"$Strain"_sscp_headers.txt $EffectorP | cut -f1 | sort | uniq -d | wc -l
echo ""
done
```

```
F.oxysporum_fsp_lactucae - AJ516
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	328
number of SSC-rich genes:	326
Number of effectors predicted by EffectorP:	375
Number of SSCPs predicted by both effectorP and this approach: 	240

F.oxysporum_fsp_lactucae - AJ520
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	299
number of SSC-rich genes:	299
Number of effectors predicted by EffectorP:	347
Number of SSCPs predicted by both effectorP and this approach: 	217
```
