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

Insert sizes of the RNA seq library were unknown until a draft alignment could
be made. To do this tophat and cufflinks were run, aligning the reads against a
single genome. The fragment length and stdev were printed to stdout while
cufflinks was running.

```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa | grep -w -e 'FON_63' -e 'Stocks4'); do
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
for StrainDir in $(ls -d ls alignment/star/*/* | grep -w -e 'FON_63' -e 'Stocks4'); do
BamFiles=$(ls $StrainDir/*/*/star_aligmentAligned.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
OutDir=$StrainDir/concatenated
mkdir -p $OutDir
samtools merge -f $OutDir/all_reads_concatenated.bam $BamFiles
done
```




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

for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa | grep -w -e 'FON_63' -e 'Stocks4'); do
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
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa | grep -w -e 'FON_63' -e 'Stocks4'); do
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
# Jobs=$(qstat | grep 'sub_cuffli' | wc -l)
# while [ $Jobs -ge 1 ]; do
# sleep 10m
# printf "."
# Jobs=$(qstat | grep 'sub_cuffli' | wc -l)
# done
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa | grep -w -e 'FON_63' -e 'Stocks4'); do
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
	# for BrakerGff in $(ls gene_pred/braker/F.*/*_braker_new/*/augustus.gff3 | grep -w -e 'Fus2'); do
for BrakerGff in $(ls gene_pred/braker/F.*/*/*/augustus.gff3 | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi'| grep -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139' | grep -e 'ncbi' | grep 'N139'); do
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
for DirPath in $(ls -d gene_pred/final_genes/F.*/*/final | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi'| grep -e '125' -e 'A23' -e 'A13' -e 'A28' -e 'CB3' -e 'PG' -e 'A8' -e 'N139'  | grep -e 'ncbi'); do
echo $DirPath;
cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
echo "";
done
```
```
gene_pred/final_genes/F.oxysporum_fsp_cepae/125_ncbi/final
17197
1546
18743

gene_pred/final_genes/F.oxysporum_fsp_cepae/A13_ncbi/final
17986
948
18934

gene_pred/final_genes/F.oxysporum_fsp_cepae/A23_ncbi/final
17006
1551
18557

gene_pred/final_genes/F.oxysporum_fsp_cepae/A28_ncbi/final
17426
1448
18874

gene_pred/final_genes/F.oxysporum_fsp_cepae/CB3_ncbi/final
16833
1345
18178

gene_pred/final_genes/F.oxysporum_fsp_cepae/PG_ncbi/final
16811
1290
18101

gene_pred/final_genes/F.oxysporum_fsp_narcissi/N139_ncbi/final
19059
1642
20701

gene_pred/final_genes/F.proliferatum/A8_ncbi/final
15421
37
15458
```


```bash
for Gff in $(ls gene_pred/final_genes/*/*/final/final_genes_appended.gff3 | grep -e 'canu_new' -e 'ncbi'); do
	Strain=$(echo $Gff | rev | cut -d '/' -f3 | rev)
	Organism=$(echo $Gff | rev | cut -d '/' -f4 | rev)
	echo "$Strain - $Organism"
	cat $Gff | grep -w 'gene' | wc -l
done
```

```
125_ncbi - F.oxysporum_fsp_cepae
18505
A13_ncbi - F.oxysporum_fsp_cepae
18790
A23_ncbi - F.oxysporum_fsp_cepae
18323
A28_ncbi - F.oxysporum_fsp_cepae
18629
CB3_ncbi - F.oxysporum_fsp_cepae
17943
Fus2_canu_new - F.oxysporum_fsp_cepae
18855
PG_ncbi - F.oxysporum_fsp_cepae
17830
N139_ncbi - F.oxysporum_fsp_narcissi
20493
A8_ncbi - F.proliferatum
15418
```

## Identification of duplicated genes in additional CodingQuary gene models

```bash
for AddGenes in $(ls gene_pred/codingquary/F.*/*/additional/additional_genes.gff | grep 'ncbi' | grep -v 'ncbi_braker' | grep 'A23'); do
Strain=$(echo $AddGenes| rev | cut -d '/' -f3 | rev)
Organism=$(echo $AddGenes | rev | cut -d '/' -f4 | rev)
OutDir=$(dirname $AddGenes)
echo "$Organism - $Strain" > $OutDir/duplicated_genes.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
$ProgDir/remove_dup_features.py --inp_gff $AddGenes >> $OutDir/duplicated_genes.txt
cat $OutDir/duplicated_genes.txt
echo ""
done
``` -->
