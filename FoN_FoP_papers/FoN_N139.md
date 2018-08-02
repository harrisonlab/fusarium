# FoN N139 publication commands


Most commands in assembly and annotation of the N139 genome are found in the
main README.md document for this repository (along with commands used in the
FoC paper). Those commands documented here are those that deviate from that
analysis.


## Removal of GPI anchor containing proteins from secreted proteins

Proteins containing GPI anchors were also removed using GPIsom


These proteins were identified through submitting the combined protein file to
the webserver at: http://gpi.unibe.ch


An output directory was made to download the file to:
"GPI anchored (C&N-term signal) (SignalP):"

```bash
  for Proteome in $(ls gene_pred/final_genes/*/*/*/final_genes_combined.pep.fasta | grep 'N139_ncbi'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/trans_mem/$Organism/$Strain/GPIsom
    mkdir -p $OutDir
  done
```

The link to the proteins results was followed embedded in the text of
"GPI anchored (C&N-term signal) (SignalP):" thes +GPI anchor proteins
were copied and pasted into:

<!-- The link to the results was followed embedded in the text of
"Seqs with C-terminal signal (GPI-SOM):" thes +GPI anchor proteins
were copied and pasted into: -->

```bash
  nano gene_pred/trans_mem/$Organism/$Strain/GPIsom/GPI_pos.fa
```

Those proteins with GPI anchors were removed from lists of Signal
peptide containing proteins

```bash
for File in $(ls gene_pred/trans_mem/*/*/GPIsom/GPI_pos.fa | grep '10300'); do
Strain=$(echo $File | rev | cut -f3 -d '/' | rev)
Organism=$(echo $File | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
TmHeaders=$(echo "$File" | sed 's/.fa/.txt/g')
cat $File | grep '>' | cut -f1 -d ' ' | sed 's/>//g' > $TmHeaders
SigP=$(ls gene_pred/combined_sigP/$Organism/$Strain/*_sp_no_trans_mem.aa)
SigPHeaders=gene_pred/combined_sigP/$Organism/$Strain/"$Strain"_sp_no_trans_mem_headers.txt
cat $SigP | grep '>' | cut -f1 | sed 's/>//g'> $SigPHeaders
GoodHeaders=$(echo "$File" | sed 's/_pos.fa/_neg.txt/g')
cat $SigPHeaders | grep -v -f $TmHeaders > $GoodHeaders
OutDir=$(dirname $SigP)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
# cat $SigP | grep -v -A1 -f $TmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem_no_GPI.aa
$ProgDir/extract_from_fasta.py --fasta $SigP --headers $GoodHeaders  > $OutDir/"$Strain"_final_sp_no_trans_mem_no_GPI.aa
echo "Number of SigP proteins:"
cat $SigP | grep '>' |wc -l
echo "Number with GPI anchors in entire proteome:"
cat $TmHeaders | wc -l
echo "Number without GPI anchors:"
cat $OutDir/"$Strain"_final_sp_no_trans_mem_no_GPI.aa | grep '>' | wc -l
echo "Number of gene models:"
cat $OutDir/"$Strain"_final_sp_no_trans_mem_no_GPI.aa | grep '>' | cut -f1 -d't' | sort | uniq |wc -l
done
```

```
  P.cactorum - 10300
  Number of SigP proteins:
  2244
  Number with GPI anchors in entire proteome:
  483
  Number without GPI anchors:
  1984
  Number of gene models:
  1973
```



## D) Secondary metabolites (Antismash and SMURF)

Antismash was run to identify clusters of secondary metabolite genes within
the genome. Antismash was run using the weserver at:
http://antismash.secondarymetabolites.org


Results of web-annotation of gene clusters within the assembly were downloaded to
the following directories:

```bash
  for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi' | grep -e 'canu_new' -e 'ncbi' | grep 'N139'); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    OutDir=analysis/secondary_metabolites/antismash/$Organism/$Strain
    mkdir -p $OutDir
  done
```

```bash
  for Zip in $(ls analysis/secondary_metabolites/antismash/*/*/*.zip); do
    OutDir=$(dirname $Zip)
    unzip -d $OutDir $Zip
  done
```

```bash
  for AntiSmash in $(ls analysis/secondary_metabolites/antismash/*/*/*/*.final.gbk); do
    Organism=$(echo $AntiSmash | rev | cut -f4 -d '/' | rev)
    Strain=$(echo $AntiSmash | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/secondary_metabolites/antismash/$Organism/$Strain
    Prefix=$OutDir/WT_antismash
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
    $ProgDir/antismash2gff.py --inp_antismash $AntiSmash --out_prefix $Prefix

    # Identify secondary metabolites within predicted clusters
    printf "Number of secondary metabolite detected:\t"
    cat "$Prefix"_secmet_clusters.gff | wc -l
    GeneGff=gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3
    bedtools intersect -u -a $GeneGff -b "$Prefix"_secmet_clusters.gff > "$Prefix"_secmet_genes.gff
    cat "$Prefix"_secmet_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' > "$Prefix"_antismash_secmet_genes.txt
    bedtools intersect -wo -a $GeneGff -b "$Prefix"_secmet_clusters.gff | grep 'mRNA' | cut -f9,10,12,18 | sed "s/ID=//g" | perl -p -i -e "s/;Parent=g\w+//g" | perl -p -i -e "s/;Notes=.*//g" > "$Prefix"_secmet_genes.tsv
    printf "Number of predicted proteins in secondary metabolite clusters:\t"
    cat "$Prefix"_secmet_genes.txt | wc -l
    printf "Number of predicted genes in secondary metabolite clusters:\t"
    cat "$Prefix"_secmet_genes.gff | grep -w 'gene' | wc -l

      # Identify cluster finder additional non-secondary metabolite clusters
      printf "Number of cluster finder non-SecMet clusters detected:\t"
      cat "$Prefix"_clusterfinder_clusters.gff | wc -l
      GeneGff=gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3
      bedtools intersect -u -a $GeneGff -b "$Prefix"_clusterfinder_clusters.gff > "$Prefix"_clusterfinder_genes.gff
      cat "$Prefix"_clusterfinder_genes.gff | grep -w 'mRNA' | cut -f9 | cut -f2 -d '=' | cut -f1 -d ';' > "$Prefix"_clusterfinder_genes.txt

      printf "Number of predicted proteins in cluster finder non-SecMet clusters:\t"
      cat "$Prefix"_clusterfinder_genes.txt | wc -l
      printf "Number of predicted genes in cluster finder non-SecMet clusters:\t"
      cat "$Prefix"_clusterfinder_genes.gff | grep -w 'gene' | wc -l
  done
```

These clusters represented the following genes. Note that these numbers just
show the number of intersected genes with gff clusters and are not confirmed by
function

```
```

SMURF was also run to identify secondary metabolite gene clusters.

Genes needed to be parsed into a specific tsv format prior to submission on the
SMURF webserver.

```bash
for Gff in $(ls gene_pred/final_genes/F.oxysporum_fsp_narcissi/N139_ncbi/final/final_genes_appended.gff3); do
	Organism=$(echo $Gff | rev | cut -f4 -d '/' | rev)
	Strain=$(echo $Gff | rev | cut -f3 -d '/' | rev)
	echo "$Organism - $Strain"
  OutDir=analysis/secondary_metabolites/smurf/$Organism/$Strain
  mkdir -p $OutDir
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
  $ProgDir/gff2smurf.py --gff $Gff > $OutDir/"$Strain"_genes_smurf.tsv
done
```

SMURF results were observed and found to be trustworthy as the highly fragmented
assembly meant that the order of gene models was not interpreted correctly.

<!--
SMURF output was received by email and downloaded to the cluster in the output
directory above.

Output files were parsed into gff format:

```bash
  OutDir=analysis/secondary_metabolites/smurf/F.venenatum/WT
  Prefix="WT"
  GeneGff=gene_pred/final/F.venenatum/WT/final/final_genes_appended_renamed.gff3
  SmurfClusters=$OutDir/Secondary-Metabolite-Clusters.txt
  SmurfBackbone=$OutDir/Backbone-genes.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/secondary_metabolites
  $ProgDir/smurf2gff.py --smurf_clusters $SmurfClusters --smurf_backbone $SmurfBackbone > $OutDir/Smurf_clusters.gff
  bedtools intersect -wo -a $GeneGff -b $OutDir/Smurf_clusters.gff | grep 'mRNA' | cut -f9,10,12,18 | sed "s/ID=//g" | perl -p -i -e "s/;Parent=g\w+//g" | perl -p -i -e "s/;Notes=.*//g" > $OutDir/"$Prefix"_smurf_secmet_genes.tsv
```

Total number of secondary metabolite clusters:

```bash
for Assembly in $(ls repeat_masked/*/*/illumina_assembly_ncbi/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep -w 'WT'); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
OutDir=analysis/secondary_metabolites/antismash/$Organism/$Strain
mkdir -p $OutDir
GeneGff=$(ls gene_pred/final/$Organism/$Strain/final/final_genes_appended_renamed.gff3)
AntismashClusters=$(ls analysis/secondary_metabolites/antismash/$Organism/$Strain/*_secmet_clusters.gff)
SmurfClusters=$(ls analysis/secondary_metabolites/smurf/$Organism/$Strain/Smurf_clusters.gff)
echo "Total number of Antismash clusters"
cat $AntismashClusters | wc -l
echo "Total number of SMURF clusters"
cat $SmurfClusters | wc -l
echo "number of Antismash clusters intersecting Smurf clusters"
bedtools intersect -a $AntismashClusters -b $SmurfClusters | wc -l
echo "number of Antismash clusters not intersecting Smurf clusters"
bedtools intersect -v -a $AntismashClusters -b $SmurfClusters | wc -l
echo "number of smurf clusters intersecting antismash clusters"
bedtools intersect -a $SmurfClusters -b $AntismashClusters | wc -l
echo "number of smurf clusters not intersecting antismash clusters"
bedtools intersect -v -a $SmurfClusters -b $AntismashClusters | wc -l
done
``` -->


# F) Genes with transcription factor annotations:


A list of PFAM domains, superfamily annotations used as part of the DBD database
and a further set of interproscan annotations listed by Shelest et al 2017 were made
http://www.transcriptionfactor.org/index.cgi?Domain+domain:all
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5415576/

```bash
  for Interpro in $(ls gene_pred/interproscan/*/*/*_interproscan.tsv | grep 'N139' | grep 'ncbi'); do
    Organism=$(echo $Interpro | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $Interpro | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=analysis/transcription_factors/$Organism/$Strain
    mkdir -p $OutDir
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transcription_factors
    $ProgDir/interpro2TFs.py --InterPro $Interpro > $OutDir/"$Strain"_TF_domains.tsv
    echo "total number of transcription factors"
    cat $OutDir/"$Strain"_TF_domains.tsv | cut -f1 | sort | uniq > $OutDir/"$Strain"_TF_gene_headers.txt
    cat $OutDir/"$Strain"_TF_gene_headers.txt | wc -l
  done
```


# Distribution of effectors through FoL contigs:

```bash
EffectorP=$(ls analysis/effectorP/F.oxysporum_fsp_narcissi/N139_ncbi/F.oxysporum_fsp_narcissi_N139_ncbi_EffectorP_secreted.gff)
CAZY=$(ls gene_pred/CAZY/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_CAZY_secreted.gff)
ContigList=$(ls less analysis/circos/F.oxysporum_fsp_narcissi/FoN_FoL_v2/FoN_FoL_linked_contigs.txt)

for FoL_contig in $(cat $ContigList | cut -f1 | sort | uniq); do
  ContigTotal="0"
  # echo $FoL_contig
  for FoN_contig in $(cat $ContigList | grep -w "$FoL_contig" | cut -f2); do
    # echo "\t"$FoN_contig
    NumEffP=$(cat $EffectorP | grep -w 'gene' | sed 's/contig/FoN_contig/g' | grep -w "$FoN_contig" | wc -l)
    # printf "\t$FoL_contig\t$FoN_contig\t$NumEffP\n"
    ContigTotal=$(($NumEffP + $ContigTotal))
  done
  printf "$FoL_contig\t$ContigTotal\n"
done > tmp.txt

for FoL_contig in $(cat $ContigList | cut -f1 | sort | uniq); do
  ContigTotal="0"
  # echo $FoL_contig
  for FoN_contig in $(cat $ContigList | grep -w "$FoL_contig" | cut -f2); do
    # echo "\t"$FoN_contig
    NumEffP=$(cat $CAZY | grep -w 'gene' | sed 's/contig/FoN_contig/g' | grep -w "$FoN_contig" | wc -l)
    # printf "\t$FoL_contig\t$FoN_contig\t$NumEffP\n"
    ContigTotal=$(($NumEffP + $ContigTotal))
  done
  printf "$FoL_contig\t$ContigTotal\n"
done > tmp2.txt

cat $EffectorP | grep -w 'gene' | wc -l
cat $EffectorP | grep -w 'gene' | cut -f1 | sort | uniq -c | sort -n -r
cat $EffectorP | grep -w 'gene' | cut -f1 | sort | uniq  | sed 's/contig/FoN_contig/g'> tmp.txt
cat $ContigList | grep -w -f tmp.txt
cat $ContigList | grep -w -f tmp.txt | sort | less -S
```







## Annotation Table


linked contigs were identified:


```bash


```

## G) Summarising annotation in annotation table

```bash
for GeneGff in $(ls gene_pred/final_genes/F.oxysporum_fsp_narcissi/N139_ncbi/final/final_genes_appended.gff3); do
Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_unmasked.fa | grep 'ncbi')
GeneConversions=$(ls genome_submission/$Organism/$Strain/gag/edited/genome_gene_conversions.tsv)
Antismash=$(ls analysis/secondary_metabolites/antismash/F.oxysporum_fsp_narcissi/N139_ncbi/*/geneclusters.txt)
# Smurf=$(ls analysis/secondary_metabolites/smurf/F.venenatum/WT/WT_smurf_secmet_genes.tsv)
TFs=$(ls analysis/transcription_factors/$Organism/$Strain/"$Strain"_TF_domains.tsv)
SigP=$(ls gene_pred/final_genes_signalp-4.1/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_final_sp.aa)
TM_out=$(ls gene_pred/trans_mem/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_TM_genes_pos.txt)
GPI_out=$(ls gene_pred/trans_mem/$Organism/$Strain/GPIsom/GPI_pos.fa)
MIMP_list=$(ls analysis/mimps/$Organism/$Strain/"$Strain"_genes_in_2kb_mimp.txt)
EffP_list=$(ls analysis/effectorP/$Organism/$Strain/"$Organism"_"$Strain"_EffectorP_headers.txt)
CAZY_list=$(ls gene_pred/CAZY/$Organism/$Strain/"$Strain"_CAZY_headers.txt)
InterPro=$(ls gene_pred/interproscan/$Organism/$Strain/*_interproscan.tsv)
SwissProt=$(ls gene_pred/swissprot/$Organism/$Strain/swissprot_vJul2016_tophit_parsed.tbl)
Orthology=$(ls analysis/orthology/orthomcl/FoN_vs_FoC_vs_FoL_vs_Fo/FoN_vs_FoC_vs_FoL_vs_Fo_orthogroups.txt)
OrthoStrainID='FoN'
OrthoStrainAll='FoN FoL FoC fo47'
OutDir=gene_pred/annotation/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes
# $ProgDir/build_annot_Fv.py --genome $Assembly --genes_gff $GeneGff --Antismash $Antismash --Smurf $Smurf --vitamins $Vitamins --TFs $TFs --InterPro $InterPro --Swissprot $SwissProt --orthogroups_PH1 $PH1_orthology --orthogroups_GR1 $GR1_orthology > $OutDir/"$Strain"_annotation_ncbi.tsv
$ProgDir/FoN_build_gene_annot_table.py --genome $Assembly --genes_gff $GeneGff --genes_renamed $GeneConversions --Antismash $Antismash --TFs $TFs --SigP $SigP --TM_list $TM_out --GPI_list $GPI_out --MIMP_list $MIMP_list --EffP_list $EffP_list --CAZY_list $CAZY_list --InterPro $InterPro --Swissprot $SwissProt --orthogroups $Orthology --strain_id $OrthoStrainID --OrthoMCL_all $OrthoStrainAll > $OutDir/"$Strain"_annotation_ncbi.tsv
done
```


Subset the table for secreted genes within 2 Kb of a mimp
Also extract fasta sequence for these genes.

```bash
for GeneGff in $(ls gene_pred/final_genes/F.oxysporum_fsp_narcissi/N139_ncbi/final/final_genes_appended.gff3); do
Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
OutDir=gene_pred/annotation/$Organism/$Strain
cat $OutDir/"$Strain"_annotation_ncbi.tsv | grep 'MIMP' > $OutDir/"$Strain"_annotation_ncbi_MIMP.tsv
cat $OutDir/"$Strain"_annotation_ncbi_MIMP.tsv | grep -v 'ncbi_gene_id' | wc -l
cat $OutDir/"$Strain"_annotation_ncbi.tsv | grep 'MIMP' | grep 'SigP' | grep -v 'GPI' | grep -v 'TM' > $OutDir/"$Strain"_annotation_ncbi_MIMP_secreted.tsv
cat $OutDir/"$Strain"_annotation_ncbi_MIMP_secreted.tsv | grep -v 'ncbi_gene_id' | wc -l
cat $OutDir/"$Strain"_annotation_ncbi_MIMP_secreted.tsv | grep -v 'ncbi_gene_id' | cut -f2 > $OutDir/"$Strain"_annotation_ncbi_MIMP_secreted_headers.txt
CDS=$(ls gene_pred/final_genes/F.oxysporum_fsp_narcissi/N139_ncbi/final/final_genes_combined.cdna.fasta)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $CDS --headers $OutDir/"$Strain"_annotation_ncbi_MIMP_secreted_headers.txt  > $OutDir/"$Strain"_annotation_ncbi_MIMP_secreted.fasta
done
```

Blast these genes vs other FoN genomes:

```bash
for RefGenome in $(ls repeat_masked/F.oxysporum_fsp_narcissi/*/*/*_contigs_unmasked.fa); do
Prefix=$(echo $RefGenome | cut -f2,3 -d '/' --output-delimiter '_')
echo $Prefix
QueryFasta=$(ls gene_pred/annotation/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_annotation_ncbi_MIMP_secreted.fasta)
OutDir=analysis/blast_homology/FoN_blast/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
rm ${Prefix}_genome.fa
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $QueryFasta dna $RefGenome $OutDir
done

for RefGenome in $(ls assembly/external_group/F.oxysporum_fsp_narcissi/Na5/GCA_002233775_1/NJCV01.1.fsa_nt); do
Prefix=$(echo $RefGenome | cut -f3,4 -d '/' --output-delimiter '_')
echo $Prefix
QueryFasta=$(ls gene_pred/annotation/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_annotation_ncbi_MIMP_secreted.fasta)
OutDir=analysis/blast_homology/FoN_blast/vs_ref_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
rm ${Prefix}_genome.fa
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $QueryFasta dna $RefGenome $OutDir
done

for RefGenome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa | grep -e 'Fus2' -e 'fo47' -e '4287'| grep -v -e 'old' -e 'broad' -e 'tgac' -e 'chromosomal' | grep '4287'); do
Prefix=$(echo $RefGenome | cut -f2,3 -d '/' --output-delimiter '_')
echo $Prefix
QueryFasta=$(ls gene_pred/annotation/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_annotation_ncbi_MIMP_secreted.fasta)
OutDir=analysis/blast_homology/FoN_blast/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
rm ${Prefix}_genome.fa
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $QueryFasta dna $RefGenome $OutDir
done

```


## Summarise blast hists

```bash
CsvFiles=$(ls analysis/blast_homology/FoN_blast/vs_*_genomes/*/*_hits.csv)
Headers=$(echo $CsvFiles | sed 's&analysis/blast_homology/FoN_blast/vs_ref_genomes/&&g' | sed 's&analysis/blast_homology/FoN_blast/vs_seq_genomes/&&g' | sed -r "s&/\w+?&&g" | sed 's/.fasta_hits.csv//g')
OutDir=analysis/blast_homology/FoN_blast/extracted
mkdir -p $OutDir
Genomes=$(ls analysis/blast_homology/FoN_blast/vs_*_genomes/*/*_genome.fa)
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/AHDB_project/blast_searches
$ProgDir/blast_parse_AHDB.py --blast_csv $CsvFiles --headers $Headers --genomes $Genomes --identity 0.50 --evalue 1e-30 --out_prefix $OutDir/FoN_mimp_secreted
```




## 5.1.B) Identifying FTF genes

Previously published FTF genes from Sanchez et al 2016 were blasted against
Fusarium genomes.

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -v -e 'HB17' | grep -v -w 'Fus2_canu_new' | grep 'N139_ncbi' | grep -v 'old'); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo $Assembly
Query=analysis/blast_homology/Fo_path_genes/FTF_cds_Sanchez_et_al_2016.fasta
OutDir=analysis/FTF/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query dna $Assembly $OutDir
done
```

BLAST hits were converted to Gff annotations and intersected with gene models:

```bash
for BlastHits in $(ls analysis/FTF/*/*/*_FTF_cds_Sanchez_et_al_2016.fasta_hits.csv | grep -v -w 'Fus2_canu_new' | grep 'N139_ncbi' | grep -v 'old'); do
Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
OutDir=analysis/FTF/$Organism/$Strain
HitsGff=$(echo $BlastHits | sed  's/.csv/.gff/g')
Column2=FTF_homolog
NumHits=1
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff

GffAppended=$(ls gene_pred/final_genes/$Organism/$Strain/final/final_genes_appended.gff3)
bedtools intersect -wao -a $HitsGff -b $GffAppended > $OutDir/"$Strain"_FTF_hits_intersected.bed
done
```

Fus2 genes g16859 and g10474 were identified as FTF1 and FTF2 homologs.

```bash
cat $OutDir/"$Strain"_FTF_hits_intersected.bed | grep -w 'gene' | cut -f18 | cut -f1 -d ';' | sort | uniq
```

Following identification of other FTF genes from other Fuarium spp. and building
a phylogeny the tree was exported from geneious and visualised using:

```r
setwd("/Users/armita/Downloads/FoN/analysis/FTF")
#===============================================================================
#       Load libraries
#===============================================================================

library(ape)
library(ggplot2)

library(ggtree)
library(phangorn)
library(treeio)

tree <- read.tree("geneious_NJ_FTF_tree.newick")

mydata <- read.csv("traits.csv", stringsAsFactors=FALSE)
rownames(mydata) <- mydata$tiplabel
mydata <- mydata[match(tree$tip.label,rownames(mydata)),]

t <- ggtree(tree, aes(linetype=nodes$support)) # Core tree
# Adjust terminal branch lengths:
# branches <- t$data
#tree$edge.length[branches$isTip] <- 1
#Tree <- branches$branch.length
#rescale_tree(t, branches$branch.length)

t <- t + geom_treescale(offset=-1.0, fontsize = 3) # Add scalebar
# t <- t + xlim(0, 0.025) # Add more space for labels



# Colouring labels by values in another df
t <- t %<+% mydata # Allow colouring of nodes by another df
#t <- t + geom_tiplab(aes(color=Source), size=3, hjust=0) +
scale_color_manual(values=c("gray39","black")) # colours as defined by col2rgb

tips <- data.frame(t$data)
tips$label <- tips$newlabel
t <- t + geom_tiplab(data=tips, aes(color=study), size=3, hjust=0, offset = 0.01) +
scale_color_manual(values=c("gray39","black")) # colours as defined by col2rgb

# Format nodes by values
nodes <- data.frame(t$data)
#nodes <- nodes[!nodes$isTip,]
nodes$label[nodes$isTip] <- ''
nodes$label[(!nodes$isTip) & (nodes$label == '')] <- 100

nodes$label <- as.numeric(nodes$label)
nodes$label <- lapply(nodes$label,round,0)

#nodes$label[nodes$label < 0.80] <- ''
nodes$support[nodes$isTip] <- 'supported'
# nodes$support[(!nodes$isTip) & (nodes$label == '100')] <- 'supported'
nodes$support[(!nodes$isTip) & (nodes$label > 80)] <- 'supported'
nodes$support[(!nodes$isTip) & (nodes$label < 80)] <- 'unsupported'
#nodes$support[(!nodes$isTip) & (nodes$label == '')] <- 'supported'
t <- t + aes(linetype=nodes$support)
nodes$label[nodes$label > 80] <- ''
t <- t + geom_nodelab(data=nodes, size=2, hjust=-0.05) # colours as defined by col2rgb


# Annotate a clade with a bar line
t <- t + geom_cladelabel(node=55, label='FTF1', align=T, colour='black', offset=0.14)
t <- t + geom_cladelabel(node=86, label='FTF2', align=T, colour='black', offset=0.14)

# Save as PDF and force a 'huge' size plot
t <- ggsave("ftf.pdf", width =20, height = 30, units = "cm", limitsize = FALSE)



```

## 5.2 Identifying PHIbase homologs
<!--
The PHIbase database was searched against the assembled genomes using tBLASTx.

```bash
for Assembly in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa| grep -v 'HB17' | grep -e 'cepae' -e 'proliferatum' -e 'narcissi' | grep -e 'Fus2_canu_new' -e 'ncbi' | grep -v 'old'); do
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/blast_pipe.sh analysis/blast_homology/PHIbase/PHI_36_accessions.fa protein $Assembly
done
```

following blasting PHIbase to the genome, the hits were filtered by effect on
virulence.

First the a tab seperated file was made in the clusters core directory containing
PHIbase. These commands were run as part of previous projects but have been
included here for completeness. -->
<!--
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
```
-->

The analysis was also performed by blasting the predicted proteins against the
PHIbase database:

The PHIbase database was searched agasinst the assembled genomes using tBLASTx.

```bash
for Proteome in $(ls gene_pred/final_genes/*/*/*/final_genes_combined.gene.fasta | grep 'N139_ncbi'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh ../../phibase/v4.4/phi_accessions.fa protein $Proteome $OutDir
# qsub /home/armita/git_repos/emr_repos/tools/pathogen/blast/blast_pipe.sh $Proteome protein ../../phibase/v4.4/phi_accessions.fa
done
```

```bash
PhibaseDir=/home/groups/harrisonlab/phibase/v4.4
PhibaseHeaders=$PhibaseDir/PHI_headers.csv
PhibaseVirulence=$PhibaseDir/PHI_virulence.csv
for BlastCSV in $(ls analysis/blast_homology/F*/*/*_phi_accessions.fa_homologs.csv | grep 'N139_ncbi'); do
Strain=$(echo $BlastCSV | rev | cut -f2 -d'/' | rev)
echo "$Strain"
OutDir=$(dirname $BlastCSV)
paste -d '\t' $PhibaseHeaders $PhibaseVirulence $BlastCSV | cut -f-3,1185- > $OutDir/"$Strain"_PHIbase_virulence.csv
cat $OutDir/"$Strain"_PHIbase_virulence.csv | cut -f2 | sort | uniq -c | tee $OutDir/"$Strain"_PHIbase_virulence.txt
done
```



9 Transcription factor genes (TF) have been identified in FoL that potentially
regulate SIX gene expression. The genes in FoL can be used to identify orthologous
FoC genes:

```bash
TF_list="FOXG_14257 FOXG_17260 FOXG_17266 FOXG_14201 FOXG_14230 FOXG_14211 FOXG_14275 FOXG_14277 FOXG_14274"
num=0
OrthoTxt=$(ls /home/groups/harrisonlab/project_files/fusarium/analysis/orthology/orthomcl/FoN_vs_FoC_vs_FoL_vs_Fo/FoN_vs_FoC_vs_FoL_vs_Fo_orthogroups.txt)
AnnotTab=$(ls gene_pred/annotation/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_annotation_ncbi.tsv)
for TF in $TF_list; do
num=$(($num +1))
echo "TF"$num
Orthogroup=$(cat $OrthoTxt | grep -w "$TF" | cut -f1 -d ':')
cat $AnnotTab | grep -w "$Orthogroup"
done > analysis/transcription_factors/FoN_TF1-9.tsv
```
