
<!--
```bash
WarwickCSV=analysis/expression/warwick/F.oxysporum_fsp_cepae/Fus2/23_06/path_vs_non_p0.001.csv
WarwickTSV=analysis/expression/warwick/F.oxysporum_fsp_cepae/Fus2/23_06/path_vs_non_p0.001.tab
cat $WarwickCSV | sed 's/,/\t/g' | tr -d '"' > $WarwickTSV

OrthogroupsTxt=$(ls analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/FoC_vs_Fo_vs_FoL_orthogroups.txt)
DEG_Orthogroups=analysis/expression/warwick/F.oxysporum_fsp_cepae/Fus2/23_06/Fus2_path_vs_non_path_orthogroups.tab

ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/orthology
$ProgDir/diff_expressed_orthogroups.py \
  --RNAseq_tab $WarwickTSV \
  --FoC_orthogroup $OrthogroupsTxt \
  --OrthoMCL_id Fus2 \
  --OrthoMCL_all Fus2 125 A23 A28 D2 PG A13 fo47 55 A1_2 CB3 HB6 4287 \
  --OrthoMCL_path Fus2 125 A23 --OrthoMCL_nonpath A28 D2 PG A13 fo47 \
  > $DEG_Orthogroups
``` -->


```bash
  # for Strain in 125 A23 Fus2_canu_new 55 A1-2 CB3 HB6 A13 A28 D2 PG; do
# for Strain in 125 A23 Fus2_canu_new CB3 A13 A28 PG; do
for Strain in 125_ncbi A23_ncbi Fus2_canu_new CB3_ncbi A13_ncbi A28_ncbi PG_ncbi; do
# for Strain in Fus2_canu_new; do
for GeneGff in $(ls gene_pred/final_genes/*/$Strain/final/final_genes_appended.gff3); do
Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Genome=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_unmasked.fa)
BlastCsv=$(ls analysis/blast_homology/$Organism/$Strain/4287_chromosomal_final_genes_combined.pep.fasta_hits.csv)
FolIntersect=$(ls analysis/blast_homology/$Organism/$Strain/4287_chromosomal_final_genes_combined_intersect.bed)
GeneGff=$(ls gene_pred/final_genes/$Organism/$Strain/final/final_genes_appended.gff3)
SigpTab=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp.tab)
TmhmmTxt=$(ls gene_pred/trans_mem/$Organism/$Strain/*_tmhmm_out.txt)
MimpTxt=$(ls analysis/mimps/$Organism/$Strain/*_prots_in_2kb_mimp.txt)
EffectorpTxt=$(ls analysis/effectorP/$Organism/$Strain/*_EffectorP.txt)
OrthogroupsTxt=$(ls analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL_publication_ncbi/FoC_vs_Fo_vs_FoL_publication_ncbi_orthogroups.txt)
InterProTsv=$(ls gene_pred/interproscan/$Organism/$Strain/*_interproscan.tsv)
SwissprotTab=$(ls gene_pred/swissprot/$Organism/$Strain/swissprot_vJul2016_tophit_parsed.tbl)
# DEG_Orthogroups=$(ls analysis/expression/warwick/F.oxysporum_fsp_cepae/Fus2/23_06/Fus2_path_vs_non_path_orthogroups.tab)

OrthoMCL_id="$Strain"
OrthoMCL_id_list="125 A23 Fus2 CB3 A13 A28 PG fo47 4287"
OrthoMCL_path_ids="125 A23 Fus2"
OrthoMCL_nonpath_ids="A13 A28 CB3 PG fo47"

if [ "$Strain" == 'Fus2_canu_new' ]; then OrthoMCL_id="Fus2"; fi
# if [ "$Strain" == 'A1-2' ]; then OrthoMCL_id="A1_2"; fi

OutDir=gene_pred/annotations/$Organism/$Strain
OutTable=$OutDir/"$Strain"_gene_annotations.tab

mkdir -p $OutDir

ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes
$ProgDir/Fo_build_gene_annot_table.py \
--blast_csv $BlastCsv \
--FoL_intersected_genes $FolIntersect \
--genome $Genome \
--FoC_genes_gff $GeneGff \
--FoC_SigP $SigpTab \
--FoC_TM_list $TmhmmTxt \
--FoC_MIMP_list $MimpTxt \
--FoC_effectorP $EffectorpTxt \
--FoC_orthogroup $OrthogroupsTxt \
--OrthoMCL_id $OrthoMCL_id \
--OrthoMCL_all $OrthoMCL_id_list \
--OrthoMCL_path $OrthoMCL_path_ids \
--OrthoMCL_nonpath $OrthoMCL_nonpath_ids \
--InterPro $InterProTsv \
--Swissprot $SwissprotTab \
> $OutTable
# --DEG_Orthogroups $DEG_Orthogroups \
# > $OutTable
done
done
```

Gene tables were made for Fp A8

```bash

for Strain in A8_ncbi; do
for GeneGff in $(ls gene_pred/final_genes/*/$Strain/final/final_genes_appended.gff3); do
Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
Strain=$(echo $GeneGff | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
Genome=$(ls repeat_masked/$Organism/$Strain/*/*_contigs_unmasked.fa)
BlastCsv=$(ls analysis/blast_homology/$Organism/$Strain/4287_chromosomal_final_genes_combined.pep.fasta_hits.csv)
FolIntersect=$(ls analysis/blast_homology/$Organism/$Strain/4287_chromosomal_final_genes_combined_intersect.bed)
GeneGff=$(ls gene_pred/final_genes/$Organism/$Strain/final/final_genes_appended.gff3)
SigpTab=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp.tab)
TmhmmTxt=$(ls gene_pred/trans_mem/$Organism/$Strain/*_tmhmm_out.txt)
MimpTxt=$(ls analysis/mimps/$Organism/$Strain/*_prots_in_2kb_mimp.txt)
EffectorpTxt=$(ls analysis/effectorP/$Organism/$Strain/*_EffectorP.txt)
OrthogroupsTxt=$(ls analysis/orthology/orthomcl/Fp_Fv_FoC_FoL_Fo/Fp_Fv_FoC_FoL_Fo_orthogroups.txt)
InterProTsv=$(ls gene_pred/interproscan/$Organism/$Strain/*_interproscan.tsv)
SwissprotTab=$(ls gene_pred/swissprot/$Organism/$Strain/swissprot_vJul2016_tophit_parsed.tbl)
# DEG_Orthogroups=$(ls analysis/expression/warwick/F.oxysporum_fsp_cepae/Fus2/23_06/Fus2_path_vs_non_path_orthogroups.tab)

OrthoMCL_id=Fp
OrthoMCL_id_list="Fp Fv Fo FoC FoL"
OrthoMCL_path_ids="Fo FoC FoL"
OrthoMCL_nonpath_ids="Fp Fv"

if [ "$Strain" == 'Fus2_canu_new' ]; then OrthoMCL_id="Fus2"; fi
# if [ "$Strain" == 'A1-2' ]; then OrthoMCL_id="A1_2"; fi

OutDir=gene_pred/annotations/$Organism/$Strain
OutTable=$OutDir/"$Strain"_gene_annotations.tab

mkdir -p $OutDir

ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes
$ProgDir/Fo_build_gene_annot_table.py \
--blast_csv $BlastCsv \
--FoL_intersected_genes $FolIntersect \
--genome $Genome \
--FoC_genes_gff $GeneGff \
--FoC_SigP $SigpTab \
--FoC_TM_list $TmhmmTxt \
--FoC_MIMP_list $MimpTxt \
--FoC_effectorP $EffectorpTxt \
--FoC_orthogroup $OrthogroupsTxt \
--OrthoMCL_id $OrthoMCL_id \
--OrthoMCL_all $OrthoMCL_id_list \
--OrthoMCL_path $OrthoMCL_path_ids \
--OrthoMCL_nonpath $OrthoMCL_nonpath_ids \
--InterPro $InterProTsv \
--Swissprot $SwissprotTab \
> $OutTable
# --DEG_Orthogroups $DEG_Orthogroups \
# > $OutTable
done
done
```

Gene tables were made for Fo fo47 and FoL 4287

```bash
for Strain in fo47 4287; do
if [ $Strain == '4287' ]; then
GeneGff=$(ls assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31_parsed.gff3)
Genome=$(ls assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.dna.chromosome_parsed.fa)
BlastCsv=$(ls analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287/4287_chromosomal_Fusox1_GeneCatalog_proteins_20110522_parsed.fa_hits.csv)
FolIntersect=$(ls analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287/4287_chromosomal_final_genes_combined_intersect.bed)
elif [ $Strain == 'fo47' ]; then
GeneGff=$(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_transcripts_parsed.gff3)
Genome=$(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs_parsed.fasta)
BlastCsv=$(ls analysis/blast_homology/F.oxysporum/fo47/4287_chromosomal_final_genes_combined.pep.fasta_hits.csv)
FolIntersect=$(ls analysis/blast_homology/F.oxysporum/fo47/4287_chromosomal_final_genes_combined_intersect.bed)
fi
Organism=$(echo $GeneGff | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
SigpTab=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain*/*_final_sp.tab)
TmhmmTxt=$(ls gene_pred/trans_mem/$Organism/$Strain*/*_tmhmm_out.txt)
MimpTxt=$(ls analysis/mimps/$Organism/$Strain*/*_genes_in_2kb_mimp.txt)
EffectorpTxt=$(ls analysis/effectorP/$Organism/$Strain*/*_EffectorP.txt)
OrthogroupsTxt=$(ls analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL_publication_ncbi/FoC_vs_Fo_vs_FoL_publication_ncbi_orthogroups.txt)
InterProTsv=$(ls gene_pred/interproscan/$Organism/$Strain*/*_interproscan.tsv)
SwissprotTab=$(ls gene_pred/swissprot/$Organism/$Strain*/swissprot_vJul2016_tophit_parsed.tbl)
# DEG_Orthogroups=$(ls analysis/expression/warwick/F.oxysporum_fsp_cepae/Fus2/23_06/Fus2_path_vs_non_path_orthogroups.tab)

OrthoMCL_id="$Strain"
OrthoMCL_id_list="125 A23 Fus2 CB3 A13 A28 PG fo47 4287"
OrthoMCL_path_ids="125 A23 Fus2"
OrthoMCL_nonpath_ids="A13 A28 CB3 PG fo47"

# if [ "$Strain" == 'Fus2_edited_v2' ]; then OrthoMCL_id="Fus2"; fi
# if [ "$Strain" == 'A1-2' ]; then OrthoMCL_id="A1_2"; fi

OutDir=gene_pred/annotations/$Organism/$Strain
OutTable=$OutDir/"$Strain"_gene_annotations.tab

mkdir -p $OutDir

ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes
$ProgDir/Fo_build_gene_annot_table.py \
--blast_csv $BlastCsv \
--FoL_intersected_genes $FolIntersect \
--genome $Genome \
--FoC_genes_gff $GeneGff \
--FoC_SigP $SigpTab \
--FoC_TM_list $TmhmmTxt \
--FoC_MIMP_list $MimpTxt \
--FoC_effectorP $EffectorpTxt \
--FoC_orthogroup $OrthogroupsTxt \
--OrthoMCL_id $OrthoMCL_id \
--OrthoMCL_all $OrthoMCL_id_list \
--OrthoMCL_path $OrthoMCL_path_ids \
--OrthoMCL_nonpath $OrthoMCL_nonpath_ids \
--InterPro $InterProTsv \
--Swissprot $SwissprotTab \
> $OutTable
# --DEG_Orthogroups $DEG_Orthogroups \
# > $OutTable
done
```




Gene tables were analysed, identifying the number of secreted proteins, expanded
gene families and genes within 2kb of a mimp.

Note - the fo47 annotation table is still a little buggy and as a result may be
producing spurious results. Data for this table should be extracted manually
from excel.

```bash
# Number of secreted proteins
  for File in $(ls gene_pred/annotations/*/*/*_gene_annotations.tab); do
    echo $(basename $File);
    cat $File | cut -f12 | grep -i 'Yes' | wc -l;
  done
  # Number of secreted proteins that are effectorP +ve
  for File in $(ls gene_pred/annotations/F*/*/*_gene_annotations.tab); do
    echo $(basename $File);
    cat $File | cut -f12,14 | grep -i "Yes.*Yes" | wc -l;
  done
  # Number of secreted proteins within 2kb of a mimp
  for File in $(ls gene_pred/annotations/F*/*/*_gene_annotations.tab); do
    echo $(basename $File);
    cat $File | cut -f12,13 | grep -i "Yes.*Yes" | wc -l;
  done
  #  Genes in expanded gene families in pathogens
  for File in $(ls gene_pred/annotations/F*/*/*_gene_annotations.tab); do
    echo $(basename $File);
    cat $File | cut -f12,13,16,18,19,32 | grep -P "\tpathogen_expanded" | wc -l
  done
  # Secreted genes in expanded gene families in pathogens
  for File in $(ls gene_pred/annotations/F*/*/*_gene_annotations.tab); do
    echo $(basename $File);
    cat $File | cut -f12,16,18,19,33 | grep -i "Yes" | grep -P "\tpathogen_expanded" | wc -l
  done
  #  Genes in expanded gene families in non-pathogens
  for File in $(ls gene_pred/annotations/F*/*/*_gene_annotations.tab | grep 'PG'); do
    echo $(basename $File);
    cat $File | cut -f12,13,16,18,19,32 | grep -w "non-pathogen_expanded"
  done | less -S
  # Secreted genes in expanded gene families in non-pathogens
  for File in $(ls gene_pred/annotations/F*/PG/*_gene_annotations.tab); do
    echo $(basename $File);
    cat $File | cut -f12,16,18,19,34 | grep -i "Yes" | grep -w 'non-pathogen_expanded' | less -S;
  done
  # Number of secreted proteins within 2kb of a mimp and in expanded orthogroups
  for File in $(ls gene_pred/annotations/F*/*/*_gene_annotations.tab); do
    echo $(basename $File);
    cat $File | cut -f12,13,19,34 | grep -i "Yes.*Yes" | grep 'expanded' | wc -l;
  done
  for File in $(ls gene_pred/annotations/F*/Fus2_canu_new/*_gene_annotations.tab); do
    echo $(basename $File);
    cat $File | cut -f1,12,13,16,18,19,20,34 | grep -i "Yes.*Yes" | grep 'expanded';
    echo ""
  done > tmp4/annot.tab
```


The number of pathogen expanded orthogroups was identified:

```bash
  cat gene_pred/annotations/*/*/*_gene_annotations.tab | cut -f 16,18,19 | grep -P "\tpathogen_expanded" | cut -f1 | sort | uniq | wc -l
  cat gene_pred/annotations/*/*/*_gene_annotations.tab | cut -f 16,17,18,19,34 | grep -P "\tpathogen_expanded" | grep -w 'path_isolates_all' | cut -f1 | sort | uniq | wc -l
  cat gene_pred/annotations/*/*/*_gene_annotations.tab | cut -f 16,17,18,19,34 | grep -P "\tpathogen_expanded" | grep -w 'path_isolates_all' | grep -v '4287(0)' | cut -f1 | sort | uniq | wc -l
```


## annotations in LS regions:

```bash
AnnotTable=gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_gene_annotations.tab
OutDir=$(dirname $AnnotTable)
cat $AnnotTable | grep -w -e 'contig_10_pilon' -e 'contig_14_pilon' -e 'contig_16_pilon' -e 'contig_19_pilon' -e 'contig_20_pilon' -e 'contig_21_pilon' -e 'contig_22_pilon' > $OutDir/Fus2_canu_new_LS._specific_regions.tab
cat $AnnotTable | grep -w -e 'contig_9_pilon' -e 'contig_11_pilon' -e 'contig_12_pilon' -e 'contig_13_pilon' -e 'contig_15_pilon' -e 'contig_17_pilon' > $OutDir/Fus2_canu_new_sp._specific_regions.tab
cat $AnnotTable | grep -w -e 'contig_1_pilon' -e 'contig_2_pilon' -e 'contig_3_pilon' -e 'contig_4_pilon' -e 'contig_5_pilon' -e 'contig_6_pilon' -e 'contig_7_pilon' -e 'contig_8_pilon' > $OutDir/Fus2_canu_new_core_regions.tab
```

## Identification of transcription factors

```bash

mkdir analysis/transcription_factors
cat gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_LS._specific_regions.tab | grep -i 'transcription' | grep -v -i -e 'transcriptional activator' | grep -i -e 'transcription factor'  | sort -k16 > analysis/transcription_factors/Fus2_canu_new_LS_specific_transcrtiption_factors.tab

cat analysis/transcription_factors/Fus2_canu_new_LS_specific_transcrtiption_factors.tab | cut -f2 | sort | uniq -c
cat analysis/transcription_factors/Fus2_canu_new_LS_specific_transcrtiption_factors.tab | cut -f1 > analysis/transcription_factors/Fus2_canu_new_LS_specific_transcrtiption_factors_headers.txt
echo "The number of transcription factors on LS contigs 10, 14, 16, 19, 20, 21 and 22 is:"
cat analysis/transcription_factors/Fus2_canu_new_LS_specific_transcrtiption_factors.tab | wc -l
echo "These are present in the following number of orthogroups:"
cat analysis/transcription_factors/Fus2_canu_new_LS_specific_transcrtiption_factors.tab | cut -f16 | sort | uniq -c | sort -r -n | wc -l
```

```
7 contig_10_pilon
8 contig_14_pilon
2 contig_16_pilon
3 contig_19_pilon
6 contig_20_pilon
2 contig_21_pilon
2 contig_22_pilon

  The number of transcription factors on LS contigs 10, 14, 16, 18, 19, 20 & 21 is:
  30
  These are present in the following number of orthogroups:
  23
```

Expression of each transcription factor was shown using:

```bash
TransFactors=$(ls analysis/transcription_factors/Fus2_canu_new_LS_specific_transcrtiption_factors.tab)
ExpressionAnnotation=$(ls analysis/expression/Fus2_expressed_genes.tsv)
for TransFactor in $(cat $TransFactors | cut -f1); do
  # echo $TransFactor
  cat $ExpressionAnnotation | grep "$TransFactor"
done > analysis/transcription_factors/Fus2_canu_new_LS_specific_transcrtiption_factors_fpkm.tab

TransFactors=$(ls analysis/transcription_factors/Fus2_canu_new_TF1-9.tsv)
ExpressionAnnotation=$(ls analysis/expression/Fus2_expressed_genes.tsv)
for TransFactor in $(cat $TransFactors | grep 'contig' | cut -f1); do
  cat $ExpressionAnnotation | grep "$TransFactor"
done > analysis/transcription_factors/Fus2_canu_new_TF1-9_fpkm.tab

```
