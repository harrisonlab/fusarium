## Preparing genome and .conf files

```bash
OutDir=analysis/circos/F.oxysporum_fsp_lactucae/R1_R4
mkdir -p $OutDir

R1_genome=$(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep 'AJ520')

ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $R1_genome --contig_prefix "R1_" > $OutDir/R1_genome.txt

R4_genome=$(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep 'AJ516')
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $R4_genome --contig_prefix "R4_" > $OutDir/R4_genome.txt

  cat $OutDir/R1_genome.txt > $OutDir/R1_R4_genome.txt
  tac $OutDir/R4_genome.txt >> $OutDir/R1_R4_genome.txt

  cat $OutDir/R1_R4_genome.txt > $OutDir/R1_R4_genome_edited.txt
```

```bash
OutDir=analysis/circos/F.oxysporum_fsp_lactucae/R1_R4
Coords=$(ls analysis/genome_alignment/mummer/F.oxysporum_fsp_lactucae/AJ520/AJ520_vs_AJ516/AJ520_vs_AJ516_coords.tsv)
ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/genome_alignment
$ProgDir/nucmer_coords2circos.py --inp_coords $Coords --queery_id R1 --ref_id R4 > $OutDir/R1_vs_R4_links.txt
cat $OutDir/R1_vs_R4_links.txt > $OutDir/R1_vs_R4_links_edited.txt
cat $OutDir/R1_vs_R4_links.txt | sed "s/R4_..//g" | sed 's/\..//g' | sort -k4,4n -k5,5n -h > $OutDir/R1_vs_R4_links_edited2.txt
```
<!--
A file showing contig orientations was made:
```bash
  cat $OutDir/R1_vs_R4_links.txt | cut -f1,4,5 | grep 'CM' | sed 's/R4_CM//g' | sed 's/\..//g' | sort -k2,2n -k3,3n -h | cut -f1 | uniq | sed "s/$/; /g" | tr -d '\n' > $OutDir/R1_contig_order.txt
  cat $OutDir/R1_vs_R4_links.txt | cut -f1,4,5 | grep -v 'CM' | sed 's/R4_CM//g' | sed 's/\..//g' | sort -k2,2n -k3,3n -h | cut -f1 | uniq | sed "s/$/; /g" | tr -d '\n' >> $OutDir/R1_contig_order.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/find_contig_orientation.py --links_file $OutDir/R1_vs_R4_links_edited2.txt > $OutDir/R1_vs_R4_contig_orientation.txt
```


linked contigs were identified:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
$ProgDir/show_linked_contigs.py --links_file $OutDir/R1_vs_R4_links_edited.txt > $OutDir/R1_R4_linked_contigs.txt
```

Chromosome break points, in reference to R4 contigs were identified:

```bash
 cat $OutDir/R1_vs_R4_links_edited2.txt | awk '!arr[$4] {arr[$4]=$0; if(prevline) print prevline; print} {prevline=$0}' | less
```

The number of bp in syntenous contigs was identified using:

```bash
  cat $OutDir/R1_vs_R4_contig_orientation.txt | tail -n3 | grep -v 'orientation' | sed 's/, /\n/g' > $OutDir/R1_syntenous_contigs.txt
  echo "Total bp in syntenous contigs:"
  # cat $OutDir/R1_genome.txt | grep -v -e "0 .... chr" -e "0 ... chr" | grep -wf $OutDir/R1_syntenous_contigs.txt | cut -f6 -d ' ' | awk '{s+=$1} END {print s}'
  cat $OutDir/R1_genome.txt | grep -wf $OutDir/R1_syntenous_contigs.txt | cut -f6 -d ' ' | awk '{s+=$1} END {print s}'
  echo "This is shown in the following number of contigs:"
  # cat $OutDir/R1_genome.txt | grep -v -e "0 .... chr" -e "0 ... chr" | grep -wf $OutDir/R1_syntenous_contigs.txt | cut -f6 -d ' ' | wc -l
  cat $OutDir/R1_genome.txt | grep -wf $OutDir/R1_syntenous_contigs.txt | cut -f6 -d ' ' | wc -l
  echo "The ramining contigs represent the following number of bp:"
  cat $OutDir/R1_genome.txt | grep -wvf $OutDir/R1_syntenous_contigs.txt | cut -f6 -d ' ' | awk '{s+=$1} END {print s}'
  cat $OutDir/R1_genome.txt | grep -wvf $OutDir/R1_syntenous_contigs.txt | cut -f3 -d ' ' > $OutDir/R1_unplaced_contigs.txt
  cat $OutDir/R1_unplaced_contigs.txt | wc -l
```


```
Total bp in syntenous contigs:
57386451
This is shown in the following number of contigs:
69
The ramining contigs represent the following number of bp:
3280976
36
```

Contig order was selected by taking the first line of that file and then also
taking the reversed order of R4 contigs using the command:

```bash
# all contigs

cat $OutDir/R1_R4_genome_edited.txt | grep 'R1' | cut -f3 -d ' '| tr -d '\n' | sed 's/R1/, R1/g' | sed 's/R4/, R4/g' > tmp3.txt
cat $OutDir/R1_R4_genome_edited.txt | grep 'R4' | cut -f3 -d ' ' | tr -d '\n' | sed 's/R4/, R4/g'  >> tmp3.txt
# Forward orientation
cat $OutDir/R1_vs_R4_contig_orientation.txt | grep -A1 'Order of all seen contigs'
cat $OutDir/R1_vs_R4_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" > tmp.txt
# reverse orientation
# cat $OutDir/R1_R4_genome_edited.txt | grep 'R1' | grep -v -e "0 .... chr" -e "0 ... chr" | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/R1/, R1/g'
cat $OutDir/R1_R4_genome_edited.txt | grep 'R1' | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/R1/, R1/g'
# reference sequences (put in reverse)
# contig order
cat $OutDir/R1_vs_R4_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1
cat $OutDir/R1_R4_genome_edited.txt | grep 'R4' | cut -f3 -d ' ' | tr -d '\n' | sed 's/R4/, R4/g'
# unplaced contigs
echo "The following contigs are unplaced:"
cat $OutDir/R1_R4_genome_edited.txt | grep -vw -f tmp.txt | grep 'R1' | cut -f3 -d ' ' | tr -d '\n' | sed 's/R1/, R1/g' > $OutDir/R1_R4_contigs_unplaced.txt

# cat $OutDir/R1_genome.txt | grep -wvf $OutDir/R1_syntenous_contigs.txt | cut -f3 -d ' ' | tr -d '\n' | sed "s/R1/, R1/g" > $OutDir/R1_R4_contigs_unplaced.txt
```

Contig orientation was used to edit the circos .conf file manually
 -->

## Preparing Effector plots

# Plot location of R1 mimps and secreted effectorP genes as a scatterplot

```bash
  GffMimp=$(ls analysis/mimps/F.oxysporum_fsp_lactucae/AJ520/AJ520_mimps.gff)
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffMimp --feature MIMP_motif --value 1 | sed -e 's/^/R1_/g' > $OutDir/R1_mimp_plot.txt
  GffCAZY=$(ls gene_pred/CAZY/F.oxysporum_fsp_lactucae/AJ520/*_CAZY_secreted.gff)
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffCAZY --feature gene --value 1 | sed -e 's/^/R1_/g' > $OutDir/R1_CAZY_plot.txt
  GffEfR1=$(ls analysis/effectorP/F.oxysporum_fsp_lactucae/AJ520/*_EffectorP_secreted.gff)
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffEfR1 --feature gene --value 1 | sed -e 's/^/R1_/g' > $OutDir/R1_effectorP_plot.txt

  BlastHits=$(ls analysis/blast_homology/F.oxysporum_fsp_lactucae/AJ520/*_six-appended_parsed.fa_homologs.gff)
  GffSix=$OutDir/R1_SIX.gff
  cat $BlastHits | grep -v -e 'MIMP' -e 'C5' -e 'CRX' | grep 'SIX' > $GffSix
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffSix --feature BLAST_hit --value 1 | sed -e 's/^/R1_/g' > $OutDir/R1_SIX_plot.txt
```


# Plot location of R4 mimps and secreted effectorP genes as a scatterplot

```bash
  GffMimp=$(ls analysis/mimps/F.oxysporum_fsp_lactucae/AJ516/AJ516_mimps.gff)
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffMimp --feature MIMP_motif --value 1 | sed -e 's/^/R4_/g' > $OutDir/R4_mimp_plot.txt
  GffCAZY=$(ls gene_pred/CAZY/F.oxysporum_fsp_lactucae/AJ516/*_CAZY_secreted.gff)
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffCAZY --feature gene --value 1 | sed -e 's/^/R4_/g' > $OutDir/R4_CAZY_plot.txt
  GffEfR4=$(ls analysis/effectorP/F.oxysporum_fsp_lactucae/AJ516/*_EffectorP_secreted.gff)
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffEfR4 --feature gene --value 1 | sed -e 's/^/R4_/g' > $OutDir/R4_effectorP_plot.txt

  BlastHits=$(ls analysis/blast_homology/F.oxysporum_fsp_lactucae/AJ516/*_six-appended_parsed.fa_homologs.gff)
  GffSix=$OutDir/R4_SIX.gff
  cat $BlastHits | grep -v -e 'MIMP' -e 'C5' -e 'CRX' | grep 'SIX' > $GffSix
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffSix --feature BLAST_hit --value 1 | sed -e 's/^/R4_/g' > $OutDir/R4_SIX_plot.txt
```

```bash
OutDir=analysis/circos/F.oxysporum_fsp_lactucae/R1_R4
Orthogroups=$(ls analysis/orthology/orthomcl/Fo_lactucae_FoL_Fo_FoC/Fo_lactucae_FoL_Fo_FoC_orthogroups.txt)
# Orthogroups=$(ls analysis/orthology/orthomcl/Fo_lactucae_FoL_Fo_FoC/formatted/Results_Oct09/Orthogroups.txt)
echo "This Number of orthogroups unique to FoLaR1:"
cat $Orthogroups | grep -e 'FoLaR1|' | grep -v -e 'FoC|' -e 'FoLy|' -e 'fo47|' -e 'FoLaR4' | wc -l
echo "This Number of orthogroups unique to FoLaR4:"
cat $Orthogroups | grep -e 'FoLaR4|' | grep -v -e 'FoC|' -e 'FoLy|' -e 'fo47|' -e 'FoLaR1' | wc -l
echo "This Number of orthogroups unique and present in both Fo fsp Lactucae:"
cat $Orthogroups | grep -e 'FoLaR1|' | grep -e 'FoLaR4|' | grep -v -e 'FoC|' -e 'FoLy|' -e 'fo47|' | wc -l
echo "Extracting these genes:"
cat $Orthogroups | grep -e 'FoLaR1|' | grep -e 'FoLaR4|' | grep -v -e 'FoC|' -e 'FoLy|' -e 'fo47|' | cut -f2- -d ' ' | sed "s/ /\n/g" | sed "s/FoLaR1|/R1_/g" | sed "s/FoLaR4|/R4_/g" > $OutDir/Fo_lactucae_unique_genes.txt
cat $Orthogroups | grep -e 'FoLaR1|' | grep -v -e 'FoLaR4|' -e 'FoC|' -e 'FoLy|' -e 'fo47|' | cut -f2- -d ' ' | sed "s/ /\n/g" | sed "s/FoLaR1|//g" > $OutDir/R1_unique_genes.txt
cat $OutDir/R1_unique_genes.txt | wc -l
cat $Orthogroups | grep -e 'FoLaR4|' | grep -v -e 'FoLaR1|' -e 'FoC|' -e 'FoLy|' -e 'fo47|' | cut -f2- -d ' ' | sed "s/ /\n/g" | sed "s/FoLaR4|//g" > $OutDir/R4_unique_genes.txt
cat $OutDir/R4_unique_genes.txt | wc -l
```

# Identification of Fo lactucae specific orthogroups:

```bash
cat $OutDir/Fo_lactucae_unique_genes.txt | grep 'R1' | sed 's/R1_//g' > $OutDir/Fo_lactucae_unique_genes_R1.txt
cat $OutDir/Fo_lactucae_unique_genes.txt | grep 'R4' | sed 's/R4_//g' > $OutDir/Fo_lactucae_unique_genes_R4.txt


Gff=$(ls gene_pred/final/F.oxysporum_fsp_lactucae/AJ520/publication/final_genes_appended_renamed_ncbi.gff3)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $OutDir/Fo_lactucae_unique_genes_R1.txt $Gff orthogroup ID > $OutDir/Fo_lactucae_unique_genes_R1.gff
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/gff2circos_scatterplot.py --gff $OutDir/Fo_lactucae_unique_genes_R1.gff --feature gene --value 1 | sed -e 's/^/R1_/g' > $OutDir/Fo_lactucae_unique_genes_R1_plot.txt

Gff=$(ls gene_pred/final/F.oxysporum_fsp_lactucae/AJ516/publication/final_genes_appended_renamed_ncbi.gff3)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl $OutDir/Fo_lactucae_unique_genes_R4.txt $Gff orthogroup ID > $OutDir/Fo_lactucae_unique_genes_R4.gff
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/gff2circos_scatterplot.py --gff $OutDir/Fo_lactucae_unique_genes_R4.gff --feature gene --value 1 | sed -e 's/^/R4_/g' > $OutDir/Fo_lactucae_unique_genes_R4_plot.txt

cat $OutDir/Fo_lactucae_unique_genes_R1_plot.txt $OutDir/Fo_lactucae_unique_genes_R4_plot.txt > $OutDir/Fo_lactucae_unique_genes_R1_R4_plot.txt
```

# Identification of unique genes for each race:

```bash
  Gff=$(ls gene_pred/final/F.oxysporum_fsp_lactucae/AJ520/publication/final_genes_appended_renamed_ncbi.gff3)
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_gff_for_sigP_hits.pl $OutDir/R1_unique_genes.txt $Gff orthogroup ID > $OutDir/R1_unique_genes.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $OutDir/R1_unique_genes.gff --feature gene --value 1 | sed -e 's/^/R1_/g' > $OutDir/R1_unique_genes_plot.txt

  Gff=$(ls gene_pred/final/F.oxysporum_fsp_lactucae/AJ516/publication/final_genes_appended_renamed_ncbi.gff3)
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_gff_for_sigP_hits.pl $OutDir/R4_unique_genes.txt $Gff orthogroup ID > $OutDir/R4_unique_genes.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $OutDir/R4_unique_genes.gff --feature gene --value 1 | sed -e 's/^/R4_/g' > $OutDir/R4_unique_genes_plot.txt

cat $OutDir/R1_unique_genes_plot.txt $OutDir/R4_unique_genes_plot.txt > $OutDir/R1_R4_unique_genes_plot.txt
```

Combine R1 and R4 effector scatterplots

```bash
cat $OutDir/R1_mimp_plot.txt $OutDir/R4_mimp_plot.txt > $OutDir/R1_R4_mimp_plot.txt
cat $OutDir/R1_CAZY_plot.txt $OutDir/R4_CAZY_plot.txt > $OutDir/R1_R4_CAZY_plot.txt
cat $OutDir/R1_effectorP_plot.txt $OutDir/R4_effectorP_plot.txt > $OutDir/R1_R4_effectorP_plot.txt
cat $OutDir/R1_SIX_plot.txt $OutDir/R4_SIX_plot.txt > $OutDir/R1_R4_SIX_plot.txt
```
<!--
```bash
OutDir=analysis/circos/F.oxysporum_fsp_lactucae/R1_R4

ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos/FoLactucae/R1_vs_R4
$ProgDir/order_Fo_contigs.py --karotype $OutDir/R1_R4_genome_edited.txt --links $OutDir/R1_vs_R4_links.txt
```
-->

## Running circos



```bash
OutDir=analysis/circos/F.oxysporum_fsp_lactucae/R1_R4
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
circos -conf $ProgDir/FoLactucae/R1_vs_R4/R1_R4_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/R1_R4_circos_full.png
mv $OutDir/circos.svg $OutDir/R1_R4_circos_full.svg
ls $PWD/$OutDir/R1_R4_circos_full.png
```

```bash
OutDir=analysis/circos/F.oxysporum_fsp_lactucae/R1_R4
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
circos -conf $ProgDir/FoLactucae/R1_vs_R4/R1_R4_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/R1_R4_circos_full_orthogroups.png
mv $OutDir/circos.svg $OutDir/R1_R4_circos_full_orthogroups.svg
ls $PWD/$OutDir/R1_R4_circos_full_orthogroups.png
```


# Further analysis of non-syntenous regions

The number of MIMPs and effectors in LS regions were identified:


```bash
echo "Number of MIMPs in syntenous regions:"
cat $OutDir/R1_mimp_plot.txt | grep -w -v -f $OutDir/R1_unplaced_contigs.txt | wc -l
cat $OutDir/R1_mimp_plot.txt | grep -w -v -f $OutDir/R1_unplaced_contigs.txt
echo "Number of MIMPs in unplaced regions:"
cat $OutDir/R1_mimp_plot.txt | grep -w -f $OutDir/R1_unplaced_contigs.txt | wc -l
```

```
Number of MIMPs in syntenous regions:
10
Number of MIMPs in unplaced regions:
197
```


The 10 MIMPs in syntenous regions were in:

```
R1_contig_16	34774	34789	1.0
R1_contig_19	394088	394103	1.0
R1_contig_19	395353	395368	1.0
R1_contig_32	188261	188276	1.0
R1_contig_172	26824	26839	1.0
R1_contig_172	62081	62096	1.0
R1_contig_172	27026	27041	1.0
R1_contig_172	62277	62292	1.0
R1_contig_196	44262	44277	1.0
R1_contig_232	11510	11525	1.0
```

```bash
cat $OutDir/R1_R4_linked_contigs.txt | grep -w -e 'R1_contig_16' -e 'R1_contig_19' -e 'R1_contig_32' -e 'R1_contig_172' -e 'R1_contig_196' -e 'R1_contig_232'
```

```
R4_CM000590.1	R1_contig_19 (chr2 x2)
R4_CM000597.1	R1_contig_196 (chr9 x1)
R4_CM000598.1	R1_contig_232 (chr10 x1)
R4_CM000599.1	R1_contig_32 (chr11 x1)
R4_CM000600.1	R1_contig_16 (chr12 x1)
R4_DS231739.1	R1_contig_172 (unplaced x4)
```
