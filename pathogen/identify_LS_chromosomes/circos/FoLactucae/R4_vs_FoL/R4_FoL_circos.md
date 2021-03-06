## Preparing genome and .conf files

```bash
OutDir=analysis/circos/F.oxysporum_fsp_lactucae/R4_FoL
mkdir -p $OutDir

R4_genome=$(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep 'AJ516')

ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $R4_genome --contig_prefix "R4_" > $OutDir/R4_genome.txt

FoL_genome=repeat_masked/F.oxysporum_fsp_lycopersici/4287_v2/fungidb_repmask/4287_v2_contigs_unmasked.fa
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $FoL_genome --contig_prefix "FoL_" > $OutDir/FoL_genome.txt

cat $OutDir/R4_genome.txt > $OutDir/R4_FoL_genome.txt
tac $OutDir/FoL_genome.txt >> $OutDir/R4_FoL_genome.txt

cat $OutDir/R4_FoL_genome.txt > $OutDir/R4_FoL_genome_edited.txt
```


```bash
OutDir=analysis/circos/F.oxysporum_fsp_lactucae/R4_FoL
Coords=$(ls analysis/genome_alignment/mummer/F.oxysporum_fsp_lactucae/AJ516/AJ516_vs_4287/AJ516_vs_4287_coords.tsv)
ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/genome_alignment
$ProgDir/nucmer_coords2circos.py --inp_coords $Coords --queery_id R4 --ref_id FoL > $OutDir/R4_vs_FoL_links.txt
cat $OutDir/R4_vs_FoL_links.txt > $OutDir/R4_vs_FoL_links_edited.txt
cat $OutDir/R4_vs_FoL_links.txt | sed "s/FoL_..//g" | sed 's/\..//g' | sort -k4,4n -k5,5n -h > $OutDir/R4_vs_FoL_links_edited2.txt
```

A file showing contig orientations was made:
```bash
  cat $OutDir/R4_vs_FoL_links.txt | cut -f1,4,5 | grep 'CM' | sed 's/FoL_CM//g' | sed 's/\..//g' | sort -k2,2n -k3,3n -h | cut -f1 | uniq | sed "s/$/; /g" | tr -d '\n' > $OutDir/R4_contig_order.txt
  cat $OutDir/R4_vs_FoL_links.txt | cut -f1,4,5 | grep -v 'CM' | sed 's/FoL_CM//g' | sed 's/\..//g' | sort -k2,2n -k3,3n -h | cut -f1 | uniq | sed "s/$/; /g" | tr -d '\n' >> $OutDir/R4_contig_order.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/find_contig_orientation.py --links_file $OutDir/R4_vs_FoL_links_edited2.txt > $OutDir/R4_vs_FoL_contig_orientation.txt
```


linked contigs were identified:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
$ProgDir/show_linked_contigs.py --links_file $OutDir/R4_vs_FoL_links_edited.txt > $OutDir/R4_FoL_linked_contigs.txt
```

Chromosome break points, in reference to FoL contigs were identified:

```bash
 cat $OutDir/R4_vs_FoL_links_edited2.txt | awk '!arr[$4] {arr[$4]=$0; if(prevline) print prevline; print} {prevline=$0}' | less
```

The number of bp in syntenous contigs was identified using:

```bash
  cat $OutDir/R4_vs_FoL_contig_orientation.txt | tail -n3 | grep -v 'orientation' | sed 's/, /\n/g' > $OutDir/R4_syntenous_contigs.txt
  echo "Total bp in syntenous contigs:"
  # cat $OutDir/R4_genome.txt | grep -v -e "0 .... chr" -e "0 ... chr" | grep -wf $OutDir/R4_syntenous_contigs.txt | cut -f6 -d ' ' | awk '{s+=$1} END {print s}'
  cat $OutDir/R4_genome.txt | grep -wf $OutDir/R4_syntenous_contigs.txt | cut -f6 -d ' ' | awk '{s+=$1} END {print s}'
  echo "This is shown in the following number of contigs:"
  # cat $OutDir/R4_genome.txt | grep -v -e "0 .... chr" -e "0 ... chr" | grep -wf $OutDir/R4_syntenous_contigs.txt | cut -f6 -d ' ' | wc -l
  cat $OutDir/R4_genome.txt | grep -wf $OutDir/R4_syntenous_contigs.txt | cut -f6 -d ' ' | wc -l
  echo "The ramining contigs represent the following number of bp:"
  cat $OutDir/R4_genome.txt | grep -wvf $OutDir/R4_syntenous_contigs.txt | cut -f6 -d ' ' | awk '{s+=$1} END {print s}'
  cat $OutDir/R4_genome.txt | grep -wvf $OutDir/R4_syntenous_contigs.txt | cut -f3 -d ' ' > $OutDir/R4_unplaced_contigs.txt
  cat $OutDir/R4_unplaced_contigs.txt | wc -l
```


```
Total bp in syntenous contigs:
58398043
This is shown in the following number of contigs:
56
The ramining contigs represent the following number of bp:
8547037
38
```

Contig order was selected by taking the first line of that file and then also
taking the reversed order of FoL contigs using the command:

```bash
# all contigs
# cat $OutDir/R4_FoL_genome_edited.txt | grep 'R4' | grep -v -e "0 .... chr" -e "0 ... chr" | cut -f3 -d ' '| tr -d '\n' | sed 's/R4/, R4/g' | sed 's/FoL/, FoL/g' > tmp3.txt
cat $OutDir/R4_FoL_genome_edited.txt | grep 'R4' | cut -f3 -d ' '| tr -d '\n' | sed 's/R4/, R4/g' | sed 's/FoL/, FoL/g' > tmp3.txt
cat $OutDir/R4_FoL_genome_edited.txt | grep 'FoL' | cut -f3 -d ' ' | tr -d '\n' | sed 's/FoL/, FoL/g'  >> tmp3.txt
# Forward orientation
cat $OutDir/R4_vs_FoL_contig_orientation.txt | grep -A1 'Order of all seen contigs'
cat $OutDir/R4_vs_FoL_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" > tmp.txt
# reverse orientation
# cat $OutDir/R4_FoL_genome_edited.txt | grep 'R4' | grep -v -e "0 .... chr" -e "0 ... chr" | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/R4/, R4/g'
cat $OutDir/R4_FoL_genome_edited.txt | grep 'R4' | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/R4/, R4/g'
# reference sequences (put in reverse)
# contig order
cat $OutDir/R4_vs_FoL_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1
cat $OutDir/R4_FoL_genome_edited.txt | grep 'FoL' | cut -f3 -d ' ' | tr -d '\n' | sed 's/FoL/, FoL/g'
# unplaced contigs
echo "The following contigs are unplaced:"
cat $OutDir/R4_FoL_genome_edited.txt | grep -vw -f tmp.txt | grep 'R4' | cut -f3 -d ' ' | tr -d '\n' | sed 's/R4/, R4/g' > $OutDir/R4_FoL_contigs_unplaced.txt

# cat $OutDir/R4_genome.txt | grep -wvf $OutDir/R4_syntenous_contigs.txt | cut -f3 -d ' ' | tr -d '\n' | sed "s/R4/, R4/g" > $OutDir/R4_FoL_contigs_unplaced.txt
```

Contig orientation was used to edit the circos .conf file manually


## Preparing Effector plots

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
  # GffAntiSmash=analysis/antismash/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_secondary_metabolite_regions.gff
  # ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  # $ProgDir/gff2circos_scatterplot.py --gff $GffAntiSmash --feature indole indole-nrps nrps nrps-t1pks other t1pks t1pks-nrps t3pks terpene --value 1 | sed -e 's/^/R4_/g' > $OutDir/R4_antismash_plot.txt

  BlastHits=$(ls analysis/blast_homology/F.oxysporum_fsp_lactucae/AJ516/*_six-appended_parsed.fa_homologs.gff)
  GffSix=$OutDir/R4_SIX.gff
  cat $BlastHits | grep -v -e 'MIMP' -e 'C5' -e 'CRX' | grep 'SIX' > $GffSix
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffSix --feature BLAST_hit --value 1 | sed -e 's/^/R4_/g' > $OutDir/R4_SIX_plot.txt
```

```bash
OutDir=analysis/circos/F.oxysporum_fsp_lactucae/R4_FoL

ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos/FoLactucae/R1_vs_FoL
$ProgDir/order_Fo_contigs.py --karotype $OutDir/R4_FoL_genome_edited.txt --links $OutDir/R4_vs_FoL_links.txt

```

## Running circos



```bash
OutDir=analysis/circos/F.oxysporum_fsp_lactucae/R4_FoL
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
circos -conf $ProgDir/FoLactucae/R4_vs_FoL/R4_FoL_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/R4_FoL_circos_full.png
mv $OutDir/circos.svg $OutDir/R4_FoL_circos_full.svg
ls $PWD/$OutDir/R4_FoL_circos_full.png
```


# Further analysis of non-syntenous regions regions

The number of MIMPs and effectors in LS regions were identified:


```bash
echo "Number of MIMPs in syntenous regions:"
cat $OutDir/R4_mimp_plot.txt | grep -w -v -f $OutDir/R4_unplaced_contigs.txt | wc -l
cat $OutDir/R4_mimp_plot.txt | grep -w -v -f $OutDir/R4_unplaced_contigs.txt
echo "Number of MIMPs in unplaced regions:"
cat $OutDir/R4_mimp_plot.txt | grep -w -f $OutDir/R4_unplaced_contigs.txt | wc -l
```

```
Number of MIMPs in syntenous regions:
10
Number of MIMPs in unplaced regions:
197
```


The 10 MIMPs in syntenous regions were in:

```
R4_contig_16	34774	34789	1.0
R4_contig_19	394088	394103	1.0
R4_contig_19	395353	395368	1.0
R4_contig_32	188261	188276	1.0
R4_contig_172	26824	26839	1.0
R4_contig_172	62081	62096	1.0
R4_contig_172	27026	27041	1.0
R4_contig_172	62277	62292	1.0
R4_contig_196	44262	44277	1.0
R4_contig_232	11510	11525	1.0
```

```bash
cat $OutDir/R4_FoL_linked_contigs.txt | grep -w -e 'R4_contig_16' -e 'R4_contig_19' -e 'R4_contig_32' -e 'R4_contig_172' -e 'R4_contig_196' -e 'R4_contig_232'
```

```
FoL_CM000590.1	R4_contig_19 (chr2 x2)
FoL_CM000597.1	R4_contig_196 (chr9 x1)
FoL_CM000598.1	R4_contig_232 (chR40 x1)
FoL_CM000599.1	R4_contig_32 (chR41 x1)
FoL_CM000600.1	R4_contig_16 (chR42 x1)
FoL_DS231739.1	R4_contig_172 (unplaced x4)
```
<!--

## Preparing genome and .conf files

```bash
OutDir=analysis/circos/F.oxysporum_fsp_narcissi/R4_FoL
mkdir -p $OutDir

R4_genome=repeat_masked/F.oxysporum_fsp_narcissi/N139_ncbi/ncbi_submission/R4_contigs_unmasked.fa
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $R4_genome --contig_prefix "R4_" > $OutDir/R4_genome.txt

FoL_genome=repeat_masked/F.oxysporum_fsp_lycopersici/4287_v2/fungidb_repmask/4287_v2_contigs_unmasked.fa
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $FoL_genome --contig_prefix "FoL_" > $OutDir/FoL_genome.txt

  cat $OutDir/R4_genome.txt > $OutDir/R4_FoL_genome.txt
  tac $OutDir/FoL_genome.txt >> $OutDir/R4_FoL_genome.txt

  # COntigs smaller than 10Kb were removed
  cat $OutDir/R4_FoL_genome.txt | grep -v 'DS231' | grep -v -e "chr2.." -e "chr3.." -e "chr4.." -e "chr5.." -e "chR491" -e "chR492" -e "chR493" -e "chR494" -e "chR495" -e "chR496" -e "chR497" -e "chR498" -e "chR499"  > $OutDir/R4_FoL_genome_edited.txt
```


```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/orthology2circos_ribbons.py --orthology analysis/orthology/orthomcl/R4_Fv_FoC_FoL_Fo/R4_Fv_FoC_FoL_Fo_orthogroups.txt --name1 R4 --gff1 gene_pred/final_genes/F.oxysporum_fsp_narcissi/N139_ncbi/final/final_genes_appended.gff3 --name2 FoL --gff2 assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_parsed.gff \
   | sort -k4,5 -V \
   > $OutDir/R4_FoL_links.txt
  # Links to FoL LS contigs 3, 6, 14 and 15 were coloured black
  # cat $OutDir/R4_FoL_links.txt \
  #   | sed '/4287_CM000591.1/ s/$/\tcolor=black/' \
  #   | sed '/4287_CM000594.1/ s/$/\tcolor=black/' \
  #   | sed '/4287_CM000602.2/ s/$/\tcolor=black/' \
  #   | sed '/4287_CM000603.1/ s/$/\tcolor=black/' \
  #   > $OutDir/R4_FoL_links_edited.txt
  cat $OutDir/R4_FoL_links.txt > $OutDir/R4_FoL_links_edited.txt
```

A file showing contig orientations was made:
```bash
cat $OutDir/R4_FoL_links_edited.txt | cut -f1 | uniq > $OutDir/R4_contig_order.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/find_contig_orientation.py --links_file $OutDir/R4_FoL_links_edited.txt > $OutDir/R4_FoL_contig_orientation.txt
```

The number of bp in syntenous contigs was identified using:
cat $OutDir/R4_FoL_links_edited.txt | cut -f1,4 | sort | uniq -c | sort -nr | grep -v ' 1 ' | sed -r "s/\s+/\t/g" | cut -f3 | sort | uniq > $OutDir/syntenous_conitg_list.txt
```bash
  cat $OutDir/R4_FoL_contig_orientation.txt | tail -n3 | grep -v 'orientation' | sed 's/, /\n/g' > $OutDir/R4_syntenous_contigs.txt
  cat $OutDir/R4_genome.txt | grep -f $OutDir/R4_syntenous_contigs.txt | cut -f6 -d ' ' | awk '{s+=$1} END {print s}'
```

Contig order was selected by taking the first line of that file and then also
taking the reversed order of FoL contigs using the command:

```bash
cat $OutDir/R4_FoL_contig_orientation.txt | grep -A1 'Order of all seen contigs'
cat $OutDir/R4_FoL_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" > tmp.txt
cat $OutDir/R4_FoL_genome_edited.txt | grep 'R4' | grep -v -e "chr2.." -e "chr3.." -e "chr4.." -e "chr5.." -e "chR491" -e "chR492" -e "chR493" -e "chR494" -e "chR495" -e "chR496" -e "chR497" -e "chR498" -e "chR499" | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/R4/, R4/g'
cat $OutDir/R4_FoL_genome_edited.txt | grep 'FoL' | cut -f3 -d ' ' | tr -d '\n' | sed 's/FoL/, FoL/g' >> tmp.txt
```

Contig orientation was used to edit the circos .conf file manually

## Preparing Effector plots

# Plot location of Fus2 mimps and secreted effectorP genes as a scatterplot

```bash
  GffMimp=analysis/mimps/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_mimps.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffMimp --feature MIMP_motif --value 1 | sed -e 's/^/R4_/g' > $OutDir/R4_mimp_plot.txt
  GffCAZY=gene_pred/CAZY/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_CAZY_secreted.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffCAZY --feature gene --value 1 | sed -e 's/^/R4_/g' > $OutDir/R4_CAZY_plot.txt
  GffEfR4=analysis/effectorP/F.oxysporum_fsp_narcissi/N139_ncbi/F.oxysporum_fsp_narcissi_N139_ncbi_EffectorP_secreted.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffEfR4 --feature gene --value 1 | sed -e 's/^/R4_/g' > $OutDir/R4_effectorP_plot.txt
  GffAntiSmash=analysis/antismash/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_secondary_metabolite_regions.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffAntiSmash --feature indole indole-nrps nrps nrps-t1pks other t1pks t1pks-nrps t3pks terpene --value 1 | sed -e 's/^/R4_/g' > $OutDir/R4_antismash_plot.txt

  BlastHits=analysis/blast_homology/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_Fo_path_genes_CRX.fa_homologs.gff
  GffSix=$OutDir/R4_SIX.gff
  cat $BlastHits | grep -v -e 'MIMP' -e 'C5' -e 'CRX' > $GffSix
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffSix --feature SIX_homolog --value 1 | sed -e 's/^/R4_/g' > $OutDir/R4_SIX_plot.txt
```

## Running circos



```bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
circos -conf $ProgDir/R4_R4/R4_FoL/R4_FoL_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/R4_FoL_circos.png
mv $OutDir/circos.svg $OutDir/R4_FoL_circos.svg
```


# Further analysis of non-syntenous regions regions

The number of MIMPs and effectors in LS regions were identified:


```bash
cat $OutDir/R4_mimp_plot.txt | grep -v -f $OutDir/R4_syntenous_contigs.txt | wc -l

```
 -->
