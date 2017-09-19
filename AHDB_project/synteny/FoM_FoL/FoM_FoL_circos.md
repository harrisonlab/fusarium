## Preparing genome and .conf files

```bash
OutDir=analysis/circos/F.oxysporum_fsp_mathioli/FoM_FoL_minion
mkdir -p $OutDir

FoM_genome=$(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep 'Stocks4')

ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $FoM_genome --contig_prefix "FoM_" > $OutDir/FoM_genome.txt

FoL_genome=repeat_masked/F.oxysporum_fsp_lycopersici/4287_v2/fungidb_repmask/4287_v2_contigs_unmasked.fa
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $FoL_genome --contig_prefix "FoL_" > $OutDir/FoL_genome.txt

cat $OutDir/FoM_genome.txt > $OutDir/FoM_FoL_genome.txt
tac $OutDir/FoL_genome.txt >> $OutDir/FoM_FoL_genome.txt

# Contigs smaller than 10Kb were removed
cat $OutDir/FoM_FoL_genome.txt | grep -v 'DS231' | grep -v -e "0 .... chr" -e "0 ... chr" > $OutDir/FoM_FoL_genome_edited.txt
```

```bash

  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/orthology2circos_ribbons.py --orthology analysis/orthology/orthomcl/Fo_FoC_FoL_FoN_FoM/Fo_FoC_FoL_FoN_FoM_orthogroups.txt --name1 FoM --gff1 gene_pred/final_genes/F.oxysporum_fsp_mathioli/Stocks4/final/final_genes_appended_renamed.gff3 --name2 FoL --gff2 assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_parsed.gff \
   | sort -k4,5 -V \
   > $OutDir/FoM_FoL_links.txt
  # Links to FoL LS contigs 3, 6, 14 and 15 were coloured black
  # cat $OutDir/FoM_FoL_links.txt \
  #   | sed '/4287_CM000591.1/ s/$/\tcolor=black/' \
  #   | sed '/4287_CM000594.1/ s/$/\tcolor=black/' \
  #   | sed '/4287_CM000602.2/ s/$/\tcolor=black/' \
  #   | sed '/4287_CM000603.1/ s/$/\tcolor=black/' \
  #   > $OutDir/FoM_FoL_links_edited.txt

  cat $OutDir/FoM_FoL_genome.txt | grep -v 'DS231' | grep -e "0 .... chr" -e "0 ... chr" | cut -f3 -d ' ' > $OutDir/FoM_FoL_excluded_contigs.txt
  cat $OutDir/FoM_FoL_links.txt | grep -v -f $OutDir/FoM_FoL_excluded_contigs.txt > $OutDir/FoM_FoL_links_edited.txt
```

A file showing contig orientations was made:
```bash
cat $OutDir/FoM_FoL_links_edited.txt | cut -f1 | uniq > $OutDir/FoM_contig_order.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/find_contig_orientation.py --links_file $OutDir/FoM_FoL_links_edited.txt > $OutDir/FoM_FoL_contig_orientation.txt
```

The number of bp in syntenous contigs was identified using:

```bash
  cat $OutDir/FoM_FoL_contig_orientation.txt | tail -n3 | grep -v 'orientation' | sed 's/, /\n/g' > $OutDir/FoM_syntenous_contigs.txt
  cat $OutDir/FoM_genome.txt | grep -v -e "0 .... chr" -e "0 ... chr" | grep -f $OutDir/FoM_syntenous_contigs.txt | cut -f6 -d ' ' | awk '{s+=$1} END {print s}'
```

```
  60310576
```

Contig order was selected by taking the first line of that file and then also
taking the reversed order of FoL contigs using the command:

```bash
# all contigs
cat $OutDir/FoM_FoL_genome_edited.txt | grep 'FoM' | grep -v -e "0 .... chr" -e "0 ... chr" | cut -f3 -d ' '| tr -d '\n' | sed 's/FoM/, FoM/g' | sed 's/FoL/, FoL/g' > tmp3.txt
# Forward orientation
cat $OutDir/FoM_FoL_contig_orientation.txt | grep -A1 'Order of all seen contigs'
cat $OutDir/FoM_FoL_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" > tmp.txt
# reverse orientation
cat $OutDir/FoM_FoL_genome_edited.txt | grep 'FoM' | grep -v -e "0 .... chr" -e "0 ... chr" | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/FoM/, FoM/g'
# reference sequences (put in reverse)
# contig order
cat $OutDir/FoM_FoL_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1
cat $OutDir/FoM_FoL_genome_edited.txt | grep 'FoL' | cut -f3 -d ' ' | tr -d '\n' | sed 's/FoL/, FoL/g'
```

Contig orientation was used to edit the circos .conf file manually
<!--
## Preparing Effector plots

# Plot location of Fus2 mimps and secreted effectorP genes as a scatterplot

```bash
  GffMimp=analysis/mimps/F.proliferatum/A8_ncbi/A8_ncbi_mimps.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffMimp --feature MIMP_motif --value 1 | sed -e 's/^/FoM_/g' > $OutDir/A8_mimp_plot.txt
  GffCAZY=gene_pred/CAZY/F.proliferatum/A8_ncbi/A8_ncbi_CAZY_secreted.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffCAZY --feature gene --value 1 | sed -e 's/^/FoM_/g' > $OutDir/A8_CAZY_plot.txt
  GffEfFoM=analysis/effectorP/F.proliferatum/A8_ncbi/F.proliferatum_A8_ncbi_EffectorP_secreted.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffEfFoM --feature gene --value 1 | sed -e 's/^/FoM_/g' > $OutDir/A8_effectorP_plot.txt
  GffAntiSmash=analysis/antismash/F.proliferatum/A8_ncbi/A8_ncbi_secondary_metabolite_regions.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffAntiSmash --feature indole indole-nrps nrps nrps-t1pks other t1pks t1pks-nrps t3pks terpene --value 1 | sed -e 's/^/FoM_/g' > $OutDir/A8_antismash_plot.txt

  BlastHits=analysis/blast_homology/F.proliferatum/A8_ncbi/A8_ncbi_Fo_path_genes_CRX.fa_homologs.gff
  GffSix=$OutDir/A8_SIX.gff
  cat $BlastHits | grep -v -e 'MIMP' -e 'C5' -e 'CRX' > $GffSix
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffSix --feature SIX_homolog --value 1 | sed -e 's/^/FoM_/g' > $OutDir/A8_SIX_plot.txt
``` -->

## Running circos


```bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/AHDB_project/synteny/FoM_FoL
circos -conf $ProgDir/FoM_FoL_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/FoM_FoL_circos.png
mv $OutDir/circos.svg $OutDir/FoM_FoL_circos.svg
```


# Further analysis of non-syntenous regions regions

The number of MIMPs and effectors in LS regions were identified:


```bash
cat $OutDir/A8_mimp_plot.txt | grep -v -f $OutDir/FoM_syntenous_contigs.txt | wc -l

```
## Preparing genome and .conf files

```bash
OutDir=analysis/circos/F.proliferatum/A8_FoL
mkdir -p $OutDir

FoM_genome=repeat_masked/F.proliferatum/A8_ncbi/ncbi_submission/A8_contigs_unmasked.fa
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $FoM_genome --contig_prefix "FoM_" > $OutDir/FoM_genome.txt

FoL_genome=repeat_masked/F.oxysporum_fsp_lycopersici/4287_v2/fungidb_repmask/4287_v2_contigs_unmasked.fa
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $FoL_genome --contig_prefix "FoL_" > $OutDir/FoL_genome.txt

  cat $OutDir/FoM_genome.txt > $OutDir/FoM_FoL_genome.txt
  tac $OutDir/FoL_genome.txt >> $OutDir/FoM_FoL_genome.txt

  # COntigs smaller than 10Kb were removed
  cat $OutDir/FoM_FoL_genome.txt | grep -v 'DS231' | grep -v -e "chr2.." -e "chr3.." -e "chr4.." -e "chr5.." -e "chr191" -e "chr192" -e "chr193" -e "chr194" -e "chr195" -e "chr196" -e "chr197" -e "chr198" -e "chr199"  > $OutDir/FoM_FoL_genome_edited.txt
```
<!--
The order of contigs was changed manually using nano
```bash
cp $OutDir/FoM_FoL_genome_edited.txt $OutDir/FoM_FoL_genome_edited2.txt
nano $OutDir/FoM_FoL_genome_edited2.txt
cp $OutDir/FoM_FoL_genome_edited2.txt $OutDir/FoM_FoL_genome_final.txt
``` -->


```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/orthology2circos_ribbons.py --orthology analysis/orthology/orthomcl/FoM_Fv_FoC_FoL_Fo/FoM_Fv_FoC_FoL_Fo_orthogroups.txt --name1 FoM --gff1 gene_pred/final_genes/F.proliferatum/A8_ncbi/final/final_genes_appended.gff3 --name2 FoL --gff2 assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_parsed.gff \
   | sort -k4,5 -V \
   > $OutDir/FoM_FoL_links.txt
  # Links to FoL LS contigs 3, 6, 14 and 15 were coloured black
  # cat $OutDir/FoM_FoL_links.txt \
  #   | sed '/4287_CM000591.1/ s/$/\tcolor=black/' \
  #   | sed '/4287_CM000594.1/ s/$/\tcolor=black/' \
  #   | sed '/4287_CM000602.2/ s/$/\tcolor=black/' \
  #   | sed '/4287_CM000603.1/ s/$/\tcolor=black/' \
  #   > $OutDir/FoM_FoL_links_edited.txt
  cat $OutDir/FoM_FoL_links.txt > $OutDir/FoM_FoL_links_edited.txt
```

A file showing contig orientations was made:
```bash
cat $OutDir/FoM_FoL_links_edited.txt | cut -f1 | uniq > $OutDir/FoM_contig_order.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/find_contig_orientation.py --links_file $OutDir/FoM_FoL_links_edited.txt > $OutDir/FoM_FoL_contig_orientation.txt
```

The number of bp in syntenous contigs was identified using:

```bash
  cat $OutDir/FoM_FoL_contig_orientation.txt | tail -n3 | grep -v 'orientation' | sed 's/, /\n/g' > $OutDir/FoM_syntenous_contigs.txt
  cat $OutDir/FoM_genome.txt | grep -v -e "chr2.." -e "chr3.." -e "chr4.." -e "chr5.." -e "chr191" -e "chr192" -e "chr193" -e "chr194" -e "chr195" -e "chr196" -e "chr197" -e "chr198" -e "chr199" | grep -f $OutDir/FoM_syntenous_contigs.txt | cut -f6 -d ' ' | awk '{s+=$1} END {print s}'
```

Contig order was selected by taking the first line of that file and then also
taking the reversed order of FoL contigs using the command:

```bash
cat $OutDir/FoM_FoL_contig_orientation.txt | grep -A1 'Order of all seen contigs'
cat $OutDir/FoM_FoL_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" > tmp.txt
cat $OutDir/FoM_FoL_genome_edited.txt | grep 'FoM' | grep -v -e "chr2.." -e "chr3.." -e "chr4.." -e "chr5.." -e "chr191" -e "chr192" -e "chr193" -e "chr194" -e "chr195" -e "chr196" -e "chr197" -e "chr198" -e "chr199" | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/FoM/, FoM/g'
cat $OutDir/FoM_FoL_genome_edited.txt | grep 'FoL' | cut -f3 -d ' ' | tr -d '\n' | sed 's/FoL/, FoL/g' >> tmp.txt
```

Contig orientation was used to edit the circos .conf file manually

## Preparing Effector plots

# Plot location of Fus2 mimps and secreted effectorP genes as a scatterplot

```bash
  GffMimp=analysis/mimps/F.proliferatum/A8_ncbi/A8_ncbi_mimps.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffMimp --feature MIMP_motif --value 1 | sed -e 's/^/FoM_/g' > $OutDir/A8_mimp_plot.txt
  GffCAZY=gene_pred/CAZY/F.proliferatum/A8_ncbi/A8_ncbi_CAZY_secreted.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffCAZY --feature gene --value 1 | sed -e 's/^/FoM_/g' > $OutDir/A8_CAZY_plot.txt
  GffEfFoM=analysis/effectorP/F.proliferatum/A8_ncbi/F.proliferatum_A8_ncbi_EffectorP_secreted.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffEfFoM --feature gene --value 1 | sed -e 's/^/FoM_/g' > $OutDir/A8_effectorP_plot.txt
  GffAntiSmash=analysis/antismash/F.proliferatum/A8_ncbi/A8_ncbi_secondary_metabolite_regions.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffAntiSmash --feature indole indole-nrps nrps nrps-t1pks other t1pks t1pks-nrps t3pks terpene --value 1 | sed -e 's/^/FoM_/g' > $OutDir/A8_antismash_plot.txt

  BlastHits=analysis/blast_homology/F.proliferatum/A8_ncbi/A8_ncbi_Fo_path_genes_CRX.fa_homologs.gff
  GffSix=$OutDir/A8_SIX.gff
  cat $BlastHits | grep -v -e 'MIMP' -e 'C5' -e 'CRX' > $GffSix
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffSix --feature SIX_homolog --value 1 | sed -e 's/^/FoM_/g' > $OutDir/A8_SIX_plot.txt
```

## Running circos



```bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/AHDB_project/synteny/FoM_FoL
circos -conf $ProgDir/FoM_FoL_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/FoM_FoL_circos.png
mv $OutDir/circos.svg $OutDir/FoM_FoL_circos.svg
```


# Further analysis of non-syntenous regions regions

The number of MIMPs and effectors in LS regions were identified:


```bash
cat $OutDir/A8_mimp_plot.txt | grep -v -f $OutDir/FoM_syntenous_contigs.txt | wc -l

```
