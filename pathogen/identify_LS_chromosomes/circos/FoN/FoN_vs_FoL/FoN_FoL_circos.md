## Preparing genome and .conf files

```bash
OutDir=analysis/circos/F.oxysporum_fsp_narcissi/FoN_FoL_v2
mkdir -p $OutDir

FoN_genome=$(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep 'N139' | grep -v 'old')

ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $FoN_genome --contig_prefix "FoN_" > $OutDir/FoN_genome.txt

FoL_genome=repeat_masked/F.oxysporum_fsp_lycopersici/4287_v2/fungidb_repmask/4287_v2_contigs_unmasked.fa
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $FoL_genome --contig_prefix "FoL_" > $OutDir/FoL_genome.txt

  cat $OutDir/FoN_genome.txt > $OutDir/FoN_FoL_genome.txt
  tac $OutDir/FoL_genome.txt >> $OutDir/FoN_FoL_genome.txt

  # Contigs smaller than 10Kb were removed
  # cat $OutDir/FoN_FoL_genome.txt | grep -v 'DS231' | grep -v -e "0 .... chr" -e "0 ... chr" > $OutDir/FoN_FoL_genome_edited.txt
    cat $OutDir/FoN_FoL_genome.txt > $OutDir/FoN_FoL_genome_edited.txt

  # cat $OutDir/FoN_FoL_genome.txt | grep -v 'DS231' > $OutDir/FoN_FoL_genome_edited.txt
```

```bash
# cat assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_parsed.gff | sed 's/_t26-1//g'

  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/orthology2circos_ribbons.py --orthology analysis/orthology/orthomcl/FoN_vs_FoC_vs_FoL_vs_Fo/FoN_vs_FoC_vs_FoL_vs_Fo_orthogroups.txt --name1 FoN --gff1 gene_pred/final_genes/F.oxysporum_fsp_narcissi/N139_ncbi/final/final_genes_appended.gff3 --name2 FoL --gff2 assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_parsed.gff \
   | sort -k4,5 -V \
   > $OutDir/FoN_FoL_links.txt
  # Links to FoL LS contigs 3, 6, 14 and 15 were coloured black
  # cat $OutDir/FoN_FoL_links.txt \
  #   | sed '/4287_CM000591.1/ s/$/\tcolor=black/' \
  #   | sed '/4287_CM000594.1/ s/$/\tcolor=black/' \
  #   | sed '/4287_CM000602.2/ s/$/\tcolor=black/' \
  #   | sed '/4287_CM000603.1/ s/$/\tcolor=black/' \
  #   > $OutDir/FoN_FoL_links_edited.txt

  # cat $OutDir/FoN_FoL_genome.txt | grep -v 'DS231' | grep -e "0 .... chr" -e "0 ... chr" | cut -f3 -d ' ' > $OutDir/FoN_FoL_excluded_contigs.txt
  echo "" > $OutDir/FoN_FoL_excluded_contigs.txt
  cat $OutDir/FoN_FoL_links.txt > $OutDir/FoN_FoL_links_edited.txt
  # cat $OutDir/FoN_FoL_links.txt | grep -v 'DS231' > $OutDir/FoN_FoL_links_edited.txt
  # cat $OutDir/FoN_FoL_links_edited.txt | wc -l
```

A file showing contig orientations was made:
```bash
  cat $OutDir/FoN_FoL_links_edited.txt | cut -f1 | uniq > $OutDir/FoN_contig_order.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/find_contig_orientation.py --links_file $OutDir/FoN_FoL_links_edited.txt > $OutDir/FoN_FoL_contig_orientation.txt

  echo "contigs with synteny are:"
  cat $OutDir/FoN_FoL_contig_orientation.txt | grep -A1 'seen contigs' | tail -n+2 | sed 's/, /\n/g' > tmp.txt
  # cat $OutDir/FoN_FoL_genome_edited.txt | grep -vw -f tmp.txt | grep 'FoN' | cut -f3 -d ' ' | sed "s/$/, /g" | tr -d '\n'
```

The number of bp in syntenous contigs was identified using:

```bash
  cat $OutDir/FoN_FoL_contig_orientation.txt | tail -n3 | grep -v 'orientation' | sed 's/, /\n/g' > $OutDir/FoN_syntenous_contigs.txt
  echo "Total bp in syntenous contigs:"
  # cat $OutDir/FoN_genome.txt | grep -v -e "0 .... chr" -e "0 ... chr" | grep -wf $OutDir/FoN_syntenous_contigs.txt | cut -f6 -d ' ' | awk '{s+=$1} END {print s}'
  cat $OutDir/FoN_genome.txt | grep -wf $OutDir/FoN_syntenous_contigs.txt | cut -f6 -d ' ' | awk '{s+=$1} END {print s}'
  echo "This is shown in the following number of contigs:"
  # cat $OutDir/FoN_genome.txt | grep -v -e "0 .... chr" -e "0 ... chr" | grep -wf $OutDir/FoN_syntenous_contigs.txt | cut -f6 -d ' ' | wc -l
  cat $OutDir/FoN_genome.txt | grep -wf $OutDir/FoN_syntenous_contigs.txt | cut -f6 -d ' ' | wc -l
  echo "The ramining contigs represent the following number of bp:"
  cat $OutDir/FoN_genome.txt | grep -wvf $OutDir/FoN_syntenous_contigs.txt | cut -f6 -d ' ' | awk '{s+=$1} END {print s}'
```

<!-- ```
Total bp in syntenous contigs:
42618748
This is shown in the following number of contigs:
267
The ramining contigs represent the following number of bp:
14898814
``` -->
```
Total bp in syntenous contigs:
42662134
This is shown in the following number of contigs:
273
The remining contigs represent the following number of bp:
14855428
```

Contig order was selected by taking the first line of that file and then also
taking the reversed order of FoL contigs using the command:

```bash
# all contigs
# cat $OutDir/FoN_FoL_genome_edited.txt | grep 'FoN' | grep -v -e "0 .... chr" -e "0 ... chr" | cut -f3 -d ' '| tr -d '\n' | sed 's/FoN/, FoN/g' | sed 's/FoL/, FoL/g' > tmp3.txt
cat $OutDir/FoN_FoL_genome_edited.txt | grep 'FoN' | cut -f3 -d ' '| tr -d '\n' | sed 's/FoN/, FoN/g' | sed 's/FoL/, FoL/g' > tmp3.txt
cat $OutDir/FoN_FoL_genome_edited.txt | grep 'FoL' | cut -f3 -d ' ' | tr -d '\n' | sed 's/FoL/, FoL/g'  >> tmp3.txt
# Forward orientation
cat $OutDir/FoN_FoL_contig_orientation.txt | grep -A1 'Order of all seen contigs'
cat $OutDir/FoN_FoL_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" > tmp.txt
# reverse orientation
# cat $OutDir/FoN_FoL_genome_edited.txt | grep 'FoN' | grep -v -e "0 .... chr" -e "0 ... chr" | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/FoN/, FoN/g'
cat $OutDir/FoN_FoL_genome_edited.txt | grep 'FoN' | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/FoN/, FoN/g'
# reference sequences (put in reverse)
# contig order
cat $OutDir/FoN_FoL_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1
cat $OutDir/FoN_FoL_genome_edited.txt | grep 'FoL' | cut -f3 -d ' ' | tr -d '\n' | sed 's/FoL/, FoL/g'
# unplaced contigs
echo "The following contigs are unplaced:"
cat $OutDir/FoN_FoL_genome_edited.txt | grep -vw -f tmp.txt | grep 'FoN' | cut -f3 -d ' ' | tr -d '\n' | sed 's/FoN/, FoN/g' > $OutDir/FoN_FoL_contigs_unplaced.txt

# cat $OutDir/FoN_genome.txt | grep -wvf $OutDir/FoN_syntenous_contigs.txt | cut -f3 -d ' ' | tr -d '\n' | sed "s/FoN/, FoN/g" > $OutDir/FoN_FoL_contigs_unplaced.txt
```

Contig orientation was used to edit the circos .conf file manually


## Preparing Effector plots

# Plot location of Fus2 mimps and secreted effectorP genes as a scatterplot

```bash
  GffMimp=$(ls analysis/mimps/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_mimps.gff)
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffMimp --feature MIMP_motif --value 1 | sed -e 's/^/FoN_/g' > $OutDir/FoN_mimp_plot.txt
  GffCAZY=gene_pred/CAZY/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_CAZY_secreted.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffCAZY --feature gene --value 1 | sed -e 's/^/FoN_/g' > $OutDir/FoN_CAZY_plot.txt
  GffEfFoN=analysis/effectorP/F.oxysporum_fsp_narcissi/N139_ncbi/F.oxysporum_fsp_narcissi_N139_ncbi_EffectorP_secreted.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffEfFoN --feature gene --value 1 | sed -e 's/^/FoN_/g' > $OutDir/FoN_effectorP_plot.txt
  # GffAntiSmash=analysis/antismash/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_secondary_metabolite_regions.gff
  # ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  # $ProgDir/gff2circos_scatterplot.py --gff $GffAntiSmash --feature indole indole-nrps nrps nrps-t1pks other t1pks t1pks-nrps t3pks terpene --value 1 | sed -e 's/^/FoN_/g' > $OutDir/FoN_antismash_plot.txt

  BlastHits=analysis/blast_homology/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_Fo_path_genes_CRX.fa_homologs.gff
  GffSix=$OutDir/FoN_SIX.gff
  cat $BlastHits | grep -v -e 'MIMP' -e 'C5' -e 'CRX' | grep 'SIX' > $GffSix
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffSix --feature SIX_homolog --value 1 | sed -e 's/^/FoN_/g' > $OutDir/FoN_SIX_plot.txt
```



## Running circos



```bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
circos -conf $ProgDir/FoN/FoN_vs_FoL/FoN_FoL_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/FoN_FoL_circos_full.png
mv $OutDir/circos.svg $OutDir/FoN_FoL_circos_full.svg
```

<!--
# Further analysis of non-syntenous regions regions

The number of MIMPs and effectors in LS regions were identified:


```bash
cat $OutDir/FoN_mimp_plot.txt | grep -v -f $OutDir/FoN_syntenous_contigs.txt | wc -l

```
## Preparing genome and .conf files

```bash
OutDir=analysis/circos/F.oxysporum_fsp_narcissi/FoN_FoL
mkdir -p $OutDir

FoN_genome=repeat_masked/F.oxysporum_fsp_narcissi/N139_ncbi/ncbi_submission/FoN_contigs_unmasked.fa
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $FoN_genome --contig_prefix "FoN_" > $OutDir/FoN_genome.txt

FoL_genome=repeat_masked/F.oxysporum_fsp_lycopersici/4287_v2/fungidb_repmask/4287_v2_contigs_unmasked.fa
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --genome $FoL_genome --contig_prefix "FoL_" > $OutDir/FoL_genome.txt

  cat $OutDir/FoN_genome.txt > $OutDir/FoN_FoL_genome.txt
  tac $OutDir/FoL_genome.txt >> $OutDir/FoN_FoL_genome.txt

  # COntigs smaller than 10Kb were removed
  cat $OutDir/FoN_FoL_genome.txt | grep -v 'DS231' | grep -v -e "chr2.." -e "chr3.." -e "chr4.." -e "chr5.." -e "chr191" -e "chr192" -e "chr193" -e "chr194" -e "chr195" -e "chr196" -e "chr197" -e "chr198" -e "chr199"  > $OutDir/FoN_FoL_genome_edited.txt
```


```bash
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/orthology2circos_ribbons.py --orthology analysis/orthology/orthomcl/FoN_Fv_FoC_FoL_Fo/FoN_Fv_FoC_FoL_Fo_orthogroups.txt --name1 FoN --gff1 gene_pred/final_genes/F.oxysporum_fsp_narcissi/N139_ncbi/final/final_genes_appended.gff3 --name2 FoL --gff2 assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_parsed.gff \
   | sort -k4,5 -V \
   > $OutDir/FoN_FoL_links.txt
  # Links to FoL LS contigs 3, 6, 14 and 15 were coloured black
  # cat $OutDir/FoN_FoL_links.txt \
  #   | sed '/4287_CM000591.1/ s/$/\tcolor=black/' \
  #   | sed '/4287_CM000594.1/ s/$/\tcolor=black/' \
  #   | sed '/4287_CM000602.2/ s/$/\tcolor=black/' \
  #   | sed '/4287_CM000603.1/ s/$/\tcolor=black/' \
  #   > $OutDir/FoN_FoL_links_edited.txt
  cat $OutDir/FoN_FoL_links.txt > $OutDir/FoN_FoL_links_edited.txt
```

A file showing contig orientations was made:
```bash
cat $OutDir/FoN_FoL_links_edited.txt | cut -f1 | uniq > $OutDir/FoN_contig_order.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/circos
  $ProgDir/find_contig_orientation.py --links_file $OutDir/FoN_FoL_links_edited.txt > $OutDir/FoN_FoL_contig_orientation.txt
```

The number of bp in syntenous contigs was identified using:
cat $OutDir/FoN_FoL_links_edited.txt | cut -f1,4 | sort | uniq -c | sort -nr | grep -v ' 1 ' | sed -r "s/\s+/\t/g" | cut -f3 | sort | uniq > $OutDir/syntenous_conitg_list.txt
```bash
  cat $OutDir/FoN_FoL_contig_orientation.txt | tail -n3 | grep -v 'orientation' | sed 's/, /\n/g' > $OutDir/FoN_syntenous_contigs.txt
  cat $OutDir/FoN_genome.txt | grep -f $OutDir/FoN_syntenous_contigs.txt | cut -f6 -d ' ' | awk '{s+=$1} END {print s}'
```

Contig order was selected by taking the first line of that file and then also
taking the reversed order of FoL contigs using the command:

```bash
cat $OutDir/FoN_FoL_contig_orientation.txt | grep -A1 'Order of all seen contigs'
cat $OutDir/FoN_FoL_contig_orientation.txt | grep -A1 'Order of all seen contigs' | tail -n1 | sed "s/, /\n/g" > tmp.txt
cat $OutDir/FoN_FoL_genome_edited.txt | grep 'FoN' | grep -v -e "chr2.." -e "chr3.." -e "chr4.." -e "chr5.." -e "chr191" -e "chr192" -e "chr193" -e "chr194" -e "chr195" -e "chr196" -e "chr197" -e "chr198" -e "chr199" | grep -w -v -f tmp.txt | cut -f3 -d ' '| tr -d '\n' | sed 's/FoN/, FoN/g'
cat $OutDir/FoN_FoL_genome_edited.txt | grep 'FoL' | cut -f3 -d ' ' | tr -d '\n' | sed 's/FoL/, FoL/g' >> tmp.txt
```

Contig orientation was used to edit the circos .conf file manually

## Preparing Effector plots

# Plot location of Fus2 mimps and secreted effectorP genes as a scatterplot

```bash
  GffMimp=analysis/mimps/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_mimps.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffMimp --feature MIMP_motif --value 1 | sed -e 's/^/FoN_/g' > $OutDir/FoN_mimp_plot.txt
  GffCAZY=gene_pred/CAZY/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_CAZY_secreted.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffCAZY --feature gene --value 1 | sed -e 's/^/FoN_/g' > $OutDir/FoN_CAZY_plot.txt
  GffEfFoN=analysis/effectorP/F.oxysporum_fsp_narcissi/N139_ncbi/F.oxysporum_fsp_narcissi_N139_ncbi_EffectorP_secreted.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffEfFoN --feature gene --value 1 | sed -e 's/^/FoN_/g' > $OutDir/FoN_effectorP_plot.txt
  GffAntiSmash=analysis/antismash/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_secondary_metabolite_regions.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffAntiSmash --feature indole indole-nrps nrps nrps-t1pks other t1pks t1pks-nrps t3pks terpene --value 1 | sed -e 's/^/FoN_/g' > $OutDir/FoN_antismash_plot.txt

  BlastHits=analysis/blast_homology/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_Fo_path_genes_CRX.fa_homologs.gff
  GffSix=$OutDir/FoN_SIX.gff
  cat $BlastHits | grep -v -e 'MIMP' -e 'C5' -e 'CRX' > $GffSix
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffSix --feature SIX_homolog --value 1 | sed -e 's/^/FoN_/g' > $OutDir/FoN_SIX_plot.txt
```

## Running circos



```bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
circos -conf $ProgDir/FoN_FoN/FoN_FoL/FoN_FoL_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/FoN_FoL_circos.png
mv $OutDir/circos.svg $OutDir/FoN_FoL_circos.svg
```


# Further analysis of non-syntenous regions regions

The number of MIMPs and effectors in LS regions were identified:


```bash
cat $OutDir/FoN_mimp_plot.txt | grep -v -f $OutDir/FoN_syntenous_contigs.txt | wc -l

```
 -->
