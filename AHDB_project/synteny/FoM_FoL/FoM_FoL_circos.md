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

cat $OutDir/FoM_FoL_genome.txt > $OutDir/FoM_FoL_genome_edited.txt

# Contigs smaller than 10Kb were removed
# cat $OutDir/FoM_FoL_genome.txt | grep -v 'DS231' | grep -v -e "0 .... chr" -e "0 ... chr" > $OutDir/FoM_FoL_genome_edited.txt
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

## Preparing Effector plots

# Plot location of mimps, secreted genes within 2Kb of a mimp a scatterplot

```bash
  GffMimpSecreted=analysis/mimps/F.oxysporum_fsp_mathioli/Stocks4/Stocks4_genes_in_2kb_mimp_secreted.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffMimpSecreted --feature gene --value 1 | sed -e 's/^/FoM_/g' > $OutDir/FoM_mimp_plot.txt
```

## Preparing Effector plots

# Plot location of FoM mimps and secreted effectorP genes as a scatterplot

```bash
Organism="F.oxysporum_fsp_mathioli"
Strain="Stocks4"
  GffMimp=$(ls analysis/mimps/${Organism}/${Strain}/${Strain}_mimps.gff)
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffMimp --feature MIMP_motif --value 1 | sed -e 's/^/FoM_/g' > $OutDir/FoM_mimp_plot.txt
  GffCAZY=$(ls gene_pred/CAZY/${Organism}/${Strain}/*_CAZY_secreted.gff)
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffCAZY --feature gene --value 1 | sed -e 's/^/FoM_/g' > $OutDir/FoM_CAZY_plot.txt
  GffEfFoM=$(ls analysis/effectorP/${Organism}/${Strain}/*_EffectorP_secreted.gff)
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffEfFoM --feature gene --value 1 | sed -e 's/^/FoM_/g' > $OutDir/FoM_effectorP_plot.txt
  # GffAntiSmash=analysis/antismash/F.oxysporum_fsp_narcissi/N139_ncbi/N139_ncbi_secondary_metabolite_regions.gff
  # ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  # $ProgDir/gff2circos_scatterplot.py --gff $GffAntiSmash --feature indole indole-nrps nrps nrps-t1pks other t1pks t1pks-nrps t3pks terpene --value 1 | sed -e 's/^/FoM_/g' > $OutDir/FoM_antismash_plot.txt

  BlastHits=$(ls analysis/blast_homology/${Organism}/${Strain}/*_six-appended_parsed.fa_homologs.gff)
  GffSix=$OutDir/FoM_SIX.gff
  cat $BlastHits | grep -v -e 'MIMP' -e 'C5' -e 'CRX' | grep 'SIX' > $GffSix
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffSix --feature BLAST_hit --value 1 | sed -e 's/^/FoM_/g' > $OutDir/FoM_SIX_plot.txt
```


Links were also drawn from whole genome alignments:

```bash
Organism="F.oxysporum_fsp_mathioli"
Strain="Stocks4"
Prefix="FoM"
OutDir=analysis/circos/F.oxysporum_fsp_mathioli/FoM_FoL_minion

Coords=$(ls analysis/genome_alignment/mummer/${Organism}/${Strain}/${Strain}_vs_4287/${Strain}_vs_4287_coords.tsv)
ProgDir=/home/armita/git_repos/emr_repos/scripts/alternaria/pathogen/genome_alignment
$ProgDir/nucmer_coords2circos.py --inp_coords $Coords --queery_id $Prefix --ref_id FoL > $OutDir/${Prefix}_vs_FoL_links_nucmer.txt
```

```bash
Organism="F.oxysporum_fsp_mathioli"
Strain="Stocks4"
Prefix="FoM"
OutDir=analysis/circos/F.oxysporum_fsp_mathioli/FoM_FoL_minion

ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos/FoLactucae/R1_vs_FoL
$ProgDir/order_Fo_contigs.py --karotype $OutDir/${Prefix}_FoL_genome_edited.txt --links $OutDir/${Prefix}_vs_FoL_links_nucmer.txt
```


## Running circos


```bash
OutDir=analysis/circos/F.oxysporum_fsp_mathioli/FoM_FoL_minion
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/AHDB_project/synteny/FoM_FoL
circos -conf $ProgDir/FoM_FoL_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/FoM_FoL_circos.png
mv $OutDir/circos.svg $OutDir/FoM_FoL_circos.svg
ls $PWD/$OutDir/FoM_FoL_circos.png
```
