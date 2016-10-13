Circos plots were generated for reassembly of the Fus2 genome

```bash

  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  Fus2_genome=assembly/canu-1.3/F.oxysporum_fsp_cepae/Fus2_canu/pilon/pilon.fasta
  OutDir=analysis/circos/F.oxysporum_fsp_cepae/Fus2_hortres
  mkdir -p $OutDir

  # Convert the Fus2 genome into circos format
  # $ProgDir/fasta2circos.py --genome $Fus2_genome --contig_prefix "" > $OutDir/Fus2_genome.txt
  # The Fus2 contigs were re-ordered in accordance with FoL chromsome names:
  cat analysis/circos/F.oxysporum_fsp_cepae/Fus2_FoL/Fus2_FoL_genome_edited.txt  | grep 'pilon' | grep -v -e 'chr23' -e 'chr24' -e 'chr25' -e 'chr26' -e 'chr27' -e 'chr28' -e 'chr29' -e 'chr30' -e 'chr31' -e 'chr32' -e 'chr33' -e 'chr34'| sed 's/Fus2_//g' > $OutDir/Fus2_genome.txt


  # Plot location of FoL gene Blast hits as a scatterplot
  for GffFile in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_canu/*_chr_*_gene_single_copy.aa_hits.gff); do
  # for GffFile in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_edited2/Fus2_edited2_chr_*_gene_single_copy.aa_hits.gff); do
    echo $GffFile
    Chr=$(echo $GffFile | rev |cut -f1 -d'/' | rev | cut -f4 -d '_')
    $ProgDir/gff2circos_scatterplot.py --gff $GffFile --feature Chr"$Chr"_gene_homolog > $OutDir/FoL_chr"$Chr"_genes.txt
  done

  # Make 100kb windows for plots
  $ProgDir/fasta2gff_windows.py --genome $Fus2_genome --size 20 > $OutDir/Fus2_20kb_windows.gff

  # Convert FoC MiSeq reads aligning in 100kb windows into coverage stats
# for ReadsBam in $(ls analysis/genome_alignment/bowtie/F.*/*/vs_Fus2_hardmasked_max1200/Fus2_canu_contigs_hardmasked_repeatmasker_TPSI_appended.fa_aligned_sorted.bam | grep -w -e 'Fus2' -e '125' -e 'A28' -e 'A13' -e 'CB3' -e 'PG' -e 'A23' -e 'A28' -e 'A8' ); do
for ReadsBam in $(ls analysis/genome_alignment/bowtie/F.*/*/vs_Fus2_unmasked_max1200/Fus2_canu_contigs_unmasked.fa_aligned_sorted.bam | grep -w -e 'Fus2' -e '125' -e 'A28' -e 'A13' -e 'CB3' -e 'PG' -e 'A23' -e 'A28' -e 'A8' ); do
Organism=$(echo $ReadsBam | rev | cut -f4 -d '/' | rev)
Strain=$(echo $ReadsBam | rev | cut -f3 -d '/' | rev)
AlignDir=$(dirname $ReadsBam)
echo "$Organism - $Strain"
bedtools coverage -abam $ReadsBam -b $OutDir/Fus2_20kb_windows.gff > $AlignDir/"$Strain"_coverage_vs_Fus2.bed
# Convert coverage bed files into circos format
$ProgDir/coverage_bed2circos.py --bed $AlignDir/"$Strain"_coverage_vs_Fus2.bed > $OutDir/"$Strain"_coverage_vs_Fus2_unmasked_scatterplot.txt
ls $OutDir/"$Strain"_coverage_vs_Fus2_unmasked_scatterplot.txt
done

  # Plot location of Fus2 mimps and secreted effectorP genes as a scatterplot
  GffMimp=analysis/mimps/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_mimps.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffMimp --feature MIMP_motif --value 1 > $OutDir/Fus2_mimp_plot.txt
  GffCAZY=gene_pred/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY_secreted.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffCAZY --feature gene --value 1 > $OutDir/Fus2_CAZY_plot.txt
  GffEffP=analysis/effectorP/F.oxysporum_fsp_cepae/Fus2_canu_new/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_secreted.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffEffP --feature gene --value 1 > $OutDir/Fus2_effectorP_plot.txt
  GffAntiSmash=analysis/antismash/Fus2_clusters.gff
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffAntiSmash --feature indole indole-nrps nrps nrps-t1pks other t1pks t1pks-nrps t3pks terpene --value 1 > $OutDir/Fus2_antismash_plot.txt

  BlastHits=analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_Fo_path_genes_CRX.fa_homologs.gff
  GffSIX=$OutDir/Fus2_SIX.gff
  cat $BlastHits | grep -v -e 'MIMP' -e 'C5' -e 'CRX' > $GffSix
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/gff2circos_scatterplot.py --gff $GffSIX --feature SIX_homolog --value 1 > $OutDir/Fus2_SIX_plot.txt


  circos -conf /home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos/Fus2/hortres_plot/Fus2_hortres_circos.conf -outputdir ./$OutDir
```
