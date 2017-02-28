```bash
  OutDir=analysis/circos/F.proliferatum/A8_FoL
  mkdir -p $OutDir

  Fp_genome=repeat_masked/F.proliferatum/A8_ncbi/ncbi_submission/A8_contigs_unmasked.fa
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/fasta2circos.py --genome $Fp_genome --contig_prefix "Fp_" > $OutDir/Fp_genome.txt

  FoL_genome=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.dna.chromosome.fa
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/fasta2circos.py --genome $FoL_genome --contig_prefix "FoL_" > $OutDir/FoL_genome.txt

```

<!--
```bash
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  Fp_genome=repeat_masked/F.proliferatum/A8_ncbi/ncbi_submission/A8_contigs_unmasked.fa
  circos -conf $ProgDir/Fus2/Fus2_FoL/Fus2_FoL_circos.conf -outputdir $OutDir
``` -->


```bash
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  Fp_genome=repeat_masked/F.proliferatum/A8_ncbi/ncbi_submission/A8_contigs_unmasked.fa

  # Convert the Fus2 genome into circos format
  $ProgDir/fasta2circos.py --genome $Fp_genome --contig_prefix "" > tmp5/Fp_genome.txt

  # Make 100kb windows for plots
  $ProgDir/fasta2gff_windows.py --genome $Fp_genome > tmp5/Fp_100kb_windows.gff

  # Extract RNAseq reads aligning in 100kb windows
  # ReadsBam=alignment/F.oxysporum_fsp_cepae/Fus2/Fus2_72hrs_rep3/accepted_hits.bam
  # bedtools coverage -abam $ReadsBam -b tmp5/Fp_100kb_windows.gff > tmp5/tmp.bed
  # Convert coverage bed files into circos format
  # $ProgDir/coverage_bed2circos.py --bed tmp5/tmp.bed > tmp5/Fus2_coverage_scatterplot.txt

  # Identify GC content in 100kb windows
  # $ProgDir/gc_content2circos.py --genome $Fp_genome --gff tmp5/Fp_100kb_windows.gff > tmp5/Fus2_GC_scatterplot.txt

  # Convert FoC MiSeq reads aligning in 100kb windows into coverage stats
  # for ReadsBam in $(ls analysis/genome_alignment/bowtie/F.*/*/vs_Fus2/Fus2_contigs_softmasked.fa_aligned_sorted.bam | grep -v 'PG' | grep 'HB6'); do
  #   Organism=$(echo $ReadsBam | rev | cut -f4 -d '/' | rev)
  #   Strain=$(echo $ReadsBam | rev | cut -f3 -d '/' | rev)
  #   AlignDir=$(dirname $ReadsBam)
  #   echo "$Organism - $Strain"
  #   bedtools coverage -abam $ReadsBam -b tmp5/Fp_100kb_windows.gff > $AlignDir/"$Strain"_coverage_vs_Fus2.bed
  #   # Convert coverage bed files into circos format
  #   $ProgDir/coverage_bed2circos.py --bed $AlignDir/"$Strain"_coverage_vs_Fus2.bed > tmp5/"$Strain"_coverage_vs_Fus2_scatterplot.txt
  # done

  # Plot location of FoL gene Blast hits as a scatterplot
  for GffFile in $(ls analysis/blast_homology/*/Fus2_pacbio_test_merged/*_chr_*_gene_single_copy.aa_hits.gff); do
    echo $GffFile
    Chr=$(echo $GffFile | rev |cut -f1 -d'/' | rev | cut -f6 -d '_')
    $ProgDir/gff2circos_scatterplot.py --gff $GffFile --feature Chr"$Chr"_gene_homolog > tmp5/FoL_chr"$Chr"_genes.txt
  done

  # Plot location of Fus2 genes in pathogen-shared orthogroups as scatterplot
  GffFile=analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/Fus2_genes/Fus2_path_shared_genes.gff
  $ProgDir/gff2circos_scatterplot.py --gff $GffFile --feature gene --value '1' > tmp5/Fus2_path_shared_genes_plot.txt
  # Plot location of Fus2 genes in pathogen & isolate 55-shared orthogroups as scatterplot
  GffFile=analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/Fus2_genes/Fus2_path_55_shared_genes.gff
  $ProgDir/gff2circos_scatterplot.py --gff $GffFile --feature gene --value '0.5' > tmp5/Fus2_path_55_shared_genes_plot.txt
  # Plot location of Fus2 genes in pathogen & intermediate-shared orthogroups as scatterplot
  GffFile=analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/Fus2_genes/Fus2_path_inter_no_55_shared_genes.gff
  $ProgDir/gff2circos_scatterplot.py --gff $GffFile --feature gene --value '0' > tmp5/Fus2_path_inter_no_55_shared_genes_plot.txt

  circos -conf /home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos/Fus2/Fus2_circos.conf -outputdir ./tmp5
```
