```bash
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  Fus2_genome=assembly/canu_spades_hybrid/F.oxysporum_fsp_cepae/Fus2/filtered_contigs/Fus2_contigs_renamed.fasta
  $ProgDir/fasta2circos.py --genome $Fus2_genome --contig_prefix "" > tmp5/Fus2_genome.txt
  for GffFile in $(ls analysis/blast_homology/*/Fus2_pacbio_test_merged/*_chr_*_gene_single_copy.aa_hits.gff); do
    echo $GffFile
    Chr=$(echo $GffFile | rev |cut -f1 -d'/' | rev | cut -f6 -d '_')
    $ProgDir/gff2circos_scatterplot.py --gff $GffFile --feature Chr"$Chr"_gene_homolog > tmp5/FoL_chr"$Chr"_genes.txt
  done
  circos -conf /home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos/Fus2/Fus2_circos.conf -outputdir ./tmp5
```


```bash
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  Fus2_genome=assembly/pacbio_test/F.oxysporum_fsp_cepae/Fus2_pacbio_merged_richards/filtered_contigs/Fus2_pacbio_merged_contigs_renamed.fasta
  $ProgDir/fasta2circos.py --genome $Fus2_genome --contig_prefix "" > tmp5/Fus2_genome.txt
  for GffFile in $(ls analysis/blast_homology/*/Fus2_pacbio_merged_richards/*_chr_*_gene_single_copy.aa_hits.gff); do
    echo $GffFile
    Chr=$(echo $GffFile | rev |cut -f1 -d'/' | rev | cut -f6 -d '_')
    $ProgDir/gff2circos_scatterplot.py --gff $GffFile --feature Chr"$Chr"_gene_homolog > tmp5/FoL_chr"$Chr"_genes.txt
  do
  circos -conf /home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos/Fus2/Fus2_circos.conf -outputdir ./tmp5
```


```bash
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
Fus2_genome=assembly/canu_spades_hybrid/F.oxysporum_fsp_cepae/Fus2/filtered_contigs/Fus2_contigs_renamed.fasta
ReadsBam=
  $ProgDir/fasta2gff_windows.py --genome $Fus2_genome > tmp.gff
  bedtools coverage -abam alignment/$Organism/$Strain/Fus2_72hrs_rep3/accepted_hits.bam -b tmp.gff > tmp.bed
```
