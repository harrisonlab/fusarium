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
  done
  circos -conf /home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos/Fus2/Fus2_circos.conf -outputdir ./tmp5
```


```bash
# Extract RNAseq reads aligning in 100kb windows
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
Fus2_genome=repeat_masked/F.oxysporum_fsp_cepae/Fus2/filtered_contigs_repmask/Fus2_contigs_unmasked.fa
ReadsBam=alignment/F.oxysporum_fsp_cepae/Fus2/Fus2_72hrs_rep3/accepted_hits.bam
$ProgDir/fasta2gff_windows.py --genome $Fus2_genome > tmp5/Fus2_100kb_windows.gff
bedtools coverage -abam $ReadsBam -b tmp5/Fus2_100kb_windows.gff > tmp5/tmp.bed
# Convert coverage bed files into circos format
$ProgDir/coverage_bed2circos.py --bed tmp5/tmp.bed > tmp5/Fus2_coverage_scatterplot.txt
# Identify GC content in 100kb windows
$ProgDir/gc_content2circos.py --genome $Fus2_genome --gff tmp5/Fus2_100kb_windows.gff > tmp5/Fus2_GC_scatterplot.txt

# Convert the Fus2 genome into circos format
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
Fus2_genome=assembly/canu_spades_hybrid/F.oxysporum_fsp_cepae/Fus2/filtered_contigs/Fus2_contigs_renamed.fasta
$ProgDir/fasta2circos.py --genome $Fus2_genome --contig_prefix "" > tmp5/Fus2_genome.txt
# Plot location of Blast hits as a scatterplot
for GffFile in $(ls analysis/blast_homology/*/Fus2_pacbio_test_merged/*_chr_*_gene_single_copy.aa_hits.gff); do
echo $GffFile
Chr=$(echo $GffFile | rev |cut -f1 -d'/' | rev | cut -f6 -d '_')
$ProgDir/gff2circos_scatterplot.py --gff $GffFile --feature Chr"$Chr"_gene_homolog > tmp5/FoL_chr"$Chr"_genes.txt
done
circos -conf /home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos/Fus2/Fus2_circos.conf -outputdir ./tmp5
```
