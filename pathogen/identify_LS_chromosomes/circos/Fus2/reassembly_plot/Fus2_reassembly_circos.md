Circos plots were generated for reassembly of the Fus2 genome

```bash

  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  Fus2_genome=assembly/canu-1.3/F.oxysporum_fsp_cepae/Fus2_canu/pilon/pilon.fasta
  OutDir=analysis/circos/F.oxysporum_fsp_cepae/Fus2_reassembly
  mkdir -p $OutDir

  # Convert the Fus2 genome into circos format
  $ProgDir/fasta2circos.py --genome $Fus2_genome --contig_prefix "" > $OutDir/Fus2_genome.txt

  # Plot location of FoL gene Blast hits as a scatterplot
  for GffFile in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_canu/*_chr_*_gene_single_copy.aa_hits.gff); do
  # for GffFile in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2_edited2/Fus2_edited2_chr_*_gene_single_copy.aa_hits.gff); do
    echo $GffFile
    Chr=$(echo $GffFile | rev |cut -f1 -d'/' | rev | cut -f4 -d '_')
    $ProgDir/gff2circos_scatterplot.py --gff $GffFile --feature Chr"$Chr"_gene_homolog > $OutDir/FoL_chr"$Chr"_genes.txt
  done

  circos -conf /home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos/Fus2/reassembly_plot/Fus2_reassembly_circos.conf -outputdir ./$OutDir
```
