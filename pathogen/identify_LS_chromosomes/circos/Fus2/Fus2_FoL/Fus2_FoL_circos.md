```bash
  OutDir=analysis/circos/F.oxysporum_fsp_cepae/Fus2_FoL
  mkdir -p $OutDir

  Fus2_genome=repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_unmasked.fa
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/fasta2circos.py --genome $Fus2_genome --contig_prefix "Fus2_" > $OutDir/Fus2_genome.txt

  FoL_genome=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_Genome_parsed.fasta
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/fasta2circos.py --genome $FoL_genome --contig_prefix "4287_" > $OutDir/FoL_genome.txt

  cat $OutDir/Fus2_genome.txt > $OutDir/Fus2_FoL_genome.txt
  tac $OutDir/FoL_genome.txt >> $OutDir/Fus2_FoL_genome.txt

  # cp $OutDir/Fus2_FoL_genome.txt $OutDir/Fus2_FoL_genome_edited.txt
  # cat $OutDir/Fus2_FoL_genome_edited.txt | grep -v -e 'chr23' -e 'chr24' -e 'chr25' -e 'chr26' -e 'chr27' -e 'chr28' -e 'chr29' -e 'chr30' -e 'chr31' -e 'chr32' -e 'chr33' -e 'chr34' > $OutDir/Fus2_FoL_genome_edited2.txt
  cat $OutDir/Fus2_FoL_genome.txt | grep -v 'DS231' | grep -v -e 'chr23' -e 'chr24' -e 'chr25' -e 'chr26' -e 'chr27' -e 'chr28' -e 'chr29' -e 'chr30' -e 'chr31' -e 'chr32' -e 'chr33' -e 'chr34' > $OutDir/Fus2_FoL_genome_edited.txt
```

The order of contigs was changed manually using nano
```bash
cp $OutDir/Fus2_FoL_genome_edited.txt $OutDir/Fus2_FoL_genome_edited2.txt
nano $OutDir/Fus2_FoL_genome_edited2.txt
cp $OutDir/Fus2_FoL_genome_edited2.txt $OutDir/Fus2_FoL_genome_final.txt
```


```bash
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/orthology2circos_ribbons.py --orthology analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL_publication_ncbi/FoC_vs_Fo_vs_FoL_publication_ncbi_orthogroups.txt --name1 Fus2 --gff1 gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3 --name2 4287 --gff2 assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_parsed.gff   > $OutDir/Fus2_FoL_links.txt
  # Links to FoL LS contigs 3, 6, 14 and 15 were coloured black
  # cat $OutDir/Fus2_FoL_links.txt \
  #   | sed '/4287_CM000591.1/ s/$/\tcolor=black/' \
  #   | sed '/4287_CM000594.1/ s/$/\tcolor=black/' \
  #   | sed '/4287_CM000602.2/ s/$/\tcolor=black/' \
  #   | sed '/4287_CM000603.1/ s/$/\tcolor=black/' \
  #   > $OutDir/Fus2_FoL_links_edited.txt
  cat $OutDir/Fus2_FoL_links.txt > $OutDir/Fus2_FoL_links_edited.txt
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  circos -conf $ProgDir/Fus2/Fus2_FoL/Fus2_FoL_circos.conf -outputdir $OutDir
  mv $OutDir/circos.png $OutDir/Fus2_FoL_circos.png
  mv $OutDir/circos.svg $OutDir/Fus2_FoL_circos.svg
```



```bash
OutDir=analysis/circos/F.oxysporum_fsp_cepae/Fus2_FoL
cat $OutDir/Fus2_FoL_genome_edited2.txt | grep -v '4287' > $OutDir/Fus2_FoL_genome_final.txt
mkdir -p $OutDir/by_FoC_chr
for Num in $(seq 1 22); do
  Chr="contig_"$Num"_pilon"
  echo "$Chr"
  OrthologyTxt=analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL_publication_ncbi/FoC_vs_Fo_vs_FoL_publication_ncbi_orthogroups.txt
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/orthology2ribbons_internal.py \
  --chr1 $Chr \
  --orthology $OrthologyTxt \
  --name1 Fus2 \
  --gff1 gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3 \
  | sort | uniq \
  > $OutDir/Fus2_FoL_links_edited.txt
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  circos -conf $ProgDir/Fus2/Fus2_FoL/Fus2_FoL_circos.conf -outputdir $OutDir
  mv $OutDir/circos.png $OutDir/by_FoC_chr/Fus2_FoL_LS_"$Chr"_circos.png
  mv $OutDir/circos.svg $OutDir/by_FoC_chr/Fus2_FoL_LS_"$Chr"_circos.svg
done
```
<!--
```bash
OutDir=analysis/circos/F.oxysporum_fsp_cepae/Fus2_FoL
cat $OutDir/Fus2_FoL_genome_edited2.txt | grep '4287' > $OutDir/Fus2_FoL_genome_final.txt
mkdir -p $OutDir/by_FoL_chr
for Num in $(seq 1 15); do
  PrevLinks=analysis/circos/F.oxysporum_fsp_cepae/Fus2_FoL/Fus2_FoL_genome_edited.txt
  Chr=$(cat $PrevLinks | grep '4287' | grep -w "chr$Num" | cut -f3 -d ' ' | sed 's/4287_//g')
  echo "$Chr"
  OrthologyTxt=analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL_publication_ncbi/FoC_vs_Fo_vs_FoL_publication_ncbi_orthogroups.txt
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/orthology2ribbons_internal.py \
  --chr1 $Chr \
  --orthology $OrthologyTxt \
  --name1 4287 \
  --gff1 assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_parsed.gff \
  | sort | uniq \
  > $OutDir/Fus2_FoL_links_edited.txt
  ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  circos -conf $ProgDir/Fus2/Fus2_FoL/Fus2_FoL_circos.conf -outputdir $OutDir
  mv $OutDir/circos.png $OutDir/by_FoL_chr/Fus2_FoL_LS_"$Chr"_circos.png
  mv $OutDir/circos.svg $OutDir/by_FoL_chr/Fus2_FoL_LS_"$Chr"_circos.svg
done
``` -->

<!--
```bash
OutDir=analysis/circos/F.oxysporum_fsp_cepae/Fus2_FoL
# cat gene_pred/annotations/F.oxysporum_fsp_lycopersici/4287/4287_gene_annotations.tab | grep -e 'transpos' | cut -f16 | sort | uniq > $OutDir/FoL_transposase_orthogroups.txt
OrthologyTxt=analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL_publication/FoC_vs_Fo_vs_FoL_publication_ncbi_orthogroups.txt
# OrthologyTxt_ed=$OutDir/FoC_vs_Fo_vs_FoL_orthogroups_ed_no_transposase.txt
# cat $OrthologyTxt | grep -v -w  -f $OutDir/FoL_transposase_orthogroups.txt > $OrthologyTxt_ed
for Num in $(seq 1 15); do
PrevLinks=analysis/circos/F.oxysporum_fsp_cepae/Fus2_FoL/Fus2_FoL_genome_edited.txt
Chr=$(cat $PrevLinks | grep '4287' | grep -w "chr$Num" | cut -f3 -d ' ')
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/orthology2ribons_circos_by_chr.py \
--chr1 $Chr \
--orthology $OrthologyTxt \
--name1 4287 \
--gff1 assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_parsed.gff \
--name2 Fus2 \
--gff2 gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3 \
> $OutDir/Fus2_FoL_LS_links.txt

ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
circos -conf $ProgDir/Fus2/Fus2_FoL/Fus2_FoL_circos.conf -outputdir $OutDir
mkdir -p $OutDir/by_FoL_chr
mv $OutDir/circos.png $OutDir/by_FoL_chr/Fus2_FoL_LS_"$Chr"_circos.png
mv $OutDir/circos.svg $OutDir/by_FoL_chr/Fus2_FoL_LS_"$Chr"_circos.svg
done
```


```bash
OutDir=analysis/circos/F.oxysporum_fsp_cepae/Fus2_FoL
# cat gene_pred/annotations/F.oxysporum_fsp_lycopersici/4287/4287_gene_annotations.tab | grep -e 'transpos' | cut -f16 | sort | uniq > $OutDir/FoL_transposase_orthogroups.txt
OrthologyTxt=analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL_publication/FoC_vs_Fo_vs_FoL_publication_orthogroups.txt
# OrthologyTxt_ed=$OutDir/FoC_vs_Fo_vs_FoL_orthogroups_ed_no_transposase.txt
# cat $OrthologyTxt | grep -v -w  -f $OutDir/FoL_transposase_orthogroups.txt > $OrthologyTxt_ed
for Chr in $(seq 1 34); do
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/orthology2ribons_circos_by_chr.py \
--chr1 $Chr \
--orthology $OrthologyTxt \
--name1 Fus2 \
--gff1 gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3 \
--name2 4287 \
--gff2 assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.gff3 \
> $OutDir/Fus2_FoL_LS_links.txt

ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
circos -conf $ProgDir/Fus2/Fus2_FoL/Fus2_FoL_circos.conf -outputdir $OutDir
mkdir -p $OutDir/by_FoC_chr
mv $OutDir/circos.png $OutDir/by_FoC_chr/Fus2_FoL_LS_"$Chr"_circos.png
mv $OutDir/circos.svg $OutDir/by_FoC_chr/Fus2_FoL_LS_"$Chr"_circos.svg
done
``` -->
