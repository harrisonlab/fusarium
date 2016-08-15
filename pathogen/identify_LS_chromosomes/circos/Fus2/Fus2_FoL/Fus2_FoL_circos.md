```bash
  OutDir=analysis/circos/F.oxysporum_fsp_cepae/Fus2_FoL
  mkdir -p $OutDir

  Fus2_genome=assembly/canu_spades_hybrid/F.oxysporum_fsp_cepae/Fus2/filtered_contigs/Fus2_contigs_renamed.fasta
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/fasta2circos.py --genome $Fus2_genome --contig_prefix "Fus2_" > $OutDir/Fus2_genome.txt

  FoL_genome=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.dna.chromosome.fa
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/fasta2circos.py --genome $FoL_genome --contig_prefix "4287_" > $OutDir/FoL_genome.txt

  cat $OutDir/FoL_genome.txt > $OutDir/Fus2_FoL_genome.txt
  tac $OutDir/Fus2_genome.txt >> $OutDir/Fus2_FoL_genome.txt

  # The file $OutDir/Fus2_FoL_genome.txt was maually edited:
  # It was called: $OutDir/Fus2_FoL_genome_edit.txt
```


```bash
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/ortholology2circos_ribbons.py --orthology analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/FoC_vs_Fo_vs_FoL_orthogroups.txt --name1 Fus2 --gff1 gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2/final/final_genes_appended.gff3 --name2 4287 --gff2 assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.gff3 > $OutDir/Fus2_FoL_links.txt
  # Links to FoL LS contigs 3, 6, 14 and 15 were coloured balck
  cat $OutDir/Fus2_FoL_links.txt \
    | sed '/4287_3/ s/$/\tcolor=black/' \
    | sed '/4287_6/ s/$/\tcolor=black/' \
    | sed '/4287_14/ s/$/\tcolor=black/' \
    | sed '/4287_15/ s/$/\tcolor=black/' \
    > $OutDir/Fus2_FoL_links_edited.txt
```

```bash
OutDir=analysis/circos/F.oxysporum_fsp_cepae/Fus2_FoL
OrthologyTxt=analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/FoC_vs_Fo_vs_FoL_orthogroups.txt
OrthologyTxt_ed=$OutDir/FoC_vs_Fo_vs_FoL_orthogroups_ed.txt
cat $OrthologyTxt | tail -n+100 > $OrthologyTxt_ed
for Chr in $(seq 1 15); do
  $ProgDir/orthology2ribons_circos_by_chr.py \
    --chr1 $Chr \
    --orthology $OrthologyTxt \
    --name1 4287 \
    --gff1 assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.gff3 \
    --name2 Fus2 \
    --gff2 gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2/final/final_genes_appended.gff3 \
    > $OutDir/Fus2_FoL_LS_links.txt

    ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
    circos -conf $ProgDir/Fus2/Fus2_FoL/Fus2_FoL_circos.conf -outputdir $OutDir
    mv $OutDir/circos.png $OutDir/Fus2_FoL_LS_"$Chr"_circos.png
  done
```

```bash
OutDir=analysis/circos/F.oxysporum_fsp_cepae/Fus2_FoL
cat gene_pred/annotations/F.oxysporum_fsp_lycopersici/4287/4287_gene_annotations.tab | grep -e 'transpos' | cut -f16 | sort | uniq > $OutDir/FoL_transposase_orthogroups.txt
OrthologyTxt=analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/FoC_vs_Fo_vs_FoL_orthogroups.txt
OrthologyTxt_ed=$OutDir/FoC_vs_Fo_vs_FoL_orthogroups_ed_no_transposase.txt
cat $OrthologyTxt | grep -v -w  -f $OutDir/FoL_transposase_orthogroups.txt > $OrthologyTxt_ed
for Chr in $(seq 1 15); do
  ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
  $ProgDir/orthology2ribons_circos_by_chr.py \
    --chr1 $Chr \
    --orthology $OrthologyTxt \
    --name1 4287 \
    --gff1 assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.gff3 \
    --name2 Fus2 \
    --gff2 gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2/final/final_genes_appended.gff3 \
    > $OutDir/Fus2_FoL_LS_links.txt

    ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
    circos -conf $ProgDir/Fus2/Fus2_FoL/Fus2_FoL_circos.conf -outputdir $OutDir
    mv $OutDir/circos.png $OutDir/Fus2_FoL_LS_"$Chr"_circos_no_transposase.png
  done
```


```bash
# for Orthogroup in $(seq 1 10); do
for Orthogroup in 68 22 36 34; do
echo "Orthogroup - $Orthogroup"
OutDir=analysis/circos/F.oxysporum_fsp_cepae/Fus2_FoL/orthogroup"$Orthogroup"
mkdir -p $OutDir
OrthologyTxt=analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/FoC_vs_Fo_vs_FoL_orthogroups.txt
OrthologyTxt_ed=$OutDir/FoC_vs_Fo_vs_FoL_orthogroups_ed_orthogroup"Orthogroup".txt
cat $OrthologyTxt | head -n $Orthogroup | tail -n1 > $OrthologyTxt_ed
for Chr in $(seq 1 15); do
echo "Chr - $Chr"
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/orthology2ribons_circos_by_chr.py \
--chr1 $Chr \
--orthology $OrthologyTxt_ed \
--name1 4287 \
--gff1 assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.gff3 \
--name2 Fus2 \
--gff2 gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2/final/final_genes_appended.gff3 \
> analysis/circos/F.oxysporum_fsp_cepae/Fus2_FoL/Fus2_FoL_LS_links.txt

ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
circos -conf $ProgDir/Fus2/Fus2_FoL/Fus2_FoL_circos.conf -outputdir $OutDir
mv $OutDir/circos.png $OutDir/Fus2_FoL_LS_"$Chr"_circos_Orthogroup"$orthogroup".png
done
done
```
