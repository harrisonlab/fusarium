


```bash

FoL_genome=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum_chromosome_and_additional_contigs.fa
# FoL_genome=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.dna.chromosome_parsed.fa
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
$ProgDir/fasta2circos.py --FoL_genome $FoL_genome > tmp4/FoL_genome.txt
CurDir=$PWD
circos -conf /home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos/circos.conf -outputdir ./tmp4
```


The position of SIX genes had been identified using the following file:

```bash
cat analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287_chromosomal/4287_chromosomal_Fo_path_genes_CRX.fa_homologs.gff | cut -f 1,4- | less -S

```


```bash
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
GffFile=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287_chromosomal/4287_chromosomal_Fo_path_genes_CRX.fa_homologs.gff
$ProgDir/gff2circosplot.py --gff $GffFile --feature SIX_homolog > tmp4/circos_graph.txt
circos -conf /home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos/circos.conf -outputdir ./tmp4
```

```bash
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
GffFile=analysis/blast_homology/F.oxysporum_fsp_lycopersici/4287_chromosomal/4287_chromosomal_Fo_path_genes_CRX.fa_homologs.gff
$ProgDir/gff2circosplot.py --gff $GffFile  --feature "SIX_homolog" --genome $FoL_genome > tmp4/circos_graph.txt
circos -conf /home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos/circos.conf -outputdir ./tmp4
```

```bash
ProgDir=~/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos
FoL_4287_gff_parsed=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31_parsed.gff3
$ProgDir/gff2circosplot.py --gff $FoL_4287_gff_parsed --feature "gene" --genome $FoL_genome > tmp4/circos_graph.txt
circos -conf /home/armita/git_repos/emr_repos/scripts/fusarium/pathogen/identify_LS_chromosomes/circos/circos.conf -outputdir ./tmp4
```
