# Commands to build reference databases for Blast searches in the AHDB fusarium project

LS effectors need to be identified that show divergence between ff.spp or are
unique to a particular ff.spp. As such, blast databases of reference genomes need
to be constructed and searched against.


Key genomes on this list are:

* Fusarium oxysporum
* Fusarium oxysporum f. sp. cepae
* Fusarium oxysporum f. sp. narcissi
* Fusarium oxysporum f. sp. mathioli
* Fusarium oxysporum f. sp. lycopersici
* Fusarium proliferatum
* Fusarium redolens
* Fusarium avenaceum
* Fusarium culmorum
* Fusarium solani
* Fusarium gramminearum
* Fusarium poae
* Fusarium tricinctum
* Fusarium equiseti
* Fusarium lactis
* Fusarium cerealis
* Fusarium sambucinum
* Fusarium sacchari
* Fusarium acuatum
* Fusarium fujikuroi
* Fusarium coercileum
* Fusarium flocciferum


An individual blast database was made for each organism on this list:

```bash
OutDir=analysis/AHDB_blast
mkdir -p $OutDir
printf \
"Fusarium oxysporum
Fusarium oxysporum f. sp. cepae
Fusarium oxysporum f. sp. narcissi
Fusarium oxysporum f. sp. mathioli
Fusarium oxysporum f. sp. lycopersici
Fusarium proliferatum
Fusarium redolens
Fusarium avenaceum
Fusarium culmorum
Fusarium solani
Fusarium gramminearum
Fusarium poae
Fusarium tricinctum
Fusarium equiseti
Fusarium lactis
Fusarium cerealis
Fusarium sambucinum
Fusarium sacchari
Fusarium acuatum
Fusarium fujikuroi
Fusarium coercileum
Fusarium flocciferum" \
| sed 's/Fusarium /F./g' | sed 's/f. sp./fsp/g' | sed 's/ /_/g' \
> $OutDir/target_organisms.txt
for Organism in $(cat $OutDir/target_organisms.txt); do
  echo "$Organism"
  ls -d assembly/external_group/$Organism/*
done
```


ncbi accession numbers for the genomes in van dam et al 2017 were downloaded

```bash
  cat analysis/metagenomics/reference_genomes/van_dam/van_dam_2017_S1.txt | tr -d '\n' | tr -d '\r' | tr -d '^M' | sed 's/Fusarium/monkeysFusarium/g' | sed 's/F\./monkeysF./g' | sed 's/monkeys/\n/g' > analysis/metagenomics/reference_genomes/van_dam/van_dam_2017_S1_parsed.txt
  for next in $(cat analysis/metagenomics/reference_genomes/van_dam/van_dam_2017_S1_parsed.txt | cut -f5 | grep -v '00000'); do
    echo $next
    ftp=$(cat analysis/metagenomics/reference_genomes/fusarium_genomes_ncbi.txt | grep $next)
    wget -P analysis/metagenomics/reference_genomes/van_dam "$ftp"/*.f*.gz;
  done
  cat analysis/metagenomics/reference_genomes/van_dam/van_dam_2017_S1_parsed.txt | cut -f5 | grep -v '00000' | wc -l
  ls analysis/metagenomics/reference_genomes/van_dam/*_genomic.fna.gz | grep -v -e '_rna_' -e '_cds_'|  wc -l

  mkdir analysis/metagenomics/reference_genomes/van_dam/fastq-dump
  for next in $(cat analysis/metagenomics/reference_genomes/van_dam/van_dam_2017_S1_parsed.txt | cut -f5 | grep '00000'); do
    echo $next
    fastq-dump -fasta -A $next -O analysis/metagenomics/reference_genomes/van_dam/fastq-dump
  done

  mkdir -p analysis/metagenomics/reference_genomes/van_dam/renamed
  for File in $(ls analysis/metagenomics/reference_genomes/van_dam/*_genomic.fna.gz | grep -v -e '_rna_' -e '_cds_'); do
    ID=$(basename $File | cut -f1,2 -d '_')
    Organism=$(cat analysis/metagenomics/reference_genomes/van_dam/van_dam_2017_S1_parsed.txt | grep "$ID" | cut -f1 | sed 's/ /_/g' | grep -v 'FOSC-3a')
    Strain=$(cat analysis/metagenomics/reference_genomes/van_dam/van_dam_2017_S1_parsed.txt | grep "$ID" | cut -f2 | grep -v 'FOSC-3a')
    echo "${Organism}_${Strain}_${ID}.fna"
    cat $File | gunzip -cf > analysis/metagenomics/reference_genomes/van_dam/renamed/"${Organism}_${Strain}_${ID}.fna"
  done

for File in $(ls analysis/metagenomics/reference_genomes/van_dam/fastq-dump/*.fasta); do
ID=$(basename $File | cut -f1 -d '.')
Organism=$(cat analysis/metagenomics/reference_genomes/van_dam/van_dam_2017_S1_parsed.txt | grep "$ID" | cut -f1 | sed 's/ /_/g' | grep -v 'FOSC-3a')
Strain=$(cat analysis/metagenomics/reference_genomes/van_dam/van_dam_2017_S1_parsed.txt | grep "$ID" | cut -f2 | grep -v 'FOSC-3a')
echo "${Organism}_${Strain}_${ID}.fna"
cp $File analysis/metagenomics/reference_genomes/van_dam/renamed/"${Organism}_${Strain}_${ID}.fna"
done
```

Create a list of genes from FoM and FoN LS regions
```bash
for File in $(ls analysis/mimps/F.oxysporum_fsp_*/*/*_genes_in_2kb_mimp_LS.fa); do
  Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
  cat $File | sed "s/>/>${Strain}|/g"
done > analysis/AHDB_blast/FoM_FoN_sec_mimps.fa
```

Blast vs the reference genome databases:

```bash
for RefGenome in $(ls analysis/metagenomics/reference_genomes/van_dam/renamed/*.fna); do
Prefix=$(basename $RefGenome .fna)
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 20s
printf "."
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
done
printf "\n"
OutDir=analysis/AHDB_blast/vs_ref_genomes/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh analysis/AHDB_blast/FoM_FoN_sec_mimps.fa dna $RefGenome $OutDir
done
```

Create a list of genes from FoC and FoL LS regions

```bash
Organism="F.oxysporum_fsp_cepae"
Strain="Fus2_canu_new"
MimpGff=$(ls analysis/mimps/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_genes_in_2kb_mimp_secreted.gff)
OutDir=$(dirname $MimpGff)
cat $MimpGff | grep -e 'contig_10' -e 'contig_16' -e 'contig_19' -e 'contig_21' | grep 'mRNA' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' > $OutDir/${Strain}_genes_in_2kb_mimp_secreted_headers.txt
Genes=$(ls gene_pred/final_genes/*/*/*/final_genes_combined.cdna.fasta | grep "$Strain")
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $OutDir/${Strain}_genes_in_2kb_mimp_secreted_headers.txt > $OutDir/${Strain}_genes_in_2kb_mimp_secreted.fa
cat $OutDir/${Strain}_genes_in_2kb_mimp_secreted.fa | grep '>' | wc -l
```

<!-- ```bash
Organism="F.oxysporum_fsp_lycopersici"
Strain="4287"
MimpGff=$(ls analysis/mimps/*/$Strain/*_genes_in_2kb_mimp_secreted.gff)
OutDir=$(dirname $MimpGff)
cat $MimpGff | grep -e 'contig_10' -e 'contig_16' -e 'contig_19' -e 'contig_21' | grep 'mRNA' | cut -f9 | cut -f1 -d ';' | cut -f2 -d '=' > $OutDir/${Strain}_genes_in_2kb_mimp_secreted_headers.txt
Genes=$(ls gene_pred/final_genes/*/*/*/final_genes_combined.cdna.fasta | grep "$Strain")
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $Genes --headers $OutDir/${Strain}_genes_in_2kb_mimp_secreted_headers.txt > $OutDir/${Strain}_genes_in_2kb_mimp_secreted.fa
cat $OutDir/${Strain}_genes_in_2kb_mimp_secreted.fa | grep '>' | wc -l
``` -->

Blast vs the reference genome databases:

```bash
cat analysis/mimps/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_genes_in_2kb_mimp_secreted.fa | sed 's/>/>FoC|/g' > analysis/AHDB_blast/FoM_FoN_FoC_sec_mimps.fa
cat analysis/AHDB_blast/FoM_FoN_sec_mimps.fa >> analysis/AHDB_blast/FoM_FoN_FoC_sec_mimps.fa
for RefGenome in $(ls analysis/metagenomics/reference_genomes/van_dam/renamed/*.fna); do
Prefix=$(basename $RefGenome .fna)
# Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
# while [ $Jobs -gt 1 ]; do
# sleep 20s
# printf "."
# Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
# done
# printf "\n"
OutDir=analysis/AHDB_blast/vs_ref_genomes/$Prefix
mkdir -p $OutDir
# CurDir=$PWD
# cd $OutDir
# cp -sf $CurDir/$RefGenome ${Prefix}_genome.fa
# cd $CurDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh analysis/AHDB_blast/FoM_FoN_FoC_sec_mimps.fa dna $RefGenome $OutDir
done
```

Blast vs sequenced genomes:

```bash
for RefGenome in $(ls repeat_masked/F.oxysporum_*/*/*/*_contigs_unmasked.fa | grep -v 'old' | grep -v '4287'); do
# Prefix=$(basename $RefGenome .fna)
Prefix=$(echo $RefGenome | cut -f2,3 -d '/' --output-delimiter '_')
# Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
# while [ $Jobs -gt 1 ]; do
# sleep 20s
# printf "."
# Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
# done
# printf "\n"
OutDir=analysis/AHDB_blast/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh analysis/AHDB_blast/FoM_FoN_FoC_sec_mimps.fa dna $RefGenome $OutDir
done
for RefGenome in $(ls repeat_masked/F.proliferatum/A8_ncbi/ncbi_submission/A8_contigs_unmasked.fa); do
# Prefix=$(basename $RefGenome .fna)
Prefix=$(echo $RefGenome | cut -f2,3 -d '/' --output-delimiter '_')
# Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
# while [ $Jobs -gt 1 ]; do
# sleep 20s
# printf "."
# Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
# done
# printf "\n"
OutDir=analysis/AHDB_blast/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh analysis/AHDB_blast/FoM_FoN_FoC_sec_mimps.fa dna $RefGenome $OutDir
done
```


## Summarise blast hists

```bash
CsvFiles=$(ls analysis/AHDB_blast/vs_*_genomes/*/*FoM_FoN_FoC_sec_mimps.fa_hits.csv | grep -e 'vs_ref_genomes' -e 'vs_seq_genomes' | grep -v 'F.oxysporum_/' | grep -v 'GCA_000599445.1')
Headers=$(echo $CsvFiles | sed 's&analysis/AHDB_blast/vs_ref_genomes/&&g' | sed 's&analysis/AHDB_blast/vs_seq_genomes/&&g' | sed -r "s&_FoM_FoN_FoC_sec_mimps.fa_hits.csv&&g" | sed -r "s&/\S*&&g"  | sed 's&/van_dam&&g')
OutDir=analysis/AHDB_blast/vs_ref_genomes/extracted
mkdir -p $OutDir
# ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
# $ProgDir/blast_parse.py --blast_csv $CsvFiles --headers $Headers --identity 0.50 --evalue 1e-30 > $OutDir/FoM_FoN_FoC_sec_mimps.tsv
Genomes=$(ls analysis/AHDB_blast/vs_*_genomes/*/*_genome.fa | grep -e 'vs_ref_genomes' -e 'vs_seq_genomes' | grep -v 'F.oxysporum_/' | grep -v 'GCA_000599445.1')
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/AHDB_project/blast_searches
rm analysis/AHDB_blast/vs_ref_genomes/extracted/*
$ProgDir/blast_parse_AHDB.py --blast_csv $CsvFiles --headers $Headers --genomes $Genomes --identity 0.50 --evalue 1e-30 --out_prefix $OutDir/FoM_FoN_FoC_sec_mimps
```


A number of interesting genes were identified whihc showed presence in many
Fo ff.spp., but were absent in other isolates. The orthogroups that these genes were contained in were downloaded and aligned for further study.


```bash
OutDir=analysis/AHDB_blast/FoM_FoN_FoC_sec_mimps
mkdir -p $OutDir
GeneList="FoC|g13457.t1 FoC|g16851.t1 FoC|g17143.t1 FoC|PGN_04464.t1 Stocks4|g16926.t1 Stocks4|g17043.t1 Stocks4|g18886.t1 Stocks4|g19207.t1 Stocks4|g19979.t1 FON_63|g16080.t1 FON_63|g20157.t1"
# Extracting proteins sequences
echo $GeneList | sed 's/Stocks4/FoM/g' | sed 's/FON_63/FoN/g' | sed "s/ /\n/g" > $OutDir/potential_marker_gene_headers.txt
OrthoFile=$(ls analysis/orthology/orthomcl/Fo_FoC_FoL_FoN_FoM/Fo_FoC_FoL_FoN_FoM_orthogroups.txt)
cat $OrthoFile | grep -f $OutDir/potential_marker_gene_headers.txt > $OutDir/potential_marker_gene_orthogroups.txt
for Orthogroup in $(cat $OutDir/potential_marker_gene_orthogroups.txt | cut -f1 -d ':'); do
  echo $Orthogroup
  cp analysis/orthology/orthomcl/Fo_FoC_FoL_FoN_FoM/orthogroups_fasta/$Orthogroup.fa $OutDir/.
done
# Extracting nucleotide sequences
echo $GeneList | sed "s/ /\n/g" | sed 's/|/_/g' > $OutDir/potential_marker_gene_headers.txt
InterestingFastas=$(ls analysis/AHDB_blast/vs_ref_genomes/extracted/FoM_FoN_FoC_sec_mimps*.fa | grep -f $OutDir/potential_marker_gene_headers.txt)
cp $InterestingFastas $OutDir/.
```


# Blast of all genes vs published genomes:

### for FoC FoL FoN FoM
```bash
  OutDir=analysis/AHDB_blast/vs_ref_genomes
  Taxon_code=FoC
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/*/final_genes_combined.cdna.fasta)
  cat $Fasta_file | sed "s/>/>$Taxon_code|/g" > $OutDir/FoC_FoN_FoM_FoL_genes.fa
  Taxon_code=FoL
  Fasta_file=$(ls assembly/external_group/F.oxysporum_fsp_lycopersici/4287_v2/fungidb/FungiDB-29_Foxysporum4287_AnnotatedTranscripts.fasta)
  cat $Fasta_file | sed "s/>/>$Taxon_code|/g" | cut -f1 -d ' ' >> $OutDir/FoC_FoN_FoM_FoL_genes.fa
  Taxon_code=FoN
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_narcissi/FON_63/final/final_genes_appended_renamed.cdna.fasta)
  cat $Fasta_file | sed "s/>/>$Taxon_code|/g" >> $OutDir/FoC_FoN_FoM_FoL_genes.fa
  Taxon_code=FoM
  Fasta_file=$(ls gene_pred/final_genes/F.oxysporum_fsp_mathioli/Stocks4/final/final_genes_appended_renamed.cdna.fasta)
  cat $Fasta_file | sed "s/>/>$Taxon_code|/g" >> $OutDir/FoC_FoN_FoM_FoL_genes.fa
```


Blast vs the reference genome databases:

```bash
for RefGenome in $(ls analysis/metagenomics/reference_genomes/van_dam/renamed/*.fna); do
Prefix=$(basename $RefGenome .fna)
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 20s
printf "."
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
done
printf "\n"
OutDir=analysis/AHDB_blast/vs_ref_genomes/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh analysis/AHDB_blast/vs_ref_genomes/FoC_FoN_FoM_FoL_genes.fa dna $RefGenome $OutDir
done
```


<!-- ## Summarise blast hists

```bash
CsvFiles=$(ls analysis/AHDB_blast/vs_ref_genomes/*/FoC_FoN_FoM_FoL_genes.fa_hits.csv | grep -v 'F.oxysporum_/')
Headers=$(echo $CsvFiles | sed 's&analysis/AHDB_blast/vs_ref_genomes/&&g' | sed 's&/van_dam_FoM_FoN_sec_mimps.fa_hits.csv&&g')
OutDir=analysis/AHDB_blast/vs_ref_genomes
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
$ProgDir/blast_parse.py --blast_csv $CsvFiles --headers $Headers --identity 0.70 --evalue 1e-30 > $OutDir/van_dam_FoM_FoN_sec_mimps.tsv
``` -->
