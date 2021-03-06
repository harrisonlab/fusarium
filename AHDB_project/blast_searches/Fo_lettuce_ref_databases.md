# Commands to build reference databases for Blast searches in the AHDB fusarium project

LS effectors need to be identified that show divergence between ff.spp or are
unique to a particular ff.spp. As such, blast databases of reference genomes need
to be constructed and searched against.

The same database was used for F. oxysporum ex. lettuce as that described in the
ref_database.md document within this git repository.

<!--
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

## Download assemblies

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


### Download additional ncbi genomes


Additional F. solani (nectria haematococca) genomes were downloaded:

```bash
  ProjectDir=/home/groups/harrisonlab/project_files/fusarium
  cd $ProjectDir
  OutDir=analysis/metagenomics/reference_genomes/additional
  Species="F.solani"
  Strain="JS-169"
  mkdir -p $OutDir/$Species/$Strain
  wget -P $OutDir/$Species/$Strain ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/MS/JJ/MSJJ02/MSJJ02.1.fsa_nt.gz
  Species="F.solani"
  Strain="IMV_00293"
  mkdir -p $OutDir/$Species/$Strain
  wget -P $OutDir/$Species/$Strain ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/NG/ZQ/NGZQ01/NGZQ01.1.fsa_nt.gz
  Species="F.solani"
  Strain="77-13-4"
  mkdir -p $OutDir/$Species/$Strain
  wget -P $OutDir/$Species/$Strain http://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/ac/ACJF01.fasta.gz
  gunzip $OutDir/*/*/*.f*.gz
  for File in $(ls $OutDir/*/*/*.f* | grep -e 'fsa' -e 'fasta'); do
    Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
    mv $File analysis/metagenomics/reference_genomes/renamed/"$Organism"_"$Strain".fna
  done
```
-->

## FoN FoM analysis

### Prepare queries FoN FoM

Create a list of secreted genes within 2Kb of a mimp from F. oxysporum ex.
lettuce.


```bash
for File in $(ls analysis/mimps/F.oxysporum_fsp_lactucae/*/*_genes_in_2kb_mimp.fa); do
  Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
  cat $File | sed "s/>/>${Strain}|/g"
done > analysis/AHDB_blast/FoM_FoN_sec_mimps.fa
```

### Perform BLAST

Blast vs the reference genome databases:

```bash
mkdir -p analysis/AHDB_Fo_lactucae_blast
for File in $(ls analysis/mimps/F.oxysporum_fsp_lactucae/*/*_genes_in_2kb_mimp_secreted.fa); do
  Strain=$(echo $File | cut -f4 -d '/')
  cat $File | sed "s/>/>${Strain}_/g"
done > analysis/AHDB_Fo_lactucae_blast/Fo_lactucae_sec_mimps.fa

for RefGenome in $(ls analysis/metagenomics/reference_genomes/van_dam/renamed/*.fna); do
Prefix=$(basename $RefGenome .fna)
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 20s
printf "."
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
done
printf "\n"
Query=$(ls analysis/AHDB_Fo_lactucae_blast/Fo_lactucae_sec_mimps.fa)
OutDir=analysis/AHDB_Fo_lactucae_blast/vs_ref_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
cp -sf $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query dna $RefGenome $OutDir
done
```

Blast vs additional ncbi genomes:

```bash
for RefGenome in $(ls analysis/metagenomics/reference_genomes/renamed/*.fna | grep 'solani' | grep -v '77-13-4'); do
Prefix=$(basename $RefGenome .fna)
echo $Prefix
OutDir=analysis/AHDB_Fo_lactucae_blast/vs_ref_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
cp -sf $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
Query=$(ls analysis/AHDB_Fo_lactucae_blast/Fo_lactucae_sec_mimps.fa)
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query dna $RefGenome $OutDir
done
```

Blast vs sequenced genomes:

```bash
for RefGenome in $(ls repeat_masked/F.oxysporum_*/*/*/*_contigs_unmasked.fa | grep -v 'old' | grep -v '4287' | grep -v -e 'FON139' -e 'FON77' -e 'FON89' -e 'FON81' -e 'FON129'); do
# Prefix=$(basename $RefGenome .fna)
Prefix=$(echo $RefGenome | cut -f2,3 -d '/' --output-delimiter '_')
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 20s
printf "."
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
done
printf "\n"
OutDir=analysis/AHDB_Fo_lactucae_blast/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
Query=$(ls analysis/AHDB_Fo_lactucae_blast/Fo_lactucae_sec_mimps.fa)
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query dna $RefGenome $OutDir
done

for RefGenome in $(ls repeat_masked/F.proliferatum/A8/*/*_contigs_unmasked.fa); do
Prefix=$(echo $RefGenome | cut -f2,3 -d '/' --output-delimiter '_')
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
# while [ $Jobs -gt 1 ]; do
# sleep 20s
# printf "."
# Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
# done
# printf "\n"
OutDir=analysis/AHDB_Fo_lactucae_blast/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
Query=$(ls analysis/AHDB_Fo_lactucae_blast/Fo_lactucae_sec_mimps.fa)
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query dna $RefGenome $OutDir
done

for RefGenome in $(ls ../fusarium_ex_narcissus/repeat_masked/F.oxysporum_*/*/*/*_contigs_unmasked.fa | grep -e 'FON139' -e 'FON77' -e 'FON89' -e 'FON81' -e 'FON129'); do
Prefix=$(echo $RefGenome | cut -f4,5 -d '/' --output-delimiter '_')
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 20s
printf "."
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
done
printf "\n"
OutDir=analysis/AHDB_Fo_lactucae_blast/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
Query=$(ls analysis/AHDB_Fo_lactucae_blast/Fo_lactucae_sec_mimps.fa)
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query dna $RefGenome $OutDir
done
```
<!--
Blast vs work in progress genomes:

```bash
for RefGenome in $(ls assembly/spades/F.*/*/filtered_contigs/contigs_min_500bp*.fasta | grep -w -e 'Straw465' -e 'F81' -e 'FOP1-EMR' -e 'R2' -e '15-074' -e 'A1-2' -e 'HB6' -e 'D2' -e 'PG8' -e 'L5' -e 'A1-2' -e 'HB6' -e 'D2' -e 'PG8' -e 'L5' -e 'A1-2' -e 'HB6'); do
Prefix=$(echo $RefGenome | cut -f3,4 -d '/' --output-delimiter '_')
echo $Prefix
OutDir=analysis/AHDB_blast/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
rm ${Prefix}_genome.fa
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh analysis/AHDB_blast/FoM_FoN_FoC_sec_mimps.fa dna $RefGenome $OutDir
done
for RefGenome in $(ls assembly/spades_pacbio/F.*/*/filtered_contigs/contigs_min_500bp.fasta | grep -w -e '55'); do
Prefix=$(echo $RefGenome | cut -f3,4 -d '/' --output-delimiter '_')
echo $Prefix
OutDir=analysis/AHDB_blast/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
rm ${Prefix}_genome.fa
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh analysis/AHDB_blast/FoM_FoN_FoC_sec_mimps.fa dna $RefGenome $OutDir
done
``` -->

## Summarise blast hists

```bash
# CsvFiles=$(ls analysis/AHDB_blast/vs_*_genomes/*/*FoM_FoN_FoC_sec_mimps.fa_hits.csv | grep -e 'vs_ref_genomes' -e 'vs_seq_genomes' | grep -v 'F.oxysporum_/' | grep -v 'GCA_000599445.1' | grep -v -e 'Fusarium_asiaticum' -e 'Fusarium_azukicol' -e 'Fusarium_brasiliense_NRRL_31757/' -e 'Fusarium_circinatum' -e 'Fusarium_cuneirostrum' -e 'Fusarium_euwallaceae' -e 'Fusarium_fujikuroi' -e 'Fusarium_graminearum' -e 'Fusarium_langsethiae' -e 'Fusarium_meridionale' -e 'Fusarium_nygamai' -e 'Fusarium_oxysporum_f._melongenae' -e 'Fusarium_oxysporum_f._sp._ciceris' -e 'Fusarium_oxysporum_f._sp._lycopersici' -e 'Fusarium_oxysporum_f._sp._medicaginis' -e 'Fusarium_oxysporum_f._sp._niveum' -e 'Fusarium_phaseoli' -e 'Fusarium_praegraminearum' -e 'Fusarium_pseudograminearum' -e 'Fusarium_temperatum' -e 'Fusarium_tucumaniae' -e 'Fusarium_udum' -e 'Fusarium_verticillioides' -e 'Fusarium_virguliforme' -e 'Fusarium_oxysporum_f._sp._conglutinans')
CsvFiles=$(ls analysis/AHDB_Fo_lactucae_blast/vs_*_genomes/*/*Fo_lactucae_sec_mimps.fa_hits.csv | grep -e 'vs_ref_genomes' -e 'vs_seq_genomes' | grep -v 'F.oxysporum_/' | grep -v -e 'GCA_000599445.1' -e '77-13-4')
# CsvFiles=$(ls analysis/AHDB_blast/vs_*_genomes/*/*FoM_FoN_FoC_sec_mimps.fa_hits.csv | grep -e 'vs_ref_genomes' | grep -v 'F.oxysporum_/' | grep -v 'GCA_000599445.1')
Headers=$(echo $CsvFiles | sed 's&analysis/AHDB_Fo_lactucae_blast/vs_ref_genomes/&&g' | sed 's&analysis/AHDB_Fo_lactucae_blast/vs_seq_genomes/&&g' | sed -r "s&_Fo_lactucae_sec_mimps.fa_hits.csv&&g" | sed -r "s&/\S*&&g"  | sed 's&/van_dam&&g')
OutDir=analysis/AHDB_Fo_lactucae_blast/vs_ref_genomes/extracted
mkdir -p $OutDir
# ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
# $ProgDir/blast_parse.py --blast_csv $CsvFiles --headers $Headers --identity 0.50 --evalue 1e-30 > $OutDir/FoM_FoN_FoC_sec_mimps.tsv
# Genomes=$(ls analysis/AHDB_blast/vs_*_genomes/*/*_genome.fa | grep -e 'vs_ref_genomes' -e 'vs_seq_genomes' | grep -v 'F.oxysporum_/' | grep -v 'GCA_000599445.1' | grep -v -e 'Fusarium_asiaticum' -e 'Fusarium_azukicol' -e 'Fusarium_brasiliense_NRRL_31757/' -e 'Fusarium_circinatum' -e 'Fusarium_cuneirostrum' -e 'Fusarium_euwallaceae' -e 'Fusarium_fujikuroi' -e 'Fusarium_graminearum' -e 'Fusarium_langsethiae' -e 'Fusarium_meridionale' -e 'Fusarium_nygamai' -e 'Fusarium_oxysporum_f._melongenae' -e 'Fusarium_oxysporum_f._sp._ciceris' -e 'Fusarium_oxysporum_f._sp._lycopersici' -e 'Fusarium_oxysporum_f._sp._medicaginis' -e 'Fusarium_oxysporum_f._sp._niveum' -e 'Fusarium_phaseoli' -e 'Fusarium_praegraminearum' -e 'Fusarium_pseudograminearum' -e 'Fusarium_temperatum' -e 'Fusarium_tucumaniae' -e 'Fusarium_udum' -e 'Fusarium_verticillioides' -e 'Fusarium_virguliforme' -e 'Fusarium_oxysporum_f._sp._conglutinans')
Genomes=$(ls analysis/AHDB_Fo_lactucae_blast/vs_*_genomes/*/*_genome.fa | grep -e 'vs_ref_genomes' -e 'vs_seq_genomes' | grep -v 'F.oxysporum_/' | grep -v 'GCA_000599445.1' | grep -e 'vs_ref_genomes' -e 'vs_seq_genomes' | grep -v 'F.oxysporum_/' | grep -v -e 'GCA_000599445.1' -e '77-13-4')
# Genomes=$(ls analysis/AHDB_blast/vs_*_genomes/*/*_genome.fa | grep -e 'vs_ref_genomes' | grep -v 'F.oxysporum_/' | grep -v 'GCA_000599445.1')
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/AHDB_project/blast_searches
# rm analysis/AHDB_blast/vs_ref_genomes/extracted/*
$ProgDir/blast_parse_AHDB.py --blast_csv $CsvFiles --headers $Headers --genomes $Genomes --identity 0.50 --evalue 1e-30 --out_prefix $OutDir/FoL_sec_mimps
```

<!--
A number of interesting genes were identified whihc showed presence in many
Fo ff.spp., but were absent in other isolates. The orthogroups that these genes were contained in were downloaded and aligned for further study.


```bash
# OutDir=analysis/AHDB_blast/FoM_FoN_FoC_sec_mimps
OutDir=analysis/AHDB_blast/FoM_FoN_FoC_sec_mimps_additional_stocks
mkdir -p $OutDir
# GeneList="FoC|g13457.t1 FoC|g16851.t1 FoC|g17143.t1 FoC|PGN_04464.t1 Stocks4|g16926.t1 Stocks4|g17043.t1 Stocks4|g18886.t1 Stocks4|g19207.t1 Stocks4|g19979.t1 FON_63|g16080.t1 FON_63|g20157.t1"
# GeneList="Stocks4|g18921.t1"
GeneList="Stocks4|g16885.t1 Stocks4|g16926.t1 Stocks4|g16947.t1 Stocks4|g16948.t1 Stocks4|g16953.t1 Stocks4|g17043.t1 Stocks4|g17165.t1 Stocks4|g17244.t1 Stocks4|g18886.t1 Stocks4|g18920.t1 Stocks4|g18921.t1 Stocks4|g19207.t1 Stocks4|g19979.t1"

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
 -->

<!-- ## Summarise blast hists

```bash
CsvFiles=$(ls analysis/AHDB_blast/vs_ref_genomes/*/FoC_FoN_FoM_FoL_genes.fa_hits.csv | grep -v 'F.oxysporum_/')
Headers=$(echo $CsvFiles | sed 's&analysis/AHDB_blast/vs_ref_genomes/&&g' | sed 's&/van_dam_FoM_FoN_sec_mimps.fa_hits.csv&&g')
OutDir=analysis/AHDB_blast/vs_ref_genomes
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
$ProgDir/blast_parse.py --blast_csv $CsvFiles --headers $Headers --identity 0.70 --evalue 1e-30 > $OutDir/van_dam_FoM_FoN_sec_mimps.tsv
``` -->



## lettuce and statice analysis

### Prepare queries lettuce and statice

Create a list of secreted genes within 2Kb of a mimp from F. oxysporum ex. lettuce and statice



### Perform BLAST

Blast vs the reference genome databases:

```bash
mkdir -p analysis/AHDB_Fo_lactucae_statice_blast2
for File in $(ls analysis/mimps/F.oxysporum_fsp_*/*/*_genes_in_2kb_mimp_secreted.fa | grep -e 'AJ516' -e 'AJ520' -e 'Stat10'); do
  Strain=$(echo $File | cut -f4 -d '/')
  cat $File | sed "s/>/>${Strain}_/g"
done > analysis/AHDB_Fo_lactucae_statice_blast/Fo_lactucae_sec_mimps.fa

for RefGenome in $(ls analysis/metagenomics/reference_genomes/van_dam/renamed/*.fna | tail -n 60); do
Prefix=$(basename $RefGenome .fna)
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 20s
printf "."
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
done
printf "\n"
Query=$(ls analysis/AHDB_Fo_lactucae_statice_blast/Fo_lactucae_sec_mimps.fa)
OutDir=analysis/AHDB_Fo_lactucae_statice_blast/vs_ref_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
cp -sf $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query dna $RefGenome $OutDir
done
```

Blast vs additional ncbi genomes:

```bash
for RefGenome in $(ls analysis/metagenomics/reference_genomes/renamed/*.fna | grep 'solani' | grep -v '77-13-4'); do
Prefix=$(basename $RefGenome .fna)
echo $Prefix
OutDir=analysis/AHDB_Fo_lactucae_statice_blast/vs_ref_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
cp -sf $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
Query=$(ls analysis/AHDB_Fo_lactucae_statice_blast/Fo_lactucae_sec_mimps.fa)
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query dna $RefGenome $OutDir
done
```

Blast vs sequenced genomes:

```bash
for RefGenome in $(ls repeat_masked/F.oxysporum_*/*/*/*_contigs_unmasked.fa | grep -v 'old' | grep -v '4287' | grep -v -e 'FON139' -e 'FON77' -e 'FON89' -e 'FON81' -e 'FON129'); do
# Prefix=$(basename $RefGenome .fna)
Prefix=$(echo $RefGenome | cut -f2,3 -d '/' --output-delimiter '_')
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 20s
printf "."
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
done
printf "\n"
OutDir=analysis/AHDB_Fo_lactucae_statice_blast/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
Query=$(ls analysis/AHDB_Fo_lactucae_statice_blast/Fo_lactucae_sec_mimps.fa)
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query dna $RefGenome $OutDir
done

for RefGenome in $(ls repeat_masked/F.proliferatum/A8/*/*_contigs_unmasked.fa); do
Prefix=$(echo $RefGenome | cut -f2,3 -d '/' --output-delimiter '_')
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
# while [ $Jobs -gt 1 ]; do
# sleep 20s
# printf "."
# Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
# done
# printf "\n"
OutDir=analysis/AHDB_Fo_lactucae_statice_blast/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
Query=$(ls analysis/AHDB_Fo_lactucae_statice_blast/Fo_lactucae_sec_mimps.fa)
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query dna $RefGenome $OutDir
done

for RefGenome in $(ls ../fusarium_ex_narcissus/repeat_masked/F.oxysporum_*/*/*/*_contigs_unmasked.fa | grep -e 'FON139' -e 'FON77' -e 'FON89' -e 'FON81' -e 'FON129'); do
Prefix=$(echo $RefGenome | cut -f4,5 -d '/' --output-delimiter '_')
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 20s
printf "."
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
done
printf "\n"
OutDir=analysis/AHDB_Fo_lactucae_statice_blast/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
Query=$(ls analysis/AHDB_Fo_lactucae_statice_blast/Fo_lactucae_sec_mimps.fa)
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query dna $RefGenome $OutDir
done
```
<!--
Blast vs work in progress genomes:

```bash
for RefGenome in $(ls assembly/spades/F.*/*/filtered_contigs/contigs_min_500bp*.fasta | grep -w -e 'Straw465' -e 'F81' -e 'FOP1-EMR' -e 'R2' -e '15-074' -e 'A1-2' -e 'HB6' -e 'D2' -e 'PG8' -e 'L5' -e 'A1-2' -e 'HB6' -e 'D2' -e 'PG8' -e 'L5' -e 'A1-2' -e 'HB6'); do
Prefix=$(echo $RefGenome | cut -f3,4 -d '/' --output-delimiter '_')
echo $Prefix
OutDir=analysis/AHDB_blast/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
rm ${Prefix}_genome.fa
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh analysis/AHDB_blast/FoM_FoN_FoC_sec_mimps.fa dna $RefGenome $OutDir
done
for RefGenome in $(ls assembly/spades_pacbio/F.*/*/filtered_contigs/contigs_min_500bp.fasta | grep -w -e '55'); do
Prefix=$(echo $RefGenome | cut -f3,4 -d '/' --output-delimiter '_')
echo $Prefix
OutDir=analysis/AHDB_blast/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
rm ${Prefix}_genome.fa
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh analysis/AHDB_blast/FoM_FoN_FoC_sec_mimps.fa dna $RefGenome $OutDir
done
``` -->

## Summarise blast hists

```bash
CsvFiles=$(ls analysis/AHDB_Fo_lactucae_statice_blast/vs_*_genomes/*/*Fo_lactucae_sec_mimps.fa_hits.csv | grep -e 'vs_ref_genomes' -e 'vs_seq_genomes' | grep -v 'F.oxysporum_/' | grep -v -e 'GCA_000599445.1' -e '77-13-4')
Headers=$(echo $CsvFiles | sed 's&analysis/AHDB_Fo_lactucae_statice_blast/vs_ref_genomes/&&g' | sed 's&analysis/AHDB_Fo_lactucae_statice_blast/vs_seq_genomes/&&g' | sed -r "s&_Fo_lactucae_sec_mimps.fa_hits.csv&&g" | sed -r "s&/\S*&&g"  | sed 's&/van_dam&&g')
OutDir=analysis/AHDB_Fo_lactucae_statice_blast/vs_ref_genomes/extracted
mkdir -p $OutDir
Genomes=$(ls analysis/AHDB_Fo_lactucae_statice_blast/vs_*_genomes/*/*_genome.fa | grep -e 'vs_ref_genomes' -e 'vs_seq_genomes' | grep -v 'F.oxysporum_/' | grep -v 'GCA_000599445.1' | grep -e 'vs_ref_genomes' -e 'vs_seq_genomes' | grep -v 'F.oxysporum_/' | grep -v -e 'GCA_000599445.1' -e '77-13-4')
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/AHDB_project/blast_searches
$ProgDir/blast_parse_AHDB.py --blast_csv $CsvFiles --headers $Headers --genomes $Genomes --identity 0.50 --evalue 1e-30 --out_prefix $OutDir/FoL_sec_mimpsAHDB
```



## lettuce analysis 2

### Prepare queries lettuce

Create a list of secreted genes within 2Kb of a mimp from F. oxysporum ex. lettuce and statice

### Perform BLAST

Blast vs the reference genome databases:

```bash
mkdir -p analysis/AHDB_Fo_lactucae_statice_blast2
for File in $(ls analysis/mimps/F.oxysporum_fsp_*/*/*_genes_in_2kb_mimp_secreted.fa | grep -e 'AJ516' -e 'AJ520'); do
  Strain=$(echo $File | cut -f4 -d '/')
  cat $File | sed "s/>/>${Strain}_/g"
done > analysis/AHDB_Fo_lactucae_statice_blast2/Fo_lactucae_sec_mimps.fa

for RefGenome in $(ls analysis/metagenomics/reference_genomes/van_dam/renamed/*.fna | tail -n 60); do
Prefix=$(basename $RefGenome .fna)
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 20s
printf "."
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
done
printf "\n"
Query=$(ls analysis/AHDB_Fo_lactucae_statice_blast2/Fo_lactucae_sec_mimps.fa)
OutDir=analysis/AHDB_Fo_lactucae_statice_blast2/vs_ref_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
cp -sf $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query dna $RefGenome $OutDir
done
```

Blast vs additional ncbi genomes:

```bash
for RefGenome in $(ls analysis/metagenomics/reference_genomes/renamed/*.fna); do
Prefix=$(basename $RefGenome .fna)
echo $Prefix
OutDir=analysis/AHDB_Fo_lactucae_statice_blast2/vs_ref_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
cp -sf $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
Query=$(ls analysis/AHDB_Fo_lactucae_statice_blast2/Fo_lactucae_sec_mimps.fa)
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query dna $RefGenome $OutDir
done
```

Blast vs sequenced genomes:

```bash
for RefGenome in $(ls repeat_masked/F.oxysporum_*/*/*/*_contigs_unmasked.fa | grep -v 'old' | grep -v 'filtered_ncbi' | grep -v '4287' | grep -v -e 'FON139' -e 'FON77' -e 'FON89' -e 'FON81' -e 'FON129'); do
# Prefix=$(basename $RefGenome .fna)
Prefix=$(echo $RefGenome | cut -f2,3 -d '/' --output-delimiter '_')
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 20s
printf "."
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
done
printf "\n"
OutDir=analysis/AHDB_Fo_lactucae_statice_blast2/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
Query=$(ls analysis/AHDB_Fo_lactucae_statice_blast2/Fo_lactucae_sec_mimps.fa)
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query dna $RefGenome $OutDir
done

for RefGenome in $(ls repeat_masked/F.proliferatum/A8/*/*_contigs_unmasked.fa); do
Prefix=$(echo $RefGenome | cut -f2,3 -d '/' --output-delimiter '_')
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
# while [ $Jobs -gt 1 ]; do
# sleep 20s
# printf "."
# Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
# done
# printf "\n"
OutDir=analysis/AHDB_Fo_lactucae_statice_blast2/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
Query=$(ls analysis/AHDB_Fo_lactucae_statice_blast2/Fo_lactucae_sec_mimps.fa)
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query dna $RefGenome $OutDir
done

for RefGenome in $(ls ../fusarium_ex_narcissus/repeat_masked/F.oxysporum_*/*/*/*_contigs_unmasked.fa | grep -e 'FON139' -e 'FON77' -e 'FON89' -e 'FON81' -e 'FON129'); do
Prefix=$(echo $RefGenome | cut -f4,5 -d '/' --output-delimiter '_')
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 20s
printf "."
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
done
printf "\n"
OutDir=analysis/AHDB_Fo_lactucae_statice_blast2/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
Query=$(ls analysis/AHDB_Fo_lactucae_statice_blast2/Fo_lactucae_sec_mimps.fa)
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh $Query dna $RefGenome $OutDir
done
```
<!--
Blast vs work in progress genomes:

```bash
for RefGenome in $(ls assembly/spades/F.*/*/filtered_contigs/contigs_min_500bp*.fasta | grep -w -e 'Straw465' -e 'F81' -e 'FOP1-EMR' -e 'R2' -e '15-074' -e 'A1-2' -e 'HB6' -e 'D2' -e 'PG8' -e 'L5' -e 'A1-2' -e 'HB6' -e 'D2' -e 'PG8' -e 'L5' -e 'A1-2' -e 'HB6'); do
Prefix=$(echo $RefGenome | cut -f3,4 -d '/' --output-delimiter '_')
echo $Prefix
OutDir=analysis/AHDB_blast/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
rm ${Prefix}_genome.fa
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh analysis/AHDB_blast/FoM_FoN_FoC_sec_mimps.fa dna $RefGenome $OutDir
done
for RefGenome in $(ls assembly/spades_pacbio/F.*/*/filtered_contigs/contigs_min_500bp.fasta | grep -w -e '55'); do
Prefix=$(echo $RefGenome | cut -f3,4 -d '/' --output-delimiter '_')
echo $Prefix
OutDir=analysis/AHDB_blast/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
rm ${Prefix}_genome.fa
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh analysis/AHDB_blast/FoM_FoN_FoC_sec_mimps.fa dna $RefGenome $OutDir
done
``` -->

## Summarise blast hists

```bash
CsvFiles=$(ls analysis/AHDB_Fo_lactucae_statice_blast2/vs_*_genomes/*/*Fo_lactucae_sec_mimps.fa_hits.csv | grep -e 'vs_ref_genomes' -e 'vs_seq_genomes' | grep -v 'F.oxysporum_/' | grep -v -e 'GCA_000599445.1' -e '77-13-4')
Headers=$(echo $CsvFiles | sed 's&analysis/AHDB_Fo_lactucae_statice_blast2/vs_ref_genomes/&&g' | sed 's&analysis/AHDB_Fo_lactucae_statice_blast2/vs_seq_genomes/&&g' | sed -r "s&_Fo_lactucae_sec_mimps.fa_hits.csv&&g" | sed -r "s&/\S*&&g"  | sed 's&/van_dam&&g')
OutDir=analysis/AHDB_Fo_lactucae_statice_blast2/vs_ref_genomes/extracted
mkdir -p $OutDir
Genomes=$(ls analysis/AHDB_Fo_lactucae_statice_blast2/vs_*_genomes/*/*_genome.fa | grep -e 'vs_ref_genomes' -e 'vs_seq_genomes' | grep -v 'F.oxysporum_/' | grep -v 'GCA_000599445.1' | grep -e 'vs_ref_genomes' -e 'vs_seq_genomes' | grep -v 'F.oxysporum_/' | grep -v -e 'GCA_000599445.1' -e '77-13-4')
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/AHDB_project/blast_searches
$ProgDir/blast_parse_AHDB.py --blast_csv $CsvFiles --headers $Headers --genomes $Genomes --identity 0.50 --evalue 1e-30 --out_prefix $OutDir/FoL_sec_mimpsAHDB

ls $OutDir/FoL_sec_mimpsAHDB.csv
```
