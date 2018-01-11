# Locus Identification

Commands used to help identify loci for metagenome sequencing of Fusarium in the field



## Downloading genome data

Fusarium genomes available on ncbi on 25/07/17 were downloaded

```bash
  ProjectDir=/home/groups/harrisonlab/project_files/fusarium
  cd $ProjectDir
  OutDir=analysis/metagenomics/reference_genomes
  mkdir -p $OutDir
  cd $OutDir
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/assembly_summary.txt
  cat assembly_summary.txt | grep 'Fusarium' | cut -f20 > fusarium_genomes_ncbi.txt
  mkdir downloads
  # for next in $(cat fusarium_genomes_ncbi.txt | grep -v 'GCA_002168265.2_ASM216826v2'); do
  for next in $(cat fusarium_genomes_ncbi.txt | awk '/GCA_002168265.2_ASM216826v2/{y=1;next}y'); do
  # wget -P downloads "$next"/*v?_genomic.fna.gz;
  wget -P downloads "$next"/*.f*.gz;
  done

  mkdir renamed
  for File in $(ls downloads/*.gz); do
    Prefix=$(echo $File | sed 's/_genomic.fna.gz//g' | cut -f2 -d '/')
    Organism=$(cat assembly_summary.txt | grep "$Prefix" | cut -f8 | sed 's/ /_/g')
    Strain=$(cat assembly_summary.txt | grep "$Prefix" | cut -f9 | sed 's/strain=//g' | sed 's/ /_/g')
    echo "$Organism - $Strain"
    mv $File renamed/"$Organism"_"$Strain".fna.gz
  done

  gunzip renamed/*.gz
```


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
    Organism=$(echo $File | rev | cut -f2 -d '/' | rev)
    Strain=$(echo $File | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    mv $File analysis/metagenomics/reference_genomes/renamed/"$Organism"_"$Strain".fna
  done
```


## TEF

Extracting TEF loci
