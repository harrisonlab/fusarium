# Stocks_Assembly
==========

This document details the commands used to assemble and annotate the Fus2 genome.

Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/fusarium

The following is a summary of the work presented in this Readme.

The following processes were applied to Fusarium genomes prior to analysis:
Data qc
Genome assembly
Repeatmasking
Gene prediction
Functional annotation


# 0. Building of directory structure
```bash
	RawDatDir=/home/miseq_data/minion/2017/FON-63_2017-05-22
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	mkdir -p $ProjectDir/raw_dna/minion/F.oxysporum_fsp_narcissi/FON_63
```

Sequence data was moved into the appropriate directories

```bash
	RawDatDir=/home/miseq_data/minion/2017/FON-63_2017-05-22
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	cp $RawDatDir/all_reads_albacore1.1.1.fastq.gz $ProjectDir/raw_dna/minion/F.oxysporum_fsp_narcissi/FON_63/.
```


## Assembly

### Removal of adapters

Splitting reads and trimming adapters using porechop
```bash
	for RawReads in $(ls  raw_dna/minion/F.oxysporum_fsp_narcissi/FON_63/all_reads_albacore1.1.1.fastq.gz); do
    Organism=$(echo $RawReads| rev | cut -f3 -d '/' | rev)
    Strain=$(echo $RawReads | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
  	OutDir=qc_dna/minion/$Organism/$Strain
  	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
  	qsub $ProgDir/sub_porechop.sh $RawReads $OutDir
  done
```

### Read correction using Canu

```bash
  for TrimReads in $(ls qc_dna/minion/F.*/*/*_trim.fastq.gz | grep 'FON_63'); do
    Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
    Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
    OutDir=assembly/canu-1.5/$Organism/"$Strain"
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
    qsub $ProgDir/sub_canu_correction.sh $TrimReads 60m $Strain $OutDir
  done
```

### Assembbly using SMARTdenovo

```bash
for CorrectedReads in $(ls assembly/canu-1.5/F.oxysporum_fsp_mathioli/Stocks4/Stocks4.trimmedReads.fasta.gz); do
Organism=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $CorrectedReads | rev | cut -f2 -d '/' | rev)
Prefix=$Strain
OutDir=assembly/SMARTdenovo/F.oxysporum_fsp_mathioli/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/SMARTdenovo
qsub $ProgDir/sub_SMARTdenovo.sh $CorrectedReads $Prefix $OutDir
done
```
