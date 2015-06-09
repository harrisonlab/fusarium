#HAPI
==========

Scripts used for the analysis of Fusarium genomes
Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/fusarium

The following is a summary of the work presented in this Readme.

The following processes were applied to Alternaria genomes prior to analysis:
Data qc
Genome assembly
Repeatmasking
Gene prediction
Functional annotation

Analyses performed on these genomes involved BLAST searching for:


ls contigs were identified using:
Alignment of raw reads to assembled genomes
Assembly of remaining reads


#Building of directory structure
```shell
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_3/ 
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	mkdir -p $ProjectDir/raw_dna/paired/F.proliferatum/A8/F
	mkdir -p $ProjectDir/raw_dna/paired/F.proliferatum/A8/R
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_narcissi/N139/F
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_narcissi/N139/R
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG3/F
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG3/R
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG18/F
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG18/R
```
Sequence data was moved into the appropriate directories

```shell
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_3
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	cp $RawDatDir/Fproliferatum_S1_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.proliferatum/A8/F/.
	cp $RawDatDir/Fproliferatum_S1_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.proliferatum/A8/R/.
	cp $RawDatDir/FoxysporumN139_S2_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_narcissi/N139/F/.
	cp $RawDatDir/FoxysporumN139_S2_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_narcissi/N139/R/.
	cp $RawDatDir/FoxysporumPG3_S3_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG3/F/.
	cp $RawDatDir/FoxysporumPG3_S3_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG3/R/.
	cp $RawDatDir/FoxysporumPG18_S4_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG18/F/.
	cp $RawDatDir/FoxysporumPG18_S4_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_pisi/PG18/R/.
```

This process was repeated for RNAseq data:

```shell
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/fusarium/rna_seq
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDB/F
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDB/R
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_CzapekDox/F
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_CzapekDox/R
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_GlucosePeptone/F
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_GlucosePeptone/R
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDA/F
	mkdir -p $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDA/R

	cp $RawDatDir/4_S1_L001_R1_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDB/F/.
	cp $RawDatDir/4_S1_L001_R2_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDB/R/.
	cp $RawDatDir/6_S2_L001_R1_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_CzapekDox/F/.
	cp $RawDatDir/6_S2_L001_R2_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_CzapekDox/R/.
	cp $RawDatDir/7_S3_L001_R1_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_GlucosePeptone/F/.
	cp $RawDatDir/7_S3_L001_R2_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_GlucosePeptone/R/.
	cp $RawDatDir/9_S4_L001_R1_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDA/F/.
	cp $RawDatDir/9_S4_L001_R2_001.fastq.gz $ProjectDir/raw_rna/paired/F.oxysporum_fsp_cepae/Fus2_PDA/R/.
```



#Data qc

programs:
  fastqc
  fastq-mcf
  kmc

Data quality was visualised using fastqc:
```shell
	for RawData in $(ls raw_dna/paired/*/*/*/*.fastq.gz | grep -v 'cepae'); do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData; 
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

Trimming was performed on data to trim adapters from 
sequences and remove poor quality data. This was done with fastq-mcf

```shell
	for StrainPath in $(ls -d raw_dna/paired/*/* | grep -v 'cepae'); do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
		IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
		ReadsF=$(ls $StrainPath/F/*.fastq*)
		ReadsR=$(ls $StrainPath/R/*.fastq*)
		echo $ReadsF
		echo $ReadsR
		qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	done
```



Data quality was visualised once again following trimming:
```shell
	for RawData in qc_dna/paired/*/*/*/*.fastq*; do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData; 
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

kmer counting was performed using kmc
This allowed estimation of sequencing depth and total genome size

```shell
	for TrimPath in qc_dna/paired/*/*; do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		TrimF=$(ls $TrimPath/F/*.fastq*)
		TrimR=$(ls $TrimPath/R/*.fastq*)
		echo $TrimF
		echo $TrimR
		qsub $ProgDir/kmc_kmer_counting.sh $TrimF $TrimR
	done
```

#Assembly

Assembly was performed using Velvet

A range of hash lengths were used and the best assembly selected for subsequent analysis

```shell
	for TrimPath in $(ls -d qc_dna/paired/*/* | grep -v 'cepae'); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/velvet
		Strain=$(printf $TrimPath | rev | cut -f1 -d '/' | rev)

		MinHash=41
		MaxHash=81
		HashStep=2
		TrimF=$(ls $TrimPath/F/*.fastq*)
		TrimR=$(ls $TrimPath/R/*.fastq*)
		GenomeSz=65

		echo $Strain
		if [ "$Strain" == 'PG18' ]; then 
			ExpCov=24
			MinCov=8
			InsLgth=600
			echo "$Strain set"
		elif [ "$Strain" == 'PG3' ]; then
			ExpCov=27
			MinCov=9
			InsLgth=600
			echo "$Strain set"
		elif [ "$Strain" == 'N139' ]; then
			ExpCov=26
			MinCov=9
			InsLgth=600
			echo "$Strain set"
		elif [ "$Strain" == 'A8' ]; then
			ExpCov=20
			MinCov=7
			InsLgth=600
			echo "$Strain set"
		fi

		qsub $ProgDir/submit_velvet_range.sh $MinHash $MaxHash $HashStep \
		$TrimF $TrimR $GenomeSz $ExpCov $MinCov $InsLgth
	done
	
```


Assemblies were summarised to allow the best assembly to be determined by eye.
Although flashed reads and additional datasets were not used in this set of analyses, there were
some results from previous assemblies in the destination folders that needed to be excluded.
In the for loop that cycled through the stats file from each assembly correctly the following 
commands was added: | grep -v 'flash' | grep -v 'combined' . This prevented assemblies
using flashed reads or additional datasets from being included in this summary. These 
two grep expressions can be excluded if copying and pasting these commands for a different project.

```shell
	for StrainPath in $(ls -d assembly/velvet/F*/* ); do
		printf "N50\tMax_contig_size\tNumber of bases in contigs\tNumber of contigs\tNumber of contigs >=1kb\tNumber of contigs in N50\tNumber of bases in contigs >=1kb\tGC Content of contigs\n" > $StrainPath/assembly_stats.csv
		for StatsFile in $(ls $StrainPath/*/stats.txt | grep -v 'flash' | grep -v 'combined'); do
		cat $StatsFile | rev | cut -f1 -d ' ' | rev | paste -d '\t' -s >> $StrainPath/assembly_stats.csv
		done
	done
	tail -n+1 assembly/velvet/F*/*/assembly_stats.csv > assembly/velvet/Fusarium_assembly_stats.csv
```



#Repeatmasking

Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking
	
```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
	BestAssPG18=assembly/velvet/F.*/PG18/F.*_PG18_53/sorted_contigs.fa
	BestAssPG3=assembly/velvet/F.*/PG3/F.*_PG3_79/sorted_contigs.fa
	BestAssN139=assembly/velvet/F.*/N139/F.*_N139_59/sorted_contigs.fa
	BestAssA8=assembly/velvet/F.*/A8/F.*_A8_51/sorted_contigs.fa
		
	qsub $ProgDir/rep_modeling.sh $BestAssPG18
	qsub $ProgDir/rep_modeling.sh $BestAssPG3
	qsub $ProgDir/rep_modeling.sh $BestAssN139
	qsub $ProgDir/rep_modeling.sh $BestAssA8

	qsub $ProgDir/transposonPSI.sh $BestAssPG18
	qsub $ProgDir/transposonPSI.sh $BestAssPG3
	qsub $ProgDir/transposonPSI.sh $BestAssN139
	qsub $ProgDir/transposonPSI.sh $BestAssA8

```

	
#Gene Prediction


Gene prediction followed three steps:
	Pre-gene prediction
		- Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified. 
	Gene model training
		- Gene models were trained for isolates 1166 and 650 using assembled RNAseq data
	Gene prediction
		- Gene models were used to predict genes in A. alternata genomes. This used RNAseq data as hints for gene models.

#Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.
```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
	cd /home/groups/harrisonlab/project_files/fusarium
	for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa | grep -v 'cepae'); do 
		echo $Genome; 
		qsub $ProgDir/sub_cegma.sh $Genome dna;
	done
```

Outputs were summarised using the commands:
```shell
	for File in $(ls gene_pred/cegma/F*/*/*_dna_cegma.completeness_report | grep -v 'cepae'); do 
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev); 
		Species=$(echo $File | rev | cut -f3 -d '/' | rev); 
		printf "$Species\t$Strain\n"; 
		cat $File | head -n18 | tail -n+4;printf "\n"; 
	done >> gene_pred/cegma/cegma_results_dna_summary.txt
	
	less gene_pred/cegma/cegma_results_dna_summary.txt
```

<!--
#Gene model training

Data quality was visualised using fastqc:
```shell
	for RawData in raw_rna/paired/*/*/*/*.fastq.gz; do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData; 
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

Trimming was performed on data to trim adapters from 
sequences and remove poor quality data. This was done with fastq-mcf

```shell
	for StrainPath in raw_rna/paired/*/*; do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
		IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/illumina_full_adapters.fa
		ReadsF=$(ls $StrainPath/F/*.fastq.gz)
		ReadsR=$(ls $StrainPath/R/*.fastq.gz)
		echo $ReadsF
		echo $ReadsR
		qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters RNA
	done
```

Data quality was visualised once again following trimming:
```shell
	for TrimData in qc_rna/paired/*/*/*/*.fastq.gz; do 
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData; 
		qsub $ProgDir/run_fastqc.sh $TrimData
	done
```	

RNAseq data was assembled into transcriptomes using Trinity
```shell
	for StrainPath in qc_rna/paired/*/*; do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/transcriptome_assembly
		ReadsF=$(ls $StrainPath/F/*.fastq.gz)
		ReadsR=$(ls $StrainPath/R/*.fastq.gz)	
		echo $ReadsF
		echo $ReadsR
		qsub $ProgDir/transcriptome_assembly_trinity.sh $ReadsF $ReadsR
	done
```	
Gene training was performed using RNAseq data. The cluster can not run this script using qlogin. As such it was run on the head node (-naughty) using screen.
Training for 650 and 1166 was performed in two instances of screen and occassionally viewed to check progress over time.
(screen is detached after opening using ctrl+a then ctrl+d. - if just ctrl+d is pressed the instance of screen is deleted. - be careful)
```shell
	screen -a
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
	Assembly650=assembly/trinity/A.alternata_ssp._gaisen/650/650_rna_contigs/Trinity.fasta 
	Genome650=repeat_masked/A.alternata_ssp._gaisen/650/A.alternata_ssp._gaisen_650_67_repmask/650_contigs_unmasked.fa
	$ProgDir/training_by_transcriptome.sh $Assembly650 $Genome650

	screen -a
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
	Assembly1166=assembly/trinity/A.alternata_ssp._tenuissima/1166/1166_rna_contigs/Trinity.fasta 
	Genome1166=repeat_masked/A.alternata_ssp._tenuissima/1166/A.alternata_ssp._tenuissima_1166_43_repmask/1166_contigs_unmasked.fa
	$ProgDir/training_by_transcriptome.sh $Assembly1166 $Genome1166
```

Quality of Trinity assemblies were assessed using Cegma to assess gene-space within the transcriptome
```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/cegma
	for Transcriptome in $(ls assembly/trinity/A.*/*/*_rna_contigs/Trinity.fasta); do  
		echo $Transcriptome;  
		qsub $ProgDir/sub_cegma.sh $Transcriptome rna; 
	done
```
Outputs were summarised using the commands:
```shell
	for File in $(ls gene_pred/cegma/A.alternata_ssp._*/*/*_rna_cegma.completeness_report); do 
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev); 
		Species=$(echo $File | rev | cut -f3 -d '/' | rev); 
		printf "$Species\t$Strain\n"; 
		cat $File | head -n18 | tail -n+4;printf "\n"; 
	done > gene_pred/cegma/cegma_results_rna_summary.txt
	
	less gene_pred/cegma/cegma_results_rna_summary.txt
```
 -->




	

#Gene prediction

Gene prediction was performed for Fusarium genomes. 
RNAseq reads were used as Hints for the location of CDS. 

A concatenated dataset of RNAseq reads from F. oxysporum fsp. cepae isolate Fus2
were used as hints for these predictions.
A gene model trained for F.oxysporum fsp. cepae was used to describe the structure of a gene.

```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/augustus
	mkdir -p qc_rna/concatenated
	RnaFiles=$(ls qc_rna/paired/._*/*/*/*.fastq.gz | paste -s -d ' ')
	RnaF=qc_rna/paired/Fus2/aligned_appended/appended_paired.1.fastq
	RnaR=qc_rna/paired/Fus2/aligned_appended/appended_paired.2.fastq
	mkdir -p qc_rna/concatenated/F.oxysporum/Fus2
	ConcatRna=qc_rna/concatenated/F.oxysporum/Fus2/Fus2_RNA_timecourse_appended.fa.gz
	cat $RnaF $RnaR | gzip -fc > $ConcatRna
	GeneModel=F.oxysporum_fsp_cepae_Fus2
	for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa | grep -v 'cepae'); do 
		qsub $ProgDir/augustus_pipe.sh $Genome $ConcatRna $GeneModel
	done
```



#Functional annotation

Interproscan was used to give gene models functional annotations. 
Annotation was run using the commands below:

Note: This is a long-running script. As such, these commands were run using 
'screen' to allow jobs to be submitted and monitored in the background. 
This allows the session to be disconnected and reconnected over time.

Screen ouput detailing the progress of submission of interporscan jobs 
was redirected to a temporary output file named interproscan_submission.log . 

```shell
	screen -a
	cd /home/groups/harrisonlab/project_files/fusarium
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan/
	for Genes in $(ls gene_pred/augustus/F.*/*/*_augustus_preds.aa | grep -v 'cepae'); do
		echo $Genes
		$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 |  tee -a interproscan_submisison.log
```

<!-- 







#Genomic analysis


The first analysis was based upon BLAST searches for genes known to be involved in toxin production

#BLAST Searches
```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast/
	Query=analysis/blast_homology/CDC_genes/A.alternata_CDC_genes.fa 
	BestAss675=assembly/velvet/A.alternata_ssp._arborescens/675/A.alternata_ssp._arborescens_675_69/sorted_contigs.fa
	BestAss970013=assembly/velvet/A.alternata_ssp._arborescens/97.0013/A.alternata_ssp._arborescens_97.0013_59/sorted_contigs.fa
	BestAss970016=assembly/velvet/A.alternata_ssp._arborescens/97.0016/A.alternata_ssp._arborescens_97.0016_77/sorted_contigs.fa
	BestAss650=assembly/velvet/A.alternata_ssp._gaisen/650/A.alternata_ssp._gaisen_650_67/sorted_contigs.fa
	BestAss1082=assembly/velvet/A.alternata_ssp._tenuissima/1082/A.alternata_ssp._tenuissima_1082_49/sorted_contigs.fa
	BestAss1164=assembly/velvet/A.alternata_ssp._tenuissima/1164/A.alternata_ssp._tenuissima_1164_67/sorted_contigs.fa
	BestAss1166=assembly/velvet/A.alternata_ssp._tenuissima/1166/A.alternata_ssp._tenuissima_1166_43/sorted_contigs.fa
	BestAss1177=assembly/velvet/A.alternata_ssp._tenuissima/1177/A.alternata_ssp._tenuissima_1177_63/sorted_contigs.fa
	BestAss24350=assembly/velvet/A.alternata_ssp._tenuissima/24350/A.alternata_ssp._tenuissima_24350_63/sorted_contigs.fa
	BestAss635=assembly/velvet/A.alternata_ssp._tenuissima/635/A.alternata_ssp._tenuissima_635_59/sorted_contigs.fa
	BestAss648=assembly/velvet/A.alternata_ssp._tenuissima/648/A.alternata_ssp._tenuissima_648_45/sorted_contigs.fa
	BestAss743=assembly/velvet/A.alternata_ssp._tenuissima/743/A.alternata_ssp._tenuissima_743_69/sorted_contigs.fa
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss675
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss970013
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss970016
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss650
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss1082
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss1164
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss1166
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss1177
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss24350
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss635
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss648
	qsub $ProgDir/blast_pipe.sh $Query dna $BestAss743
```
BLAST search results were summarised into a presence/absence table

The presence /absence table determines presence if a hit is present and the alignment represents >50% of the query sequence.
This thresholding means that some hits have not been summarised including AMT11, AMT15 and ALT1.
```shell
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast/
	InFiles=$(ls analysis/blast_homology/A.alternata_ssp._*/*/*_A.alternata_CDC_genes.fa_homologs.csv | paste -s -d ' ')
	echo $InFiles
	$ProgDir/blast_differentials.pl $InFiles
	mv *.csv analysis/blast_homology/CDC_genes/.
``` 






#CDC Assembly

Raw reads were aligned against assembled genomes to identify contigs that were unique to a isolate or clade
```shell
	for Pathz in $(ls -d qc_dna/paired/A.alternata_ssp._*/*); do  
		Strain=$(echo $Pathz | cut -d '/' -f4)
		echo "using reads for $Strain"
		ProgPath=/home/armita/git_repos/emr_repos/tools/pathogen/lineage_specific_regions
		F_IN=$(ls $Pathz/F/*.fastq.gz)
		R_IN=$(ls $Pathz/R/*.fastq.gz)
		for Assemblyz in $(ls repeat_masked/A.alternata_ssp._*/*/*/*_contigs_unmasked.fa); do
			basename $Assemblyz
			qsub "$ProgPath"/bowtie2_alignment_pipe.sh $F_IN $R_IN $Assemblyz
		done
	done
```
A summary file was made from the alignment logs.
The percentage of reads aligning to each set of assembled contigs was determined.

```shell
	SummaryFile=analysis/ls_contigs/alignment_summaries.txt
	printf "" > "$SummaryFile"
	for OUTPUT in $(ls bowtie2_alignment_pipe.sh.e*); do 
		ID=$(echo $OUTPUT | rev | cut -d 'e' -f1 | rev | less); 
		cat bowtie2_alignment_pipe.sh.o"$ID" | grep -E "Trimmed .* reads .*/F/|Output files: " | sed -e 's/.*\/F\///g' | cut -f1 -d ')' | cut -f2 -d '"' >> "$SummaryFile"; 
		cat $OUTPUT >> "$SummaryFile"; 
		printf "\n" >> "$SummaryFile"; 
	done
	cat $SummaryFile | grep ' overall alignment rate' | cut -f1 -d ' ' | less
```
These are the reads percentages reported by bowtie but do not actually reflect the percentage unaligned reads where neither pair matched.

These were identified using SAM FLAGS to extract unaligned pairs.
```shell
	for Pathz in $(ls assembly/ls_contigs/A.alternata_ssp._*/*/vs_*/*_sorted.bam); do  
		OutFileF=$(echo $Pathz | sed 's/.bam/_unaligned_F.txt/g')
		OutFileR=$(echo $Pathz | sed 's/.bam/_unaligned_R.txt/g')
		OutFileSum=$(echo $Pathz | sed 's/.bam/_sum.txt/g')
		samtools view -f 77 "$Pathz" | cut -f1 > $OutFileF
		samtools view -f 141 "$Pathz" | cut -f1 > $OutFileR
		NoReads=$(samtools view -f 1 "$Pathz" | wc -l)
		NoReadsF=$(cat "$OutFileF" | wc -l)
		NoReadsR=$(cat "$OutFileR" | wc -l)
		printf "File\tNo. paired reads\tF reads unaligned\tR reads unaligned\n" > "$OutFileSum"
		printf "$Pathz\t$NoReads\t$NoReadsF\t$NoReadsR\n" >> "$OutFileSum"
	done

	for File in $(ls assembly/ls_contigs/A.alternata_ssp._*/*/vs_*/*_sum.txt); do 
		cat $File | tail -n+2; 
	done > analysis/ls_contigs/assembly_summaries2.txt 
```

The number of reads aligning per bp of assembly was determined. Typical alignment values were 0.20 reads per bp. Contigs were detemined as unique to that alignment if they contained an average of 0 reads per bp.
The number of bp unique to each assembly were identified.

```shell
	mkdir -p analysis/ls_contigs
	Outfile=analysis/ls_contigs/ls_contig_size.csv
	printf "Reads" > $Outfile
	for Genome in $(ls -d assembly/ls_contigs/A.alternata_ssp._arborescens/675/vs_A.alternata_ssp._*); do
		NameG=$(basename "$Genome" | sed 's/vs_//g' | sed 's/_repmask_contigs//g')
		printf "\t$NameG" >> $Outfile
	done
	printf "\n" >> $Outfile
	for Reads in $(ls -d assembly/ls_contigs/A.alternata_ssp._*/*); do
		NameR=$(basename "$Reads")
		printf "$NameR" >> $Outfile
		for Genome in $(ls -d $Reads/*); do
			printf "\t" >> $Outfile
			cat $Genome/*_sorted_indexstats_coverage.csv | head -n-1 | grep -E -w "0.00$" | cut -f2 | awk '{s+=$1} END {printf s}' >> $Outfile
		done
		printf "\n" >> $Outfile
	done
```
This did not give clear results.

Contigs that had no reads align to them were identified.
 -->
