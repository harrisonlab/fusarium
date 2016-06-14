#Fus2_Assembly
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
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq/fusarium/HAPI_seq_3/
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/F/
	mkdir -p $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/R/
```

Sequence data was moved into the appropriate directories

```bash
	RawDatDir=/home/groups/harrisonlab/raw_data/raw_seq
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	cp $RawDatDir/fusarium/HAPI_seq_1/FUS2_S2_L001_R1_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/F/.
	cp $RawDatDir/fusarium/HAPI_seq_1/FUS2_S2_L001_R2_001.fastq.gz $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/R/.
	cp $RawDatDir/raw_data/raw_seq/fusarium/warwick_seqs/fus2/s_6_1_sequence.txt $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/F/.
	cp $RawDatDir/raw_data/raw_seq/fusarium/warwick_seqs/fus2/s_6_2_sequence.txt $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/R/.
```

One set of sequence reads were in .txt format and unzipped. These were renamed and zipped.
```bash
	mv $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/F/s_6_1_sequence.txt $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/F/s_6_1_sequence.fastq
	mv $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/R/s_6_2_sequence.txt $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/R/s_6_2_sequence.fastq
	gzip $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/F/s_6_1_sequence.fastq
	gzip $ProjectDir/raw_dna/paired/F.oxysporum_fsp_cepae/Fus2/R/s_6_2_sequence.fast
```

# 1. Data qc

programs:
  fastqc
  fastq-mcf
  kmc

Data quality was visualised using fastqc:
```bash
	for RawData in $(ls raw_dna/paired/*/*/Fus2/*.fastq.gz); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```


# 2. Assembly

Commands used to assemble the Fus2 genome can be found in:
assembly/Fus2_assembly_commands.sh


# 3. Repeatmasking

Commands used to repeatmask the Fus2 genome can be found in:
assembly/Fus2_assembly_commands.sh


# 4. Gene prediction

Gene prediction was performed using Braker1. The commands that were used for
gene prediction can be found in:
gene_prediction/braker/braker_predictions.md

# 5. Functional annotation

# Functional annotation


## A) Interproscan

Interproscan was used to give gene models functional annotations.
Annotation was run using the commands below:

Note: This is a long-running script. As such, these commands were run using
'screen' to allow jobs to be submitted and monitored in the background.
This allows the session to be disconnected and reconnected over time.

Screen ouput detailing the progress of submission of interporscan jobs
was redirected to a temporary output file named interproscan_submission.log .

```bash
	screen -a
	cd /home/groups/harrisonlab/project_files/fusarium
	getAnnoFasta.pl gene_pred/braker/F.oxysporum_fsp_cepae/Fus2/F.oxysporum_fsp_cepae_Fus2/augustus.gff
	cat gene_pred/braker/F.oxysporum_fsp_cepae/Fus2/F.oxysporum_fsp_cepae_Fus2/augustus.gff | grep -v '#' > gene_pred/braker/F.oxysporum_fsp_cepae/Fus2/F.oxysporum_fsp_cepae_Fus2/augustus_extracted.gff
	Proteome=gene_pred/braker/F.oxysporum_fsp_cepae/Fus2/F.oxysporum_fsp_cepae_Fus2/augustus.aa
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	$ProgDir/sub_interproscan.sh $Proteome
```

Following interproscan annotation split files were combined using the following
commands:

```bash
	PredGenes=gene_pred/braker/F.oxysporum_fsp_cepae/Fus2/F.oxysporum_fsp_cepae_Fus2/augustus.aa
	InterProRaw=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2/raw
	OutDir=gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2
	Organism=F.oxysporum_fsp_cepae
	Strain=Fus2
	mkdir -p $OutDir
	printf "" > $OutDir/"$Strain"_interproscan.tsv
	printf "" > $OutDir/"$Strain"_interproscan.xml
	printf "" > $OutDir/interpro_features.gff
	printf "" > $OutDir/"$Strain"_interpro.gff3
	for File in $(ls -v $InterProRaw/*_split_*.tsv); do
		cat $File >> $OutDir/"$Strain"_interproscan.tsv
	done
	for File in $(ls -v $InterProRaw/*_split_*.xml); do
		cat $File >> $OutDir/"$Strain"_interproscan.xml
	done
	for File in $(ls -v $InterProRaw/*_split_*.gff3); do
		FastaStart=$(cat $File | grep -E "^##FASTA" -n | cut -d ':' -f1)
		cat $File | head -n "$FastaStart" | grep -v -E "^#" >> $OutDir/interpro_features.gff
	done
	cat $OutDir/interpro_features.gff $PredGenes >> $OutDir/"$Strain"_interpro.gff
	rm $OutDir/interpro_features.gff
```


## B) SwissProt

```bash
  screen -a
  qlogin
  ProjDir=/home/groups/harrisonlab/project_files/fusarium
  cd $ProjDir
  for Proteome in $(ls gene_pred/braker/F.oxysporum_fsp_cepae/Fus2/*/augustus.aa); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    OutDir=$ProjDir/gene_pred/swissprot/$Species/$Strain/
    mkdir -p $OutDir
    blastp \
    -db /home/groups/harrisonlab/uniprot/swissprot/uniprot_sprot \
    -query $ProjDir/$Proteome \
    -out $OutDir/swissprot_v2015_10_hits.tbl \
    -evalue 1e-100 \
    -outfmt 6 \
    -num_threads 16 \
    -num_alignments 10
  done
```



# 5. Genomic analysis


## 5.1  Identifying SIX gene homologs

## 5.1.a) Performing BLAST searches

BLast searches were performed against the Fus2 genome:

```bash
  Fus2_ass=repeat_masked/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_49/Fus2_combined_49_contigs_unmasked.fa
  for Assembly in $Fus2_ass; do
    ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
    Query=analysis/blast_homology/six_genes/six-appended_parsed.fa
    qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
  done
```

## 5.1.b) Converting BLAST results to gff annotations

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
	for BlastHits in $(ls analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2/*_six-appended_parsed.fa_homologs.csv); do
		Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)  
		Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
		ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
		HitsGff=analysis/blast_homology/$Organism/$Strain/"$Strain"_six-appended_parsed.fa_homologs.gff
		Column2=BLAST_hit
		NumHits=1
		$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
	done
```



### 4.3 Characterising BLAST hits

Genes overlapped by BLAST hits were identified:

```bash
	Fus2AugGff=$(ls gene_pred/braker/F.oxysporum_fsp_cepae/Fus2/*/augustus_extracted.gff)
	for Proteins in $Fus2AugGff; do
		Species=$(echo "$Proteins" | rev | cut -f4 -d '/' | rev)
		Strain=$(echo "$Proteins" | rev | cut -f3 -d '/' | rev)
		HitsGff=$(ls analysis/blast_homology/$Species/$Strain/"$Strain"_six-appended_parsed.fa_homologs.gff)
		ORFs=$(ls gene_pred/ORF_finder/$Species/$Strain/"$Strain"_ORF_mod.gff)
		OutDir=analysis/blast_homology/$Species/$Strain/"$Strain"_six-appended_parsed_analysis
		mkdir -p $OutDir
		BrakerIntersect=$OutDir/"$Strain"_six-appended_parsed_BrakerIntersect.gff
		BrakerNoIntersect=$OutDir/"$Strain"_six-appended_parsed_BrakerNoIntersect.gff
		ORFIntersect=$OutDir/"$Strain"_six-appended_parsed_ORFIntersect.gff
		Fus2Proteins=$OutDir/"$Strain"_six-appended_parsed_BrakerIntersect.txt
		Fus2ORFs=$OutDir/"$Strain"_six-appended_parsed_ORFIntersect.txt

		echo "The following number of Avr genes had blast hits:"
		cat $HitsGff | wc -l
		echo "The following number of blast hits intersected Braker gene models:"
		bedtools intersect -wa -u -a $HitsGff -b $Proteins > $BrakerIntersect
		cat $BrakerIntersect | wc -l
		echo "The following Blast hits did not intersect Braker gene models:"
		bedtools intersect -v -a $HitsGff -b $Proteins > $BrakerNoIntersect
		cat $BrakerNoIntersect | wc -l
		echo "The following number of blast hits intersected ORF gene models:"
		bedtools intersect -wa -u -a $BrakerNoIntersect -b $ORFs > $ORFIntersect
		cat $ORFIntersect | wc -l
		bedtools intersect -wao -a $HitsGff -b $Proteins | grep -w 'transcript' | cut -f9,18 | tr -d '"' | tr -d ';' | sed 's/ID=//g' > $Fus2Proteins
		bedtools intersect -wao -a $BrakerNoIntersect -b $ORFs | grep -w 'transcript' | cut -f9,18 | cut -d ';' -f1,4 | tr -d '"' | sed 's/;/\t/g' | sed 's/ID=//g' | sed 's/Name=//g' > $Fus2ORFs
	done
```

```
	The following number of Avr genes had blast hits:
	14
	The following number of blast hits intersected Braker gene models:
	12
	The following Blast hits did not intersect Braker gene models:
	2
	The following number of blast hits intersected ORF gene models:
	2
```

Of the SIX genes searched for, 14 had blast hits within the
FoC Fus2 genome. 12 of these intersected Braker gene models and the remaining
2 intersected ORF gene models.

The genes relating to each Fus2 BLAST hit were extracted.

```bash
	Fus2_pep=gene_pred/braker/F.oxysporum_fsp_cepae/Fus2/F.oxysporum_fsp_cepae_Fus2/augustus.aa
	for Proteins in $Fus2_pep; do
		Species=$(echo "$Proteins" | rev | cut -f4 -d '/' | rev)
		Strain=$(echo "$Proteins" | rev | cut -f3 -d '/' | rev)
		echo "$Species - $Strain"
		OutDir=analysis/blast_homology/$Species/$Strain/"$Strain"_six-appended_parsed_analysis
		BlastProteins=$(ls $OutDir/"$Strain"_six-appended_parsed_BrakerIntersect.txt)
		while read Line; do
			BlastName=$(echo $Line | cut -f1 -d ' ' )
			HitName=$(echo $Line | cut -f2 -d ' ')
			echo "$BlastName - $HitName"
			echo "$HitName" > tmp.txt
			OutDir=analysis/blast_homology/$Species/$Strain/"$Strain"_six-appended_parsed_analysis/fasta
			mkdir -p $OutDir
			OutFasta=$(echo ""$Species"_"$Strain"_""$BlastName""_homolog.fa" | sed 's/_BlastHit_1//g')
			ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
			$ProgDir/extract_from_fasta.py --fasta $Proteins --headers tmp.txt > $OutDir/$OutFasta
		done < $BlastProteins
		OutDir=analysis/blast_homology/$Species/$Strain/"$Strain"_six-appended_parsed_analysis
		echo ""
		ORFs=$(ls gene_pred/ORF_finder/$Species/$Strain/"$Strain".aa_cat.fa)
		BlastORFs=$(ls $OutDir/"$Strain"_six-appended_parsed_ORFIntersect.txt)
		for BlastName in $(cat $BlastORFs | cut -f1 | sort | uniq); do
			cat $BlastORFs | grep "$BlastName" | cut -f2 > tmp.txt
			echo "$BlastName - ORF hits"
			OutDir=analysis/blast_homology/$Species/$Strain/"$Strain"_six-appended_parsed_analysis/fasta
			mkdir -p $OutDir
			OutFasta=$(echo ""$Species"_"$Strain"_""$BlastName""_homolog.fa" | sed 's/_BlastHit_1//g')
			ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
			$ProgDir/extract_from_fasta.py --fasta $ORFs --headers tmp.txt > $OutDir/$OutFasta
			sed -E "s/^>/>$Species|/g" -i $OutDir/$OutFasta
		done
		echo ""
	done
```

```
	F.oxysporum_fsp_cepae - Fus2
	Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six10_(SIX10)_mRNA,_complete_cds_BlastHit_1 - g9363.t1
	Fusarium_oxysporum_f._sp._lycopersici_isolate_FOL-MM10_secreted_in_xylem_3_(SIX3)_gene,_complete_cds_BlastHit_1 - g11703.t1
	Fusarium_oxysporum_f._sp._lycopersici_isolate_IPO3_secreted_in_xylem_3_(SIX3)_gene,_complete_cds_BlastHit_1 - g11703.t1
	Fusarium_oxysporum_f._sp._lycopersici_isolate_14844_secreted_in_xylem_3_(SIX3)_gene,_complete_cds_BlastHit_1 - g11703.t1
	Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_3_(SIX3)_gene,_complete_cds_BlastHit_1 - g11703.t1
	Fusarium_oxysporum_f._sp._lycopersici_SIX3_gene_for_Secreted_in_xylem_3_protein_BlastHit_1 - g11703.t1
	Fusarium_oxysporum_f._sp._lycopersici_Six5_mRNA,_complete_cds_BlastHit_1 - g11815.t1
	Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_5_(SIX5)_gene,_partial_cds_BlastHit_1 - g11815.t1
	Fusarium_oxysporum_f._sp._lycopersici_secreted_in_xylem_Six7_(SIX7)_mRNA,_complete_cds_BlastHit_1 - g9364.t1
	Fusarium_oxysporum_f._sp._lilii_isolate_NRRL_28395_secreted_in_xylem_7-like_protein_(SIX7)_gene,_complete_cds_BlastHit_1 - g9364.t1
	Fusarium_oxysporum_f._sp._lycopersici_isolate_BFOL-51_secreted_in_xylem_7_(SIX7)_gene,_complete_cds_BlastHit_1 - g9364.t1
	Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six9_(SIX9)_mRNA,_complete_cds_BlastHit_1 - g11984.t1

	Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six12_(SIX12)_mRNA,_complete_cds_BlastHit_1 - ORF hits
	Fusarium_oxysporum_f._sp._lycopersici_strain_Fol007_Six14_(SIX14)_mRNA,_complete_cds_BlastHit_1 - ORF hits
```
