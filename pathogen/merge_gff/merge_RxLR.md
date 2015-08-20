

These commands were performed to merge gff evidence of Putative RxLR effectors
from fusarium spp. isolates.

# Merge Fus2 ORF and Augustus RxLR motif effectors, WY hmm model effectors and RxLR hmm model effectors

## Convert ORF predictions into correct gff3 format

```bash
  ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
  ORF_Gff=gene_pred/ORF_finder/F.oxysporum_fsp_cepae/Fus2/Fus2_ORF.gff
  ORF_Gff_mod=gene_pred/ORF_finder/F.oxysporum_fsp_cepae/Fus2/Fus2_ORF_mod.gff
  $ProgDir/gff_corrector.pl $ORF_Gff > $ORF_Gff_mod
```

## Extract names of effectors from each source

### Augustus genes identified as putative effectors

Extracting RxLR Regex genes

```bash
  Aug_Regex_RxLR_FA=analysis/RxLR_effectors/RxLR_EER_regex_finder/F.oxysporum_fsp_cepae/Fus2/Fus2_Aug_RxLR_EER_regex.fa
  Aug_Regex_RxLR=analysis/RxLR_effectors/RxLR_EER_regex_finder/F.oxysporum_fsp_cepae/Fus2/Fus2_Aug_RxLR_regex_names.txt
  Aug_Regex_RxLR_EER=analysis/RxLR_effectors/RxLR_EER_regex_finder/F.oxysporum_fsp_cepae/Fus2/Fus2_Aug_RxLR_EER_regex_names.txt
  cat $Aug_Regex_RxLR_FA | grep '>' | sed 's/>//g' | cut -f1 | sed 's/ //g' > $Aug_Regex_RxLR
  cat $Aug_Regex_RxLR_FA | grep '>' | grep 'EER_motif' | sed 's/>//g' | cut -f1 | sed 's/ //g' > $Aug_Regex_RxLR_EER
```

Extracting WY hmm domain containing genes

```bash
  Aug_hmm_WY_FA=analysis/RxLR_effectors/hmmer_WY/F.oxysporum_fsp_cepae/Fus2/Fus2_Aug_WY_hmmer.fa
  Aug_hmm_WY=analysis/RxLR_effectors/hmmer_WY/F.oxysporum_fsp_cepae/Fus2/Fus2_Aug_WY_hmmer_names.txt
  cat $Aug_hmm_WY_FA | grep '>' | sed 's/>//g' | cut -f1 | sed 's/ //g' > $Aug_hmm_WY
```

Extracting RxLR hmm domain containing genes

```
  Aug_hmm_RxLR_FA=analysis/RxLR_effectors/hmmer_RxLR/F.oxysporum_fsp_cepae/Fus2/Fus2__Aug_RxLR_hmmer.fa
  Aug_hmm_RxLR=analysis/RxLR_effectors/hmmer_RxLR/F.oxysporum_fsp_cepae/Fus2/Fus2_Aug_RxLR_hmmer_names.txt
  cat $Aug_hmm_RxLR_FA | grep '>' | sed 's/>//g' | cut -f1 | sed 's/ //g' > $Aug_hmm_RxLR
```

### ORF fragments identified as putative effectors

Extracting RxLR containing ORFs

```bash
  ORF_Regex_RxLR_FA=analysis/RxLR_effectors/RxLR_EER_regex_finder/F.oxysporum_fsp_cepae/Fus2/Fus2_ORF_RxLR_EER_regex.fa
  ORF_Regex_RxLR=analysis/RxLR_effectors/RxLR_EER_regex_finder/F.oxysporum_fsp_cepae/Fus2/Fus2_ORF_RxLR_regex_names.txt
  ORF_Regex_RxLR_EER=analysis/RxLR_effectors/RxLR_EER_regex_finder/F.oxysporum_fsp_cepae/Fus2/Fus2_ORF_RxLR_EER_regex_names.txt
  cat $ORF_Regex_RxLR_FA | grep '>' | sed 's/>//g' | cut -f1 | sed 's/ //g' > $ORF_Regex_RxLR
  cat $ORF_Regex_RxLR_FA | grep '>' | grep 'EER_motif' | sed 's/>//g' | cut -f1 | sed 's/ //g' > $ORF_Regex_RxLR_EER
```

Extracting WY hmm domain containing ORFs

```bash
  ORF_hmm_WY_FA=analysis/RxLR_effectors/hmmer_WY/F.oxysporum_fsp_cepae/Fus2/Fus2_ORF_WY_hmmer.fa
  ORF_hmm_WY=analysis/RxLR_effectors/hmmer_WY/F.oxysporum_fsp_cepae/Fus2/Fus2_ORF_WY_hmmer_names.txt
  cat $ORF_hmm_WY_FA | grep '>' | sed 's/>//g' | cut -f1 | sed 's/ //g' > $ORF_hmm_WY
```

Extracting RxLR hmm domain containing ORFs

```
  ORF_hmm_RxLR_FA=analysis/RxLR_effectors/hmmer_RxLR/F.oxysporum_fsp_cepae/Fus2/Fus2_ORF_RxLR_hmmer.fa
  ORF_hmm_RxLR=analysis/RxLR_effectors/hmmer_RxLR/F.oxysporum_fsp_cepae/Fus2/Fus2_ORF_RxLR_hmmer_names.txt
  cat $ORF_hmm_RxLR_FA | grep '>' | sed 's/>//g' | cut -f1 | sed 's/ //g' > $ORF_hmm_RxLR
```



## Set variables
```bash
  ProgDir=~/git_repos/emr_repos/scripts/phytophthora/pathogen/merge_gff

	# Infiles
	Aug_Gff=gene_pred/augustus/F.oxysporum_fsp_cepae/Fus2/Fus2_augustus_preds.gtf
  ORF_Gff_mod=gene_pred/ORF_finder/F.oxysporum_fsp_cepae/Fus2/Fus2_ORF_mod.gff

  Aug_Regex_RxLR=analysis/RxLR_effectors/RxLR_EER_regex_finder/F.oxysporum_fsp_cepae/Fus2/Fus2_Aug_RxLR_regex_names.txt
  Aug_Regex_RxLR_EER=analysis/RxLR_effectors/RxLR_EER_regex_finder/F.oxysporum_fsp_cepae/Fus2/Fus2_Aug_RxLR_EER_regex_names.txt
  Aug_hmm_WY=analysis/RxLR_effectors/hmmer_WY/F.oxysporum_fsp_cepae/Fus2/Fus2_Aug_WY_hmmer_names.txt
  Aug_hmm_RxLR=analysis/RxLR_effectors/hmmer_RxLR/F.oxysporum_fsp_cepae/Fus2/Fus2_Aug_RxLR_hmmer_names.txt
  Aug_Mimp_1500=analysis/mimps/F.oxysporum_fsp_cepae/Fus2/Fus2_mimps_intersected_Aug_genes_names.txt

  ORF_Regex_RxLR=analysis/RxLR_effectors/RxLR_EER_regex_finder/F.oxysporum_fsp_cepae/Fus2/Fus2_ORF_RxLR_regex_names.txt
  ORF_Regex_RxLR_EER=analysis/RxLR_effectors/RxLR_EER_regex_finder/F.oxysporum_fsp_cepae/Fus2/Fus2_ORF_RxLR_EER_regex_names.txt
  ORF_hmm_WY=analysis/RxLR_effectors/hmmer_WY/F.oxysporum_fsp_cepae/Fus2/Fus2_ORF_WY_hmmer_names.txt
  ORF_hmm_RxLR=analysis/RxLR_effectors/hmmer_RxLR/F.oxysporum_fsp_cepae/Fus2/Fus2_ORF_RxLR_hmmer_names.txt
  ORF_Mimp_1500=analysis/mimps/F.oxysporum_fsp_cepae/Fus2/Fus2_mimps_intersected_ORF_genes_names.txt

  # Outfiles
  OutDir=analysis/database/F.oxysporum_fsp_cepae/Fus2
	AugDB=$OutDir/Fus2_Aug.db
  OrfDB=$OutDir/Fus2_ORF.db

	Aug_RxLR_DB=$OutDir/Aug_RxLR_rxlr.db
  Aug_RxLR_EER_DB=$OutDir/Aug_RxLR_EER_rxlr.db
  Aug_hmm_WY_DB=$OutDir/Aug_RxLR_EER_WY_rxlr.db
  Aug_hmm_RxLR_DB=$OutDir/Aug_RxLR_EER_WY_RxLR_rxlr.db
  Aug_mimp_DB=$OutDir/Aug_RxLR_EER_WY_RxLR_mimps_rxlr.db
  Aug_named_parents_DB=$OutDir/Fus2_Aug_annotated_parents.db

  ORF_RxLR_DB=$OutDir/ORF_RxLR_rxlr.db
  ORF_RxLR_DB_EER=$OutDir/ORF_RxLR_EER_rxlr.db
  ORF_hmm_WY_DB=$OutDir/ORF_RxLR_EER_WY_rxlr.db
  ORF_hmm_RxLR_DB=$OutDir/ORF_RxLR_EER_WY_RxLR_rxlr.db
  ORF_mimp_DB=$OutDir/ORF_RxLR_EER_WY_RxLR_mimps_rxlr.db
  ORF_named_parents_DB=$OutDir/Fus2_ORF_annotated_parents.db

  ORF_Out_Gff=$OutDir/Fus2_ORF_effectors.gff

  Orf_combined_DB=$OutDir/Fus2_ORF_combined.db
  Orf_merged_DB=$OutDir/Fus2_ORF_merged.db

  Combined_DB=$OutDir/Fus2_Aug_ORF_combined.db
  Merged_DB=$OutDir/Fus2_Aug_ORF_merged.db
  Merged_Gff=$OutDir/Fus2_Aug_ORF_merged.gff
  Effectors_Gff=$OutDir/Fus2_effectors.gff

  # WyDB=414_WY.db
	# WyID=WY_id.txt
	# WyDB_mod=414_WY_note.db
	#
	# RxlrID=rxlr_id.txt
	# RxlrDB_mod=414_rxlr_note.db
	# Rxlr_Wy_DB=414_rxlr_WY.db
  #
	# OrfMerged=414_rxlr_WY_merged.db
	# MergedDB=414_Aug_ORF_merged.db
	# FinalDB=414_Aug_ORF.db
	# FinalGff=414_Aug_ORF.gff
```

## Make a db of aug genes and effector ORFs
```bash
	$ProgDir/make_gff_database.py --inp $Aug_Gff --db $AugDB
  $ProgDir/make_gff_database.py --inp $ORF_Gff_mod --db $OrfDB
	# $ProgDir/make_gff_database.py --inp $ORF_WYs --db $WyDB
	# $ProgDir/make_gff_database.py --inp $ORF_RxLRs --db $RxlrDB
```



## Get the IDs of all the genes in the RxLR and WY ORF databases and add notes

Add notes to Augustus genes with effector evidence

```bash
	# $ProgDir/get_db_id.py --db $WyDB --type gene --out $WyID

  $ProgDir/note2db.py --in_db $AugDB --out_db $Aug_RxLR_DB --id_file $Aug_Regex_RxLR --str Aug_RxLR_motif --attribute ID
  $ProgDir/note2db.py --in_db $Aug_RxLR_DB --out_db $Aug_RxLR_EER_DB --id_file $Aug_Regex_RxLR_EER --str Aug_RxLR_EER_motif --attribute ID
  $ProgDir/note2db.py --in_db $Aug_RxLR_EER_DB --out_db $Aug_hmm_WY_DB --id_file $Aug_hmm_WY --str Aug_WY_hmm --attribute ID
  $ProgDir/note2db.py --in_db $Aug_hmm_WY_DB --out_db $Aug_hmm_RxLR_DB --id_file $Aug_hmm_RxLR --str Aug_RxLR_hmm --attribute ID
  $ProgDir/note2db.py --in_db $Aug_hmm_RxLR_DB --out_db $Aug_mimp_DB --id_file $Aug_Mimp_1500 --str Aug_mimp_intersect --attribute ID
  $ProgDir/notes2parents.py --in_db $Aug_mimp_DB --out_db $Aug_named_parents_DB
```

Add notes to ORF fragments with effector evidence

```bash

  $ProgDir/note2db.py --in_db $OrfDB --out_db $ORF_RxLR_DB --id_file $ORF_Regex_RxLR --str ORF_RxLR_motif --attribute Name
  $ProgDir/note2db.py --in_db $ORF_RxLR_DB --out_db $ORF_RxLR_DB_EER --id_file $ORF_Regex_RxLR_EER --str ORF_RxLR_EER_motif --attribute Name
  $ProgDir/note2db.py --in_db $ORF_RxLR_DB_EER --out_db $ORF_hmm_WY_DB --id_file $ORF_hmm_WY --str ORF_WY_hmm --attribute Name
  $ProgDir/note2db.py --in_db $ORF_hmm_WY_DB --out_db $ORF_hmm_RxLR_DB --id_file $ORF_hmm_RxLR --str ORF_RxLR_hmm --attribute Name
  $ProgDir/note2db.py --in_db $ORF_hmm_RxLR_DB --out_db $ORF_mimp_DB --id_file $ORF_Mimp_1500 --str ORF_mimp_intersect --attribute ID
  $ProgDir/notes2parents.py --in_db $ORF_mimp_DB --out_db $ORF_named_parents_DB
	# $ProgDir/note2db.py --in_db $WyDB --out_db $WyDB_mod --id_file $WyID --str ORF_WY_hmm --attribute ID
	# $ProgDir/get_db_id.py --db $RxlrDB --type gene --out $RxlrID
	# $ProgDir/note2db.py --in_db $RxlrDB --out_db $RxlrDB_mod --id_file $RxlrID --str ORF_RxLR_atg --attribute ID
```

$ProgDir/extract_by_note.py --db $Aug_named_parents_DB --str Aug_RxLR_motif Aug_RxLR_EER_motif Aug_WY_hmm Aug_RxLR_hmm Aug_mimp_intersect ORF_RxLR_motif ORF_RxLR_EER_motif ORF_WY_hmm ORF_RxLR_hmm ORF_mimp_intersect --out tmp.gff --type gene transcript

```bash
  $ProgDir/extract_by_note.py --db $ORF_named_parents_DB --str ORF_RxLR_motif ORF_RxLR_EER_motif ORF_WY_hmm ORF_RxLR_hmm ORF_mimp_intersect --out $ORF_Out_Gff --type gene transcript
  $ProgDir/make_gff_database.py --inp $ORF_Out_Gff --db $Orf_combined_DB
```

## Merge all features in the ORF database
```bash
	$ProgDir/merge_db_features.py --inp $Orf_combined_DB --id ORF_finder --source ORF_finder --out $Orf_merged_DB
```

$ProgDir/extract_by_note.py --db $ORF_hmm_RxLR_DB --str ORF_RxLR_motif ORF_RxLR_EER_motif ORF_WY_hmm ORF_RxLR_hmm ORF_mimp_intersect --out $ORF_Out_Gff --type gene transcript

## Merge the effector ORF database and the augustus database together
```bash
	$ProgDir/merge_db.py --inp $Aug_named_parents_DB $Orf_merged_DB --db $Combined_DB
```

## Identify all effector ORFs contained within Augustus genes
```bash
	$ProgDir/contained_features.py --inp $Combined_DB --out_db $Merged_DB --A AUGUSTUS --B ORF_finder --out_gff $Merged_Gff
```

The total number of Augustus genes are: 11055
The total number of atg genes are:      732
Of these, this many were merged:        662
Into this many features:        326
And this many remain unmerged:  11125
The final dataset contains the following number of features:    11451

The total number of Augustus genes are:	11055
The total number of atg genes are:	732
Of these, this many were merged:	662
Into this many features:	326
And this many remain unmerged:	11125
The final dataset contains the following number of features:	11451

```bash
  $ProgDir/extract_by_note.py --db $Merged_DB --str Aug_RxLR_motif Aug_RxLR_EER_motif Aug_WY_hmm Aug_RxLR_hmm Aug_mimp_intersect ORF_RxLR_motif ORF_RxLR_EER_motif ORF_WY_hmm ORF_RxLR_hmm ORF_mimp_intersect --out $Effectors_Gff --type gene transcript
  for i in 1; do
  printf "Aug_RxLR_motif:\t"
  cat $Effectors_Gff | grep 'Aug_RxLR_motif' | wc -l
  printf "Aug_RxLR_EER_motif:\t"
  cat $Effectors_Gff | grep 'Aug_RxLR_EER_motif' | wc -l
  printf "Aug_WY_hmm:\t"
  cat $Effectors_Gff | grep 'Aug_WY_hmm' | wc -l
  printf "Aug_RxLR_hmm:\t"
  cat $Effectors_Gff | grep 'Aug_RxLR_hmm' | wc -l
  printf "Aug_mimp_intersect:\t"
  cat $Effectors_Gff | grep 'Aug_mimp_intersect' | wc -l
  printf "ORF_RxLR_motif:\t"
  cat $Effectors_Gff | grep 'ORF_RxLR_motif' | wc -l
  printf "ORF_RxLR_EER_motif:\t"
  cat $Effectors_Gff | grep 'ORF_RxLR_EER_motif' | wc -l
  printf "ORF_WY_hmm:\t"
  cat $Effectors_Gff | grep 'ORF_WY_hmm' | wc -l
  printf "ORF_WY_hmm:\t"
  cat $Effectors_Gff | grep 'ORF_RxLR_hmm' | wc -l
  printf "ORF_mimp_intersect:\t"
  cat $Effectors_Gff | grep 'ORF_mimp_intersect' | wc -l
  done
```

Aug_RxLR_motif:	17
Aug_RxLR_EER_motif:	1
Aug_WY_hmm:	0
Aug_RxLR_hmm:	2
Aug_mimp_intersect:	2
ORF_RxLR_motif:	564
ORF_RxLR_EER_motif:	16
ORF_WY_hmm:	4
ORF_WY_hmm:	30
ORF_mimp_intersect:	124

Higher confidence putative effectors were extracted using:
```bash
  $ProgDir/extract_by_note.py --db $Merged_DB --str Aug_RxLR_motif Aug_RxLR_EER_motif Aug_WY_hmm Aug_RxLR_hmm Aug_mimp_intersect ORF_RxLR_EER_motif ORF_WY_hmm ORF_RxLR_hmm ORF_mimp_intersect --out $Effectors_Gff --type gene transcript
  for i in 1; do
  printf "Aug_RxLR_motif:\t"
  cat $Effectors_Gff | grep 'Aug_RxLR_motif' | wc -l
  printf "Aug_RxLR_EER_motif:\t"
  cat $Effectors_Gff | grep 'Aug_RxLR_EER_motif' | wc -l
  printf "Aug_WY_hmm:\t"
  cat $Effectors_Gff | grep 'Aug_WY_hmm' | wc -l
  printf "Aug_RxLR_hmm:\t"
  cat $Effectors_Gff | grep 'Aug_RxLR_hmm' | wc -l
  printf "Aug_mimp_intersect:\t"
  cat $Effectors_Gff | grep 'Aug_mimp_intersect' | wc -l
  printf "ORF_RxLR_motif:\t"
  cat $Effectors_Gff | grep 'ORF_RxLR_motif' | wc -l
  printf "ORF_RxLR_EER_motif:\t"
  cat $Effectors_Gff | grep 'ORF_RxLR_EER_motif' | wc -l
  printf "ORF_WY_hmm:\t"
  cat $Effectors_Gff | grep 'ORF_WY_hmm' | wc -l
  printf "ORF_WY_hmm:\t"
  cat $Effectors_Gff | grep 'ORF_RxLR_hmm' | wc -l
  printf "ORF_mimp_intersect:\t"
  cat $Effectors_Gff | grep 'ORF_mimp_intersect' | wc -l
  printf "Total number of putative effectors:\t"
  cat $Effectors_Gff | grep 'gene' | wc -l
  done
```

Aug_RxLR_motif:	17
Aug_RxLR_EER_motif:	1
Aug_WY_hmm:	0
Aug_RxLR_hmm:	2
Aug_mimp_intersect:	2
ORF_RxLR_motif:	24
ORF_RxLR_EER_motif:	16
ORF_WY_hmm:	4
ORF_WY_hmm:	30
Total number of putative effectors:	188
