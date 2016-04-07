# Fus2Preds.md

This file documents commands used to run coding quary for gene prediction of the
FoC isolate Fus2.


## Trial Run

As this was a trial run of codingquary a qlogin sesison was started in screen.

```bash
  screen -a
  qlogin -pe smp 16 -l virtual_free=2G
  cd /home/groups/harrisonlab/project_files/fusarium
  ProjDir=$PWD
   WorkDir=/tmp/CodingQuarryPM
  mkdir -p $WorkDir
  cd $WorkDir
```

```bash
  Assembly=$ProjDir/repeat_masked/F.oxysporum_fsp_cepae/Fus2/filtered_contigs_repmask/Fus2_contigs_softmasked.fa
  CufflinksGTF=$ProjDir/gene_pred/cufflinks/F.oxysporum_fsp_cepae/Fus2/concatenated/transcripts.gtf
  CufflinksGFF3=$WorkDir/transcripts.gff3
  CufflinksGTF_to_CodingQuarryGFF3.py $CufflinksGTF > $CufflinksGFF3
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
  $ProgDir/run_CQ-PM_unstranded_EMR.sh $Assembly $CufflinksGFF3 2>&1 | tee codingquaryPM_log.txt
```

```bash
  OutDir=$ProjDir/gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2
  mkdir -p $OutDir/out
  mv codingquaryPM_log.txt $OutDir
  mv out/* $OutDir/out/.
  rm -r $WorkDir
```

## Through qsub

```bash
  cd /home/groups/harrisonlab/project_files/fusarium
  Assembly=repeat_masked/F.oxysporum_fsp_cepae/Fus2/filtered_contigs_repmask/Fus2_contigs_softmasked.fa
  CufflinksGTF=gene_pred/cufflinks/F.oxysporum_fsp_cepae/Fus2/concatenated/transcripts.gtf
  OutDir=gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
  qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
```


```bash
  cd /home/groups/harrisonlab/project_files/fusarium
  Assembly=repeat_masked/F.oxysporum_fsp_cepae/Fus2/filtered_contigs_repmask/Fus2_contigs_softmasked.fa
  CufflinksGTF=timecourse/2016_genes/Fus2/72hrs/cufflinks/transcripts.gtf
  OutDir=gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_72hr
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
  qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
```


# Supplimenting BRaker gene models with CodingQuary genes

```bash
  BrakerGff=gene_pred/braker/F.oxysporum_fsp_cepae/Fus2/F.oxysporum_fsp_cepae_Fus2_braker/augustus_extracted.gff
  CodingQuaryGff=gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2/out/PredictedPass.gff3
  PGNGff=gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2/out/PGN_predictedPass.gff3
  OutDir=gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2/additional
  AddGenesList=$OutDir/additional_genes.txt
  AddGenesGff=$OutDir/additional_genes.gff
  FinalGff=$OutDir/combined_genes.gff
  mkdir -p $OutDir

  bedtools intersect -v -s -a $CodingQuaryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
  bedtools intersect -v -s -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
  ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
  $ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
  $ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
  cat $BrakerGff $AddGenesGff | bedtools sort > $FinalGff
```



# Mining intergenic regions of Braker predictions for Path genes:

CodingQuary can be used to search between gene features:

```bash
  screen -a
  qlogin -pe smp 16 -l virtual_free=2G
  cd /home/groups/harrisonlab/project_files/fusarium
  ProjDir=$PWD
   WorkDir=/tmp/CodingQuarryPM
  mkdir -p $WorkDir
  cd $WorkDir
```

```bash
  Assembly=$ProjDir/repeat_masked/F.oxysporum_fsp_cepae/Fus2/filtered_contigs_repmask/Fus2_contigs_softmasked.fa
  BrakerPreds=$ProjDir/gene_pred/braker/F.oxysporum_fsp_cepae/Fus2_fungi_softmasked/F.oxysporum_fsp_cepae_Fus2_braker_softmasked/augustus.gtf
  CufflinksGTF=$ProjDir/gene_pred/cufflinks/F.oxysporum_fsp_cepae/Fus2/concatenated/transcripts.gtf
  CufflinksGFF3=$WorkDir/transcripts.gff3
  CufflinksGTF_to_CodingQuarryGFF3.py $CufflinksGTF > $CufflinksGFF3
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
  $ProgDir/run_CQ-PM_unstranded_EMR.sh $Assembly $CufflinksGFF3 2>&1 | tee codingquaryPM_log.txt
```
