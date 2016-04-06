# Fus2Preds.md

This file documents commands used to run coding quary for gene prediction of the
FoC isolate Fus2.


## login to session

As this was a trial run of codingquary a qlogin sesison was started in screen.

```bash
  screen -a
  qlogin -pe smp 16 -l virtual_free=2G
  cd /home/groups/harrisonlab/project_files/fusarium
  # WorkDir=$TMPDIR/CodingQuarry
  ProjDir=$PWD
  # WorkDir=/tmp/CodingQuarry
   WorkDir=/tmp/CodingQuarryPM
  mkdir -p $WorkDir
  cd $WorkDir
```

```bash
  Assembly=$ProjDir/repeat_masked/F.oxysporum_fsp_cepae/Fus2/filtered_contigs_repmask/Fus2_contigs_softmasked.fa
  CufflinksGTF=$ProjDir/gene_pred/cufflinks/F.oxysporum_fsp_cepae/Fus2/concatenated/transcripts.gtf
  CufflinksGFF3=$WorkDir/transcripts.gff3
  CufflinksGTF_to_CodingQuarryGFF3.py $CufflinksGTF > $CufflinksGFF3
  # CodingQuarry -f $Assembly -t $CufflinksGFF3 -p 16 -d 2>&1 | tee codingquary_log.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
  $ProgDir/run_CQ-PM_unstranded_EMR.sh $Assembly $CufflinksGFF3 2>&1 | tee codingquaryPM_log.txt
```
