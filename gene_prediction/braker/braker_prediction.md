
Braker_prediction.md

This document details commands used for an initial run of Braker to train and
predict gene models.

Braker is a pipeline based upon Augustus and GeneMark-est.

set variables
```bash
qlogin
WorkDir=/tmp/braker
ProjDir=/home/groups/harrisonlab/project_files/fusarium
Assembly=$ProjDir/repeat_masked/F.oxysporum_fsp_cepae/Fus2/Fus2_combined_49/Fus2_combined_49_contigs_unmasked.fa
OutDir=$ProjDir/gene_pred/braker/F.oxysporum_fsp_cepae/Fus2
```

move to working directory
```
mkdir -p $WorkDir
cd $WorkDir
```

merge .bam files into a single alignment file
```bash
  ls $ProjDir/timecourse/v2_genes/*/*/accepted_hits.bam > bamlist.txt
  # BamFiles=$(ls $ProjDir/timecourse/v2_genes/*/*/accepted_hits.bam | tr '\n' ' ')


  bamtools merge -list bamlist.txt -out timecourse/v2_genes/Fus2/merged/fus2_accepted_hits_merged.bam
```

<!-- Run Braker using examples:
```bash
braker.pl \
--cores 1 \
--fungus \
--SAMTOOLS_PATH=/home/armita/prog/samtools-0.1.18 \
--GENEMARK_PATH=/home/armita/prog/genemark/gm_et_linux_64/gmes_petap \
--BAMTOOLS_PATH=/home/armita/prog/bamtools/bamtools/bin \
--AUGUSTUS_CONFIG_PATH=/home/armita/prog/augustus-3.1/config \
--genome=/home/armita/prog/braker1/example.fa \
--species=new_species \
--bam=/home/armita/prog/braker1/example.bam
``` -->




```bash
  braker.pl \
    --cores 16 \
    --fungus \
    --genome=$Assembly \
    --GENEMARK_PATH=/home/armita/prog/genemark/gm_et_linux_64/gmes_petap \
    --BAMTOOLS_PATH=/home/armita/prog/bamtools/bamtools/bin \
    --species=F.oxysporum_fsp_cepae_Fus2 \
    --bam=$ProjDir/timecourse/v2_genes/Fus2/merged/fus2_accepted_hits_merged.bam
```

```bash
  mkdir -p $OutDir
  cp -r braker/* $OutDir/.

  rm -r $WorkDir
```

# Extract gff and amino acid sequences

```bash
  for File in $(ls gene_pred/braker/F.*/*/*/augustus.gff); do
    getAnnoFasta.pl $File
    OutDir=$(dirname $File)
    echo "##gff-version 3" > $OutDir/augustus_extracted.gff
    cat $File | grep -v '#' >> $OutDir/augustus_extracted.gff
  done
```
