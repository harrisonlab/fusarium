

```bash
  InBam=alignment/F.oxysporum_fsp_cepae/Fus2/Fus2_72hrs_rep1/accepted_hits.bam
  OutFq=tmp.fastq
  bamtools convert -in $InBam -out $OutFq -format fastq
  mkdir tmp2
  cd tmp2
  Trinity.pl --seqType fq --JM 4G --single ../$OutFq --run_as_paired --CPU 1 --output tmp_trinity
```



```bash
  cd /home/groups/harrisonlab/project_files/fusarium
  mkdir tmp3
  cd tmp3
  AcceptedHitsList=accepted_hits_list.txt
  OutFq=Fus2_aligned_reads.fq
  ls ../alignment/F.oxysporum_fsp_cepae/Fus2/*Fus2*/accepted_hits.bam > $AcceptedHitsList
  bamtools convert -list $AcceptedHitsList -out $OutFq -format fastq
  cd tmp3
  Trinity.pl --seqType fq --JM 64G --single $OutFq --run_as_paired --CPU 16 --output tmp_trinity
```

```bash
  gmap_build -d Fus2_genome -D . -k 13 ../repeat_masked/F.oxysporum_fsp_cepae/Fus2_edited_v2/filtered_contigs_repmask/Fus2_edited_v2_contigs_unmasked.fa
  gmap -n 0 -D . -d Fus2_genome tmp_trinity/Trinity.fasta -f gff3_gene > trinity_gmap.gff3
  gmap -n 0 -D . -d Fus2_genome tmp_trinity/Trinity.fasta -f gff3_match_cdna > trinity_gmap_cdna.gff3
  gmap -n 0 -D . -d Fus2_genome tmp_trinity/Trinity.fasta -f gff3_match_est > trinity_gmap_est.gff3
```
