
# Characterising_LS_contigs.md
---

### 4.1.a Assembly statistics

The length of core and LS regions was calculated from contigs sizes present in:
```bash
  cat analysis/circos/F.oxysporum_fsp_cepae/Fus2_reassembly/Fus2_genome.txt
```

### 4.1.b Repetative content

The % repeatmasking was identified for each contig:

```bash
OutDir=analysis/LS_contigs/F.oxysporum_fsp_cepae/Fus2_canu_new
mkdir -p $OutDir
Assembly=repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_hardmasked_repeatmasker_TPSI_appended.fa
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
$ProgDir/count_Ns.py --inp_fasta $Assembly --out_txt $OutDir/N_content.txt
```

### 4.1.c Gene density and % content

```bash
/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/feature_density.py --inp_fasta repeat_masked/F.oxysporum_fsp_cepae/Fus2_canu_new/edited_contigs_repmask/Fus2_canu_contigs_hardmasked_repeatmasker_TPSI_appended.fa --inp_gff gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3 --features gene
```


### 4.1.c
The number of putative effectors was identified in Core, lineage specific and species specific contigs:

	```bash
	echo "Number of genes in LS, spp. specific and core regions:"
	Gff=gene_pred/final_genes/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3
	cat $Gff | grep -w 'gene' | grep -w -e 'contig_10_pilon' -e 'contig_14_pilon' -e 'contig_16_pilon' -e 'contig_19_pilon' -e 'contig_20_pilon' -e 'contig_21_pilon' -e 'contig_22_pilon' | wc -l
	cat $Gff | grep -w 'gene' | grep -w -e 'contig_9_pilon' -e 'contig_11_pilon' -e 'contig_12_pilon' -e 'contig_13_pilon' -e 'contig_15_pilon' -e 'contig_17_pilon' -e 'contig_18_pilon' | wc -l
	cat $Gff | grep -w 'gene' | grep -w -e 'contig_1_pilon' -e 'contig_2_pilon' -e 'contig_3_pilon' -e 'contig_4_pilon' -e 'contig_5_pilon' -e 'contig_6_pilon' -e 'contig_7_pilon' -e 'contig_8_pilon' | wc -l
	```
	1927, 3731, 13126

	```bash
	echo "Number of Secreted CAZY genes in LS, spp. specific and core regions:"
	CazyGff=gene_pred/CAZY/F.oxysporum_fsp_cepae/Fus2_canu_new/Fus2_canu_new_CAZY_secreted.gff
	cat $CazyGff | grep -w 'gene' | grep -w -e 'contig_10_pilon' -e 'contig_14_pilon' -e 'contig_16_pilon' -e 'contig_19_pilon' -e 'contig_20_pilon' -e 'contig_21_pilon' -e 'contig_22_pilon' | wc -l
	cat $CazyGff | grep -w 'gene' | grep -w -e 'contig_9_pilon' -e 'contig_11_pilon' -e 'contig_12_pilon' -e 'contig_13_pilon' -e 'contig_15_pilon' -e 'contig_17_pilon' -e 'contig_18_pilon' | wc -l
	cat $CazyGff | grep -w 'gene' | grep -w -e 'contig_1_pilon' -e 'contig_2_pilon' -e 'contig_3_pilon' -e 'contig_4_pilon' -e 'contig_5_pilon' -e 'contig_6_pilon' -e 'contig_7_pilon' -e 'contig_8_pilon' | wc -l
	```
	18, 157, 211

	```bash
	echo "Number of Secreted effectorP genes in LS, spp. specific and core regions:"
	EffectorPGff=analysis/effectorP/F.oxysporum_fsp_cepae/Fus2_canu_new/F.oxysporum_fsp_cepae_Fus2_canu_new_EffectorP_secreted.gff
	cat $EffectorPGff | grep -w 'gene' | grep -w -e 'contig_10_pilon' -e 'contig_14_pilon' -e 'contig_16_pilon' -e 'contig_19_pilon' -e 'contig_20_pilon' -e 'contig_21_pilon' -e 'contig_22_pilon' | wc -l
	cat $EffectorPGff | grep -w 'gene' | grep -w -e 'contig_9_pilon' -e 'contig_11_pilon' -e 'contig_12_pilon' -e 'contig_13_pilon' -e 'contig_15_pilon' -e 'contig_17_pilon' -e 'contig_18_pilon' | wc -l
	cat $EffectorPGff | grep -w 'gene' | grep -w -e 'contig_1_pilon' -e 'contig_2_pilon' -e 'contig_3_pilon' -e 'contig_4_pilon' -e 'contig_5_pilon' -e 'contig_6_pilon' -e 'contig_7_pilon' -e 'contig_8_pilon' | wc -l
	```

	40, 111, 104
