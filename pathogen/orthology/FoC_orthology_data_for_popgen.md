This work was performed in the following project directory:

```bash
ProjectDirectory=/home/groups/harrisonlab/project_files/fusarium
```

Genome assemblies can be found in:

```bash
ls repeat_masked/F.*/*/filtered_contigs_repmask/*_contigs_unmasked.fa| grep -e 'F.oxysporum_fsp_cepae' | grep -w -v 'Fus2'
```

```
repeat_masked/F.oxysporum_fsp_cepae/125/filtered_contigs_repmask/125_contigs_unmasked.fa
repeat_masked/F.oxysporum_fsp_cepae/55/filtered_contigs_repmask/55_contigs_unmasked.fa
repeat_masked/F.oxysporum_fsp_cepae/A1-2/filtered_contigs_repmask/A1-2_contigs_unmasked.fa
repeat_masked/F.oxysporum_fsp_cepae/A13/filtered_contigs_repmask/A13_contigs_unmasked.fa
repeat_masked/F.oxysporum_fsp_cepae/A23/filtered_contigs_repmask/A23_contigs_unmasked.fa
repeat_masked/F.oxysporum_fsp_cepae/A28/filtered_contigs_repmask/A28_contigs_unmasked.fa
repeat_masked/F.oxysporum_fsp_cepae/CB3/filtered_contigs_repmask/CB3_contigs_unmasked.fa
repeat_masked/F.oxysporum_fsp_cepae/D2/filtered_contigs_repmask/D2_contigs_unmasked.fa
repeat_masked/F.oxysporum_fsp_cepae/Fus2_edited_v2/filtered_contigs_repmask/Fus2_edited_v2_contigs_unmasked.fa
repeat_masked/F.oxysporum_fsp_cepae/HB17/filtered_contigs_repmask/HB17_contigs_unmasked.fa
repeat_masked/F.oxysporum_fsp_cepae/HB6/filtered_contigs_repmask/HB6_contigs_unmasked.fa
repeat_masked/F.oxysporum_fsp_cepae/PG/filtered_contigs_repmask/PG_contigs_unmasked.fa
```


Protein sequences for gene models can be found in:

```bash
  ls gene_pred/codingquary/F*/*/final/final_genes_combined.pep.fasta | grep -e 'F.oxysporum_fsp_cepae' | grep -w -v 'Fus2'
```

```
gene_pred/codingquary/F.oxysporum_fsp_cepae/125/final/final_genes_combined.pep.fasta
gene_pred/codingquary/F.oxysporum_fsp_cepae/55/final/final_genes_combined.pep.fasta
gene_pred/codingquary/F.oxysporum_fsp_cepae/A1-2/final/final_genes_combined.pep.fasta
gene_pred/codingquary/F.oxysporum_fsp_cepae/A13/final/final_genes_combined.pep.fasta
gene_pred/codingquary/F.oxysporum_fsp_cepae/A23/final/final_genes_combined.pep.fasta
gene_pred/codingquary/F.oxysporum_fsp_cepae/A28/final/final_genes_combined.pep.fasta
gene_pred/codingquary/F.oxysporum_fsp_cepae/CB3/final/final_genes_combined.pep.fasta
gene_pred/codingquary/F.oxysporum_fsp_cepae/D2/final/final_genes_combined.pep.fasta
gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/final/final_genes_combined.pep.fasta
gene_pred/codingquary/F.oxysporum_fsp_cepae/HB17/final/final_genes_combined.pep.fasta
gene_pred/codingquary/F.oxysporum_fsp_cepae/HB6/final/final_genes_combined.pep.fasta
gene_pred/codingquary/F.oxysporum_fsp_cepae/PG/final/final_genes_combined.pep.fasta
```

Nucleotide sequences for gene models can be foudn in:

```bash
  ls gene_pred/codingquary/F*/*/final/final_genes_combined.gene.fasta | grep -e 'F.oxysporum_fsp_cepae' | grep -w -v 'Fus2'
```

```
  gene_pred/codingquary/F.oxysporum_fsp_cepae/125/final/final_genes_combined.gene.fasta
  gene_pred/codingquary/F.oxysporum_fsp_cepae/55/final/final_genes_combined.gene.fasta
  gene_pred/codingquary/F.oxysporum_fsp_cepae/A1-2/final/final_genes_combined.gene.fasta
  gene_pred/codingquary/F.oxysporum_fsp_cepae/A13/final/final_genes_combined.gene.fasta
  gene_pred/codingquary/F.oxysporum_fsp_cepae/A23/final/final_genes_combined.gene.fasta
  gene_pred/codingquary/F.oxysporum_fsp_cepae/A28/final/final_genes_combined.gene.fasta
  gene_pred/codingquary/F.oxysporum_fsp_cepae/CB3/final/final_genes_combined.gene.fasta
  gene_pred/codingquary/F.oxysporum_fsp_cepae/D2/final/final_genes_combined.gene.fasta
  gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/final/final_genes_combined.gene.fasta
  gene_pred/codingquary/F.oxysporum_fsp_cepae/HB17/final/final_genes_combined.gene.fasta
  gene_pred/codingquary/F.oxysporum_fsp_cepae/HB6/final/final_genes_combined.gene.fasta
  gene_pred/codingquary/F.oxysporum_fsp_cepae/PG/final/final_genes_combined.gene.fasta
```

gff annotations of gene models can be foudn at:

```bash
ls gene_pred/codingquary/F*/*/final/final_genes_appended.gff3  | grep -e 'F.oxysporum_fsp_cepae' | grep -w -v 'Fus2'
```

```
  gene_pred/codingquary/F.oxysporum_fsp_cepae/125/final/final_genes_appended.gff3
  gene_pred/codingquary/F.oxysporum_fsp_cepae/55/final/final_genes_appended.gff3
  gene_pred/codingquary/F.oxysporum_fsp_cepae/A1-2/final/final_genes_appended.gff3
  gene_pred/codingquary/F.oxysporum_fsp_cepae/A13/final/final_genes_appended.gff3
  gene_pred/codingquary/F.oxysporum_fsp_cepae/A23/final/final_genes_appended.gff3
  gene_pred/codingquary/F.oxysporum_fsp_cepae/A28/final/final_genes_appended.gff3
  gene_pred/codingquary/F.oxysporum_fsp_cepae/CB3/final/final_genes_appended.gff3
  gene_pred/codingquary/F.oxysporum_fsp_cepae/D2/final/final_genes_appended.gff3
  gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/final/final_genes_appended.gff3
  gene_pred/codingquary/F.oxysporum_fsp_cepae/HB17/final/final_genes_appended.gff3
  gene_pred/codingquary/F.oxysporum_fsp_cepae/HB6/final/final_genes_appended.gff3
  gene_pred/codingquary/F.oxysporum_fsp_cepae/PG/final/final_genes_appended.gff3
```


downloaded genomes fo47 (non-pathogen) and 4287 (FoL pathogen) were also
included in the following analyses. Downloaded data for these isolates can be
found:

fo47:

```
  assembly/external_group/F.oxysporum/fo47/broad/
```

with the following proteins being used in the orthology analysis:

```
  assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_proteins.fasta
```

4287:

```
  assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1
```

with the following proteins being used in the orthology analysis:

```
  assembly/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522.aa.fasta
```

## Orthology results


Othology analysis was performed on these isolates.
The commands were documented in the following files under Methology 3:

https://github.com/harrisonlab/fusarium/blob/master/pathogen/orthology/F.oxysporum_fsp.cepae_pathogen_vs_non-pathogen_orthology.md

The total set of proteins used in the analysis can be found at:
- note from this point onwards each gene name is preceeded by an isolate identifier in the fasta files.

```
  analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/goodProteins/goodProteins.fasta
```

The list of genes present in each ortholog group is present at:

```bash
  analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/FoC_vs_Fo_vs_FoL_orthogroups.txt
```

Fasta files of orthogroups
can be found within the directory:

```bash
  analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/fasta/all_orthogroups
```

A venn diagram showing orthogroups common to all pathogens, common to all non-pathogens
or shared between both of these groups can be found at:

```bash
  analysis/orthology/orthomcl/FoC_vs_Fo_vs_FoL/FoC_vs_Fo_vs_FoL_orthogroups.pdf
```
