# Files available for population genetic work


## location

This work was performed in the following project directory:

```bash
ProjectDirectory=/home/groups/harrisonlab/project_files/fusarium
```

## Assemblies

Genome assemblies can be found in:

```bash
ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa| grep -e 'F.oxysporum_fsp_cepae' -e 'narcissi' -e 'F.proliferatum' | grep -v -e 'edited' -e 'Fus2_contigs' -e 'HB17'
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
repeat_masked/F.oxysporum_fsp_cepae/HB6/filtered_contigs_repmask/HB6_contigs_unmasked.fa
repeat_masked/F.oxysporum_fsp_cepae/PG/filtered_contigs_repmask/PG_contigs_unmasked.fa
repeat_masked/F.oxysporum_fsp_narcissi/N139/filtered_contigs_repmask/N139_contigs_unmasked.fa
repeat_masked/F.proliferatum/A8/filtered_contigs_repmask/A8_contigs_unmasked.fa
```

## Gene models

Protein sequences for gene models can be found in:

```bash
  ls gene_pred/codingquary/F*/*/final/final_genes_combined.pep.fasta | grep -e 'F.oxysporum_fsp_cepae' -e 'narcissi' -e 'F.proliferatum' | grep -v -e 'edited' -e 'Fus2_contigs' -e 'HB17'
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
gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_combined.pep.fasta
gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2/final/final_genes_combined.pep.fasta
gene_pred/codingquary/F.oxysporum_fsp_cepae/HB6/final/final_genes_combined.pep.fasta
gene_pred/codingquary/F.oxysporum_fsp_cepae/PG/final/final_genes_combined.pep.fasta
gene_pred/codingquary/F.oxysporum_fsp_narcissi/N139/final/final_genes_combined.pep.fasta
gene_pred/codingquary/F.proliferatum/A8/final/final_genes_combined.pep.fasta
```

Nucleotide sequences for gene models can be foudn in:

```bash
  ls gene_pred/codingquary/F*/*/final/final_genes_combined.gene.fasta | grep -e 'F.oxysporum_fsp_cepae' -e 'narcissi' -e 'F.proliferatum' | grep -v -e 'edited' -e 'Fus2_contigs' -e 'HB17'
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
  gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_combined.gene.fasta
  gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2/final/final_genes_combined.gene.fasta
  gene_pred/codingquary/F.oxysporum_fsp_cepae/HB6/final/final_genes_combined.gene.fasta
  gene_pred/codingquary/F.oxysporum_fsp_cepae/PG/final/final_genes_combined.gene.fasta
  gene_pred/codingquary/F.oxysporum_fsp_narcissi/N139/final/final_genes_combined.gene.fasta
  gene_pred/codingquary/F.proliferatum/A8/final/final_genes_combined.gene.fasta
```

gff annotations of gene models can be found at:

```bash
ls gene_pred/codingquary/F*/*/final/final_genes_appended.gff3  | grep -e 'F.oxysporum_fsp_cepae' -e 'narcissi' -e 'F.proliferatum' | grep -v -e 'edited' -e 'Fus2_contigs' -e 'HB17'
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
gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_canu_new/final/final_genes_appended.gff3
gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2/final/final_genes_appended.gff3
gene_pred/codingquary/F.oxysporum_fsp_cepae/HB6/final/final_genes_appended.gff3
gene_pred/codingquary/F.oxysporum_fsp_cepae/PG/final/final_genes_appended.gff3
gene_pred/codingquary/F.oxysporum_fsp_narcissi/N139/final/final_genes_appended.gff3
gene_pred/codingquary/F.proliferatum/A8/final/final_genes_appended.gff3
```

## Downloaded genomes

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

## CEGMA

Searches for core Eukaryotic genes was performed with CEGMA. Output files can be
found at:

```bash
ls gene_pred/cegma/*/*/*_dna_cegma.cegma.gff | grep -e 'F.oxysporum_fsp_cepae' -e 'narcissi' -e 'F.proliferatum' | grep -v -e 'edited' -e 'Fus2_contigs' -e 'HB17'
```

```
  gene_pred/cegma/F.oxysporum_fsp_cepae/125/125_dna_cegma.cegma.gff
  gene_pred/cegma/F.oxysporum_fsp_cepae/55/55_dna_cegma.cegma.gff
  gene_pred/cegma/F.oxysporum_fsp_cepae/A1-2/A1-2_dna_cegma.cegma.gff
  gene_pred/cegma/F.oxysporum_fsp_cepae/A13/A13_dna_cegma.cegma.gff
  gene_pred/cegma/F.oxysporum_fsp_cepae/A23/A23_dna_cegma.cegma.gff
  gene_pred/cegma/F.oxysporum_fsp_cepae/A28/A28_dna_cegma.cegma.gff
  gene_pred/cegma/F.oxysporum_fsp_cepae/CB3/CB3_dna_cegma.cegma.gff
  gene_pred/cegma/F.oxysporum_fsp_cepae/D2/D2_dna_cegma.cegma.gff
  gene_pred/cegma/F.oxysporum_fsp_cepae/Fus2/Fus2_dna_cegma.cegma.gff
  gene_pred/cegma/F.oxysporum_fsp_cepae/HB17/HB17_dna_cegma.cegma.gff
  gene_pred/cegma/F.oxysporum_fsp_cepae/HB6/HB6_dna_cegma.cegma.gff
  gene_pred/cegma/F.oxysporum_fsp_cepae/PG/PG_dna_cegma.cegma.gff
```

## Annotation

Interproscan was used to annotate predicted proteins:

```bash
  ls gene_pred/interproscan/F.*/*/*_interproscan.tsv  | grep -e 'F.oxysporum_fsp_cepae' -e 'narcissi' -e 'F.proliferatum' | grep -v -e 'edited' -e 'Fus2_contigs' -e 'HB17'
```

```
  gene_pred/interproscan/F.oxysporum/fo47/fo47_interproscan.tsv
  gene_pred/interproscan/F.oxysporum_fsp_cepae/125/125_interproscan.tsv
  gene_pred/interproscan/F.oxysporum_fsp_cepae/55/55_interproscan.tsv
  gene_pred/interproscan/F.oxysporum_fsp_cepae/A1-2/A1-2_interproscan.tsv
  gene_pred/interproscan/F.oxysporum_fsp_cepae/A13/A13_interproscan.tsv
  gene_pred/interproscan/F.oxysporum_fsp_cepae/A23/A23_interproscan.tsv
  gene_pred/interproscan/F.oxysporum_fsp_cepae/A28/A28_interproscan.tsv
  gene_pred/interproscan/F.oxysporum_fsp_cepae/CB3/CB3_interproscan.tsv
  gene_pred/interproscan/F.oxysporum_fsp_cepae/D2/D2_interproscan.tsv
  gene_pred/interproscan/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_interproscan.tsv
  gene_pred/interproscan/F.oxysporum_fsp_cepae/HB17/HB17_interproscan.tsv
  gene_pred/interproscan/F.oxysporum_fsp_cepae/HB6/HB6_interproscan.tsv
  gene_pred/interproscan/F.oxysporum_fsp_cepae/PG/PG_interproscan.tsv
  gene_pred/interproscan/F.oxysporum_fsp_lycopersici/4287/4287_interproscan.tsv
```

as was swissprot:

```bash
ls gene_pred/swissprot/*/*/swissprot_v2015_tophit_parsed.tbl  | grep -e 'F.oxysporum_fsp_cepae' -e 'fo47' -e '4287' | grep -w -v 'Fus2'
```

```
gene_pred/swissprot/F.oxysporum/fo47/swissprot_v2015_tophit_parsed.tbl
gene_pred/swissprot/F.oxysporum_fsp_cepae/125/swissprot_v2015_tophit_parsed.tbl
gene_pred/swissprot/F.oxysporum_fsp_cepae/55/swissprot_v2015_tophit_parsed.tbl
gene_pred/swissprot/F.oxysporum_fsp_cepae/A1-2/swissprot_v2015_tophit_parsed.tbl
gene_pred/swissprot/F.oxysporum_fsp_cepae/A13/swissprot_v2015_tophit_parsed.tbl
gene_pred/swissprot/F.oxysporum_fsp_cepae/A23/swissprot_v2015_tophit_parsed.tbl
gene_pred/swissprot/F.oxysporum_fsp_cepae/A28/swissprot_v2015_tophit_parsed.tbl
gene_pred/swissprot/F.oxysporum_fsp_cepae/CB3/swissprot_v2015_tophit_parsed.tbl
gene_pred/swissprot/F.oxysporum_fsp_cepae/D2/swissprot_v2015_tophit_parsed.tbl
gene_pred/swissprot/F.oxysporum_fsp_cepae/Fus2_edited_v2/swissprot_v2015_tophit_parsed.tbl
gene_pred/swissprot/F.oxysporum_fsp_cepae/HB17/swissprot_v2015_tophit_parsed.tbl
gene_pred/swissprot/F.oxysporum_fsp_cepae/HB6/swissprot_v2015_tophit_parsed.tbl
gene_pred/swissprot/F.oxysporum_fsp_cepae/PG/swissprot_v2015_tophit_parsed.tbl
gene_pred/swissprot/F.oxysporum_fsp_lycopersici/4287/swissprot_v2015_tophit_parsed.tbl
```

## Final annotation tables:

Annotation tables summarising gene locations, orthology and function and
expression were generated:

```bash
ls gene_pred/annotations/*/*/*_gene_annotations.tab | grep -e 'F.oxysporum_fsp_cepae' -e 'fo47' -e '4287' | grep -w -v 'Fus2'
```

```
gene_pred/annotations/F.oxysporum/fo47/fo47_gene_annotations.tab
gene_pred/annotations/F.oxysporum_fsp_cepae/125/125_gene_annotations.tab
gene_pred/annotations/F.oxysporum_fsp_cepae/55/55_gene_annotations.tab
gene_pred/annotations/F.oxysporum_fsp_cepae/A1-2/A1-2_gene_annotations.tab
gene_pred/annotations/F.oxysporum_fsp_cepae/A13/A13_gene_annotations.tab
gene_pred/annotations/F.oxysporum_fsp_cepae/A23/A23_gene_annotations.tab
gene_pred/annotations/F.oxysporum_fsp_cepae/A28/A28_gene_annotations.tab
gene_pred/annotations/F.oxysporum_fsp_cepae/CB3/CB3_gene_annotations.tab
gene_pred/annotations/F.oxysporum_fsp_cepae/D2/D2_gene_annotations.tab
gene_pred/annotations/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_gene_annotations.tab
gene_pred/annotations/F.oxysporum_fsp_cepae/HB6/HB6_gene_annotations.tab
gene_pred/annotations/F.oxysporum_fsp_cepae/PG/PG_gene_annotations.tab
gene_pred/annotations/F.oxysporum_fsp_lycopersici/4287/4287_gene_annotations.tab
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
