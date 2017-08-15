

```bash
WorkDir=/home/groups/harrisonlab/project_files/fusarium
cat $WorkDir/genome_submission/F.oxysporum_fsp_cepae/Fus2/gag/edited/genome2_gene_conversions.tsv | sed 's/vAg/g/g' > $WorkDir/genome_submission/F.oxysporum_fsp_cepae/Fus2/gag/edited/genome2_gene_conversions_parsed.tsv

cat  5\ Supplementary\ Table\ 5.txt | sed 's/^M/\n/g' | tail -n+2 > 5_Supplementary_Table_5.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/rename_genes
$ProgDir/rename_genes.py --new_names $WorkDir/genome_submission/F.oxysporum_fsp_cepae/Fus2/gag/edited/genome2_gene_conversions_parsed.tsv --input_table 5_Supplementary_Table_5.txt --ignore_header True --id_column 2 > 5_Supplementary_Table_5_renamed.txt

# every time the ^M character is shown then you have to edit with ctrl+V ctrl+M manually
cat  helpful\ files/Supp.\ Tab.\ 2.txt  | sed 's/^M/\n/g' | tail -n+2 | head -n-3 > Supplementary_Table_2.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/rename_genes
$ProgDir/rename_genes.py --new_names $WorkDir/genome_submission/F.oxysporum_fsp_cepae/Fus2/gag/edited/genome2_gene_conversions_parsed.tsv --input_table Supplementary_Table_2.txt --ignore_header True --id_column 3 > Supplementary_Table_2_renamed.txt


cat  9\ Supplementary\ Table\ 9.txt  | sed 's/^M/\n/g' | tail -n+2 > Supplementary_Table_9.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/rename_genes
$ProgDir/rename_genes.py --new_names $WorkDir/genome_submission/F.oxysporum_fsp_cepae/Fus2/gag/edited/genome2_gene_conversions_parsed.tsv --input_table Supplementary_Table_9.txt --ignore_header True --id_column 4 > 9_Supplementary_Table_9_renamed.txt

```
