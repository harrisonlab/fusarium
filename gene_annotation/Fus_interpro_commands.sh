#!/usr/bin/bash
# Fusarium interproscan commands

SCRIPT_DIR=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan

for strain_path in $(ls -d gene_pred/augustus/F.oxysporum_fsp_cepae/*); do 
	echo $strain_path
	STRAIN=$(basename $strain_path)
	echo $STRAIN
	$SCRIPT_DIR/sub_interproscan.sh $strain_path/"$STRAIN"_aug_out.aa
done

for strain_path in $(ls -d gene_pred/interproscan/*/*); do 
	STRAIN=$(basename $strain_path)
	rm $strain_path/*_split_*
	printf "" > $strain_path/"$STRAIN"_interproscan.tsv
	printf "" > $strain_path/"$STRAIN"_interproscan.gff3
	printf "" > $strain_path/"$STRAIN"_interproscan.xml
	for file in $(ls -v $strain_path/split/"$STRAIN"_aug_out.aa_split_*.tsv); do 
		cat $file >> $strain_path/"$STRAIN"_interproscan.tsv
	done
	for file in $(ls -v $strain_path/split/"$STRAIN"_aug_out.aa_split_*.gff3); do 
		cat $file >> $strain_path/"$STRAIN"_interproscan.gff3
	done
	for file in $(ls -v $strain_path/split/"$STRAIN"_aug_out.aa_split_*.xml); do 
		cat $file >> $strain_path/"$STRAIN"_interproscan.xml
	done
done

cat gene_pred/interproscan/Fus2/split/Fus2_interproscan.gff3 | grep -E "^g" | sort -t g -k2 -n > gene_pred/interproscan/Fus2/split/Fus2_interproscan5.gff3

#-------------------------------------------
# Make the .gff3 file readable by geneious.
#-------------------------------------------
## Split each .gff3 file into gff features
## and corresponding fasta sequence files
## append all the features into one file
## append all the sequence into one file
## append the features and sequence file.
# 
# mkdir -p tmp3
# printf "" > tmp3/interpro_features.gff
# printf "" > tmp3/interpro_seqs.fa
# printf "" > tmp3/interpro.gff3
# for FILE in $(ls -v gene_pred/interproscan/Fus2/split/raw/Fus2_aug_out.aa_split_*.gff3); do 
# 	FILE_NAME=$(basename $FILE)
# #FILE=gene_pred/interproscan/Fus2/split/raw/Fus2_aug_out.aa_split_500.fa.gff3
# 	FASTA_START=$(cat $FILE | grep -E "^##FASTA" -n | cut -d ':' -f1)
# 	cat $FILE | head -n "$FASTA_START" | grep -v -E "^#" >> tmp3/interpro_features.gff
# 	FASTA_ENDS=$(cat $FILE | grep -E "^>match" -m 1 -n -B1 | head -n1 | cut -d '-' -f1)
# 	cat $FILE | head -n "$FASTA_ENDS" | tail -n +"$FASTA_START" | grep -v -E "^#" >> tmp3/interpro_seqs.fa
# done
# cat tmp3/interpro_features.gff tmp3/interpro_seqs.fa > tmp3/interpro.gff3

# This means that genes without features 
# don't appear in the .gff3 file.
# Need to use all predicted genes instead.

mkdir -p tmp3
printf "" > tmp3/interpro_features.gff
printf "" > tmp3/interpro_seqs.fa
printf "" > tmp3/interpro.gff3
for FILE in $(ls -v gene_pred/interproscan/Fus2/split/raw/Fus2_aug_out.aa_split_*.gff3); do 
	FILE_NAME=$(basename $FILE)
	FASTA_START=$(cat $FILE | grep -E "^##FASTA" -n | cut -d ':' -f1)
	cat $FILE | head -n "$FASTA_START" | grep -v -E "^#" >> tmp3/interpro_features.gff
done
cat tmp3/interpro_features.gff gene_pred/augustus/F.oxysporum_fsp_cepae/Fus2/Fus2_aug_out.aa > tmp3/interpro.gff3

