# sixgenes2sixcontigs.sh


for STRAIN in assembly/velvet/F.oxysporum_fsp_cepae/*/; do
	STRAIN=$(echo $STRAIN | cut -d '/' -f4)
	mkdir -p analysis/blast_homology/'CDCs'/F.oxysporum_fsp_cepae/"$STRAIN"/
done

cat analysis/blast_homology/F.oxysporum_fsp_cepae/125/125_six-appended_parsed.fa_homologs.csv | tail -n +2 | cut -f 9 | sort | uniq > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/125/125_six_contigs.txt 
printf '' > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/125/125_six_contigs.fa
printf '' > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/125/125_six_contigs.gff; 
while read line; do 
	grep -A1 -w ">$line" assembly/velvet/F.oxysporum_fsp_cepae/125/125_assembly.71/sorted_contigs.fa >> analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/125/125_six_contigs.fa; 
	grep -A1 -w "$line" gene_pred/augustus/F.oxysporum_fsp_cepae/125/125_aug.gff >> analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/125/125_six_contigs.gff; 
done<analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/125/125_six_contigs.txt

cat analysis/blast_homology/F.oxysporum_fsp_cepae/55/55_six-appended_parsed.fa_homologs.csv | tail -n +2 | cut -f 9 | sort | uniq > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/55/55_six_contigs.txt 
printf '' > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/55/55_six_contigs.fa
printf '' > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/55/55_six_contigs.gff; 
while read line; do 
	grep -A1 -w ">$line" assembly/velvet/F.oxysporum_fsp_cepae/55/55_assembly.61/sorted_contigs.fa >> analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/55/55_six_contigs.fa; 
	grep -A1 -w "$line" gene_pred/augustus/F.oxysporum_fsp_cepae/55/55_aug.gff >> analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/55/55_six_contigs.gff; 
done<analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/55/55_six_contigs.txt

cat analysis/blast_homology/F.oxysporum_fsp_cepae/A23/A23_six-appended_parsed.fa_homologs.csv | tail -n +2 | cut -f 9 | sort | uniq > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/A23/A23_six_contigs.txt 
printf '' > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/A23/A23_six_contigs.fa
printf '' > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/A23/A23_six_contigs.gff; 
while read line; do 
	grep -A1 -w ">$line" assembly/velvet/F.oxysporum_fsp_cepae/A23/A23_assembly.61/sorted_contigs.fa >> analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/A23/A23_six_contigs.fa; 
	grep -A1 -w "$line" gene_pred/augustus/F.oxysporum_fsp_cepae/A23/A23_aug.gff >> analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/A23/A23_six_contigs.gff; 
done<analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/A23/A23_six_contigs.txt

cat analysis/blast_homology/F.oxysporum_fsp_cepae/A28/A28_six-appended_parsed.fa_homologs.csv | tail -n +2 | cut -f 9 | sort | uniq > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/A28/A28_six_contigs.txt 
printf '' > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/A28/A28_six_contigs.fa
printf '' > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/A28/A28_six_contigs.gff; 
while read line; do 
	grep -A1 -w ">$line" assembly/velvet/F.oxysporum_fsp_cepae/A28/A28_assembly.61/sorted_contigs.fa >> analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/A28/A28_six_contigs.fa; 
	grep -A1 -w "$line" gene_pred/augustus/F.oxysporum_fsp_cepae/A28/A28_aug.gff >> analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/A28/A28_six_contigs.gff; 
done<analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/A28/A28_six_contigs.txt

cat analysis/blast_homology/F.oxysporum_fsp_cepae/D2/D2_six-appended_parsed.fa_homologs.csv | tail -n +2 | cut -f 9 | sort | uniq > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/D2/D2_six_contigs.txt 
printf '' > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/D2/D2_six_contigs.fa
printf '' > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/D2/D2_six_contigs.gff; 
while read line; do 
	grep -A1 -w ">$line" assembly/velvet/F.oxysporum_fsp_cepae/D2/D2_assembly.41/sorted_contigs.fa >> analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/D2/D2_six_contigs.fa; 
	grep -A1 -w "$line" gene_pred/augustus/F.oxysporum_fsp_cepae/D2/D2_aug.gff >> analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/D2/D2_six_contigs.gff; 
done<analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/D2/D2_six_contigs.txt

cat analysis/blast_homology/F.oxysporum_fsp_cepae/Fus2/Fus2_six-appended_parsed.fa_homologs.csv | tail -n +2 | cut -f 9 | sort | uniq > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/Fus2/Fus2_six_contigs.txt 
printf '' > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/Fus2/Fus2_six_contigs.fa
printf '' > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/Fus2/Fus2_six_contigs.gff; 
while read line; do 
	grep -A1 -w ">$line" assembly/velvet/F.oxysporum_fsp_cepae/Fus2/Fus2_assembly.81/sorted_contigs.fa >> analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/Fus2/Fus2_six_contigs.fa; 
	grep -A1 -w "$line" gene_pred/augustus/F.oxysporum_fsp_cepae/Fus2/Fus2_aug.gff >> analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/Fus2/Fus2_six_contigs.gff; 
done<analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/Fus2/Fus2_six_contigs.txt

cat analysis/blast_homology/F.oxysporum_fsp_cepae/HB17/HB17_six-appended_parsed.fa_homologs.csv | tail -n +2 | cut -f 9 | sort | uniq > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/HB17/HB17_six_contigs.txt 
printf '' > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/HB17/HB17_six_contigs.fa
printf '' > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/HB17/HB17_six_contigs.gff; 
while read line; do 
	grep -A1 -w ">$line" assembly/velvet/F.oxysporum_fsp_cepae/HB17/HB17_assembly.51/sorted_contigs.fa >> analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/HB17/HB17_six_contigs.fa; 
	grep -A1 -w "$line" gene_pred/augustus/F.oxysporum_fsp_cepae/HB17/HB17_aug.gff >> analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/HB17/HB17_six_contigs.gff; 
done<analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/HB17/HB17_six_contigs.txt

cat analysis/blast_homology/F.oxysporum_fsp_cepae/PG/PG_six-appended_parsed.fa_homologs.csv | tail -n +2 | cut -f 9 | sort | uniq > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/PG/PG_six_contigs.txt 
printf '' > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/PG/PG_six_contigs.fa
printf '' > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/PG/PG_six_contigs.gff; 
while read line; do 
	grep -A1 -w ">$line" assembly/velvet/F.oxysporum_fsp_cepae/PG/PG_assembly.81/sorted_contigs.fa >> analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/PG/PG_six_contigs.fa; 
	grep -A1 -w "$line" gene_pred/augustus/F.oxysporum_fsp_cepae/PG/PG_aug.gff >> analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/PG/PG_six_contigs.gff; 
done<analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/PG/PG_six_contigs.txt

for STRAIN in analysis/blast_homology/'CDCs'/F.oxysporum_fsp_cepae/*/; do
	STRAIN=$(echo $STRAIN | cut -d '/' -f5)
	/home/armita/git_repos/emr_repos/scripts/fusarium/tools/rename_fasta.pl analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/$STRAIN/"$STRAIN"_six_contigs.fa F.oxysporum_fsp_cepae $STRAIN > analysis/blast_homology/CDCs/F.oxysporum_fsp_cepae/$STRAIN/"$STRAIN"_six_contigs_renamed.fa
done

#	mv analysis/blast_homology/CDCs/"$STRAIN"_six_contigs_renamed.fa analysis/blast_homology/CDCs/"$STRAIN"_six_contigs.fa

