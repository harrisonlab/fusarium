set -- 125 55 A23 A28 D2 Fus2 HB17 PG
for a; do 
	shift
	for b; do
	printf "%s - %s\n" "$a" "$b"
	STRAIN1="$a"
	STRAIN2="$b"
	FEATURES2=../../../gene_pred/augustus/F.oxysporum_fsp_cepae/"$b"/"$b"_aug.gff
	cd "$a"-"$b"
	FEATURENAME2="$STRAIN2"_uniq_vs_"$STRAIN1"
	cut -f4 table."$STRAIN1"-"$STRAIN2" | sed "s/ $STRAIN2/\n$STRAIN2/g" | cut -d' ' -f1 | tail -n+2 | cat - "$STRAIN2"_seqs.txt | sort | uniq -u > "$STRAIN2"_uniq_vs_"$STRAIN1".txt
	cat "$STRAIN2"_uniq_vs_"$STRAIN1".txt | cut -d '|' -f2 | xargs -I{} grep -w {} $FEATURES2 | sed "s/AUGUSTUS/$FEATURENAME2/" > "$STRAIN2"_uniq_vs_"$STRAIN1".gff
	/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/filter_gff_StartStop.pl "$STRAIN2"_uniq_vs_"$STRAIN1".gff > "$STRAIN2"_uniq_vs_"$STRAIN1"_filtered.gff
	cd ../	
	done 
done

grep 'stop_codon' -c *-*/*_uniq_vs_*_filtered.gff


set -- 125 55 A23 A28 D2 Fus2 HB17 PG
for a; do 
	shift
	for b; do
	printf "%s - %s\n" "$a" "$b"
	STRAIN1="$a"
	STRAIN2="$b"
	cd "$a"-"$b"
	grep "1.000 $a" table."$a"-"$b" > "$a"_duplications.txt
	grep "1.000 $b" table."$a"-"$b" > "$b"_duplications.txt
	cat "$a"_duplications.txt "$b"_duplications.txt | sort | uniq > "$a"-"$b"_duplications.txt
	cd ../	
	done 
done

