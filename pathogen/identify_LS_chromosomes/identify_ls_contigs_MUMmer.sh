nucmer --prefix=ref_qry analysis/identify_LS_contigs/F.oxysporum_fsp.lycopersici/4287/F.oxysporum_fsp.lycopersici_4287_LS_contigs.fasta assembly/velvet/F.oxysporum_fsp_cepae/125/125_assembly.71/sorted_contigs.fa 

show-coords -rcl ref_qry.delta > ref_qry.coords

mv ref* analysis/identify_LS_contigs/F.oxysporum_fsp_cepae/125/.
mkdir -p analysis/identify_LS_contigs/F.oxysporum_fsp_cepae/125/aligns/


show-tiling  analysis/identify_LS_contigs/F.oxysporum_fsp_cepae/125/ref_qry.delta > analysis/identify_LS_contigs/F.oxysporum_fsp_cepae/125/ref_qry.tiling

# .tiling file format
#	start	end	gap_to_next	length	%_query	%_id	orientation	contig_ID	
#	108102  134159  72266   26058   96.68   96.47   +       NODE_41_length_25988_cov_33.290710

while read line; do
	if [ $(echo $line | head -c 1) = ">" ]; then
    	REF_SEQ=$(echo $line | cut -d ' ' -f1 | sed 's/>//')
		echo "ref_seq set to: $REF_SEQ"
	else
		QRY_SEQ=$(echo $line | cut -d ' ' -f 8)
		echo "qry_seq set to: $QRY_SEQ"
		show-aligns analysis/identify_LS_contigs/F.oxysporum_fsp_cepae/125/ref_qry.delta $REF_SEQ $QRY_SEQ > analysis/identify_LS_contigs/F.oxysporum_fsp_cepae/125/aligns/ref_qry_"$REF_SEQ"_"$QRY_SEQ".aligns
	fi
done<analysis/identify_LS_contigs/F.oxysporum_fsp_cepae/125/ref_qry.tiling

# for PAIR in $(tail -n +6 analysis/identify_LS_contigs/F.oxysporum_fsp_cepae/125/ref_qry.coords | cut -d '|' -f 7 | sed 's/\s/:/g' | sort | uniq); do
# 	REF_SEQ=$(printf $PAIR | cut -d ':' -f 2)
# 	QRY_SEQ=$(printf $PAIR | cut -d ':' -f 3)
# #	echo "$REF_SEQ dont take no shit from $QRY_SEQ"
# done