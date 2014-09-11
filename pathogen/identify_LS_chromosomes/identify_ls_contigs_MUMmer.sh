#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l virtual_free=0.9G

USAGE="identify_ls_contigs_MUMmer.sh <reference_contigs.fa> <assembled_contigs.fa>"

CUR_PATH=$PWD

REF_CONTIGS=$1
QRY_CONTIGS=$2

IN_REF=$(echo $REF_CONTIGS | rev | cut -d "/" -f1 | rev)

ORGANISM=$(echo $QRY_CONTIGS | rev | cut -d "/" -f4 | rev)
STRAIN=$(echo $QRY_CONTIGS | rev | cut -d "/" -f3 | rev)
IN_QRY=$(echo $QRY_CONTIGS | rev | cut -d "/" -f1 | rev)
OUTFILE=$(echo $QRY_CONTIGS | sed s/.aa//)

WORK_DIR=$TMPDIR/"$STRAIN"_MUMmer
mkdir -p $WORK_DIR
cd $WORK_DIR

cp $CUR_PATH/$REF_CONTIGS .
cp $CUR_PATH/$QRY_CONTIGS .

nucmer --prefix="$STRAIN"_MUMmer $IN_REF $IN_QRY 
show-coords -rcl "$STRAIN"_MUMmer.delta > "$STRAIN"_MUMmer.coords
show-tiling "$STRAIN"_MUMmer.delta > "$STRAIN"_MUMmer.tiling

mkdir MUMmer_alignments

while read line; do
	if [ $(echo $line | head -c 1) = ">" ]; then
    	REF_SEQ=$(echo $line | cut -d ' ' -f1 | sed 's/>//')
		echo "ref_seq set to: $REF_SEQ"
	else
		QRY_SEQ=$(echo $line | cut -d ' ' -f 8)
		echo "qry_seq set to: $QRY_SEQ"
		show-aligns "$STRAIN"_MUMmer.delta $REF_SEQ $QRY_SEQ > MUMmer_alignments/"$REF_SEQ"_"$QRY_SEQ".aligns
	fi
done<"$STRAIN"_MUMmer.tiling

mkdir -p $CUR_PATH/analysis/identify_LS_contigs/$ORGANISM/$STRAIN/
rm $IN_REF $IN_QRY
cp -r $WORK_DIR/ $CUR_PATH/analysis/identify_LS_contigs/$ORGANISM/$STRAIN/.
rm -r $TMPDIR

exit