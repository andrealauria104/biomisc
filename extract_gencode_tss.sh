#!/bin/bash
ANNOT=$1
BED=${2:-true}

echo ""
echo " [+] Extract TSS, annotation = ${ANNOT}"
echo ""

if [[ $BED == true ]]; then
	OUTFILE=${ANNOT%.gtf}.tss.bed
	echo "  -- output BED format (0-based), file = ${OUTFILE}"
	echo ""
	awk 'BEGIN{OFS="\t"};$3 == "gene" {
		if($7 == "+") tss = $4;
		else if($7 == "-") tss = $5;
		gsub("\"|;","");
		print $1,tss-1,tss,$7,$10}
		' < $ANNOT | sort -k1,1 -k2,2n -V > $OUTFILE
else
	OUTFILE=${ANNOT%.gtf}.tss.txt
	echo "  -- output file = ${OUTFILE}"
	echo ""
	awk 'BEGIN{OFS="\t"};$3 == "gene" {
		if($7 == "+") tss = $4;
		else if($7 == "-") tss = $5;
		gsub("\"|;","");
		print $1,tss,$7,$10}
		' < $ANNOT > $OUTFILE
fi


