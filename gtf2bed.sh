#!/bin/bash
INFILE=$1
OUTFILE="${INFILE%.gtf}.bed"

echo "[*] Convert gtf to bed [*]"
echo ""
echo " - input file: $INFILE"
echo " - output file: $OUTFILE"

awk 'OFS="\t" {print $1,$4-1,$5,$10,$16,$7}' $INFILE | tr -d '";' > $OUTFILE

echo ""
echo "done."
