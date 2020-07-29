#!/bin/bash
#
# Convert SAM to sorted and indexed BAM.
#
# Function definition ---
fsam2bam () {
	tmpbam=${1%.*}.tmp.bam
	outbam=${1%.*}.bam

	echo "[*] Converting SAM to BAM ..."

	if [ ! -f $1 ]; then
		echo "$1 does not exist."
		return 1
	fi

	samtools view -bS $1 > $tmpbam
	if [ $? -ne 0 ]; then
		echo "	[!] SAM2BAM conversion not sucessful."
		echo "	[!] $1 remains unchanged."
		rm $tmpbam
		return 1
	fi

	echo "[*] Sorting BAM ..."
	samtools sort -o $outbam $tmpbam 
	if [ $? -ne 0 ]; then
		echo "	[!] BAM file sorting not sucessful."
		echo "	[!] $outbam is in unsorted BAM format".
		mv $tmpbam $outbam
		return 1
	fi
	rm $tmpbam

	echo "[*] Indexing BAM ..."
	samtools index $outbam
	if [ $? -ne 0 ]; then
		echo "	[!] BAM file indexing not sucessful."
		return 1
	fi

	return 0
}

# execution ---
INPUT=$1
BAMDIR=${2:-BAM}

echo "============================================="
echo "=== Convert SAM to sorted and indexed BAM ==="
echo "============================================="
echo ""

if [ -f "${INPUT}" ]; then
	fsam2bam $INPUT
	echo "[*] done."
elif [ -d "${INPUT}" ]; then
	SAMDIR=$INPUT
	mkdir $BAMDIR
	for SAM in $SAMDIR/*.sam
	do
		echo " -- processing file: $SAM"
		fsam2bam $SAM
	done
	echo " -- moving files"
	mv $SAMDIR/*.bam* $BAMDIR
	echo "[*] done."
fi
