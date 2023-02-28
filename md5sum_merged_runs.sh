#!/usr/bin/env bash
# Variables
CONFIGFILE=$1
DRY=${2:-false}

# Check variables
if [ ! -e $CONFIGFILE ] || [ -z $CONFIGFILE ]; then
    echo ""
    echo "[!] Error: please provide valid path to merge-config file."
	echo ""
    exit 1
fi
echo "================================="
echo "= Check merged runs with md5sum ="
echo "================================="
echo ""
echo "Config file: ${CONFIGFILE}"

for i in $(awk 'NR>1{print $2}' $CONFIGFILE | sort -u);do
	echo -e "\nprocessing file: ${i}"
	if [ $DRY = true ]; then
		echo "zcat ${i} | md5sum"
		echo "zcat $(grep $i $CONFIGFILE | cut -f1) | md5sum"
	else
		echo "zcat ${i} | md5sum"
		OUT_DEST=$(zcat ${i} | md5sum | cut -d' ' -f1)
		echo "zcat $(grep $i $CONFIGFILE | cut -f1) | md5sum"
		OUT_SRC=$(zcat $(grep $i $CONFIGFILE | cut -f1) | md5sum | cut -d' ' -f1)
		
		echo ""
		echo $OUT_DEST
		echo $OUT_SRC
		
		if [ $OUT_SRC = $OUT_DEST ];then
			echo -e "\n[+] Merging is correct for ${i}"
		else
			echo -e "\n[!] Merging is wrong for ${i}"
			echo ${i} >> wrong_merging_list.txt
		fi
	fi
done