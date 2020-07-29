#!/bin/bash
METHRATIO=$1
PATTERN=${2:-"txt"}

wrap_awk () {
	awk 'BEGIN {OFS ="\t";print "sample","n_C","n_CT","ratio"} \
	NR>1 && $$4!="CG" {sum_C += $7; sum_CT += $8} \
	END {gsub(".*\\/|.txt","",FILENAME); \
	print FILENAME,sum_C,sum_CT,1-sum_C/sum_CT}'
}

calcratio () {	
	echo " -- processing file ${1}"
	if [[ $1 == "*.gz?" ]]; then
		gzip -cd $1 | awk 'BEGIN {OFS ="\t";print "sample","n_C","n_CT","ratio"} \
        			NR>1 && $$4!="CG" {sum_C += $7; sum_CT += $8} \
        			END {gsub(".*\\/|.txt","",FILENAME); \
        			print FILENAME,sum_C,sum_CT,1-sum_C/sum_CT}'  
	else
		awk 'BEGIN {OFS ="\t";print "sample","n_C","n_CT","ratio"} \
		NR>1 && $$4!="CG" {sum_C += $7; sum_CT += $8} \
        	END {gsub(".*\\/|.txt","",FILENAME); \
        	print FILENAME,sum_C,sum_CT,1-sum_C/sum_CT}' $1
	fi		
}

if [[ -d $METHRATIO ]]; then
	FILES=$(ls $METHRATIO | grep $PATTERN)
	for i in $FILES; do
		calcratio $METHRATIO/$i
	done
else
	calcratio $METHRATIO
fi

