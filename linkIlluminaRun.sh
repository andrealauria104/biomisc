#!/bin/bash
RUNDIR=$1
PATTERN=$2

if [[ -z "$RUNDIR" || -z "$PATTERN" ]]
then
        echo "[!] Missing reguired arguments: linkIlluminaRun <path_to_run> <pattern_in_fastq>"
else
        for i in $RUNDIR/$PATTERN
        do      
                echo " -- symlink fastq file: ${i}"
                ln -s $i $( basename $i | sed -e 's/_S[[:digit:]]\+\|_001//g')
        done
fi
