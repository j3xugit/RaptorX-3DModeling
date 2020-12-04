#!/bin/bash

if [ $# -lt 1 ]; then 
	echo $0 "seqFile"
	echo "	This script estimates the amount of GPU memory needed to predict distance matrix of a protein"
	exit 1
fi

seqFile=$1
seqLen=`tail -1 $seqFile | wc -c | cut -f1 -d' ' `
neededRAM=`expr 367729730 + 93485 \* $seqLen + 3990 \* $seqLen \* $seqLen + 500000000`
echo $neededRAM
