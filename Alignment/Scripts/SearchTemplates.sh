#!/bin/bash

if [ $# -lt 1 ]; then
        echo $0 "query [resultDir] where query shall ends with .hhm or .a3m"
        exit 1
fi

query=$1

ResultDir=`pwd`
if [ $# -gt 1 ]; then
        ResultDir=$2
fi

if [ ! -d $ResultDir ]; then
	mkdir -p $ResultDir
fi

seqName=`basename $query .hhm `

PDB70DB=/mnt/data/RaptorXCommon/HHblits/DB/PDB70_29Apr2020/pdb70
hhsearch -i $query -d $PDB70DB -mact 0.1 -o $ResultDir/$seqName.pdb70.hhr -z 500 -b 500
