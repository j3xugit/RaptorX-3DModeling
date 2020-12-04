#!/bin/bash

if [ $# -lt 3 ]; then
	echo $0 proteinList inDir outDir
	exit 1
fi

proteinList=$1
inDir=$2
outDir=$3

program=$ModelingHome/Common/PostProcessDistMatrix.py

for i in `cat $proteinList`
do
	for f in $inDir/${i}.*pkl
	do
		python $program $f $outDir
	done
done


