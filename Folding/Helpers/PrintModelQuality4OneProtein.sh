#!/bin/bash

if [ $# -lt 2 ]; then
	echo $0 FoldingResultsDir MyDMName
	exit 1
fi

resultDir=$1
target=$2

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

for i in $resultDir/${target}*Results
do
	for j in $i/T*quality*sorted
	do
		if [ -s $j ]; then
			python $cmdDir/FindBest.py $j
		fi
	done
done
