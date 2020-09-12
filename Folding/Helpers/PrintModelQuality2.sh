#!/bin/sh

if [ $# -lt 1 ]; then
	echo $0 FoldingResultsDir [domainListFile]
	exit 1
fi

resultDir=$1
if [ ! -d $resultDir ]; then
	echo "ERROR: invalid folder for quality.txt.sorted files $resultDir"
	exit 1
fi

listFile=CASP13/CASP13TBMHard.list
if [ $# -ge 2 ]; then
	listFile=$2
fi

if [ ! -f $listFile ]; then
	echo "ERROR: invalid protein list file $listFile"
	exit 1
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

for i in `cat $listFile`
do
	for j in $resultDir/${i}*quality.txt.sorted
	do
		if [ -s $j ]; then
			python $cmdDir/FindBest.py $j
		fi
	done
done
