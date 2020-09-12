#!/bin/sh

if [ $# -lt 1 ]; then
	echo $0 FoldingResultsDir [MyDMListFile]
	exit 1
fi

resultDir=$1

listFile=''
if [ $# -ge 2 ]; then
	listFile=$2
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

if [ -z "$listFile" ]; then
	for i in $resultDir/*Results
	do
		for j in $i/T*quality*sorted
		do
			if [ -s $j ]; then
				python $cmdDir/FindBest.py $j
			fi
		done
	done

	exit
fi

for i in `cat $listFile`
do
	for j in $resultDir/${i}*Results
	do
		for k in $j/T*quality*sorted
		do 
			if [ -s $k ]; then
				python $cmdDir/FindBest.py $k
			fi
		done
	done
done

