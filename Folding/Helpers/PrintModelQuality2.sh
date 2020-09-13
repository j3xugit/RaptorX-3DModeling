#!/bin/sh

if [ $# -lt 1 ]; then
	echo $0 QualityResultsDir [domainListFile]
	echo "  This script prints out the best quality in a set of -quality.txt.sorted files; each contains a list of quality for the 3D models of one protein"
        echo "  when MyDMListFile is not provided, this script will check out all the -quality.txt.sorted files in ResultsDir/"
        echo "  Otherwise, it will check out only those files with name like target*-quality.txt.sorted where target is a protein name in MyDMListFile"

	exit 1
fi

resultDir=$1
if [ ! -d $resultDir ]; then
	echo "ERROR: invalid folder for the -quality.txt.sorted files $resultDir"
	exit 1
fi

listFile=''
if [ $# -ge 2 ]; then
	listFile=$2
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

if [ -z "$listFile" ]; then
	for j in $resultDir/*-quality.txt.sorted
        do
                if [ -s $j ]; then
                        python $cmdDir/FindBest.py $j
                fi
        done
	exit
fi

for i in `cat $listFile`
do
	for j in $resultDir/${i}*-quality.txt.sorted
	do
		if [ -s $j ]; then
			python $cmdDir/FindBest.py $j
		fi
	done
done
