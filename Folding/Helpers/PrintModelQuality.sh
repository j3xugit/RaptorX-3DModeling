#!/bin/bash

if [ $# -lt 1 ]; then
	echo $0 ResultDir [MyDMListFile]	
	echo "	This script summarizes the quality score in a set of .sorted files; each .sorted file contains a list of quality for the 3D models of one protein"
	echo "	ResultDir: the meta folder containg a list of subfolders, each shall contain a set of 3D models and a .sorted file"
	echo "	MyDMListFile: a list of protein domains defined by users"
	echo "		when it is not provided, this script will check out all the subfolders *-*Results in ResultDir/"
	echo "		Otherwise, it will check out only those subfolders with name like target*Results where target is a protein name in MyDMListFile"
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
		for j in $i/*-quality.txt.sorted
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
		for k in $j/*-quality.txt.sorted
		do 
			if [ -s $k ]; then
				python $cmdDir/FindBest.py $k
			fi
		done
	done
done

