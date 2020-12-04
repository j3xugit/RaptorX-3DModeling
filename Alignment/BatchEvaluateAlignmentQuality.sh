#!/bin/bash

if [ $# -lt 4 ]; then
	echo "$0 alignmentListFile alignDir queryPDBDir templatePDBDir [ResultDir]"
	echo "	Each line in alignmentListFile is a string with format queryName-templateName"
	exit 1
fi

alnList=$1
alnDir=$2
queryDir=$3
tplDir=$4

ResultDir=`pwd`
if [ $# -ge 5 ]; then
	ResultDir=$5
	mkdir -p $ResultDir
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`
program=$cmdDir/EvaluateAlignmentQuality.py

for i in `cat $alnList`
do
	queryName=`echo $i | cut -f1 -d'-'`
	tplName=`echo $i | cut -f2 -d'-'`
	alnFile=$alnDir/${tplName}-${queryName}.fasta
	queryPDB=$queryDir/$queryName.pdb
	tplPDB=$tplDir/$tplName.pdb
	python $program $alnFile $queryPDB $tplPDB $ResultDir 
done
