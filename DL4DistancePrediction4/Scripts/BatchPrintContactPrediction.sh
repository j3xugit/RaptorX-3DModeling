#!/bin/bash

if [ $# -lt 2 ]; then
	echo "$0 proteinListFile inputFolder [ResDir]"
	echo "	This script prints out predicted contact matrices from predicted distance/orientation file for a list of proteins"
	echo "	proteinListFile: a file for a list of protein names, each in a row"
	echo "	inputFolder: a folder containing predicted distance/orientation information for proteins in the list file"
	echo "	ResDir: the folder for result saving, default current work directory"
	echo "	For each protein, two text files will be generated and saved to ResDir: XXX.CASP.rr and XXX.CM.txt where XXX is the proteinName"
	exit 1
fi

proteins=`cat $1`
inputfolder=$2
if [ ! -d $inputfolder ]; then
	echo "ERROR: invalid folder for predicted distance information: $inputfolder"
	exit 1
fi

ResDir=`pwd`
if [ $# -ge 3 ]; then
	ResDir=$3
	if [ ! -d $ResDir ]; then
		mkdir -p $ResDir
	fi
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

program=$cmdDir/PrintContactPrediction.sh

for target in $proteins
do
	predDistMatrix=$inputfolder/${target}.predictedDistMatrix.pkl
	if [ ! -f $predDistMatrix ]; then
		echo "WARNING: cannot find the predicted dist matrix file for $target: $predDistMatrix "
		continue
	fi
	$program $predDistMatrix $ResDir
done
