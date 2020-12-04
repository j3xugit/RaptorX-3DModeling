#!/bin/bash

if [[ -z "${DL4DistancePredHome}" ]]; then
        echo "ERROR: please set environmental variable DL4DistancePredHome, e.g., $ModelingHome/DL4DistancePrediction4"
        exit 1
fi

if [ $# -lt 1 ]; then
	echo "$0 predictDistMatrixPKLfile [ResDir]"
	echo "	This script prints predicted contact matix from predicted distance/orientation matrix file for one protein"
	echo "	predictedDistMatrixPKLfile: a file for predicted distance/orientation information of one protein"
	echo "	ResDir: the folder for result saving, default current work directory"
	echo "	Two text files will be generated and saved to ResDir: XXX.CASP.rr and XXX.CM.txt where XXX is the protein name"
	exit 1
fi

inputfile=$1

if [ ! -f $inputfile ]; then
	echo "ERROR: invalid dist/ori matrix file $inputfile"
	exit 1
fi

ResDir=`pwd`
if [ $# -ge 2 ]; then
	ResDir=$2
fi
if [ ! -d $ResDir ]; then
	mkdir -p $ResDir
fi

program=$DL4DistancePredHome/PrintContactPrediction.py
python $program -d $ResDir $inputfile

target=`basename $inputfile | cut -f1 -d'.'`

currDir=`pwd`
cd $ResDir
python ${DL4DistancePredHome}/Utils/PlotContactMapByMatrix.py ${target}.CM.txt
cd $currDir
