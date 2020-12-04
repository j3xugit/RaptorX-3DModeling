#!/bin/bash

if [[ -z "${DL4DistancePredHome}" ]]; then
        echo "ERROR: Please set environmental variable DL4DistancePredHome to the installation folder of DL4DistancePrediction4 "
        exit 1
fi

if [ $# -lt 2 ]; then
	echo $0 proteinListFile inputFolder
	echo "	This script reads input features for a list of proteins and save each into an file inputFeatures.pkl"
	echo "	The features are read one-by-one protein in serial mode"
	echo "	proteinListFile: the file containing a list of proteins"
	echo "	inputFolder shall contain a list of subfolders with name feat_proteinName_contact "
	exit 1
fi

proteinFile=$1
inputFolder=$2

ResDir=`basename $inputFolder`_PKL
if [ ! -d $ResDir ]; then
	mkdir -p $ResDir
fi

program=$DL4DistancePredHome/ReadSingleInputFeature.py

for i in `cat $proteinFile`
do
	python $program $i $inputFolder/feat_${i}_contact/
	mv $i.inputFeatures.pkl ${ResDir}/
done
