#!/bin/bash

## this script generates a tgt file from an .hhm or .a3m file
## the tgt file is named after proteinName.tgt.pkl, which is supposed to be compatible with the old TGT file

if [ $# -lt 1 ]; then
	echo $0 MSAfile [ResultDir]
	echo "	This script generates a TGT file from an .hhm or .a3m file"
	echo "	The resultant file is named after proteinName.tgt.pkl and saved in ResultDir (if provided) or current work directory"
	exit 1
fi

if [[ -z "$DL4PropertyPredHome" ]]; then
	echo "ERROR: please set environmental variable DL4PropertyPredHome to the instllation directory of DL4PropertyPrediction, e.g., $$HOME/3DModeling/DL4PropertyPrediction"
	exit 1
fi

if [[ -z "$DistFeatureHome" ]]; then
	echo "ERROR: please set environmental variable DistFeatureHome to the instllation directory of DL4PropertyPrediction, e.g., $$HOME/3DModeling/BuildFeatures"
	exit 1
fi

MSAfile=$1

ResultDir=`pwd`
if [ $# -ge 2 ]; then
	ResultDir=$2
	mkdir -p $ResultDir
fi

fulnam=`basename $MSAfile`
target=${fulnam%.*}

hhmfile=${ResultDir}/${target}.hhm

if [[ "$fulnam" == *.hhm ]]; then
	cp $MSAfile $hhmfile
else
	$DistFeatureHome/util/hhmake -i $MSAfile -o $hhmfile
fi

## convert hhmfile to a feature file for property prediction
python $DL4PropertyPredHome/GenPropertyFeaturesFromMultiHHMs.py $target $hhmfile $ResultDir

if [ ! -f $ResultDir/$target.propertyFeatures.pkl ]; then
	echo "ERROR: failed to generate $ResultDir/$target.propertyFeatures.pkl "
	exit 1
fi

## predict properties
#echo $target $ResultDir/$target.propertyFeatures.pkl $ResultDir
$DL4PropertyPredHome/Scripts/PredictProperty4OneInputPKL.sh $target $ResultDir/$target.propertyFeatures.pkl $ResultDir

if [ ! -f $ResultDir/$target.predictedProperties.pkl ]; then
	echo "ERROR: failed to generate $ResultDir/$target.predictedProperties.pkl "
	exit 1
fi

cmd=`readlink -f $0 `
cmdDir=`dirname $cmd`

## generate a tgt.pkl file from both hhmfile and predictedProperties.pkl
python $cmdDir/HHM2TGT.py $hhmfile $ResultDir/$target.predictedProperties.pkl $ResultDir

if [ ! -f $ResultDir/$target.tgt.pkl ]; then
	echo "ERROR: failed to generate $ResultDir/$target.tgt.pkl "
	exit 1
fi

