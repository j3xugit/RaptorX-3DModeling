#!/bin/bash

if [[ -z "${DL4DistancePredHome}" ]]; then
	echo "ERROR: please set environmental variable DL4DistancePredHome to installation folder of DL4DistancePrediction4"
	exit 1
fi

function Usage 
{
	echo $0 "predFolder truthFile [targetName]"
	echo "	This script evaluates all predicted contacts for one protein"
	echo "	predFolder: the folder for all predicted distance files of one protein, ending with .predicteDistMatrix.pkl"
	echo "	truthFile: the ground truth file ending with .native.pkl"
}

if [ $# -lt 2 ]; then
	Usage
	exit 1
fi

predFolder=$1
if [ ! -d $predFolder ]; then
	echo "ERROR: invalid folder for predicted contact matrix $predFolder"
	exit 1
fi

truth=$2
if [ ! -f $truth ]; then
	echo "ERROR: invalid ground truth file $truth"
	exit 1
fi

if [ $# -ge 3 ]; then
	target=$3
else
	target=`basename $truth .native.pkl`
fi

program=$DL4DistancePredHome/EvaluateContactAccuracy.py

for predFile in $predFolder/${target}*.predictedDistMatrix.pkl
do
	if [ ! -f $predFile ]; then
		continue
	fi
	echo contact pred accuracy of $predFile
	python $program $predFile $truth $target
done
