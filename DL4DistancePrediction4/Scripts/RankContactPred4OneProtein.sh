#!/bin/bash

if [[ -z "${DL4DistancePredHome}" ]]; then
	echo "ERROR: please set environmental variable DL4DistancePredHome to installation folder of DL4DistancePrediction4, e.g., RaptorX-3DModeling/DL4DistancePrediction4/"
	exit 1
fi

function Usage 
{
	echo $0 proteinName predFolder
	echo "	This script ranks predictions of one protein by top L long-range contact probability sum"
	echo "	proteinName: the protein name"
	echo "	predFolder: the folder for predicted dist/orientation matrix files"
	echo "	this script prints the ranking list onto screen. The left column is the predicted dist matrix files"
}

if [ $# -lt 2 ]; then
	Usage
	exit 1
fi

protein=$1
inputFolder=$2
if [ ! -d $inputFolder ]; then
	echo "ERROR: invalid folder for predicted dist/ori matrix files: $inputFolder"
	exit 1
fi

program=$DL4DistancePredHome/RankPredictionsByTopProb.py
if [ ! -f $program ]; then
	echo "ERROR: unavailable script program $program"
	exit 1
fi

predFiles=`ls $inputFolder/${protein}*predictedDistMatrix.pkl | grep -v ${protein}D | grep -v ${protein}d `
if [ ! -z "$predFiles" ]; then
	python $program $predFiles
fi
