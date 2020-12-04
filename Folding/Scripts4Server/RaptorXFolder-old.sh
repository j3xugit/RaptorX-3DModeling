#!/bin/bash

if [ $# -lt 4 ]; then
	echo $0 "seqFile outDir No3Dmodel A3Minput [gpu (default -1)] "
	echo "   outDir is the result direcotry for all the results, e.g., ./ for current work directory"
	echo "   if gpu is set to -1, then automatically select GPU with the maximum amount of free memory"
	exit 1
fi

seqFile=$1
target=`basename $seqFile .fasta`
outDir=$2

No3DModel=$3
A3Minput=$4

GPU=-1
if [ $# -ge 5 ]; then
	GPU=$5
fi

if [ -z "${DistFeatureHome}" ]; then
        echo "Please set the environmental variable DistFeatureHome to the installation folder of the BuildFeatures module, e.g., $HOME/3DModleing/BuildFeatures"
        exit 1
fi

if [ -z "${DL4DistancePredHome}" ]; then
        echo "Please set the environmental variable DL4DistancePredHome to the installation folder of the DL4DistancePrediction module, e.g., $HOME/3DModleing/DL4DistancePrediction2"
        exit 1
fi

if [ -z "${DL4PropertyPredHome}" ]; then
        echo "Please set the environmental variable DL4PropertyPredHome to the installation folder of the DL4PropertyPrediction module, e.g., $HOME/3DModleing/DL4PropertyPrediction"
        exit 1
fi

if [ -z "${DistanceFoldingHome}" ]; then
        echo "Please set the environmental variable DistanceFoldingHome to the installation folder of the DL4PropertyPrediction module, e.g., $HOME/3DModleing/Folding"
        exit 1
fi

## buildFeatures always try to use a GPU with the maximum amount of free memory
$DistFeatureHome/buildFeatures.sh $seqFile $outDir $A3Minput $GPU
if [ $? -ne 0 ]; then
	echo "Failed in generating feature files for $target"
	exit 1
fi

$DL4PropertyPredHome/Utils/PredictProperty4Server.sh $target $outDir/${target}_OUT/ $GPU
if [ $? -ne 0 ]; then
	echo "Failed in predicting property for $target"
	exit 1
fi

$DL4DistancePredHome/Utils/PredictDistance4Server.sh $target $outDir/${target}_OUT/ $GPU
if [ $? -ne 0 ]; then
	echo "Failed in predicting distance information for $target"
	exit 1
fi

if [ $No3DModel -eq 0 ]; then 
	## using 10 to generate 200 models and 20 to generate 400 models
	$DistanceFoldingHome/scripts/FoldProtein4Server.sh $target $outDir/${target}_OUT/ 10
	if [ $? -ne 0 ]; then
		echo "Failed in folding $target"
		exit 1
	fi

fi

echo "Finished running contact/distance prediction (and ab initio folding) for $target. Please check out results in $outDir/${target}_OUT/ "
