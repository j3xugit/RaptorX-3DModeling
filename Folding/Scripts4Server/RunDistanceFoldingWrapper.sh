#!/bin/bash

if [ $# -lt 5 ]
then
        echo "$0 resultLocation seqFile No3Dmodel A3Minput email "
        echo "       resultLocation shall be something like RaptorX@raptorx:resultDir, saving the input and output files "
        echo "       seqFile looks like XXXX.seq or XXX.fasta and is available at resultLocation "
	echo "	     No3Dmodel=1 indicates that only do contact/distance prediction, but not folding "
	echo "	     A3Minput=1 indicates that seqFile actually is a single multiple sequence alignment with the first sequence being the query protein"
        exit 1
fi

#export ModelingHome=/mnt/data/RaptorXCommon/RaptorX-3DModeling/
#export DistFeatureHome=$ModelingHome/BuildFeatures/
#export DL4DistancePredHome=$ModelingHome/DL4DistancePrediction2/
#export DL4PropertyPredHome=$ModelingHome/DL4PropertyPrediction/
#export DistanceFoldingHome=$ModelingHome/Folding/


if [ -z "${DistanceFoldingHome}" ]; then
	echo "Please set the environmental variable DistanceFoldingHome to the installation folder of the Folding module, e.g., $HOME/3DModeling/Folding"
	exit 1
fi

resultLocation=$1
seqFile=$2
No3Dmodel=$3
A3Minput=$4
email=$5

## automatic selection of a GPU
GPU=-1

Server=`echo $resultLocation | cut -f1 -d':' | cut -f2 -d'@' | cut -f1 -d'.' `

fullName=`basename $seqFile `
target=${fullName%.*}

##We do all contact/distance prediction jobs under Work4ContactServer
WorkDir=Work4Server-${Server}
mkdir -p $WorkDir
currDir=`pwd`
cd $WorkDir

##fetch the sequence file by scp
scp ${resultLocation}/${seqFile} ./$target.fasta
rcode=$?
if [ $rcode -ne 0 ]
then
        echo "Failed to scp the query sequence file from remote location ${resultLocation}"
        exit 1
fi

## clean up old files if they exist due to crash
if [ -d ${target}_OUT ]; then
	rm -rf ${target}_OUT/
fi

$DistanceFoldingHome/scripts/RaptorXFolder.sh $target.fasta ./ $No3Dmodel $A3Minput $GPU
rcode=$?
if [ $rcode -ne 0 ]
then
        echo "Failed to run contact/distance prediction and ab initio folding for $target"
        exit 1
fi

## collect results and check if required files have been generated or not
## the result files are saved to a new folder ${target}_Results
$DistanceFoldingHome/scripts/CollectResults4Server.sh $target ${target}_OUT

##copy the results back to the resultLocation
chmod -R a+rX ${target}_Results/
scp -r ${target}_Results/* ${resultLocation}/
rcode=$?
if [ $rcode -ne 0 ]
then
        echo "Failed to scp the prediction results back to ${resultLocation}"
        exit 1
fi

##clean up
rm -rf ${target}_OUT ${target}_Results ${target}.fasta

cd $currDir
