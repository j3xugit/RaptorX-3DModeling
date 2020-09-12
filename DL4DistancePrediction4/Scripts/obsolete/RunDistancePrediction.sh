#!/bin/sh

GPU=cuda0

if [ $# -lt 2 ]; then
	echo "RunDistancePrediction.sh ModelName DataName [gpu] "
	echo "   See DistanceModelFiles.txt for ModelName and DistanceDataFiles.txt for DataName"
	echo "   gpu can be cuda0 to cuda3, default is cuda0"
	exit 
fi

ModelName=$1
DataName=$2

if [ $# -ge 3 ]; then
	GPU=$3
fi

#echo ModelName=$ModelName
#echo GPU=$GPU

## you may want to modify the installation directory of the package
export DL4DistancePredHome=$HOME/3DModeling/DL4DistancePrediction4/

program=$DL4DistancePredHome/RunDistancePredictor.py
if [ ! -f $program ]; then
	echo ERROR: the main program $program does not exist
	exit -1
fi

## load data and model file names
. result4HBBeta/DistanceV3ModelFiles.txt
. data4HBBeta/DistanceDataFiles.txt

## the input feature for prediction
predFile=`eval echo '$'${DataName}`
#echo predFile=$predFile
if [ $predFile == "" ]; then
	echo ERROR: predFiles is empty
	exit -1
fi

NativeDir=`eval echo '$'${DataName}Native`
#echo NativeDir=$NativeDir

ModelFiles=`eval echo '$'${ModelName}`
#echo ModelFiles=$ModelFiles

if [ $ModelFiles == "" ]; then
	echo ERROR: ModelFiles is empty
	exit -1
fi

ResultDir=LocalResults/Dist_${DataName}_${ModelName}/
#echo ResultDir=$ResultDir
mkdir -p $ResultDir

if [[ $NativeDir != "" && -d "$NativeDir" ]]; then
	THEANO_FLAGS=blas.ldflags=,device=$GPU,floatX=float32 python $program  -p $predFile -m $ModelFiles -g $NativeDir -d $ResultDir 
else
	THEANO_FLAGS=blas.ldflags=,device=$GPU,floatX=float32 python $program  -p $predFile -m $ModelFiles -d $ResultDir 
fi
