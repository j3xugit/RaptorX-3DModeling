#!/bin/sh

GPU=cuda0

if [ $# -lt 2 ]; then
	echo "RunDistPred4OneInputFile.sh ModelName InputFeature_PKL [gpu] "
	echo "   See DistanceModelFiles.txt for ModelName "
	echo "   InputFeature is a PKL file, containing a list of proteins and their features"
	echo "   gpu can be cuda0 to cuda3, default is cuda0"
	exit 
fi

ModelName=$1
predFile=$2
DataName=`basename $predFile .pkl | cut -f1 -d'.'`

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

## load model file names
. result4HBBeta/DistanceV3ModelFiles.txt

if [ $predFile == "" ]; then
	echo ERROR: predFile is empty
	exit -1
fi

if [ ! -f $predFile ]; then
	echo ERROR: input feature file does not exist: $predFile
	exit -1
fi

ModelFiles=`eval echo '$'${ModelName}`
#echo ModelFiles=$ModelFiles

if [ $ModelFiles == "" ]; then
	echo ERROR: ModelFiles is empty
	exit -1
fi

ResultDir=LocalResults/Dist_${DataName}_${ModelName}/
#echo ResultDir=$ResultDir
mkdir -p $ResultDir

#THEANO_FLAGS=blas.ldflags=,device=$GPU,floatX=float32,dnn.include_path=$HOME/CUDNN/cuda/include,dnn.library_path=$HOME/CUDNN/cuda/lib64 python $program  -p $predFile -m $ModelFiles -d $ResultDir 
THEANO_FLAGS=blas.ldflags=,device=$GPU,floatX=float32 python $program  -p $predFile -m $ModelFiles -d $ResultDir 
