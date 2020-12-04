#!/bin/bash


if [[ -z "${DL4PropertyPredHome}" ]]; then
	echo "ERROR: please set environmental variable DL4PropertyPredHome, e.g., $HOME/3DModeling/DL4PropertyPrediction"
	exit 1
fi

if [[ -z "${CUDA_ROOT}" ]]; then
        echo "ERROR: please set environmental variable CUDA_ROOT"
        exit 1
fi

DeepModelFile=$DL4PropertyPredHome/params/ModelFile4PropertyPred.txt
if [ ! -f $DeepModelFile ]; then
        echo "ERROR: cannot find the file for deep models: $DeepModelFile"
        exit 1
fi

if [ $# -lt 3 ]; then
        echo "$0 ModelName inputPKLfile ResultDir [gpu] "
	echo "	This script predicts properties from one input PKL file"
	echo "	ModelName: a name for a set of deep learning models, e.g., PhiPsiSet10820Models. See $DeepModelFile for more names"
	echo "	inputPKLfile: an input feautre file in PKL format"
	echo "	ResultDir: the folder for result saving"
        echo "	gpu: cuda0 to cuda3, default cuda0"
        exit 1
fi

ModelName=$1

## the input feature for prediction
predFile=$2
if [ ! -f $predFile ]; then
        echo "ERROR: the input PKL file does not exist: $predFile"
        exit 1
fi
DataName=`basename $predFile | cut -f1 -d'.'`

#ResultDir=LocalResults/Property_${DataName}_${ModelName}/
ResultDir=$3
if [ ! -d $ResultDir ]; then
	mkdir -p $ResultDir
fi

GPU=cuda0
if [ $# -ge 4 ]; then
        GPU=$4
fi

program=$DL4PropertyPredHome/RunPropertyPredictor.py
if [ ! -f $program ]; then
        echo "ERROR: the main program $program does not exist"
        exit 1
fi

. $DeepModelFile

ModelFiles=`eval echo '$'${ModelName}`
if [ $ModelFiles == "" ]; then
        echo ERROR: ModelFiles is empty
        exit 1
fi


if [[ $nativeDir == "" || ! -d "$nativeDir" ]]; then
	#echo "The folder for native properties does not exist"
	THEANO_FLAGS=blas.ldflags=,device=$GPU,floatX=float32,dnn.include_path=${CUDA_ROOT}/include,dnn.library_path=${CUDA_ROOT}/lib64 python $program -p $predFile -m $ModelFiles -d $ResultDir 
	ret=$?
else
	THEANO_FLAGS=blas.ldflags=,device=$GPU,floatX=float32,dnn.include_path=${CUDA_ROOT}/include,dnn.library_path=${CUDA_ROOT}/lib64 python $program -p $predFile -m $ModelFiles -g $nativeDir -d $ResultDir 
	ret=$?
fi

if [ $ret -ne 0 ]; then
	echo "ERROR: failed to predict properties for $predFile using $ModelName"
	exit 1
fi
