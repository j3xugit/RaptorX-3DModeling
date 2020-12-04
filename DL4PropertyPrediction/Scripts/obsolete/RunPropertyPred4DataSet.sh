#!/bin/bash

GPU=cuda0

if [[ -z "${DL4PropertyPredHome}" ]]; then
        echo "ERROR: please set environmental variable DL4PropertyPredHome, e.g., $HOME/3DModeling/DL4PropertyPrediction"
        exit 1
fi

if [[ -z "${CUDA_ROOT}" ]]; then
        echo "ERROR: please set environmental variable CUDA_ROOT"
        exit 1
fi

#DeepModelFile=$DL4PropertyPredHome/params/ModelFile4PropertyPred.txt
DeepModelFile=/mnt/data/RaptorXCommon/TrainTestData/ProteinProperty_Project/Jinbo_Folder/result4property/PropertyModelFiles.txt
if [ ! -f $DeepModelFile ]; then
        echo "ERROR: cannot find the file for deep models: $DeepModelFile"
        exit 1
fi

DataFile=/mnt/data/RaptorXCommon/TrainTestData/ProteinProperty_Project/Jinbo_Folder/data4property/PropertyDataFiles.txt
if [ ! -f $DataFile ]; then
        echo "ERROR: cannot find the data file: $DataFile"
        exit 1
fi

if [ $# -lt 2 ]; then
        echo "$0 ModelName DataName [gpu] "
	echo "	This script predicts properties for one dataset with name DataName"
	echo "	ModelName: a name for a set of deep learning models, e.g., PhiPsiSet10820Models. See $DeepModelFile for more names"
        echo "	gpu: cuda0 to cuda3, default cuda0"
	echo "	The result will be saved to LocalResults/PhiPsiSS8_DataName_ModelName/ "
        exit
fi

ModelName=$1
#echo ModelName=$ModelName

DataName=$2

if [ $# -ge 3 ]; then
        GPU=$3
fi

program=$DL4PropertyPredHome/RunPropertyPredictor.py
if [ ! -f $program ]; then
        echo ERROR: the main program $program does not exist
        exit 1
fi

## load data and model file names
. $DeepModelFile
. $DataFile

## the input feature for prediction
predFile=`eval echo '$'${DataName}`
#echo predFile=$predFile
if [ $predFile == "" ]; then
        echo ERROR: predFile is empty
        exit 1
fi

nativeDir=`eval echo '$'${DataName}Native`
##echo nativeDir=${nativeDir}

ModelFiles=`eval echo '$'${ModelName}`
#echo ModelFiles=$ModelFiles

if [ $ModelFiles == "" ]; then
        echo ERROR: ModelFiles is empty
        exit 1
fi

ResultDir=LocalResults/PhiPsiSS8_${DataName}_${ModelName}/
#echo ResultDir=$ResultDir
mkdir -p $ResultDir

if [[ $nativeDir == "" || ! -d "$nativeDir" ]]; then
	echo "The folder for native properties does not exist"
	THEANO_FLAGS=blas.ldflags=,device=$GPU,floatX=float32,dnn.include_path=${CUDA_ROOT}/include,dnn.library_path=${CUDA_ROOT}/lib64 python $program -p $predFile -m $ModelFiles -d $ResultDir 
	ret=$?
else
	THEANO_FLAGS=blas.ldflags=,device=$GPU,floatX=float32,dnn.include_path=${CUDA_ROOT}/include,dnn.library_path=${CUDA_ROOT}/lib64 python $program -p $predFile -m $ModelFiles -g $nativeDir -d $ResultDir 
	ret=$?
fi

if [ $ret -ne 0 ]; then
	echo "ERROR: failed to run property prediction for $DataName using $ModelName"
	exit 1
fi

