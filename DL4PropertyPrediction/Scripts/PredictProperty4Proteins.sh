#!/bin/bash

if [[ -z "${DL4PropertyPredHome}" ]]; then
        echo "ERROR: please set environmental variable DL4PropertyPredHome to the installation folder of DL4PropertyPrediction"
        exit 1
fi

if [[ -z "${CUDA_ROOT}" ]]; then
        echo "ERROR: please set environmental variable CUDA_ROOT"
        exit 1
fi

DeepModelFile=$DL4PropertyPredHome/params/ModelFile4PropertyPred.txt
#DeepModelFile=/mnt/data/RaptorXCommon/TrainTestData/ProteinProperty_Project/Jinbo_Folder/result4property/PropertyModelFiles.txt
ModelName=PhiPsiSet10820Models
GPU=cuda0
ResultDir=`pwd`

function Usage {
        echo $0 "[ -f DeepModelFile | -m ModelName | -d ResultDir | -g gpu ] proteinListFile inputFolder"
	echo "	This script predicts properties for a set of proteins using a local machine"
	echo "	proteinListFile: a file for a list of proteins, each in one row"
	echo "	inputFolder: a meta folder containing a list of subfolders, each subfolder has name XXX_OUT where XXX is the protein name in proteinListFile"
	echo " "
	echo "	DeepModelFile: a file containing a set of deep model names, default $DeepModelFile"
	echo "	ModelName: a model name defined in DeepModelFile representing a set of deep learning models, default $ModelName"
	echo "	ResultDir: the folder for result saving, default current work directory"
        echo "	gpu: cuda0, cuda1, cuda2 and cuda3, default $GPU"
}

while getopts ":f:m:d:g:" opt; do
	case ${opt} in
		f )
		  DeepModelFile=$OPTARG
		  ;;
		m )
		  ModelName=$OPTARG
		  ;;
		d )
		  ResultDir=$OPTARG
		  ;;
		g )
		  GPU=$OPTARG
		  ;;
  		\? )
                  echo "Invalid Option: -$OPTARG" 1>&2
                  exit 1
                  ;;
                : )
                  echo "Invalid Option: -$OPTARG requires an argument" 1>&2
                  exit 1
                  ;;
        esac
done
shift $((OPTIND -1))

if [ $# -ne 2 ]; then
        Usage
        exit 1
fi


proteinListFile=$1
if [ ! -f $proteinListFile ]; then
	echo "ERROR: invalid protein list file $proteinListFile"
	exit 1
fi

inputFolder=$2
if [ ! -d $inputFolder ]; then
	echo "ERROR: invalid folder for input feature files: $inputFolder"
	exit 1
fi

if [ ! -d $ResultDir ]; then
	mkdir -p $ResultDir
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`
## collect features
${cmdDir}/BatchCollectPropertyFeatures.sh $proteinListFile $inputFolder $ResultDir
if [ $? -ne 0 ]; then
	echo "ERROR: failed to run ${cmdDir}/BatchCollectPropertyFeatures.sh $proteinListFile $inputFolder $ResultDir"
	exit 1
fi

program=$DL4PropertyPredHome/RunBatchPropertyPredictor.py
if [ ! -f $program ]; then
        echo ERROR: the main program $program does not exist
        exit 1
fi

## load data and model file names
if [ ! -f $DeepModelFile ]; then
        echo "ERROR: cannot find the file for deep models: $DeepModelFile"
        exit 1
fi
. $DeepModelFile

ModelFiles=`eval echo '$'${ModelName}`
#echo ModelFiles=$ModelFiles
if [ $ModelFiles == "" ]; then
        echo "ERROR: ModelFiles for $ModelName is empty"
        exit 1
fi

#THEANO_FLAGS=blas.ldflags=,device=$GPU,floatX=float32,dnn.include_path=${CUDA_ROOT}/include,dnn.library_path=${CUDA_ROOT}/lib64 python $program -p $proteinListFile -i $inputFolder -m $ModelFiles -d $ResultDir 
THEANO_FLAGS=blas.ldflags=,device=$GPU,floatX=float32,dnn.include_path=${CUDA_ROOT}/include,dnn.library_path=${CUDA_ROOT}/lib64 python $program -p $proteinListFile -i $ResultDir -m $ModelFiles -d $ResultDir 
if [ $? -ne 0 ]; then
	echo "ERROR: failed to run property prediction for $proteinListFile with $inputFolder using $ModelName"
	exit 1
fi
