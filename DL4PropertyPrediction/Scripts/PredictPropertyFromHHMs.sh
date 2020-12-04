#!/bin/bash

if [[ -z "${DL4PropertyPredHome}" ]]; then
        echo "ERROR: please set environmental variable DL4PropertyPredHome to the install folder of DL4PropertyPrediction"
        exit 1
fi

if [[ -z "${CUDA_ROOT}" ]]; then
        echo "ERROR: please set environmental variable CUDA_ROOT"
        exit 1
fi

DeepModelFile=$DL4PropertyPredHome/params/ModelFile4PropertyPred.txt
ModelName=PhiPsiSet10820Models
GPU=cuda0
ResultDir=`pwd`

Usage () {
        echo $0 "[ -f DeepModelFile | -m ModelName | -d ResultDir | -g gpu ] proteinListFile hhmFolder"
	echo "	This script predicts local property for a list of proteins"
	echo "	proteinListFile: a file for a list of protein names, each in one row"
	echo "	hhmFolder: the folder for .hhm files, generated by hhmake from MSA files"
	echo " "
        echo "	DeepModelFile: a file containing a set of deep model names, default $DeepModelFile"
        echo "	ModelName: a model name (defined in DeepModelFile) representing a set of deep learning models, default $ModelName"
        echo "	ResultDir: the folder for result saving, default current work directory"
        echo "	gpu: cuda0 (default), cuda1, cuda2, cuda3"
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

if [ ! -f $DeepModelFile ]; then
        echo "ERROR: cannot find the file for deep models: $DeepModelFile"
        exit 1
fi

if [ ! -d $ResultDir ]; then
	mkdir -p $ResultDir
fi

proteinListFile=$1
hhmFolder=$2

##create a temporary directory for feature file
#tmpdir=$(mktemp -d -t tmpWork4PropertyPrediction-XXXXXXXXXX)
bname=`basename $proteinListFile .list`
bname=`basename $bname .txt`
predFile=$ResultDir/$bname.propertyFeatures.pkl

python $DL4PropertyPredHome/GenPropertyFeaturesForMultiProteins.py $proteinListFile $hhmFolder $predFile

program=$DL4PropertyPredHome/RunPropertyPredictor.py
if [ ! -f $program ]; then
        echo ERROR: the main program $program does not exist
        exit 1
fi

. $DeepModelFile

ModelFiles=`eval echo '$'${ModelName}`
if [ $ModelFiles == "" ]; then
        echo "ERROR: ModelFiles for $ModelName is empty"
        exit 1
fi

THEANO_FLAGS=blas.ldflags=,device=$GPU,floatX=float32,dnn.include_path=${CUDA_ROOT}/include,dnn.library_path=${CUDA_ROOT}/lib64 python $program -p $predFile -m $ModelFiles -d $ResultDir 
if [ $? -ne 0 ]; then
	echo "ERROR: failed to run python $program -p $predFile -m $ModelFiles -d $ResultDir on $GPU"
	exit 1
fi

rm -rf $tmpdir
