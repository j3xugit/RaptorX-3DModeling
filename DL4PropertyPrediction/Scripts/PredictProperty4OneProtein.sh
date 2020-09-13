#!/bin/bash

if [ -z "${DL4PropertyPredHome}" ]; then
        echo "ERROR: Please set the environmental variable DL4PropertyPredHome to the install folder of DL4PropertyPrediction"
        exit 1
fi

DeepModelFile=$DL4PropertyPredHome/params/ModelFile4PropertyPred.txt
ModelName=PhiPsiSet10820Models
GPU=-1
ResultDir=`pwd`

function Usage {
        echo $0 "[ -f DeepModelFile | -m ModelName | -d ResultDir | -g gpu ] proteinName inputFolder"
        echo "	This script predicts properties for a set of proteins using a local machine with GPUs"
 	echo "	inputFolder: a folder containing files for property prediction, e.g., 1pazA_OUT/. It shall contain subfolders proteinName_contact/ and proteinName_thread/, which in turn shall contain .a3m or .hhm files"
	echo " "
        echo "	DeepModelFile: a file containing a set of deep model names, default $DeepModelFile"
        echo "	ModelName: a model name defined in DeepModelFile representing a set of deep learning models, default $ModelName"
        echo "	ResultDir: the folder for result saving, default current work directory"
	echo "	gpu: -1 (default) or 0-3. If set to -1, automatically select a GPU"
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

if [ $# -lt 2 ]; then
        Usage
        exit 1
fi

target=$1

rootDir=$2
if [ ! -d $rootDir ]; then
	echo "ERROR: invalid folder for .a3m or .hhm files: $rootDir"
	exit 1
fi

if [ ! -d $ResultDir ]; then
        mkdir -p $ResultDir
fi

inputFeature=$ResultDir/${target}.propertyFeatures.pkl
$DL4PropertyPredHome/Scripts/CollectPropertyFeatures.sh $target $rootDir $ResultDir
if [ $? -ne 0 ]; then
	echo "ERROR: failed to collect input property feature for $target"
	exit 1
fi

if [ ! -f $inputFeature ]; then
	echo "ERROR: cannot find input property feature file $inputFeature"
        exit 1
fi

if [ ! -f $DeepModelFile ]; then
        echo "ERROR: cannot find the file for deep models: $DeepModelFile"
        exit 1
fi

#echo "Running $DL4PropertyPredHome/Scripts/PredictPropertyLocal.sh -f $DeepModelFile -m $ModelName -d $ResultDir -g $GPU $inputFeature"

$DL4PropertyPredHome/Scripts/PredictPropertyLocal.sh -f $DeepModelFile -m $ModelName -d $ResultDir -g $GPU $inputFeature
if [ $? -ne 0 ]; then
	echo "ERROR: failed to predict property for $target using information in $rootDir"
	exit 1
fi

rm -f $inputFeature
