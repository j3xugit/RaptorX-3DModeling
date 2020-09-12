#!/bin/sh

GPU=cuda0

if [ $# -lt 2 ]; then
	echo "RunDistancePred4OneProtein.sh ModelName InputFeature_PKL [gpu] "
	echo "   See DistanceModelFiles.txt for ModelName "
	echo "   InputFeature is a PKL file, which has a list of features"
	echo "   gpu can be cuda0 to cuda3, default is cuda0"
	exit 1
fi

ModelName=$1
predFile=$2
proteinName=`basename $predFile .distanceFeatures.pkl`
ResultDir=`dirname $predFile`

if [ $# -ge 3 ]; then
	GPU=$3
fi

echo ModelName=$ModelName
echo GPU=$GPU

if [[ -z "${DL4DistancePredHome}" ]]; then
	echo "Please set the environmental variable DL4DistancePredHome to the installation dir of the  DL4DistancePrediction2 module."
	exit 1
fi


if [ ! -f $program ]; then
	echo ERROR: the python code for distance prediction $program does not exist
	exit 1
fi

##correctly we use hard code here. Need to change them later.
ModelPoolRootDir=/mnt/data/RaptorXCommon/TrainTestData/Distance_Contact_TrainData/Jinbo_Folder/
ModelPoolFile=$ModelPoolRootDir/result4HBBeta/DistanceModelFiles.txt

## load model file names
. $ModelPoolFile

if [ $predFile == "" ]; then
	echo ERROR: predFile is empty
	exit 1
fi

if [ ! -f $predFile ]; then
	echo ERROR: input feature file does not exist: $predFile
	exit 1
fi

ModelFiles=`eval echo '$'${ModelName}`
#echo ModelFiles=$ModelFiles

if [ $ModelFiles == "" ]; then
	echo ERROR: ModelFiles is empty
	exit 1
fi

ResultDir=$ResultDir/Dist_Server_${ModelName}/
#echo ResultDir=$ResultDir
mkdir -p $ResultDir

program=$DL4DistancePredHome/RunDistancePredictor2.py

CUDNNDIR=/usr/local/
if [ -d $HOME/CUDNN ]; then
	CUDNNDIR=$HOME/CUDNN
fi

THEANO_FLAGS=blas.ldflags=,device=$GPU,floatX=float32,dnn.include_path=$CUDNNDIR/cuda/include,dnn.library_path=$CUDNNDIR/cuda/lib64 python $program  -p $predFile -m $ModelFiles -d $ResultDir 
#THEANO_FLAGS=blas.ldflags=,device=$GPU,floatX=float32,dnn.conv.algo_fwd=none python $program -p $predFile -m $ModelFiles -d $ResultDir 

if [ $? -ne 0 ]; then
	echo "Failed to predict distance info using $ModelName from $predFile"
	exit 1
fi

currDir=`pwd`
cd ${ResultDir}

python $DL4DistancePredHome/EstimateAtomDistBounds.py -i ${proteinName}.predictedDistMatrix.pkl -c txt -b -p
ret1=$?

python $DL4DistancePredHome/MergeDistanceBins.py ${proteinName}.correctedDistMatrix.pkl 12C
ret2=$?

if [ $ret1 -ne 0 -o $ret2 -ne 0 ]; then
	echo "Failed to postprocess predicted distance from $predFile"
	exit 1
fi

## visualize the predicted contact map
python $DL4DistancePredHome/Scripts/PlotContactMapByMatrix.py ${proteinName}.gcnn 
if [ $? -ne 0 ]; then
	echo "Failed to visualize the predicted contact matrix for ${proteinName}. "
	exit 1
fi

## visualize the distance matrix
python $DL4DistancePredHome/Scripts/PlotDistanceMatrix.py ${proteinName}.bound.txt
if [ $? -ne 0 ]; then
	echo "Failed to visualize the predicted distance matrix for ${proteinName}. "
	exit 1
fi

cd $currDir
