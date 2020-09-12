#!/bin/sh

if [ $# -lt 4 ]; then
	echo $0 proteinListFile SeqDir Folder2PredictedPairInfo Folder2PredictedPhiPsi [savefolder] [numModels]
	echo "	this script folds proteins in the list file using the predicted pairwise and property information"
	echo "	SeqDir: the folder for sequence file, which shall be in FASTA format and file name ends with .fasta"
	echo "	Folder2PredictedPairInfo: folder for predicted pairwise information, which shall be in PKL format and the file name ends with .predictedDistMatrix.pkl"
	echo "	Folder2PredictedPhiPsi: folder for predicted Phi/Psi, which shall be in PKL format and the file name ends with .predictedProperties.pkl"
	echo "	numModels: the number of 3D models generated for each protein in the list, default 10"
	exit 1
fi

proteins=$1
seqDir=$2
pairFolder=$3
PhiPsiFolder=$4

savefolder=`pwd`
if [ $# -ge 5 ]; then
	savefolder=$5
	if [ ! -d $savefolder ]; then
		mkdir -p $savefolder
	fi
fi

numModels=10
if [ $# -ge 6 ]; then
	numModels=$6
fi

cmdDir=$DistanceFoldingHome/Scripts4Rosetta/
for i in `cat $proteins`
do
	$cmdDir/GenPotentialNFold.sh $seqDir/${i}.fasta $pairFolder/${i}.predictedDistMatrix.pkl ${PhiPsiFolder}/${i}.predictedProperties.pkl $savefolder $numModels
done

