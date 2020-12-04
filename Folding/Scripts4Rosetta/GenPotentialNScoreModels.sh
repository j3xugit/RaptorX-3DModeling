#!/bin/bash

if [ $# -lt 4 ]; then
	echo $0 targetName inFolder predictedMatrixFile predictedProperty [savefolder]
	echo "	this script generates Rosetta constraints from predicted inter-residue relationship and PhiPsi angles and then score the 3D models in inFolder"
	echo "	inFolder: the folder for a set of PDB files"
	echo "	predictedMatrixFile: a PKL file (ending with .predictedDistMatrix.pkl) for all predicted inter-residue information including distance and orientation"
	echo "	predictedProperty: a PKL file (ending with .predictedPropertyFile) for predicted Phi Psi angles."
	echo "		If this entry is a special string cst, then predictedMatrixFile is interpreted as a Rosetta constraint file instead of a PKL file"
	echo "	savefolder: the folder for result saving, default current work directory"
	exit 1
fi

if [[ -z "${DistanceFoldingHome}" ]]; then
        echo "ERROR: Please set the environmental variable DistanceFoldingHome, e.g., $HOME/3DModeling/Folding"
        exit 1
fi

target=$1
inFolder=$2
pairMatrixFile=$3
propertyFile=$4

savefolder=`pwd`
if [ $# -ge 5 ]; then
	savefolder=$5
	mkdir -p $savefolder
fi

if [ ! -f $pairMatrixFile ]; then
	echo "ERROR: the file for predicted distance/orientation info does not exist: $pairMatrixFile"
	exit 1
fi

if [ "$propertyFile" == "cst" ]; then
        ## treat pairMatrixFile as a constraint file instead of a PKL file 
        cstfile=$pairMatrixFile
else
	if [ ! -f $propertyFile ]; then
		echo "ERROR: the file for predicted property info does not exist: $propertyFile"
		exit 1
	fi
        ## create a temporary cstfolder
        cstfolder=$(mktemp -d -t tmpCstDir4${target}-XXXXXXXXXX)
        cstfolder=`readlink -f $cstfolder`
        $DistanceFoldingHome/Scripts4Rosetta/GenRosettaPotential.sh -d $cstfolder $pairMatrixFile $propertyFile
        cstfile=$cstfolder/$target.pairPotential4Rosetta.SPLINE.txt
fi

program=$DistanceFoldingHome/Scripts4Rosetta/ScoreModels.py

savefile=$target.score
python $program $inFolder $cstfile $savefile

if [ "$propertyFile" != "cst" ]; then
        rm -rf $cstfolder
fi

date
