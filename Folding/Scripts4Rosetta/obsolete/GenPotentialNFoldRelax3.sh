#!/bin/sh

if [ $# -lt 3 ]; then
	echo "$0 inFile predictedPairInfo predictedPropertyInfo [savefolder [numModels]]"
	echo "	this script generates Rosetta constraints from predicted inter-residue relationship and PhiPsi angles and then fold the protein"
	echo "	inFile: the primary sequence file in FASTA format"
	echo "	predictedPairInfo: a PKL file for all predicted inter-residue information including distance and orientation"
	echo "	predictedPropertyInfo: a PKL file for predicted Phi Psi angles"
	echo "	savefolder: the folder for result saving, default current work directory"
	echi "	numModels: the number of models to be generated, default 1"
	exit 1
fi

if [[ -z "${DistanceFoldingHome}" ]]; then
        echo "ERROR: Please set the environmental variable DistanceFoldingHome, e.g., $HOME/3DModeling/Folding"
        exit 1
fi

inFile=$1
target=`basename $inFile `
target=`echo $target | cut -f1 -d'.'`

pairMatrixFile=$2
propertyFile=$3

savefolder=`pwd`
if [ $# -ge 4 ]; then
	savefolder=$4
	mkdir -p $savefolder
fi

numModels=1
if [ $# -ge 5 ]; then
	numModels=$5
fi

if [ ! -f $pairMatrixFile ]; then
	echo "ERROR: the file for predicted distance/orientation info does not exist: $pairMatrixFile"
	exit 1
fi

date

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
        $DistanceFoldingHome/Scripts4Rosetta/GenRosettaPotential.sh -s 3 -d $cstfolder $pairMatrixFile $propertyFile
        cstfile=$cstfolder/$target.pairPotential4Rosetta.SPLINE.txt
fi

ResDir=$savefolder/${target}-FoldRelaxResults
if [ ! -d $ResDir ]; then
	mkdir -p ${ResDir}
fi

program=$DistanceFoldingHome/Scripts4Rosetta/FoldNRelax3.py

i=0
while [ $i -lt $numModels ];
do
        python $program $inFile $cstfile -s $ResDir/ -p
        i=`expr $i + 1`
done

if [ "$propertyFile" != "cst" ]; then
        rm -rf $cstfolder
fi

date
