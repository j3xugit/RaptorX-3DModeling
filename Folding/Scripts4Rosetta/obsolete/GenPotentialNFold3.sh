#!/bin/sh

if [ $# -lt 3 ]; then
	echo $0 inFile predictedMatrixFile predictedProperty [savefolder [numModels] ]
	echo "	this script generates Rosetta constraints from predicted inter-residue relationship and PhiPsi angles and then fold the protein"
	echo "	inFile: the primary sequence file in FASTA format or an initial PDB file"
	echo "		this file shall have name ending with .seq, .fasta or .pdb"
	echo "	predictedMatrixFile: a PKL file (ending with .predictedDistMatrix.pkl) for all predicted inter-residue information including distance and orientation"
	echo "	predictedProperty: a PKL file (ending with .predictedPropertyFile) for predicted Phi Psi angles."
	echo "		If this entry is a special string cst, then predictedMatrixFile is interpreted as a Rosetta constraint file instead of a PKL file"
	echo "	savefolder: the folder for result saving, default current work directory"
	echo "	numModels: the number of models to be generated, default 1"
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

ResDir=$savefolder/${target}-FoldResults
if [ ! -d $ResDir ]; then
	mkdir -p ${ResDir}
fi

program=$DistanceFoldingHome/Scripts4Rosetta/FoldNRelax3.py

i=0
while [ $i -lt $numModels ];
do
        python $program $inFile $cstfile -s $ResDir/ -q -p
        i=`expr $i + 1`
done

if [ "$propertyFile" != "cst" ]; then
        rm -rf $cstfolder
fi

date
