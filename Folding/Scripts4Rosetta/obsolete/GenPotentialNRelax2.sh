#!/bin/sh

if [ $# -lt 3 ]; then
	echo $0 initModelFile predictedPairInfo predictedPropertyInfo [savefolder]
	echo "	this script generates Rosetta constraints from predicted inter-residue relationship and PhiPsi angles and then fold the protein"
	echo "	initModelFile: the pdb file for an initial model"
	echo "	predictedPairInfo: a PKL file for all predicted inter-residue information including distance and orientation"
	echo "	predictedPropertyInfo: a PKL file for predicted Phi Psi angles"
	echo "	savefolder: the folder for result saving, default current work directory"
	exit 1
fi

if [[ -z "${DistanceFoldingHome}" ]]; then
        echo "ERROR: Please set the environmental variable DistanceFoldingHome, e.g., $HOME/3DModeling/Folding"
        exit 1
fi

initModelFile=$1
pairMatrixFile=$2
propertyFile=$3

bname=`basename $initModelFile .pdb`
target=`echo $bname | cut -f1 -d'.'`

savefolder=`pwd`
if [ $# -ge 4 ]; then
	savefolder=$4
	mkdir -p $savefolder
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

ResDir=$savefolder/${target}-RelaxResults
if [ ! -d $ResDir ]; then
	mkdir -p ${ResDir}
fi

program=$DistanceFoldingHome/Scripts4Rosetta/Relax2.py
#python $program $initModelFile $cstfile -e 1.5 -d 0.2 -w 0.2 -a 0.2 -s $ResDir/ 
python $program $initModelFile $cstfile -r -e 1.5 -s $ResDir/ 

if [ "$propertyFile" != "cst" ]; then
        rm -rf $cstfolder
fi

date
