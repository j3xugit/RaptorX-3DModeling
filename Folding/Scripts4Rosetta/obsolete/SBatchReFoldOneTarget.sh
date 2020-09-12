#!/bin/sh

if [ $# -lt 3 ]; then
	echo $0 folder4models predictedPairInfo predictedPropertyInfo [savefolder]
	echo "	this script generates Rosetta constraints from predicted inter-residue relationship and PhiPsi angles and then fold the protein"
	echo "	folder4models: the folder for initial 3D models, each model file has PDB format "
	echo "	predictedPairInfo: a PKL file for all predicted inter-residue information including distance and orientation"
	echo "	predictedPropertyInfo: a PKL file for predicted Phi Psi angles"
	echo "	savefolder: the folder for result saving, default current work directory"
	exit 1
fi

if [[ -z "${DistanceFoldingHome}" ]]; then
        echo "ERROR: Please set the environmental variable DistanceFoldingHome, e.g., $HOME/3DModeling/Folding"
        exit 1
fi

folder4initModels=$1
pairMatrixFile=$2
propertyFile=$3

target=`basename $folder4initModels`
target=`echo $target | cut -f1 -d'.' `

QUEUE=contrib-cpu

savefolder=`pwd`
if [ $# -ge 4 ]; then
	savefolder=$4
	mkdir -p $savefolder
fi

if [ ! -f $pairMatrixFile ]; then
	echo "ERROR: the file for predicted distance/orientation info does not exist: $pairMatrixFile"
	exit 1
fi

if [ "$propertyFile" != "cst" ]; then
	if [ ! -f $propertyFile ]; then
		echo "ERROR: the file for predicted property info does not exist: $propertyFile"
		exit 1
	fi
fi

program=${DistanceFoldingHome}/Scripts4Rosetta/GenPotentialNFold2.sh

for model in $folder4initModels/*.pdb
do
        sbatch -p $QUEUE -J ${target}ReFold -o ${target}.fold.out -e ${target}.fold.out -C avx $program $model $pairMatrixFile $propertyFile $savefolder 1
done
