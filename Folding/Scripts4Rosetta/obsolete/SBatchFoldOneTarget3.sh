#!/bin/sh

if [ $# -lt 3 ]; then
	echo $0 proteinSeqFile predictedPairInfo predictedPropertyInfo [savefolder [numModels] ]
	echo "	this script generates Rosetta constraints from predicted inter-residue relationship and PhiPsi angles and then fold the protein"
	echo "	inFile: the primary sequence file in FASTA format"
	echo "	predictedPairInfo: a PKL file for all predicted inter-residue information including distance and orientation"
	echo "	predictedPropertyInfo: a PKL file for predicted Phi Psi angles"
	echo "	savefolder: the folder for result saving, default current work directory"
	echo "	numModels: the number of models to be generated, default 10"
	exit 1
fi

if [[ -z "${DistanceFoldingHome}" ]]; then
        echo "ERROR: Please set the environmental variable DistanceFoldingHome, e.g., $HOME/3DModeling/Folding"
        exit 1
fi

inFile=$1
pairMatrixFile=$2
propertyFile=$3

target=`basename $inFile`
target=`echo $target | cut -f1 -d'.' `

seqLen=`tail -1 $seqFile | wc -c`
QUEUE=contrib-cpu
if [ $seqLen -ge 380 ]; then
        QUEUE=contrib-cpu-long
fi

savefolder=`pwd`
if [ $# -ge 4 ]; then
	savefolder=$4
	mkdir -p $savefolder
fi

numModels=10
if [ $# -ge 5 ]; then
	numModels=$5
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

program=${DistanceFoldingHome}/Scripts4Rosetta/GenPotentialNFold3.sh
i=0
while [ $i -lt $numModels ];
do
        sbatch -p $QUEUE -J ${target}Fold -o ${target}.fold.out -e ${target}.fold.out -C avx $program $inFile $pairMatrixFile $propertyFile $savefolder 1
        i=`expr $i + 1 `
done
