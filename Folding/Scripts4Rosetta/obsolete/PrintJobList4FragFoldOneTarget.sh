#!/bin/sh

if [ $# -lt 4 ]; then
	echo "Usage: $0 seqFile FragDir predictedMatrixFile predictedPropertyFile [savefolder [numModels] ]"
	echo "	This scripts fold a protein sequence using fragments and predicted distance/orientation/angle information"
	echo "	seqFile: the sequence file in FASTA format"
	echo "	FragDir: the folder for fragments"
	echo "	predictedMatrixFile: the PKL file for predicted distance/orientaton information"
	echo "	predictedPropertyFile: the PKL file for predicted Phi/Psi info"
	echo "	savefolder: the folder for result saving"
	echo "	numModels: the number of 3D models to be generated, default 150"
	exit 1
fi

seqfile=$1
target=`basename $seqfile .fasta`
target=`basename $target .seq`

FragDir=$2
predMatrixFile=$3
predPropertyFile=$4

savefolder=`pwd`
if [ $# -ge 5 ]; then
	savefolder=$5
	mkdir -p $savefolder
fi

numModels=150
if [ $# -ge 6 ]; then
	numModels=$6
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`
program=$cmdDir/GenPotentialNFragFold.sh
FlagTemplateFile=$cmdDir/FlagTemplates/flags4FragFold.txt

numModelsPerJob=1
i=0
while [ $i -lt $numModels ];
do
	echo "$program $seqfile $FragDir $predMatrixFile $predPropertyFile $FlagTemplateFile $savefolder $numModelsPerJob "
	i=` expr $i + $numModelsPerJob `
done
