#!/bin/sh

if [ $# -lt 4 ]; then
	echo $0 proteinListFile inFileDir cstFolder savefolder [numModels]
	echo "	this script generates Rosetta constraints from predicted inter-residue relationship and PhiPsi angles and then fold the protein"
	echo "	inFileDir: the folder for input sequence file in FASTA format or a PDB file for an initial model"
	echo "	cstFolder: the folder for all Rosetta constraints including those for inter-residue distance and orientation and Phi/Psi angles"
	echo "	savefolder: the folder for result saving, default current work directory"
	echo "	numModels: the number of models to be generated, default 150"
	exit 1
fi

proteinList=$1
inDir=$2
cstFolder=$3
savefolder=$4

if [ ! -d $savefolder ]; then
	mkdir -p $savefolder
fi

numModels=150
if [ $# -ge 5 ]; then
	numModels=$5
fi

if [ ! -d $cstFolder ]; then
	echo "ERROR: the folder for predicted distance/orientation info does not exist: $cstFolder"
	exit 1
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

for target in `cat $proteinList`
do
	$cmdDir/PrintJobList4FoldOneTarget3.sh $inDir/${target}.fasta $cstFolder/${target}.pairPotential4Rosetta.SPLINE.txt cst $savefolder $numModels
done
