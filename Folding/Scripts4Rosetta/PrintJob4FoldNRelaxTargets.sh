#!/bin/bash

savefolder=`pwd`
numModels=300
runningMode=0
UsePerturbation=false

alpha=1.61

function Usage 
{
	echo $0 "[ -d savefolder | -n numModels | -r runningMode| -a alpha | -p ] proteinListFile inFileDir cstFolder [propertyFolder]"
	echo "	This script prints out folding/relaxation jobs for a list of proteins"
	echo "	proteinListFile: a file for a list of proteins, each in one row"
	echo "	inFileDir: the folder for input sequence file in FASTA format or PDB file for an initial model"
	echo "	cstFolder: the folder for Rosetta constraint files, each having name targetName.pairPotential4Rosetta.SPLINE.txt"
	echo "	propertyFolder: if not specified, cstFolder shall be a folder for Rosetta constraints. Otherwise cstFolder and propertyFolder are folders for predicted distance and property files"
	echo "	-d: the folder for saving the resultant 3D models, default current work directory"
        echo "	-n: the number of models to be generated, default $numModels"
        echo "	-r: running mode: 0 (fold only) or 1 (fold+relax), default $runningMode"
	echo "	-a: alpha used in DFIRE potential, default $alpha; if >20, a random value between 1.57 and 1.63 will be used"
	echo "	-p: use perturbation at folding stage, default No"
}

while getopts ":n:d:r:a:p" opt; do
        case ${opt} in
                n )
                  numModels=$OPTARG
                  ;;
                d )
                  savefolder=$OPTARG
                  ;;
                r )
                  runningMode=$OPTARG
                  ;;
		a )
		  alpha=$OPTARG
		  ;;
		p )
		  UsePerturbation=true
		  ;;
                \? )
                  echo "Invalid Option: -$OPTARG" 1>&2
                  exit 1
                  ;;
                : )
                  echo "Invalid Option: -$OPTARG requires an argument" 1>&2
                  exit 1
                  ;;
        esac
done
shift $((OPTIND -1))

if [ $runningMode -ne 0 -a $runningMode -ne 1 ]; then
	echo "ERROR: invalid running mode"
	exit 1
fi

if [ $# -lt 3 -o $# -gt 4 ]; then
	Usage
	exit 1
fi

proteinList=$1
if [ ! -f $proteinList ]; then
	echo "ERROR: invalid protein list file $proteinListFile"
	exit 1
fi

inDir=$2
if [ ! -d $inDir ]; then
	echo "ERROR: invalid folder for input files $inDir"
	exit 1
fi

cstFolder=$3
if [ ! -d $cstFolder ]; then
	echo "ERROR: invalid folder for predicted distance/orientation PKL files or Rosetta constraints $cstFolder"
	exit 1
fi

propertyFolder='cst'
if [ $# -eq 4 ]; then
	propertyFolder=$4
fi
if [ "$propertyFolder" != "cst" ]; then
	if [ ! -d $propertyFolder ]; then
		echo "ERROR: invalid folder for predicted property files: $propertyFolder"
		exit 1
	fi
	
fi

if [ ! -d $savefolder ]; then
	mkdir -p $savefolder
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

for target in `cat $proteinList`
do
	if [ $runningMode -eq 0 ]; then
		ResDir=$savefolder/${target}-FoldResults
	else
		ResDir=$savefolder/${target}-RelaxResults
	fi
	#mkdir -p $ResDir
	
	seqFile=$inDir/${target}.fasta
        if [ ! -f $seqFile ]; then
                seqName=`echo $target | cut -f1 -d'-' `
                seqFile=$inDir/${seqName}.fasta
        fi

	if [ "$propertyFolder" == "cst" ]; then
		cstFile=$cstFolder/${target}.pairPotential4Rosetta.SPLINE.txt
		propertyFile=""
	else
		cstFile=$cstFolder/${target}.predictedDistMatrix.pkl
		propertyFile=${propertyFolder}/${target}.predictedProperties.pkl
		if [ ! -f $propertyFile ]; then
			seqName=`echo $target | cut -f1 -d'-' `
			propertyFile=${propertyFolder}/${seqName}.predictedProperties.pkl
		fi
	fi

	if $UsePerturbation; then
		$cmdDir/PrintJob4FoldNRelaxOneTarget.sh -d $ResDir -n $numModels -r $runningMode -a $alpha -p $seqFile $cstFile $propertyFile
	else
		$cmdDir/PrintJob4FoldNRelaxOneTarget.sh -d $ResDir -n $numModels -r $runningMode -a $alpha $seqFile $cstFile $propertyFile
	fi
done
