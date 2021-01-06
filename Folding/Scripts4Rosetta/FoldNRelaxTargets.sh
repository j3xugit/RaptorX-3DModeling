#!/bin/sh

savefolder=`pwd`
numModels=300
numCPUs=18
runningMode=0
UsePerturbation=false

function Usage 
{
	echo $0 "[ -d savefolder | -n numModels | -c numCPUs | -r runningMode |  -p ] proteinListFile inFileDir cstFolder [propertyFolder]"
	echo "	This script folds/relaxes jobs for a list of proteins"
	echo "	proteinListFile: a file for a list of proteins, each in one row"
	echo "	inFileDir: folder for input sequence file in FASTA format or PDB file for an initial model"
	echo "	cstFolder: folder for predicted distance/orientation files (in PKL format) or folder for Rosetta constraint files with suffix .pairPotential4Rosetta.SPLINE.txt"
	echo "	propertyFolder: if not specified, cstFolder shall be a folder for Rosetta constraints. Otherwise cstFolder and propertyFolder are folders for predicted distance and property files, respectively"
	echo "	-d: the folder for result saving, default current work directory"
        echo "	-n: the number of models to be generated, default $numModels"
	echo "	-c: the number of CPUs to be used simultaneously, default $numCPUs"
        echo "	-r: running mode: 0 (fold only) or 1 (fold+relax), default $runningMode"
	echo "	-p: use perturbation at folding stage, default No"
}

while getopts ":n:c:d:r:p" opt; do
        case ${opt} in
                n )
                  numModels=$OPTARG
                  ;;
		c )
		  numCPUs=$OPTARG
		  ;;
                d )
                  savefolder=$OPTARG
                  ;;
                r )
                  runningMode=$OPTARG
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
	echo "ERROR: invalid folder for Rosetta constraints $cstFolder"
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
		$cmdDir/FoldNRelaxOneTarget.sh -d $ResDir -n $numModels -c $numCPUs -r $runningMode -p $seqFile $cstFile $propertyFile
	else
		$cmdDir/FoldNRelaxOneTarget.sh -d $ResDir -n $numModels -c $numCPUs -r $runningMode $seqFile $cstFile $propertyFile
	fi
done
