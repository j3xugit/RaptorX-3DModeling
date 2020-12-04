#!/bin/bash

savefolder=`pwd`
numModels=120
alpha=1.61

runningMode=0
UsePerturbation=false

machineType=0

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

function Usage 
{
	echo $0 "[ -d savefolder | -n numModels | -r runningMode | -a alpha | -t MachineType | -p ] seqFile predictedPairInfo [predictedPropertyInfo]"
	echo "	This script folds a protein using local CPUs"
	echo "	seqFile: the primary sequence file in FASTA format"
	echo "	predictedPairInfo: a Rosetta constraint file or a PKL file for predicted distance/orientation"
	echo "	predictedPropertyInfo: could be empty, a string 'cst' or a PKL file for predicted Phi/Psi angles"
	echo "		when empty or 'cst', it indicates that predictedPairInfo is a Rosetta constraint text file instead of a PKL file containing predicted distance/orientation probability"
	echo "	-n: the number of models to be generated, default $numModels"
	echo "	-d: the folder for result saving, default current work directory"
	echo "	-r: 0 (default) or 1,  0 for fold only and 1 for fold+relax"
	echo "	-p: use perturbation at folding stage, default No"
	echo "	-q: the job queue of slurm cluster, default $QUEUE"
	echo "	-a: alpha value for DFIRE potential, default 1.61. If > 20, a random value between 1.57 and 1.63 will be used"
	echo "	-t: the type of machine for folding jobs: 0 (default) for self-determination (which may not work for your machine), 1 for a mulit-CPU Linux computer with GNU parallel installed,"
        echo "          2 for a slurm cluster with similar nodes, 3 for a slurm cluster with hetergenous nodes and 4 for a multi-CPU Linux computer without GNU parallel installed"
}

while getopts ":n:d:r:a:t:p" opt; do
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
		p )
		  UsePerturbation=true
		  ;;
		a )
		  alpha=$OPTARG
		  ;;
		t )
		  machineType=$OPTARG
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
	Usage
	exit 1
fi

if [ $# -lt 2 -o $# -gt 3 ]; then
	Usage
	exit 1
fi

inFile=$1
if [ ! -f $inFile ]; then
	echo "ERROR: invalid input file for folding $inFile"
	exit 1
fi
target=`basename $inFile`
target=`echo $target | cut -f1 -d'.' `

pairMatrixFile=$2
if [ ! -f $pairMatrixFile ]; then
	echo "ERROR: invalid file for Rosetta constraints or predicted distance/orientation info: $pairMatrixFile"
	exit 1
fi

propertyFile="cst"
if [ $# -eq 3 ]; then
	propertyFile=$3
fi

if [ "$propertyFile" != "cst" ]; then
	if [ ! -f $propertyFile ]; then
		echo "ERROR: invalid file for predicted property: $propertyFile"
		exit 1
	fi
fi

host=`hostname`

if [ $machineType -eq 0 ]; then
	program=Scripts4Rosetta/FoldNRelaxOneTarget.sh

	if [[ "$host" == raptorx10.uchicago.edu ]]; then
		program=ParallelFoldNRelaxOneTarget.sh

	elif [[ "$host" == raptorx3.uchicago.edu ]]; then
		program=SRunFoldNRelaxOneTarget.sh

	elif [[ "$host" == slurm.ttic.edu ]]; then
		program=SlurmFoldNRelaxOneTarget.sh

	elif [[ "$host" == raptorx[456789].uchicago.edu ]]; then
		program=ParallelFoldNRelaxOneTarget.sh
	fi

elif [ $machineType -eq 1 ]; then
	program=ParallelFoldNRelaxOneTarget.sh

elif [ $machineType -eq 2 ]; then
	program=SRunFoldNRelaxOneTarget.sh

elif [ $machineType -eq 3 ]; then
	program=SlurmFoldNRelaxOneTarget.sh

elif [ $machineType -eq 4 ]; then
	program=Scripts4Rosetta/FoldNRelaxOneTarget.sh
else
	echo "ERROR: unsupported machine type: $machineType"
	exit 1
fi

###### you may edit the below setting for your own machine
parallelOptions=""
if [ "$program" == "ParallelFoldNRelaxOneTarget.sh" ]; then
	if [ "$host" == "raptorx10.uchicago.edu" ]; then
		parallelOptions=" --memfree 120G --load 92% "
	elif [ "$host" == "raptorx7.uchicago.edu" -o "$host" == "raptorx8.uchicago.edu" -o "$host" == raptorx9.uchicago.edu ]; then 
		parallelOptions=" --memfree 20G --load 98% "
	fi
fi
############################################

if [ -z "$parallelOptions" ]; then
	$cmdDir/$program -d $savefolder -n $numModels -r $runningMode $inFile $pairMatrixFile $propertyFile
	rcode=$?
else
	$cmdDir/$program -o "${parallelOptions}" -d $savefolder -n $numModels -r $runningMode $inFile $pairMatrixFile $propertyFile
	rcode=$?
fi

if [ $rcode -ne 0 ]; then
        echo "ERROR: failed to execute $cmdDir/$program -d $savefolder -n $numModels -r $runningMode $inFile $pairMatrixFile $propertyFile"
        exit 1
fi
