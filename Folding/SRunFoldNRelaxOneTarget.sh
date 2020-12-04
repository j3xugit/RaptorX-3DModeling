#!/bin/bash

savefolder=`pwd`
numModels=120
alpha=1.61

runningMode=0
UsePerturbation=false

propertyFile="cst"

host=`hostname`

####  You may adjust the below setting for your own slurm cluster
QUEUE=cpu
SPEC=" -c20 --exclude=cpu3 --nodes=1-1 "
if [ "$host" == "slurm.ttic.edu" ]; then
	SPEC=" -C avx --nodes=1-1 "
	QUEUE=contrib-cpu
fi
###############################

if [[ -z "${DistanceFoldingHome}" ]]; then
        echo "ERROR: Please set the environmental variable DistanceFoldingHome to installation folder of Folding"
        exit 1
fi

function Usage 
{
	echo $0 "[-d savefolder | -n numModels | -r runningMode | -q partition | -a alpha | -p ] seqFile predictedPairInfo [predictedPropertyInfo]"
	echo "	This script folds a protein using predicted distance/orientation/angles on a slurm cluster with nodes of similar configurations"
	echo "		GNU parallel shall be installed for all nodes since it will be called to run jobs on each node"
	echo "	seqFile: the primary sequence file in FASTA format"
	echo "	predictedPairInfo: a Rosetta constraint file or a PKL file for predicted distance/orientation"
	echo "	predictedPropertyInfo: could be empty, a string 'cst' or a PKL file for predicted Phi/Psi angles"
	echo "		when empty or 'cst', it indicates that predictedPairInfo is a Rosetta constraint text file instead of a PKL file"
	echo "	-n: the number of models to be generated, default $numModels"
	echo "	-d: the folder for result saving, default current work directory"
	echo "	-r: 0 (default) or 1,  0 for fold only and 1 for fold+relax"
	echo "	-p: use perturbation at folding stage, default No"
	echo "	-q: the job queue of slurm cluster, default $QUEUE"
	echo "	-a: alpha value for DFIRE potential, default 1.61. If >20, a random value will be used"
}

while getopts ":n:d:r:q:a:p" opt; do
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
		q )
		  QUEUE=$OPTARG
		  ;;
		a )
		  alpha=$OPTARG
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
seqLen=`tail -n +1 $inFile | wc -c`

if [ "$host" == "slurm.ttic.edu" ]; then
	if [ $seqLen -ge 450 ]; then
        	QUEUE=contrib-cpu-long
	fi
fi

pairMatrixFile=$2
if [ ! -f $pairMatrixFile ]; then
	echo "ERROR: invalid file for Rosetta constraints or predicted distance/orientation info: $pairMatrixFile"
	exit 1
fi

if [ $# -eq 3 ]; then
	propertyFile=$3
fi

if [ "$propertyFile" != "cst" ]; then
	if [ ! -f $propertyFile ]; then
		echo "ERROR: invalid file for predicted property: $propertyFile"
		exit 1
	fi
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`
#program=${DistanceFoldingHome}/ParallelFoldNRelaxOneTarget.sh
program=${cmdDir}/ParallelFoldNRelaxOneTarget.sh
if [ ! -x $program ]; then
	echo "ERROR: invalid or non-exectuable program $program"
	exit 1
fi

###### You may revise the below setting depending on the amount of CPU memory at a node
## numTasks: the number of tasks at each batch
numTasks=20
if [ $seqLen -gt 500 ]; then
	numTasks=`expr 10000 / $seqLen`
fi
if [ $numTasks -lt 10 ]; then
	numTasks=10
fi
#############################################################

numModels=`expr $numModels + $numTasks - 1 `
numRuns=`expr $numModels / $numTasks `

SRUN="srun -p $QUEUE -J ${target} -o ${target}.fold.out -e ${target}.fold.out $SPEC"

command=" $program -d $savefolder -n $numTasks -r $runningMode -a $alpha "
if $UsePerturbation; then
	command=$command" -p "
fi

parallelOptions=" --memfree 20G "

i=0
while [ $i -lt $numRuns ];
do
	$SRUN $command -o "${parallelOptions}" $inFile $pairMatrixFile $propertyFile &
        i=`expr $i + 1`
	sleep 1
done

wait
