#!/bin/bash

savefolder=`pwd`
numModels=120
alpha=1.61

runningMode=0
UsePerturbation=false

host=`hostname`

SPEC='-C avx'
QUEUE=contrib-cpu
if [ "$host" == "raptorx3.uchicago.edu" ]; then
	QUEUE=cpu
	SPEC=''
fi

if [[ -z "${DistanceFoldingHome}" ]]; then
        echo "ERROR: Please set the environmental variable DistanceFoldingHome to installation folder of Folding"
        exit 1
fi

function Usage 
{
	echo $0 "[-d savefolder | -n numModels | -r runningMode | -q partition | -a alpha | -p ] seqFile predictedPairInfo predictedPropertyInfo"
	echo "	This script folds a protein on a slurm cluster using predicted distance/orientation/angles"
	echo "	seqFile: the primary sequence file in FASTA format"
	echo "	predictedPairInfo: a Rosetta constraint file or a PKL file for predicted distance/orientation"
	echo "	predictedPropertyInfo: a string cst or a PKL file for predicted Phi/Psi angles"
	echo "		when it is cst, it indicates that predictedPairInfo is a Rosetta constraint text file instead of a PKL file"
	echo "	-n: the number of models to be generated, default $numModels"
	echo "	-d: the folder for result saving, default current work directory"
	echo "	-r: 0 (default) or 1,  0 for fold only and 1 for fold+relax"
	echo "	-p: use perturbation at folding stage, default No"
	echo "	-q: the job queue of slurm cluster, default $QUEUE"
	echo "  -a: alpha value for DFIRE potential, default 1.61. If > 20, then a random alpha will be used"
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

if [ $# -lt 3 ]; then
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
	if [ $seqLen -ge 400 ]; then
        	QUEUE=contrib-cpu-long
	fi
fi

if [ "$host" == "raptorx3.uchicago.edu" ]; then
	if [ $seqLen -ge 400 ]; then
		SPEC=$SPEC" -c2 "
	fi
fi

pairMatrixFile=$2
if [ ! -f $pairMatrixFile ]; then
	echo "ERROR: invalid file for Rosetta constraints or predicted distance/orientation info: $pairMatrixFile"
	exit 1
fi

propertyFile=$3
if [ "$propertyFile" != "cst" ]; then
	if [ ! -f $propertyFile ]; then
		echo "ERROR: invalid file for predicted property: $propertyFile"
		exit 1
	fi
fi

program=${DistanceFoldingHome}/Scripts4Rosetta/GenPotentialNFoldRelax.sh
if [ ! -x $program ]; then
	echo "ERROR: invalid or not exectuable program $program"
	exit 1
fi

i=0
while [ $i -lt $numModels ];
do
	if $UsePerturbation; then
        	sbatch -p $QUEUE -J ${target} -o ${target}.fold.out -e ${target}.fold.out $SPEC $program -d $savefolder -r $runningMode -a $alpha -p $inFile $pairMatrixFile $propertyFile
	else
        	sbatch -p $QUEUE -J ${target} -o ${target}.fold.out -e ${target}.fold.out $SPEC $program -d $savefolder -r $runningMode -a $alpha $inFile $pairMatrixFile $propertyFile
	fi
        i=`expr $i + 1`
done
