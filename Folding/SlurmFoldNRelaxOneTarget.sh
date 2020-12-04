#!/bin/bash

#if [[ -z "${ModelingHome}" ]]; then
#        echo "ERROR: Please set the environmental variable ModelingHome to installation folder of RaptorX-3DModeling"
#        exit 1
#fi

#if [[ -z "${DistanceFoldingHome}" ]]; then
#        echo "ERROR: Please set the environmental variable DistanceFoldingHome to installation folder of Folding"
#        exit 1
#fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

host=`hostname`
####### You may have to revise the below setting to make this script work on your own cluster
SPEC="-C avx"
QUEUE=contrib-cpu
if [[ "$host" == raptorx3.uchicago.edu ]]; then
	SPEC=""
	QUEUE=cpu
fi
###############################################################

savefolder=`pwd`
numModels=120
alpha=1.61
runningMode=0
UsePerturbation=false

propertyFile="cst"

NotBlock=False

function Usage 
{
	echo $0 "[-d savefolder | -n numModels | -r runningMode | -q partition | -a alpha | -p | -b ] seqFile predictedPairInfo [predictedPropertyInfo]"
	echo "	This script runs folding jobs mainly on a slurm cluster with hetergeneous nodes which may have very different configurations"
	echo "		you may have to slighlty edit this script to make it work on your own slurm cluster"
	echo "	seqFile: the primary sequence file in FASTA format"
	echo "	predictedPairInfo: a Rosetta constraint file or a PKL file for predicted distance/orientation"
	echo "	predictedPropertyInfo: could be empty, a string 'cst' or a PKL file for predicted Phi/Psi angles"
	echo "		when empty or 'cst', it indicates that predictedPairInfo is a Rosetta constraint text file instead of a PKL file"
	echo "	-n: the number of models to be generated, default $numModels"
	echo "	-d: the folder for result saving, default current work directory"
	echo "	-r: 0 (default) or 1,  0 for fold only and 1 for fold+relax"
	echo "	-p: use perturbation at folding stage, default No"
	echo "	-q: the job queue of slurm cluster, default $QUEUE"
	echo "  -a: alpha value for DFIRE potential, default 1.61. If >20, a random value will be used"
	echo "	-b: if specified, do not wait for the completition of the jobs. Otherwise wait"
}

while getopts ":n:d:r:q:a:pb" opt; do
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
		b ) 
		  NotBlock=True
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

#####################you may have to revise the below setting ################
if [ $seqLen -ge 500 -a "$host" == slurm.ttic.edu ]; then
        QUEUE=contrib-cpu-long
fi
############################################################################

####you may revise the below setting depending on how much CPU memory is avaiable on a node
if [ $seqLen -ge 1200 ]; then
	SPEC=$SPEC" -c3 "
elif [ $seqLen -ge 600 ]; then
	SPEC=$SPEC" -c2 "
else
	SPEC=$SPEC" -c1 "
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

#program=${DistanceFoldingHome}/Scripts4Rosetta/PrintJob4FoldNRelaxOneTarget.sh
program=${cmdDir}/Scripts4Rosetta/PrintJob4FoldNRelaxOneTarget.sh
if [ ! -x $program ]; then
	echo "ERROR: invalid or not exectuable program $program"
	exit 1
fi


mkdir -p $savefolder
## aim to generate slightly more models in case some jobs crash
numModels2=`expr $numModels + 10 `
JobPrintCmd="${program} -d $savefolder -n $numModels2 -r $runningMode -a $alpha"
if $UsePerturbation; then
	JobPrintCmd=$JobPrintCmd" -p "
fi
JobPrintCmd=$JobPrintCmd" $inFile $pairMatrixFile $propertyFile "
JobFile=${target}.FR.jobs
$JobPrintCmd > $JobFile

ArrayJobGenerator=${cmdDir}/Helpers/ArrayJobs.py
echo "Running ${ArrayJobGenerator} -J ${target} $SPEC $JobFile $QUEUE"
${ArrayJobGenerator} -J ${target} $SPEC $JobFile $QUEUE
if [ ! -f ${target}-sbatch-script-last.txt -o ! -f ${target}-batch-commands-last.txt ]; then
	echo "ERROR: failed to generate two needed files for sbatch"
	exit 1
fi

## submit jobs
retStr=$(sbatch -x cpu23 ${target}-sbatch-script-last.txt)
#echo $retStr
jobID=`echo $retStr | cut -f4 -d' '`
#echo $jobID

if (( $NotBlock )); then
	exit 0
fi

## block until job is done or enough protein models have been generated
while true
do
	## get job status
	jobStatusLineNumber=`squeue -j $jobID | wc -l `
	if [ $jobStatusLineNumber -eq 1 ]; then
		## job is done
		echo "The folding job $jobID is done for $target"
		break
	fi
	
	numDone=`ls ${savefolder}/ | wc -l `
	if [ $numDone -ge $numModels ]; then
		## done, then exit
		echo "Enough 3D models have been generated for $target"
		#scancel $jobID
		break
	fi

	## sleep to wait
	sleep 2m
done
