#!/bin/bash

savefolder=`pwd`
numModels=120
alpha=1.61

runningMode=0
UsePerturbation=false

propertyFile="cst"
machineType=0

RemoteAccountInfo=""
RemoteAccount=""
RemoteWorkDirBase=""

function Usage 
{
	echo $0 "[-R remoteAccount | -d savefolder | -n numModels | -r runningMode | -a alpha | -t machineType | -p ] seqFile predictedPairInfo [predictedPropertyInfo]"
	echo "	This script folds a protein using a remote account"
	echo "		GNU parallel shall be installed for all nodes since it will be called to run jobs on each node"
	echo "	seqFile: the primary sequence file in FASTA format"
	echo "	predictedPairInfo: a Rosetta constraint file or a PKL file for predicted distance/orientation"
	echo "	predictedPropertyInfo: could be empty, a string 'cst' or a PKL file for predicted Phi/Psi angles"
	echo "		when empty or 'cst', it indicates that predictedPairInfo is a Rosetta constraint text file instead of a PKL file"
	echo "	-R: the remote account, e.g., raptorx@raptorx3.uchicago.edu:RemoteWorkDir or raptorx10.uchicago.edu:RemoteWorkDir"
	echo "		if RemoteWorkDir is empty, then remote job will be run directly under the home directory of the remote account"
	echo "		if -R is not used, then run folding jobs locally (default)"
	echo "		By default, the remote machine is assumed to be a slurm cluster such as raptorx@slurm.ttic.edu. If not, you may need to modify some code in this script to make it work"
	echo "	-n: the number of models to be generated, default $numModels"
	echo "	-d: the folder for result saving, default current work directory"
	echo "	-r: 0 (default) or 1, 0 for fold only and 1 for fold+relax"
	echo "	-p: use perturbation at folding stage, default No"
	echo "	-q: the job queue of slurm cluster, default $QUEUE"
	echo "	-a: alpha value for DFIRE potential, default 1.61. If > 20, then a random alpha will be used"
	echo "	-t: the type of machine for folding jobs: 0 (default) for self-determination (which may not work for your machine), 1 for a mulit-CPU Linux computer with GNU parallel installed,"
        echo "          2 for a slurm cluster with similar nodes, 3 for a slurm cluster with hetergenous nodes and 4 for a multi-CPU Linux computer without GNU parallel installed"
}

while getopts ":R:n:d:r:a:t:p" opt; do
        case ${opt} in
		R )
		  RemoteAccountInfo=$OPTARG
		  RemoteAccount=`echo $OPTARG | cut -f1 -d':'`
		  RemoteWorkDirBase=`echo $OPTARG | cut -f2 -d':'`
		  ;;
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

if [ $# -eq 3 ]; then
	propertyFile=$3
fi

if [ "$propertyFile" != "cst" ]; then
	if [ ! -f $propertyFile ]; then
		echo "ERROR: invalid file for predicted property: $propertyFile"
		exit 1
	fi
fi

## create a remote work directory
localMachine=`hostname | cut -f1 -d'.'`
if [ -z "${RemoteWorkDirBase}" ]; then
	RemoteWorkDir=tmpWorkDir4RemoteDistFolding-${target}-$localMachine-$$
else
	RemoteWorkDir=$RemoteWorkDirBase/tmpWorkDir4RemoteDistFolding-${target}-$localMachine-$$
fi

ssh -o StrictHostKeyChecking=no $RemoteAccount "mkdir -p $RemoteWorkDir"
if [ $? -ne 0 ]; then
	echo "ERROR: failed to create $RemoteWorkDir at $RemoteAccount!"
        exit 1
fi

## transfer three files to remote account
scp $inFile $pairMatrixFile $propertyFile $RemoteAccount:$RemoteWorkDir/
if [ $? -ne 0 ]; then
        echo "ERROR: failed to transfer three files to $RemoteWorkDir at $RemoteAccount"
        exit 1
fi

## execute the folding jobs at the remote account
inFileBase=`basename $inFile`
pairMatrixFileBase=`basename $pairMatrixFile`
propertyFileBase=`basename $propertyFile`

program=LocalFoldNRelaxOneTarget.sh
ssh -o StrictHostKeyChecking=no $RemoteAccount "cd $RemoteWorkDir; \$DistanceFoldingHome/$program -t $machineType -d ${target}-RelaxResults -n $numModels -r $runningMode $inFileBase $pairMatrixFileBase $propertyFileBase"
if [ $? -ne 0 ]; then
        echo "ERROR: failed to execute $program -t $machineType -d ${target}-RelaxResults -n $numModels -r $runningMode $inFileBase $pairMatrixFileBase $propertyFileBase at $RemoteWorkDir of $RemoteAccount"
        exit 1
fi

## transfer the results back
mkdir -p $savefolder
rsync -av $RemoteAccount:$RemoteWorkDir/${target}-RelaxResults/ $savefolder/
if [ $? -ne 0 ]; then
        echo "ERROR: failed to run rsync -av $RemoteAccount:$RemoteWorkDir/${target}-RelaxResults/ $savefolder/"
        exit 1
fi

## remove the remote work dir
ssh -o StrictHostKeyChecking=no $RemoteAccount "rm -rf $RemoteWorkDir"
if [ $? -ne 0 ]; then
        echo "WARNING: failed to delete $RemoteWorkDir at $RemoteAccount!"
fi
