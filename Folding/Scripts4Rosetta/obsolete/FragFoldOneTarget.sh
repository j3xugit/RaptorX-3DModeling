#!/bin/sh

if [ $# -lt 4 ]; then
	echo "Usage: $0 seqFile FragDir predictedMatrixFile predictedPropertyFile [savefolder [numModels] ]"
	echo "	This scripts fold a protein sequence using fragments and predicted distance/orientation/angle information"
	echo "	seqFile: the sequence file in FASTA format"
	echo "	FragDir: the folder for fragments"
	echo "	predictedMatrixFile: the PKL file for predicted distance/orientaton information"
	echo "	predictedPropertyFile: the PKL file for predicted Phi/Psi info"
	echo "	savefolder: the folder for result saving"
	echo "	numModels: the number of 3D models to be generated"
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
keywords=GenPotentialNFragFold

## at most this number of jobs to be run simultaneously
numAllowedJobs=15
numModelsPerJob=3
i=0
while [ $i -lt $numModels ];
do
	while true
	do
		sleep 1
		a=`expr $RANDOM % 5 `
                sleep $a

                ##check the number of running jobs
                numRunningJobs=`ps -x | grep ${keywords} | wc -l  `

                if [ $numRunningJobs -lt `expr $numAllowedJobs + 1 ` ]; then
			$program $seqfile $FragDir $predMatrixFile $predPropertyFile $savefolder $numModelsPerJob &
			break
		fi
	done
	i=`expr $i + $numModelsPerJob `
done

wait

