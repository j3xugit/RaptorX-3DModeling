#!/bin/sh

if [ $# -lt 2 ]; then
	echo "Usage: $0 seqFile numModels cstFolder "
	echo "		seqFile: the sequence file in FASTA format"
	echo "		numModels: the number of 3D models to be generated"
	echo "		cstFolder is CASP13MyDMUCE3URE3_EC34CL99S35DFIRE17.5-s5/ by default"
	echo "		MyDM2DMfile is CASP13/CASP13FM2MyDM.txt by default"
	exit 1
fi

seqfile=$1
target=`basename $seqfile .fasta`
target=`basename $target .seq`

numModels=$2

cstDir=CASP13MyDMUCE3URE3_EC34CL99S35DFIRE17.5-s5/
if [ $# -ge 3 ]; then
	cstDir=$3
fi

potSuffix=SPLINE
cstfile=${cstDir}/${target}.pairPotential4Rosetta.${potSuffix}.txt

base=`basename $cstDir`
ResDir=${target}-Fold-${base}
mkdir -p ${ResDir}

keywords=FoldNRelax.py
myself=`basename $0 `
dname=`dirname "$(readlink -f "$0")"`

## at most this number of jobs to be run simultaneously
numAllowedJobs=15

i=0
while [ $i -lt $numModels ];
do
	while true
	do
		sleep 1
		a=`expr $RANDOM % 10 `
                sleep $a

                ##check the number of running jobs
                #echo ${keywords}
                numRunningJobs=`ps -x | grep ${keywords} | grep -v ${myself} | wc -l  `
                #echo $numRunningJobs

                if [ $numRunningJobs -lt `expr $numAllowedJobs + 1 ` ]; then
			python $dname/FoldNRelax.py $seqfile $cstfile -n 5000 -t 0.00001 -v $ResDir/ -q &
			break
		fi
	done
	i=`expr $i + 1`
done

wait

