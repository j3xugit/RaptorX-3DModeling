#!/bin/sh

if [ $# -lt 2 ]; then
	echo "Usage: $0 target initModelFolder [cstFolder [MyDM2DMfile] ]"
	echo "		cstFolder is CASP13MyDMUCE3URE3_EC34CL99S35DFIRE17.5-s5/ by default"
	echo "		MyDM2DMfile is CASP13/CASP13FM2MyDM.txt by default"
	exit 1
fi

target=$1

initFolder=$2
base=`basename $initFolder`
ResDir=${target}-Relax-${base}
mkdir -p ${ResDir}

cstDir=CASP13MyDMUCE3URE3_EC34CL99S35DFIRE17.5-s5/
if [ $# -ge 3 ]; then
	cstDir=$3
fi
potSuffix=SPLINE
cstfile=${cstDir}/${target}.pairPotential4Rosetta.${potSuffix}.txt

keywords=Relax.py
myself=`basename $0 | cut -f1 -d'.'`
dname=`dirname "$(readlink -f "$0")"`

## at most this number of jobs to be run simultaneously
numAllowedJobs=15

## refine all the init models
for model in $initFolder/${target}.*.pdb
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
			python $dname/Relax.py $model $cstfile -v $ResDir &
			break
		fi
	done
done

wait

MyDM2DMfile=CASP13/CASP13FM2MyDM.txt
if [ $# -ge 4 ]; then
	MyDM2DMfile=$4
fi
## caculate quality of the refined models
DMs=`grep -w :${target} ${MyDM2DMfile} | cut -f1 -d':' `
for dm in $DMs
do
	$dname/CalcQualityMyDM.sh $dm ${target} $ResDir
done
