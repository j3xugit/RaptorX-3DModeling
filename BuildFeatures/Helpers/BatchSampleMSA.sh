#!/bin/bash

if [ $# -lt 3 ]; then
	echo $0 proteinListFile inDir ratio [numSamples [ outDir ] ]
	echo "	inDir: the folder for all input a3m files"
	echo "	outDir: the folder for result files"
	exit 1
fi

list=$1
MSADir=$2
ratio=$3

numSamples=8
if [ $# -ge 4 ]; then
	numSamples=$4
fi

ResDir=SampledMSAs
if [ $# -ge 5 ]; then
	ResDir=$5
fi
mkdir -p $ResDir

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`
program=$cmdDir/SampleA3M.py

numAllowedJobs=15
keywords=`basename $program`
myself=`basename $0 `


for target in `cat $list `
do
	while true
        do
                ## check the number of running jobs
                numRunningJobs=`ps -x | grep ${keywords} | grep -v ${myself} | wc -l`
                if [ $numRunningJobs -lt `expr $numAllowedJobs + 1 ` ]; then
			(python $program $MSADir/${target}.a3m $ratio $numSamples; mv ${target}.a3m_S? $ResDir/) &
                        break
                else
                        a=`expr $RANDOM % 4 `
                        sleep $a
                fi

        done

done

wait

