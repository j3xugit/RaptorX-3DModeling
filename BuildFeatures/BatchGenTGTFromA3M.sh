#!/bin/bash

if [ $# -lt 2 ]; then
	echo "$0 proteinList A3MDir [ResDir]"
	echo "	This script generates .tgt files from a3m files for a list of proteins in parallel"
	echo "	proteinList: a file containing a list of proteins, each in one row"
	echo "	A3MDir: the folder for the MSA in a3m format"
	echo "	ResDir: the folder for result files, default current work directory"
	exit 1
fi

proteinList=$1
A3MDir=$2

ResDir=`pwd`
if [ $# -ge 3 ]; then
	ResDir=$3
	if [ ! -d $ResDir ]; then
		mkdir -p $ResDir
	fi
fi

machine=`hostname`
if [[ "$machine" == "raptorx10.uchicago.edu" ]]; then
	numAllowedJobs=250
else
	numAllowedJobs=19
fi

if [ $# -ge 4 ]; then
	numAllowedJobs=$4
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

program=$cmdDir/GenTGTFromA3M.sh
keywords=`basename $program`
myself=`basename $0 `

for target in `cat $proteinList `
do
	while true
        do
         	## check the number of running jobs
                numRunningJobs=`ps -x | grep ${keywords} | grep -v ${myself} | wc -l`
                if [ $numRunningJobs -lt $numAllowedJobs ]; then
			$program $A3MDir/$target.a3m $ResDir &
                        break
                else
                        a=`expr $RANDOM % 3 `
                        sleep $a
                fi

        done

done

wait
