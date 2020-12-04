#!/bin/bash

if [ $# -lt 2 ]; then
	echo "$0 proteinListFile TGTDir [ResDir]"
	exit 1
fi

list=$1
TGTDir=$2

RESDIR=TGTPKL
if [ $# -ge 3 ]; then
	RESDIR=$3
fi
mkdir -p $RESDIR

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

program=$cmdDir/TGT2PKL.py

numAllowedJobs=15
if [ $# -ge 4 ]; then
	numAllowedJobs=$4
fi

keywords=`basename $program`
myself=`basename $0`

for i in `cat $list `
do
	while true
	do
		numRunningJobs=`ps -x | grep ${keywords} | grep -v ${myself} | wc -l`
                if [ $numRunningJobs -lt `expr $numAllowedJobs + 1 ` ]; then
        		python $program $TGTDir/$i.tgt $RESDIR &
			break
		else
			a=`expr $RANDOM % 3 `
                        sleep $a
                fi
        done

done

wait

