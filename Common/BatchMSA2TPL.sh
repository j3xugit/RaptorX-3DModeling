#!/bin/bash

if [ $# -lt 1 ]; then
	echo "$0 proteinListFile [MSADir [PDBDir [ResDir [numCPUs] ] ] ]"
	echo " numCPUs: the number of parallel jobs to be run, default 20"
	exit 1
fi

list=$1

MSADIR=MSA_BC100
if [ $# -ge 2 ]; then
	MSADIR=$2
fi

PDBDIR=pdb_BC100
if [ $# -ge 3 ]; then
	PDBDIR=$3
fi

RESDIR=TPLPKL_BC100
if [ $# -ge 4 ]; then
	RESDIR=$4
fi
mkdir -p $RESDIR

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

program=$cmdDir/MSA2TPL.sh

numAllowedJobs=20
if [ $# -ge 5 ]; then
	numAllowedJobs=$5
fi

keywords=`basename $program`
myself=`basename $0 `


for i in `cat $list `
do
	while true
	do
		numRunningJobs=`ps -x | grep ${keywords} | grep -v ${myself} | wc -l`
                if [ $numRunningJobs -lt `expr $numAllowedJobs + 1 ` ]; then
			for structurefile in ${PDBDIR}/${i}.*
			do
        			$program $MSADIR/$i.a3m $structurefile $RESDIR &
			done
			break
		else
			a=`expr $RANDOM % 3 `
                        sleep $a
                fi
        done

done

