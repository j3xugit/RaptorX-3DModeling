#!/bin/bash

if [ $# -lt 3 ]; then
	echo "$0 proteinListFile TPLDir PDBDir [ResDir [numCPUs] ]"
	echo " numCPUs: the number of parallel jobs to be run, default 20"
	exit 1
fi

list=$1
TPLDIR=$2
PDBDIR=$3

RESDIR=TPLPKL_BC40
if [ $# -ge 4 ]; then
	RESDIR=$4
fi
mkdir -p $RESDIR

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

program=$cmdDir/TPL2PKL2.py

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
			pdbfile=$PDBDIR/$i.pdb
        		python $program $TPLDIR/$i.tpl ${pdbfile} $RESDIR &
			break
		else
			a=`expr $RANDOM % 3 `
                        sleep $a
                fi
        done

done

wait
