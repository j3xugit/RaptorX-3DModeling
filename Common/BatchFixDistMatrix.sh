#!/bin/bash

if [ $# -lt 3 ]; then
	echo $0 proteinListFile OriginalGTDir PDBDir [ResDir]
	exit 1
fi

list=$1
OriginalGTDir=$2
PDBDIR=$3

RESDIR=`basename $OriginalGTDir`-revised
if [ $# -ge 4 ]; then
	RESDIR=$4
fi
if [ ! -d $RESDIR ]; then
	mkdir -p $RESDIR
fi

program=$HOME/3DModeling/Common/FixDistMatrix.py

numAllowedJobs=15
keywords=`basename $program`
myself=`basename $0`
	

for i in `cat $list `
do
	while true
	do
		## check the number of running jobs
                numRunningJobs=`ps -x | grep ${keywords} | grep -v ${myself} | wc -l`
                if [ $numRunningJobs -lt `expr $numAllowedJobs + 1 ` ]; then
			pdbfile=${PDBDIR}/${i}.pdb
                        ciffile=${PDBDIR}/${i}.cif
			if [ -f $pdbfile ]; then
				structfile=$pdbfile
			else
				structfile=$ciffile
			fi
			gtfile=$OriginalGTDir/${i}.native.pkl
			tplfile=$OriginalGTDir/${i}.tpl.pkl
			if [ -f $gtfile ]; then
				orgfile=$gtfile
			else
				orgfile=$tplfile
			fi

        		python $program $orgfile $structfile $RESDIR &
                        break
                else
                        a=`expr $RANDOM % 3 `
                        sleep $a
                fi
	done

done

