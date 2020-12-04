#!/bin/bash

if [ $# -lt 3 ]; then
	echo "$0 proteinListFile HHMDir PDBDir [ResDir]"
	exit 1
fi

list=$1
if [ ! -f $list ]; then
	echo "ERROR: invalid protein list file $list"
	exit 1
fi

HHMDIR=$2
if [ ! -d $HHMDIR ]; then
	echo "ERROR: invalid folder for .hhm files $HHMDIR"
	exit 1
fi

PDBDIR=$3
if [ ! -d $PDBDIR ]; then
	echo "ERROR: invalid folder for .pdb/.cif files $PDBDIR"
	exit 1
fi

RESDIR=`pwd`
if [ $# -ge 4 ]; then
	RESDIR=$4
fi
mkdir -p $RESDIR

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

program=$cmdDir/HHM2TPL.py

for i in `cat $list `
do
	pdbfile=${PDBDIR}/${i}.pdb
	if [ ! -f $pdbfile ]; then
		pdbfile=${PDBDIR}/${i}.cif
	fi
	if [ ! -f $pdbfile ]; then
		echo "WARNING: invalid structure file for $i"
		continue
	fi
        python $program $HHMDIR/$i.hhm ${pdbfile} $RESDIR
done

