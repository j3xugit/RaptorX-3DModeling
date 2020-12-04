#!/bin/bash

if [ $# -lt 2 ]; then
	echo "$0 proteinListFile PDBDir [ResDir]"
	exit 1
fi

list=$1
if [ ! -f $list ]; then
	echo "ERROR: invalid protein list file $list"
	exit 1
fi

PDBDIR=$2
if [ ! -d $PDBDIR ]; then
	echo "ERROR: invalid folder for raw protein structure files"
	exit 1
fi

RESDIR=`pwd`
if [ $# -ge 3 ]; then
	RESDIR=$3
fi
if [ ! -d $RESDIR ]; then
	mkdir -p $RESDIR
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

program=$cmdDir/SelectProteinChain.sh

for i in `cat $list `
do
  	$program $i ${PDBDIR} $RESDIR
done
