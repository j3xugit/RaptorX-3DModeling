#!/bin/bash

if [ $# -lt 2 ]; then
	echo "$0 proteinList A3MDir [ResDir]"
	echo "	This script generates TGT files for a list of protein in serial but not in parallel"
	echo "	proteinList: a file containing a list of proteins, each in one row"
	echo "	A3MDir: the folder for MSAs in a3m format"
	echo "	ResDir: the folder for result files, default current work directory"
	exit 1
fi

proteinList=$1
if [ ! -f $proteinList ]; then
	echo "ERROR: invalid protein list file $proteinList"
	exit 1
fi

A3MDir=$2
if [ ! -d $A3MDir ]; then
	echo "ERROR: invalid folder for MSAs $A3MDir "
	exit 1
fi


cmd=`readlink -f $0`
cmdDir=`dirname $cmd`
program=$cmdDir/GenTGTFromA3M.sh
if [ ! -x $program ]; then
	echo "ERROR: invalid executable program $program"
	exit 1
fi

ResDir=`pwd`
if [ $# -ge 3 ]; then
	ResDir=$3
fi

if [ ! -d $ResDir ]; then
	mkdir -p $ResDir
fi

for target in `cat $proteinList `
do
	$program $A3MDir/$target.a3m $ResDir
done
