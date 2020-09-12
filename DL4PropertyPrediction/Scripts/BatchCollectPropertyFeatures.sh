#!/bin/bash

if [ $# -lt 2 ]; then
	echo $0 "proteinListFile MetaDir [resultDir [MSA_mode] ]"
	echo "	proteinListFile: a file for a list of proteins, each in one row"
	echo "	MetaDir: a folder containing a set of XXX_OUT/ where XXX is a protein name. XXX_OUT shall contain subfolders XXX_contact/ and XXX_thread/, which in turn shall contain .a3m or .hhm files"
	echo "	resultDir: a folder for the resultant feature file named as XXX.propertyFeatures.pkl, default current work directory"
	echo "	MSA_mode: which MSAs to be used, all (default), uce3, ure5 or user?"
	exit 1
fi

proteinListFile=$1
if [ ! -f $proteinListFile ]; then
	echo "ERROR: invalid protein list file $proteinListFile"
	exit 1
fi

RootDir=$2
if [ ! -d $RootDir ]; then
	echo "ERROR: invalid meta folder for collection of feature files: $RootDir"
	exit 1
fi

resDir=`pwd`
if [ $# -ge 3 ]; then
	resDir=$3
fi
if [ ! -d $resDir ]; then
	mkdir -p $resDir
fi

MSA="all"
if [ $# -ge 4 ]; then
	MSA=$4
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`
program=$cmdDir/CollectPropertyFeatures.sh

for target in `cat $proteinListFile`
do
	$program $target $RootDir/${target}_OUT $resDir $MSA
done

