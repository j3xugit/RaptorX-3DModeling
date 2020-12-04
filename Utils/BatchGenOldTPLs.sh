#!/bin/bash

## this script generates old tpl files for a list of proteins
if [ $# -lt 3 ]; then
	echo "$0 proteinListFile TGT_folder/A3M_folder StructureFolder [ResDir]"
	echo "	proteinListFile: a file for a list of proteins, each in one row"
	echo "	TGT_folder/A3M_folder: a folder for MSA files or TGT files"
	echo "	StructureFolder: a folder for protein structure files ending with .pdb or .cif"
	echo "	ResDir: the folder for result saving, default current work directory"
	exit 1
fi

proteinListFile=$1
MSAFolder=$2
StructFolder=$3
ResDir=`pwd`

if [ $# -ge 4 ]; then
	ResDir=$4
	mkdir -p $ResDir
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`

for target in `cat $proteinListFile`
do
	msafile=$MSAFolder/${target}.tgt
	if [ ! -f $msafile ]; then
		msafile=$MSAFolder/${target}.a3m
	fi
	if [ ! -f $msafile ]; then
		echo "ERROR: invalid input file for $target in $MSAFolder"
		exit 1
	fi

	structfile=$StructFolder/${target}.pdb
	if [ ! -f $structfile ]; then
		structfile=$StructFolder/${target}.cif
	fi
	if [ ! -f $structfile ]; then
		echo "ERROR: invalid structure file for $target in $StructFolder"
		exit 1
	fi

	$cmdDir/GenerateOldTPL.sh $msafile $structfile $ResDir
done
