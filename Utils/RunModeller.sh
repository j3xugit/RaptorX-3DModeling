#!/bin/bash

if [ $# -lt 2 ]; then
	echo $0 aliFile PDBDir [ResDir]
	echo "	This script builds a 3D model from an aliFile in FASTA format"
	echo "	aliFile: an alignment file in FASTA format. The last entry is supposed to be the query sequence"
	echo "	PDBDir: the folder for PDB files of a template chain, e.g., pdb_BC100. It is better not to use the raw PDB files"
	echo "	The resultant 3D model will have name XXX.pdb where XXX is the basename of aliFile and saved to ResDir"
	exit 1
fi

aliFile=`readlink -f $1`
PDBDir=`readlink -f $2`

ResDir=`pwd`
if [ $# -ge 3 ]; then
	ResDir=$3
	mkdir -p $ResDir
fi

cmd=`readlink -f $0`
cmdDir=`dirname $cmd`
program=$cmdDir/ModellerUtils.py

if [ ! -f $program ]; then
	echo "ERROR: invalid program $program"
	exit 1
fi

bname=`basename $aliFile .fasta`
workDir=$(mktemp -d -t tmpModellerDir4${bname}-XXXXXXXXXX)
currDir=`pwd`

cd $workDir
python $program $aliFile $PDBDir
target=`echo $bname | cut -f1 -d'-' `
if [ ! -f $target.B99990001.pdb ]; then
	echo "ERROR: failed to generate a 3D model from $aliFile"
	exit 1
fi

pid=$$
mv $target.B99990001.pdb $bname.$pid.pdb
cd $currDir
mv $workDir/$bname.$pid.pdb $ResDir
rm -rf $workDir
