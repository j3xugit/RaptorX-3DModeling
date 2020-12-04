#!/bin/bash

## this script generates a tpl file from an .hhm or .a3m file plus a pdb file
## the tpl file is named after proteinName.tpL.pkl, which is supposed to be compatible with the old TPL file

if [[ -z "$DistFeatureHome" ]]; then
        echo "ERROR: please set environmental variable DistFeatureHome to the instllation directory of BuildFeatures"
        exit 1
fi

if [ $# -lt 2 ]; then
	echo $0 "MSAfile pdbFile [ResultDir]"
	echo "	This script generates a .tpl.pkl file from an MSA and a structure file"
	echo "	MSAfile: an .hhm or .a3m file"
	echo "	pdbFile: a protein structure file"
	echo "	ResultDir: the folder for result saving, default current work directory"
	echo "	The result file is named after proteinName.tpl.pkl"
	exit 1
fi


MSAfile=$1
if [ ! -f $MSAfile ]; then
	echo 'ERROR: invalid MSA or HHM file: ' $MSAfile
	exit 1
fi

pdbFile=$2
if [ ! -f $pdbFile ]; then
	echo "ERROR: invalid protein structure file $pdbFile"
	exit 1
fi

ResultDir=`pwd`
if [ $# -ge 3 ]; then
	ResultDir=$3
	mkdir -p $ResultDir
fi

fulnam=`basename $MSAfile`
if [[ "$fulnam" == *.hhm ]]; then
	target=${fulnam%.hhm}
elif [[ "$fulnam" == *.a3m ]]; then
	target=${fulnam%.a3m}
else
	echo "ERROR: unsupproted file suffix for multiple sequence alignment: " $MSAfile
	exit 1
fi

hhmfile=${ResultDir}/${target}.hhm
if [[ "$fulnam" == *.hhm ]]; then
	cp $MSAfile $hhmfile
else
	$DistFeatureHome/util/hhmake -i $MSAfile -o $hhmfile
fi

cmd=`readlink -f $0 `
cmdDir=`dirname $cmd`

python $cmdDir/HHM2TPL.py $hhmfile $pdbFile $ResultDir

if [ ! -f $ResultDir/$target.tpl.pkl ]; then
	echo ERROR: failed to generate $ResultDir/$target.tpl.pkl
	exit 1
fi
